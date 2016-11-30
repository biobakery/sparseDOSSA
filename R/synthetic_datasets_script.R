sparseDOSSA = function(strNormalizedFileName = "SyntheticMicrobiome.pcl",
                       strCountFileName = "SyntheticMicrobiome-Counts.pcl",
                       parameter_filename = "SyntheticMicrobiomeParameterFile.txt",
                       bugs_to_spike = 0,
                       spikeFile = NA,
                       calibrate = NA,
                       datasetCount = 1,
                       read_depth = 8030,
                       number_features = 300,
                       bugBugCorr =  "0.5",
                       spikeCount = "1",
                       lefse_file = NULL,
                       percent_spiked = 0.03,
                       minLevelPercent = 0.1,
                       number_samples = 50,
                       max_percent_outliers = 0.05,
                       number_metadata = 5,
                       spikeStrength = "1.0",
                       seed =  NA,
                       percent_outlier_spikins = 0.05,
                       minOccurence =  0,
                       verbose =  TRUE,
                       minSample =  0,
                       scalePercentZeros = 1,
                       association_type =  "linear",
                       noZeroInflate =  FALSE,
                       noRunMetadata = FALSE,
                       runBugBug =  FALSE) {
  int_base_metadata_number = number_metadata
  if (int_base_metadata_number < 1)
    stop("Please provide the base number for metadata generation as 1 or greater.")
  
  int_number_features = number_features
  if (int_number_features < 1)
    stop("Please provide a number of features of atleast 1")
  
  strAssociationType = association_type
  if (!strAssociationType %in% c("linear", "rounded_linear"))
    stop("Only linear associations supported at this point.")
  
  int_number_samples = number_samples
  if (int_number_samples < 1)
    stop("Please provide a number of samples of atleast 1")
  
  dPercentOutliers = max_percent_outliers
  if ((dPercentOutliers > 1) |
      (dPercentOutliers < 0))
    stop("Please provide a percent outliers in the range of 0 to 1")
  
  dPercentOutlierSpikins = percent_outlier_spikins
  if ((dPercentOutlierSpikins > 1) |
      (dPercentOutlierSpikins < 0))
    stop("Please provide a percent spikins in the range of 0 to 1")
  
  iReadDepth = read_depth
  if (iReadDepth < max(int_number_features, int_number_samples))
    stop(
      "Please provide a read depth of atleast equal to feature size or sample size (which ever is larger)"
    )
  
  iNumAssociations = bugs_to_spike
  if (iNumAssociations < 0)
    stop(
      "Please provide a number of associations (bug-bug correlation) greater than or equal to 0"
    )
  
  iNumberDatasets = datasetCount
  if (iNumberDatasets < 1)
    stop("Please provide a number of datasets which is at least 1.")
  
  dPercentMultSpiked = percent_spiked
  if ((dPercentMultSpiked > 1) |
      (dPercentMultSpiked < 0))
    stop("Please provide a percent multivariate spike in the range of 0 to 1")
  
  strCalibrationFile = calibrate
  
  dMinLevelCountPercent = minLevelPercent
  if ((dMinLevelCountPercent > 1) |
      (dMinLevelCountPercent < 0))
    stop("Please provide a min level percent in the range of 0 to 1")
  
  dMinOccurenceCount = minOccurence
  if (dMinOccurenceCount < 0)
    stop("Please provide a min occurence greater than or equal to 0")
  
  dMinOccurenceSample = minSample
  if (dMinOccurenceSample < 0)
    stop("Please provide a min sample greater than or equal to 0")
  if (dMinOccurenceSample > int_number_samples)
  {
    dMinOccurenceSample = int_number_samples
    message(
      paste(
        "The min sample (for QC) was larger than the actual sample size, reset the min sample to the sample size, minSample is now equal to number_samples which is ",
        int_number_samples
      )
    )
  }
  
  if (scalePercentZeros < 0) {
    stop("Please provide a scale percent zero greater than 0.")
  }
  
  vdSpikeCount = as.integer(unlist(strsplit(spikeCount, ",")))
  
  vdSpikeStrength = as.double(unlist(strsplit(spikeStrength, ",")))
  
  if (length(vdSpikeCount) != length(vdSpikeStrength))
  {
    if (length(vdSpikeCount) == 1)
    {
      vdSpikeCount = rep(vdSpikeCount, length(vdSpikeStrength))
    } else if (length(vdSpikeStrength) == 1)
    {
      vdSpikeStrength = rep(vdSpikeStrength, length(vdSpikeCount))
    } else {
      stop(
        "Please make sure to either provide the same length of spike-in counts and spike-n strengths or provide one that is 1 value. These values will be paired in order or one will be repeated to pair with the other values."
      )
    }
  }
  
  strBugBugCorr = bugBugCorr
  if (!is.null(strBugBugCorr)) {
    vecBugBugCorr = as.numeric(strsplit(strBugBugCorr, ',')[[1]])
    if (any(abs(vecBugBugCorr) >= 1)) {
      stop("Correlation values must be between -1 and 1.")
    }
  } else {
    vecBugBugCorr <- 0
  }
  
  strSpikeFile <- spikeFile
  if (is.na(strSpikeFile) && iNumAssociations == 0) {
    warning(
      paste0(
        "number of associations = 0, and no spike file ",
        "specified; no bug-bug spike-ins will be done."
      )
    )
  }
  if (!is.na(strSpikeFile) && iNumAssociations > 0) {
    warning(
      paste0(
        "number of associations > 0 but spike file ",
        "specified; bug-bug spike-ins will be done according ",
        "to spike file."
      )
    )
  }
  
  
  mtrxSpikeConfig = cbind(vdSpikeCount, vdSpikeStrength)
  
  fVerbose     = verbose
  fZeroInflate = !noZeroInflate
  fRunMetadata = !noRunMetadata
  fRunBugBug   = runBugBug
  
  # seed the random number generator
  if (!is.na(seed))
  {
    set.seed(seed)
  }
  
  # List of associations and bugs
  vParametersAssociations = vector(length = 1e5)
  iPA = 1
  
  list_of_bugs = vector('list', length = 1e4)
  ilb = 1
  
  # Holds key words for the feature names of the microbiomes
  lsMicrobiomeKeys = vector(length = 1e4)
  
  # Default number of metadata
  number_metadata = 0
  
  if (fRunMetadata) {
    # generate the metadata
    lsMetadataInfo = func_generate_metadata(int_base_metadata_number,
                                            int_number_samples,
                                            dMinLevelCountPercent)
    mat_metadata =  lsMetadataInfo[["mat_metadata"]]
    metadata_parameters = lsMetadataInfo[["mtrxParameters"]]
    
    vParametersAssociations[iPA:(iPA + length(lsMetadataInfo[["mtrxParameters"]]) -
                                   1)] = lsMetadataInfo[["mtrxParameters"]]
    iPA = iPA + length(lsMetadataInfo[["mtrxParameters"]])
    
    # generate plain random lognormal bugs
    # Get the fitted values for calibrating rlnorm
    vdExp = NA
    vdMu = NA
    vdSD = NA
    vdPercentZero = NA
    #dSDBeta = c_dSDBeta
    #dBetaZero = c_dBetaZero
    #dGrandBeta = c_dBetaGrandSD
    
    message("Parameters BEFORE Calibration File")
    message(
      paste(
        "Length exp",
        NA,
        "Length vdMu",
        NA,
        "length vdSD",
        NA,
        "length vdPercentZero",
        NA,
        "Read depth",
        iReadDepth
      )
    )
    
    if (!is.na(strCalibrationFile) & (strCalibrationFile != "NA"))
    {
      # Get the fit for the data
      message("Calibrating...")
      lsFit = funcCalibrateRLNormToMicrobiome(strCalibrationFile, fVerbose)
      vdExp = lsFit[["exp"]]
      vdMu = lsFit[["mu"]]
      vdSD = lsFit[["sd"]]
      vdPercentZero = lsFit[["percentZero"]]
      #iReadDepth = lsFit[["dAverageReadDepth"]]
      int_number_features = lsFit[["iFeatureCount"]]
    }
    message(
      "Parameters AFTER Calibration File (if no calibration file is used, defaults are shown)"
    )
    message(
      paste(
        "Length exp",
        length(vdExp),
        "Length vdMu",
        length(vdMu),
        "length vdSD",
        length(vdSD),
        "length vdPercentZero",
        length(vdPercentZero),
        "Read depth",
        iReadDepth,
        "Feature Count",
        int_number_features
      )
    )
    
    mat_random_lognormal_bugs = func_generate_random_lognormal_matrix(
      int_number_features = int_number_features,
      int_number_samples  = int_number_samples,
      iMinNumberCounts    = dMinOccurenceCount,
      iMinNumberSamples   = dMinOccurenceSample,
      iReadDepth          = iReadDepth,
      vdExp               = vdExp,
      vdMu                = vdMu,
      vdPercentZero       = vdPercentZero,
      vdSD                = vdSD,
      fZeroInflate        = fZeroInflate,
      lSDRel              = list(
        BetaSD = c_d$SDBeta,
        InterceptSD = c_d$SDIntercept
      ),
      lPercentZeroRel     = list(
        InterceptZero = c_d$InterceptZero,
        BetaZero      = c_d$BetaZero,
        Scale         = scalePercentZeros
      ),
      dBetaGrandSD        = c_d$BetaGrandSD,
      fVerbose            = fVerbose
    )
    
    
    vParametersAssociations[iPA:(iPA - 1 + length(mat_random_lognormal_bugs[["mtrxParameters"]]))] = mat_random_lognormal_bugs[["mtrxParameters"]]
    iPA = iPA + length(mat_random_lognormal_bugs[["mtrxParameters"]])
    
    list_of_bugs[[ilb]] = mat_random_lognormal_bugs[["mat_bugs"]]
    ilb = ilb + 1
    
    # generate random lognormal with outliers
    mat_random_lognormal_outliers_bugs = func_generate_random_lognormal_with_outliers(
      int_number_features = int_number_features,
      int_number_samples  = int_number_samples,
      dMaxPercentOutliers = dPercentOutliers,
      dPercentSamples     = dPercentOutlierSpikins,
      mtrxBugs            = mat_random_lognormal_bugs[["mat_bugs"]],
      fVerbose            = fVerbose
    )
    vParametersAssociations[iPA:(iPA - 1 + length(mat_random_lognormal_outliers_bugs[["mtrxParameters"]]))] = mat_random_lognormal_outliers_bugs[["mtrxParameters"]]
    iPA = iPA + length(mat_random_lognormal_outliers_bugs[["mtrxParameters"]])
    
    list_of_bugs[[ilb]] = mat_random_lognormal_outliers_bugs[["mat_bugs"]]
    ilb = ilb + 1
    
    # Holds key words for the feature names of the microbiomes
    lsMicrobiomeKeys[[1]] = c_str$Random
    lsMicrobiomeKeys[[2]] = c_str$Outlier
    ilM = 3
    
    # There are 4 groups of metadata (2 continuous, binary, and quarternery)
    number_metadata = c_i$CountTypesOfMetadata * int_base_metadata_number
    
    # Hold the info to freeze the levels in the data
    llsLevels = vector('list', length = nrow(mtrxSpikeConfig))
    lliMetadata = vector('list', length = nrow(mtrxSpikeConfig))
    lliData = vector('list', length = nrow(mtrxSpikeConfig))
    
    # Generate random lognormal with varying amounts of spikes
    for (iIndex in seq_len(nrow(mtrxSpikeConfig)))
    {
      vMultConfig = mtrxSpikeConfig[iIndex,]
      iSpikeCount = vMultConfig[1]
      sKey = as.character(iSpikeCount)
      iSpikeStrength = vMultConfig[2]
      lsLevels = NULL
      liData = NULL
      lviMetadata = NULL
      
      if (sKey %in% names(llsLevels))
      {
        lsLevels = llsLevels[[sKey]]
        liData = lliData[[sKey]]
        lviMetadata = lliMetadata[[sKey]]
      }
      
      mat_random_lognormal_multivariate_spikes = func_generate_random_lognormal_with_multivariate_spikes(
        int_number_features     = int_number_features,
        int_number_samples      = int_number_samples,
        percent_spikes          = dPercentMultSpiked,
        multiplier              = iSpikeStrength,
        metadata_matrix         = mat_metadata,
        multivariate_parameter  = iSpikeCount,
        dMinLevelCountPercent   = dMinLevelCountPercent,
        mtrxBugs                = mat_random_lognormal_bugs[["mat_bugs"]],
        fZeroInflated           = fZeroInflate,
        lviFrozeMetadataIndices = lviMetadata,
        liFrozeDataIndicies     = liData,
        lsFrozeLevels           = lsLevels,
        fVerbose                = fVerbose
      )
      mat_random_lognormal_multivariate_spikes_bugs = mat_random_lognormal_multivariate_spikes[["mat_bugs"]]
      
      lliMetadata[[sKey]] = mat_random_lognormal_multivariate_spikes[["MetadataIndices"]]
      lliData[[sKey]] = mat_random_lognormal_multivariate_spikes[["DataIndices"]]
      llsLevels[[sKey]] = mat_random_lognormal_multivariate_spikes[["Levels"]]
      
      # generate known associations for random lognormal with spikes
      vParametersAssociations[iPA:(iPA - 1 + length(mat_random_lognormal_multivariate_spikes[["m_parameter_rec"]]))] = mat_random_lognormal_multivariate_spikes[["m_parameter_rec"]]
      iPA = iPA + length(mat_random_lognormal_multivariate_spikes[["m_parameter_rec"]])
      
      list_of_bugs[[ilb]] = mat_random_lognormal_multivariate_spikes_bugs
      ilb = ilb + 1
      
      lsMicrobiomeKeys[[ilM]] = paste(c_str$Spike,
                                      "n",
                                      iSpikeCount,
                                      "m",
                                      sub(".", "_", paste(iSpikeStrength), fixed = TRUE),
                                      sep = "_")
      ilM = ilM + 1
    }
  }
  
  if (fRunBugBug) {
    v_log_corr_values = vecBugBugCorr
    
    # Get the fitted values for calibrating rlnorm
    vdExp = NA
    vdMu = NA
    vdSD = NA
    vdPercentZero = NA
    
    message("Parameters for Bug-Bug spikes BEFORE Calibration File")
    message(
      paste(
        "Length exp",
        NA,
        "Length vdMu",
        NA,
        "length vdSD",
        NA,
        "length vdPercentZero",
        NA,
        "Read depth",
        iReadDepth
      )
    )
    
    if (!is.na(strCalibrationFile) & (strCalibrationFile != "NA"))
    {
      # Get the fit for the data
      message("Calibrating...")
      lsFit = funcCalibrateRLNormToMicrobiome(strCalibrationFile, fVerbose)
      vdExp = lsFit[["exp"]]
      vdMu = lsFit[["mu"]]
      vdSD = lsFit[["sd"]]
      vdPercentZero = lsFit[["percentZero"]]
      iReadDepth = lsFit[["dAverageReadDepth"]]
      int_number_features = lsFit[["iFeatureCount"]]
    }
    message(
      "Parameters for Bug-Bug spikes AFTER Calibration File (if no calibration file is used, defaults are shown)"
    )
    message(
      paste(
        "Length exp",
        length(vdExp),
        "Length vdMu",
        length(vdMu),
        "length vdSD",
        length(vdSD),
        "length vdPercentZero",
        length(vdPercentZero),
        "Read depth",
        iReadDepth,
        "Feature Count",
        int_number_features
      )
    )
    
    
    # Get the initial mu vector for generating features so that the datasets are iid
    lsInitialDistribution = funcGenerateFeatureParameters(
      int_number_features = int_number_features,
      int_number_samples  = int_number_samples,
      iMinNumberSamples   = dMinOccurenceSample,
      iReadDepth          = iReadDepth,
      vdExp               = vdExp,
      vdMu                = vdMu,
      vdSD                = vdSD,
      vdPercentZero       = vdPercentZero,
      lSDRel              = list(
        BetaSD = c_d$SDBeta,
        InterceptSD = c_d$SDIntercept
      ),
      lPercentZeroRel     = list(
        InterceptZero = c_d$InterceptZero,
        BetaZero      = c_d$BetaZero,
        Scale         = scalePercentZeros
      ),
      dBetaGrandSD        = c_d$BetaGrandSD,
      fVerbose            = fVerbose
    )
    
    
    
    # Update the Mu, SD and Percent zero bugs and report on distributions
    vdMu = lsInitialDistribution[["mu"]]
    vdSD = lsInitialDistribution[["sd"]]
    vdPercentZero = lsInitialDistribution[["PercentZero"]]
    vdExp = lsInitialDistribution[["exp"]]
    
    mtrxDistributionParameters = matrix(NA, nrow = 5, ncol = 1)
    mtrxDistributionParameters[1, 1] = paste(c_str$SyntheticMicrobiome,
                                             c_str$DistributionParameters,
                                             sep = "")
    mtrxDistributionParameters[2, 1] = paste(c_str$MuVector, toString(round(vdMu, 2)))
    mtrxDistributionParameters[3, 1] = paste(c_str$SDVector, toString(round(vdSD, 2)))
    mtrxDistributionParameters[4, 1] = paste(c_str$PercentZeroVector, toString(round(vdPercentZero, 2)))
    mtrxDistributionParameters[5, 1] = paste(c_str$ExpVector, toString(round(vdExp, 2)))
    vParametersAssociations[iPA:(iPA + 4)] = c(mtrxDistributionParameters)
    iPA = iPA + 5
    
    
    # Get the indices for the associations
    if (is.na(strSpikeFile)) {
      lsAssociations = func_get_log_corr_mat_from_num(
        num_features = int_number_features,
        num_spikes = iNumAssociations,
        log_corr_values = v_log_corr_values
      )
    } else {
      lsAssociations = func_get_log_corr_mat_from_file(num_features = int_number_features,
                                                       file_name = strSpikeFile)
    }
    
    # Flag to add the parameters only once, since the datasets are iid
    fAddedParameters    = FALSE
    fAddBasisParameters = FALSE
    
    for (iDataset in 1:iNumberDatasets) {
      # generate plain random lognormal bugs
      mat_random_lognormal_bugs = func_generate_random_lognormal_matrix(
        int_number_features = int_number_features,
        int_number_samples  = int_number_samples,
        iMinNumberCounts    = dMinOccurenceCount,
        iMinNumberSamples   = dMinOccurenceSample,
        iReadDepth          = iReadDepth,
        vdExp               = vdExp,
        vdMu                = vdMu,
        vdPercentZero       = vdPercentZero,
        vdSD                = vdSD,
        mdLogCorr           = lsAssociations[["mdLogCorr"]],
        fZeroInflate        = fZeroInflate,
        lSDRel              = list(
          BetaSD = c_d$SDBeta,
          InterceptSD = c_d$SDIntercept
        ),
        lPercentZeroRel     = list(
          InterceptZero = c_d$InterceptZero,
          BetaZero      = c_d$BetaZero,
          Scale         = scalePercentZeros
        ),
        dBetaGrandSD        = c_d$BetaGrandSD,
        fVerbose            = fVerbose
      )
      
      if (!fAddedParameters) {
        vParametersAssociations[iPA:(iPA - 1 + length(mat_random_lognormal_bugs[["mtrxBasisParameters"]]))] = mat_random_lognormal_bugs[["mtrxBasisParameters"]]
        vParametersAssociations[iPA + length(mat_random_lognormal_bugs[["mtrxBasisParameters"]])] = paste(c_str$NumberDatasets, iNumberDatasets)
        iPA = iPA + length(mat_random_lognormal_bugs[["mtrxBasisParameters"]]) + 1
        
        parameter_mat <- matrix(NA, nrow = 9, ncol = 1)
        
        parameter_mat[1, 1]  <- paste0(c_str$SyntheticMicrobiome,
                                       c_str$BugBugAssociations)
        parameter_mat[2, 1]  <- paste(c_str$NumberOfFeatures,
                                      int_number_features)
        parameter_mat[3, 1]  <- paste(c_str$NumberOfSamples,
                                      int_number_samples)
        parameter_mat[4, 1]  <-
          paste(c_str$NumberOfAssociations,
                lsAssociations[["strNumberOfAssociations"]])
        parameter_mat[5, 1]  <-
          paste(c_str$CorrDomainBugs,
                lsAssociations[["strNumberCorrDomainBugs"]])
        parameter_mat[6, 1]  <-
          paste(c_str$CorrRangeBugsIdx,
                lsAssociations[["strIdxCorrRangeBugs"]])
        parameter_mat[7, 1]  <-
          paste(c_str$CorrDomainBugsIdx,
                lsAssociations[["strIdxCorrDomainBugs"]])
        parameter_mat[8, 1] <-
          paste(c_str$DirOfAssociations,
                lsAssociations[["strDirAssociations"]])
        parameter_mat[9, 1] <-
          paste(c_str$LogCorrValues,
                lsAssociations[["strLogCorrValues"]])
        
        vParametersAssociations[iPA:(iPA + 8)] = c(parameter_mat)
        iPA = iPA + 9
        
        fAddedParameters = TRUE
      }
      list_of_bugs[[ilb]] = mat_random_lognormal_bugs[["mat_basis"]]
      ilb  = ilb + 1
      
      lsMicrobiomeKeys[[ilM]] = paste(c_str$BugBugAssociations,
                                      "a",
                                      iNumAssociations,
                                      "d",
                                      iDataset,
                                      sep = "_")
      ilM = ilM + 1
    }
  }
  
  # preallocate final pcl  matrix
  length(list_of_bugs) = ilb - 1
  length(lsMicrobiomeKeys) = ilM - 1
  
  final_matrix = matrix(
    data = NA,
    nrow = (number_metadata + int_number_features * length(list_of_bugs)) +
      1,
    ncol = (int_number_samples + 1)
  )
  final_matrix[1, 1] = '#SampleID'
  final_matrix[1, 2:(int_number_samples + 1)] = paste('Sample', 1:int_number_samples, sep =
                                                        '')
  if (number_metadata > 0) {
    final_matrix[2:(number_metadata + 1), 1] = paste(c_str$Metadata, 1:number_metadata, sep =
                                                       '')
    vdDim = dim(mat_metadata)
    mat_metadata[(floor(vdDim[1] / 2) + 1):vdDim[1], ] = paste("Group_", mat_metadata[(floor(vdDim[1] /
                                                                                               2) + 1):vdDim[1], ], sep = "")
    final_matrix[2:(number_metadata + 1), 2:(int_number_samples + 1)] = mat_metadata
  }
  
  # Make a matrix for counts (use the other for normalized)
  mtrxFinalCounts = final_matrix
  
  start = 2 + number_metadata
  end = (2 + number_metadata) + (int_number_features - 1)
  iFirstData = start
  
  for (iMatIndex in seq_along(list_of_bugs))
  {
    final_matrix[start:end, 1] = paste(paste(c_str$Feature, lsMicrobiomeKeys[iMatIndex], sep =
                                               "_"),
                                       1:int_number_features,
                                       sep = '_')
    # Normalize each column by thier sum and add to output
    final_matrix[start:end, 2:(int_number_samples + 1)] = funcNormalizeMicrobiome(list_of_bugs[[iMatIndex]])
    mtrxFinalCounts[start:end, 1] = paste(paste(c_str$Feature, lsMicrobiomeKeys[iMatIndex], sep =
                                                  "_"),
                                          1:int_number_features,
                                          sep = '_')
    mtrxFinalCounts[start:end, 2:(int_number_samples + 1)] = list_of_bugs[[iMatIndex]]
    start = start + int_number_features
    end = end + int_number_features
  }
  
  # Write the table as normalized counts
  write.table(
    final_matrix,
    file = strNormalizedFileName,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = '\t'
  )
  
  # Write the tables as counts
  mtrxCounts = final_matrix
  write.table(
    mtrxFinalCounts,
    file = strCountFileName,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = '\t'
  )
  
  # Write truth tables
  length(vParametersAssociations) = iPA - 1
  write.table(
    as.matrix(vParametersAssociations),
    file = parameter_filename,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = '\t'
  )
  
  return(
    output.files = list(
      pcl_counts = strCountFileName,
      pcl_normalized = strNormalizedFileName,
      truth_file = parameter_filename
    )
  )
}