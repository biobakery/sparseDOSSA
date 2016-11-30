# This is the core script for generating synthetic datasets. This will call each function in turn to generate the metadata, and the bugs, and combine them into the appropriate pcl files.
################### Methods

funcCalibrateRLNormToMicrobiome = function(
  ### Given a TSV file
  ### Parameters for distribution generatio are given
  ### To be used in matrix generation
  ### Estimated parameters include
  ### SD excluding zeros
  ### Mus excluding zeros
  ### Percent zeros
  ### The beta for estimating the SD given the mu
  ### The beta for estimating the percent zero given the mu
  ### All Mus, SD, grand MU, and grand SD are ready for rlnorm (have been measured from a logged (rlnorm)
  sCalibrationFile,
  ### File to be read in and used to calibrate constansts and relationships in the underlying data.
  fVerbose = FALSE
  ### Flag to turn on logging and pdf creation
){
  message( "start funcCalibrateRLNormToMicrobiome" )
  # Read in file
  message( "Reading file." )
  dfData = read.table( sCalibrationFile )
  row.names( dfData ) = dfData[[1]]
  dfData = dfData[ -1, -1 ]
  
  # Get read depth of the samples (From a tsv file, samples = rows)
  ldReadDepths = rowSums( dfData )
  #sapply( seq_len(nrow( dfData )), function( x ) sum( as.numeric( as.matrix( dfData )[ x, ] ) ) )
  
  # Get the vector of Exps, Mus and SD ignoring 0s
  # Logged will be used to make features, they can be directly used in the rlnorm function
  # Not logged will be used to later be logged and estimate a grand mean and grand SD for the initial distribution,
  # Given that every point in the larger vector is the mu of the feature vectors.
  # Which could then be used by rlnorm to generate a each feature if needed.
  # This is not needed here but the calculation of this value is useful. It is used as an initial value for synthetic creation
  # of the initial vector of feature mus.
  
  # Mean of the nonzero data
  vdExp = vector(length=ncol(dfData))
  # Mean of the logged nonzero data
  vdMu = vector(length=ncol(dfData))
  # Standard deviation of the logged nonzero data
  vdLogSD = vector(length=ncol(dfData))
  # Percent zeros in data
  vdPercentZero = vector(length=ncol(dfData))
  
  # Calculate parameters for each feature (column)
  for( iIndex in seq_len(ncol( dfData )) ){
    # Get the percent zero before removing zeros for other measurements
    vdCur = as.numeric( as.vector( as.matrix( dfData[ iIndex ] ) ) )
    vdPercentZero[iIndex] = mean( vdCur == 0 )
    
    # Measure expectation of the feature with zeros
    vdExp[iIndex] = mean( vdCur )
    
    # Remove zeros from data
    vdCur = vdCur[ which( vdCur != 0 ) ]
    
    #### Note
    #### rlnorm needs a mean and sd from a logged rlnorm distribution which would match this
    #### without further manipulation. The "mean" in the formula is actually not the expectation
    #### The expectation is e^mu+.5*sd^2 so this is always more than mu.
    
    # Log nonzero data
    vdLogCur = log( vdCur )
    vdLogSD[iIndex] = sd( vdLogCur )
    vdMu[iIndex] = funcGetMu( vdExp[iIndex], exp(sd( vdLogCur ) ) )
  }
  
  # Estimate the distribution parameters from the expectation vector
  # Includes the relationship between the grand mu and grand sd
  # The grand mu, grand expectation (logged) and the grand sd
  lParams = funcGenerateExpVectorParameters( vdExp, TRUE )
  
  ##### Get relationship between logSD and log(Exp)
  # Log to make a linear relationship
  # This is both values logged as SD is based on vdLogCur
  vdExpLog = log( vdExp )
  
  lmod = lm( vdLogSD ~ vdExpLog )
  dBetaSD = coef( lmod )[ "vdExpLog" ]
  dInterceptSD = coef(lmod)[ "(Intercept)" ]
  
  #### Percent Zero and Exp
  ### Estimated with a polynomial to the second degree
  #dBeta2Zero = NA
  dBetaZero = NA
  dInterceptZero = NA
  if( sum( vdPercentZero ) > 0 ){
    # Fit with logisitic regression
    # modified by bor
    lmodGlm = glm( vdPercentZero ~ vdExpLog, family = binomial )
    dInterceptZero = lmodGlm$coefficients[1]
    dBetaZero = lmodGlm$coefficients[2]
  }
  if( fVerbose ){
    # Indicate the relationships found
    message("The following relationships were found.")
    message("***Grand distribution***")
    message(paste("Expectation of feature expectations vs SD of feature expectations:",lParams$GrandLogSDBeta))
    message("***Feature distributions***")
    message("Log Exp (with zeros) vs Log SD (without zeros):")
    message(paste("Intercept=", dInterceptSD, ", Beta=", dBetaSD))
    message("Log Exp (with zeros) vs Percent Zeros:")
    message(paste("Intercept=", dInterceptZero, ", Beta=", dBetaZero, sep=""))
    
    # Percent zero
    if( !is.na( dBetaZero ) ){
      message(c(dBetaZero,dInterceptZero))
      #plot( vdExpLog, vdPercentZero, main = "Estimating Percent Zero no intercept" )
      #points( x = vdExpLog, y = funcEstimatePercentZero( vdExpLog, dBetaZero, dBeta2Zero, 0 ), col = 'violet' ) 
    }
  }
  
  message( "stop funcCalibrateRLNormToMicrobiome" )
  return( list( exp = vdExp, mu = vdMu, sd = exp( vdLogSD ), percentZero = vdPercentZero,
                dAverageReadDepth = mean( ldReadDepths ), iFeatureCount = ncol( dfData ) ) )
  
  ### When returning the grand Mu remember that you are returning the Mu that gives the expectation for the Mus
  ### given the rlnorm function so this is different than the mus measured in the logged distribution (rlnorm)
  ### exp: Not Logged expectation of distribution ( mean(x) )
  ### mu: The exponentiated logMu (calculated not measured)
  ### sd: The exponentiated logSD (based on the exp calculated without zeros)
  ### dSDBeta: Relationship between logSD and log(mean(x))
  ### percentZero: Percent of values which are zero
  ### dZeroBeta: Relationship between percent zero and log(mean(x))
  ### dGrandBeta: Relationship between the logSD and Exp
  ### dAverageReadDepth: The average read depth per samples
  ### iFeatureCount: The number of features measured
}


# 5 Tests 9/4/2013
funcEstimateFeatureSD = function(
  ### Estimate the SD given the Log Exp and parameters modeling the relationship between the Log Exp and the SD
  vdExpLog,
  ### The measured mean of feature values
  dBetaSD,
  ### The beta for the relationship between the mu and the SD
  dInterceptSD = 0
){
  return( dInterceptSD  + ( vdExpLog * dBetaSD ) )
  ### Returns the LogSD for a feature
}


# 5 Tests 10/22/2013
funcEstimateGrandSD = function(
  ### Estimate the grand SD given the grand Mu and parameters modeling the relationship between the grand mu and grand sd
  dExp,
  dBetaGrandSD
){
  return(dExp * dBetaGrandSD)
  ### Returns the LogSD for the expecation vector
}


# 6 Tests 10/22/2013
# Change to logistic regression method
# modified by bor
funcEstimatePercentZero = function(
  ### Estimate the percent zero given the logged expectation and parameters modeling the relationship between the log Exp and the percent zero
  vdExpLog,
  ### The measured mean of feature values
  dBetaZero,
  ### The beta for the relationship between the logged exp and the percent zero (for degree 1)
  # No need to use this as now using logisitic regression
  # dBeta2Zero,
  ### The beta for the relationship between the logged exp and the percent zero (for degree 2)
  dInterceptZero,
  ### The intercept for the polynomial relationship, must have a value
  dScale = 1
  ### A scale for the percent zero after being generated by the relationship with exp log.
){
  message("***Scale***")
  message(dScale)
  vdScales = exp( dInterceptZero  + ( vdExpLog * dBetaZero ) )/( 1 + exp( dInterceptZero  + ( vdExpLog * dBetaZero ) ) ) * dScale
  return( vdScales )
}


# 9 Test 10/22/2013
funcForceMinCountsInMinSamples = function(
  ### Force a feature to pass the requirement of having a certain minimal count in a minimal number of samples
  vdFeature,
  ### Vector of signal (integers).
  iMinNumberCounts = 0,
  ### Minimum number of counts for a sample to pass the filter
  iMinNumberSamples = 0
  ### Min number of samples to have the minimum number of counts
){
  # No need to perform function if there is nothing to adjust
  if( ( iMinNumberCounts + iMinNumberSamples ) == 0 ){ return( vdFeature ) }
  
  # Check to make sure there are enough non-zeros to add signal to for the min number of samples
  iNonZeroSamplesNeeded = min( 0, iMinNumberSamples - length( which( vdFeature == 0 ) ) )
  if(iNonZeroSamplesNeeded > 0)
  {
    # If more samples are needed, add them back in as the mean
    dSignalMean = round( mean( vdFeature[ which( vdFeature > 0 ) ] ) )
    viUpdate = funcSample( which( vdFeature == 0 ), iNonZeroSamplesNeeded )
    vdFeature[ viUpdate ] = dSignalMean
  }
  
  # Look to see how many measurements are above the min.
  # If there are not enough counts more then the min number of counts
  # then inflate the distribution up until there is enough counts.
  iNeededExtraValues = iMinNumberSamples - length( which( vdFeature >= iMinNumberCounts ) )
  
  # Probability for each sample to get a count is based on it's current percentage of samples (multinomial distribution)
  vdProbabilities = vdFeature / sum( vdFeature )
  
  # Add and remove counts to indices that are not zero
  viAddIndices = which( vdFeature > 0 )
  
  # While we need to add counts
  while( iNeededExtraValues > 0 )
  {
    # Index to add to, if multiple possibilities, select using the current percentage counts.
    iIndexAdd = viAddIndices
    if( length( iIndexAdd ) > 1 )
    {
      iIndexAdd = funcSample( viAddIndices, 1, prob = vdProbabilities[ viAddIndices ] )
    }
    
    vdFeature[ iIndexAdd ] = vdFeature[ iIndexAdd ] + 1
    iNeededExtraValues = iMinNumberSamples - length( which( vdFeature >= iMinNumberCounts ) )
  }
  return( vdFeature )
}


# 2 Tests 10/22/2013
funcGenerateExpVectorParameters = function(
  ### Get the point estimate for the relationship between the mu and sd
  vdExpectations,
  ### These are the mus of each of the features from the calibration file.
  ### These are untransformed (not logged) and the expectation with zeros
  ### This will not be a sparse/zero-inflated distribution given that it is the expectation of
  ### each feature.
  fVerbose = FALSE
  ### If pdfs should be made
){
  message("Start funcGenerateExpVectorParameters")
  
  #Point estimates of the distributions
  #The exp you should get if you use these SD and Mu directly in the rlnorm
  # Not logged
  dExp = mean(vdExpectations)
  #This gives the logSD value
  dLogSD = sd(log(vdExpectations))
  
  # Get the Mu that would generate the Exp given the SD
  # Already logged (logMu)
  dLogMu = log(funcGetMu(dExp, exp(dLogSD)))
  
  # The relationship between the Exp and the LogSD
  dGrandBeta = dLogSD / dExp
  
  return(list(GrandLogMu=dLogMu, GrandExp=dExp, GrandLogSD=dLogSD, dReadDepth=sum(vdExpectations), GrandLogSDBeta=dGrandBeta))
  ### GrandMu: The logMus of the distributions
  ### GrandExp: The expectation of the distributions mean(x)
  ### GrandSD: The logSDs of the distributions
  ### dReadDepth: The total average read depth of the samples
  ### GrandSDBeta: The relationship between the logSD and the Exp
}

# 31 Tests 9/4/2013
funcGenerateFeatureParameters = function(
  ### Generate the feature mu paramater vector and associated SD parameter, expectation, and percent zeros if needed
  # If a mu parameter vector is given of the size of the int_number_samples then pass through
  # Otherwise sample to that size with replacement.
  int_number_features,
  ### Number of features
  int_number_samples,
  ### Number of samples
  iMinNumberSamples,
  ### Minimum number of samples not zero
  iReadDepth,
  ### Simulated read depth
  vdExp = NA,
  ### Vector of expectation for the features (should add up to read depth)
  vdMu = NA,
  ### Vector of Mu parameters for the original expectation distribution (means of features) if not supplied, one will be generated by rlnorm
  vdSD = NA,
  ### The vector of SD parameters matching the vdMus.
  vdPercentZero = NA,
  ### The vector of percent zeros matching the vdMus
  lSDRel = NA,
  ### If vdSD is not given, SD will be generated by a relationship with Exp parameters using the parameters given here
  lPercentZeroRel = NA,
  ### If vdPercentZero is not given, vdPercentZero will be generated by a relationship with Exp using the parameters given here
  dBetaGrandSD = c_d$BetaGrandSD,
  ### The beta for the relationship between the mu of mus (of feature distributions) and the SD of mus (of feature distributions)
  fVerbose = FALSE
  ### Controls the plotting of graphic pdf logging (Default FALSE, TRUE indicates logging occurs)
){
  message("start funcGenerateFeatureParameters")
  
  if(all(is.na(vdExp)))
  {
    message("funcGenerateFeatureParameters: Generating vdExp Vector.")
    
    # Draw a vector of expectations
    # This will be the template distribution all bugs within a sample will be based on.
    # This allows to have structure within the data set of bugs more or less prevalent with a level of consistency
    # The mean of the distribution of bugs is derived from the max number of bug counts for a sample, making the
    # total number of bugs per sample the same within a small random margin
    lsParams = NULL
    if(c_f$FreezeSDGrandMu)
    {
      message("The Grand SD is frozen")
      lsParams = list(dBestExp=iReadDepth/int_number_features, dBestSD=1, dBestMu=iReadDepth/int_number_features, dBestDepth=iReadDepth, dTargetDepth=iReadDepth)
    } else {
      lsParams = funcGetParamsForReadDepth(dReadDepthPerFeature=iReadDepth/int_number_features, dBeta=dBetaGrandSD)
    }
    
    # Remember this is a vector of expectation not of log mu parameters.
    message( paste( "LogMu", lsParams$dLogMu, "LogSD", lsParams$dLogSD, "Threshold",(c_i$TimesSDIsOutlier*exp(lsParams$dLogSD))+exp(lsParams$dLogMu)))
    vdExp = funcTruncatedRLNorm( iNumberMeasurements=int_number_features,
                                 vdLogMean=lsParams$dLogMu,
                                 vdLogSD=lsParams$dLogSD,
                                 viThreshold=( c_i$TimesSDIsOutlier*exp( lsParams$dLogSD ) ) + exp( lsParams$dLogMu ) )
    ### Update the distribution to the sum (read depth) requested.
    ### Depending on how many features are requested this is more or less needed
    ### This is not needed at the limit with many features.
    vdExp = funcUpdateDistributionToSum( vdExp, iReadDepth )
    
    # Mus, SD, and percent zero depend on expectations. Since new expectations have been created these are reset to NA so they will be regenerated.
    vdMu = NA
    vdSD = NA
    vdPercentZero = NA
  }
  
  # Make sure there are no zeros in the expectation vector
  vdExp[ which( vdExp == 0 ) ] = c_d$ALittleMoreThanZero
  
  # Modified by bor
  if(all(is.na(vdSD)))
  {
    # Generate vector of SD based on mu since it is not known
    message("funcGenerateFeatureParameters: Generating vdSD Vector.")
    vdSD = exp(funcEstimateFeatureSD(log(vdExp), lSDRel$BetaSD, lSDRel$InterceptSD))
    
    # If there is only 1 feature then the feature can not vary because its expectation would
    # Also be the read depth that needs to be preserved. Given there is only 1 feature no other
    # Features will be available to compensate for fluctations.
    if(length(vdExp)==1){vdSD = 0}
    
    # Floor to close to 1 because less can not be logged
    viInvalid = which( vdSD < 1 )
    if( length( viInvalid ) > 0 )
    {
      message( paste( "funcGenerateFeatureParameters: Changing low SDs to a little more than 0. # occurences = ", length( viInvalid ) ) )
      vdSD[ viInvalid ] = 1+c_d$ALittleMoreThanZero
    }
    
    if(fVerbose)
    {
      vdVisExp = vdExp
      vdVisExp[which(vdVisExp==0)] = c_d$ALittleMoreThanZero
      vdVisSD = vdSD
      vdVisSD[which(vdVisSD==0)] = c_d$ALittleMoreThanZero
    }
  }
  # modified by bor
  if(all(is.na(vdMu)))
  {
    message("funcGenerateFeatureParameters: Generating vdMu Vector.")
    # We know the vdExp for each sample
    # We know the SD for each sample
    vdMu = sapply(seq_along(vdExp), function(x) funcGetMu(vdExp[x],vdSD[x]))
  }
  # modified by bor
  if(all(is.na(vdPercentZero)))
  {
    # Generate vector of percent zero based on exp since it is not known
    message("funcGenerateFeatureParameters: Generating vdPercentZero Vector.")
    vdPercentZero = funcEstimatePercentZero( log(vdExp), lPercentZeroRel$BetaZero, lPercentZeroRel$InterceptZero, lPercentZeroRel$Scale ) 
    message("***vdPercentZero***")
    message(summary(vdPercentZero))
    viLessThanZero = which(vdPercentZero < 0)
    if(length(viLessThanZero>0))
    {
      message(paste("funcGenerateFeatureParameters: Changing negative Percent Zeros to 0. # occurences = ",length(viLessThanZero)))
      vdPercentZero[vdPercentZero < 0] = 0
    }
  }
  
  message("stop funcGenerateFeatureParameters")
  
  # QC and contraints for percent zero
  # Make sure the percent zero passes the max
  # If there are not enough nonzeros, there is no signal to use.
  # This should be the max. Given there is a certain number of samples that have to have signal
  # The percentage of zeros must allow for those samples not to be zero and so restricts
  # the max the percent zero can be.
  dMaxPercent = 1 - ( iMinNumberSamples / int_number_samples )
  vdPercentZero[ which( vdPercentZero > dMaxPercent ) ] = dMaxPercent
  return( list( exp = vdExp, mu = vdMu, sd = vdSD, PercentZero = vdPercentZero ) )
  
  ### exp Vector of expectations (the expectation vector)
  ### mu Vector of mu (not logMu) associated with the vdExp (by index)
  ### sd Vector of sd (not logSD) associated with the vdExp (by index)
  ### PercentZero Vector of percent zeros (0-1) associated with the vdExp (by index)
}


func_generate_metadata = function(
  int_base_metadata_number,
  int_number_samples,
  ### Number of samples
  dMinLevelPercent
  ### The minimum percent of samples a level can have
){
  # Preallocate matrix
  mat_metadata = matrix(data=NA,nrow=(int_base_metadata_number*c_i$CountTypesOfMetadata),ncol=int_number_samples)
  
  # Used to report on metadata
  mtrxParameters = vector(length=1e4)
  mtrxParameters[[1]] = c_str$MetadataDetails
  imP = 2
  
  # Continous metadata, mean = 0 for all
  # Generating and padding the list of potential mean values
  li_mean_value_list = as.list(rep(0, int_base_metadata_number*2))
  # li_mean_value_list = list(runif(1, c_d$RunifMin, c_d$RunifMax),1,100)
  # if(length(li_mean_value_list) < (int_base_metadata_number*2))
  # {
  #  for (k in (length(li_mean_value_list)+1):(int_base_metadata_number*2))
  #  {
  #  li_mean_value_list = c(li_mean_value_list,runif(1, c_d$RunifMin, c_d$RunifMax))
  #  }
  # } else {
  #  li_mean_value_list = li_mean_value_list[1:(int_base_metadata_number*2)]
  # }
  
  # # generating the continuous metadata
  # # Should be random normal (is not bugs)
  for ( i in seq_along(li_mean_value_list) )
  {
    mat_metadata[i,] = rnorm(int_number_samples,mean=li_mean_value_list[[i]],sd=1)
    mtrxParameters[[imP]] = paste(c_str$Metadata, i, " ", c_str$Continuous, sep ="")
    imP = imP + 1
  }
  
  # generating the binary variables
  # Set up what the distributions for the different binary metadata will be
  prob_values_tmp = c( 0.5, runif( int_base_metadata_number-1, c_d$MinBinary, c_d$MaxBinary ) )
  li_probability_values_list = lapply( prob_values_tmp, function(x) c(x,1-x) )
  
  # if(length(li_probability_values_list)+1 < int_base_metadata_number)
  # {
  #   for (k in (length(li_probability_values_list)+1):int_base_metadata_number)
  #   {
  #     probability = runif(1, c_d$MinBinary, c_d$MaxBinary)
  #     li_probability_values_list = c(li_probability_values_list,list(c(probability, 1-probability)))
  #   }
  # } else {
  #   li_probability_values_list = li_probability_values_list[1: int_base_metadata_number]
  # }
  
  # Check to make sure minimal number of occurences is possible given the number of levels
  # for example, if the min percentage of occurences was originally 60%, it would not be possible to have 60% of samples
  # in both binary grouping, the max in this case would be 50%, the max for quarternery would be 25%
  # In this case the min value will be set to the max number of samples in a level in an even distribution * the percentage
  # So an original 60% samples for a binary case would be ceiling(60%*(#samples/2)*#samples) and a quarternery case 
  # would be ceiling(60%*(#samples/4)*#samples)
  iMinLevelBinaryCount = ceiling(int_number_samples*dMinLevelPercent)
  iMinLevelQuarterneryCount = iMinLevelBinaryCount
  if((dMinLevelPercent)>.25)
  {
    iMinLevelQuarterneryCount = ceiling(dMinLevelPercent*(int_number_samples/4)*int_number_samples)
    if((dMinLevelPercent)>.5)
    {
      iMinLevelBinaryCount = ceiling(dMinLevelPercent*(int_number_samples/2)*int_number_samples)
    }
  }
  
  # Draw from the binary values in binary_names given the previously generated binary distributions
  binary_names = c(1,2)
  for(i in seq_along(li_probability_values_list))
  {
    vsCurMetadata = funcSample( binary_names, size=int_number_samples, prob=li_probability_values_list[[i]], replace = TRUE)
    fMetadataFailed = !funcIsFactorMetadataValid(vsCurMetadata, iMinLevelBinaryCount)
    iBinaryMetadataLoop = 1
    while(fMetadataFailed)
    {
      vsCurMetadata = funcSample( binary_names, size=int_number_samples, prob=li_probability_values_list[[i]], replace = TRUE)
      fMetadataFailed = !funcIsFactorMetadataValid(vsCurMetadata, iMinLevelBinaryCount)
      iBinaryMetadataLoop = iBinaryMetadataLoop + 1
      if(iBinaryMetadataLoop > c_i$LoopingControlIncrement)
      {
        fMetadataFailed = FALSE
        message(paste("Suboptional metadata was created, did not pass quality control, is too imbalanced. Minimum level is preferred to be ", iMinLevelBinaryCount,"."))
      }
    }
    mat_metadata[i+int_base_metadata_number*2,] = vsCurMetadata
    mtrxParameters[[imP]] = paste(c_str$Metadata, i+int_base_metadata_number*2, " ", c_str$Factor, " ", paste(levels(as.factor(mat_metadata[i+int_base_metadata_number*2,])), collapse = " "), sep="")
    imP = imP + 1
  }
  
  # # generating the quarternary metadata with simple distributions and at least one uniform
  # # names of the feature choices
  li_list_of_distributions = lapply( 1:int_base_metadata_number, function(x) as.vector( rdirichlet(1, alpha=c(1,1,1,1)) ) )
  
  # if(length(li_list_of_distributions) < int_base_metadata_number)
  # {
  #  for (k in (length(li_list_of_distributions)+1):int_base_metadata_number)
  #  {
  # li_list_of_distributions = c(li_list_of_distributions, list(c(.25,.25,.25,.25)))
  #  }
  # } else {
  #  li_list_of_distributions = li_list_of_distributions[1:int_base_metadata_number]
  # }
  
  # # handpicked n-nomial distributions (metadata values)
  quarternary_names = c(1,2,3,4)
  
  for (i in seq_along(li_list_of_distributions))
  {
    vsCurMetadata = funcSample( quarternary_names, size=int_number_samples, prob=li_list_of_distributions[[i]], replace = TRUE)
    fMetadataFailed = !funcIsFactorMetadataValid(vsCurMetadata, iMinLevelQuarterneryCount)
    iBinaryMetadataLoop = 1
    while(fMetadataFailed)
    {
      vsCurMetadata = funcSample( quarternary_names, size=int_number_samples, prob=li_list_of_distributions[[i]], replace = TRUE)
      fMetadataFailed = !funcIsFactorMetadataValid(vsCurMetadata, iMinLevelQuarterneryCount)
      iBinaryMetadataLoop = iBinaryMetadataLoop + 1
      if(iBinaryMetadataLoop > c_i$LoopingControlIncrement)
      {
        fMetadataFailed = FALSE
        #message(paste("Suboptional metadata was created, did not pass quality control, is too imbalanced. Minimum level is preferred to be ", iMinLevelQuarterneryCount,"."))
        #message(vcCurMetadata)
      }
    }
    mat_metadata[i+3*int_base_metadata_number,] = vsCurMetadata
    
    mtrxParameters[[imP]] = paste(c_str$Metadata, i+3*int_base_metadata_number, " ", c_str$Factor, " ", paste(levels(as.factor(mat_metadata[i+3*int_base_metadata_number,])), collapse = " "), sep="")
    imP = imP + 1
  }
  length(mtrxParameters) = imP - 1
  return(list(mat_metadata=mat_metadata, mtrxParameters=mtrxParameters))
}

### Modified by ehs
func_generate_random_lognormal_matrix = function(
  int_number_features,
  ### Number of features
  int_number_samples,
  ### Number of samples,
  iMinNumberCounts,
  ### Minimum number of counts for a feature to be considered in a sample for the QC filtering
  iMinNumberSamples,
  ### Created bugs must have a minimum number of samples (iMinNumberSamples) that have a minimum number of counts (iMinNumberSamples)
  iReadDepth,
  ### Simulated read depth for sample creation
  vdExp = NA,
  ### The vector of expectations for each feature. If not provided one will be generated and vdMu, vdSD, and vdPercentZero will be reset
  vdMu = NA,
  ### Vector of Mu parameters for the original exp distribution (means of features) if not supplied, one will be generated
  vdPercentZero = NA,
  ### Vector of percent zero parameters for the original exp distribution (means of features) if not supplied, one will be generated
  vdSD = NA,
  ### Vector of SD parameters for the original exp distribution (means of features) if not supplied, one will be generated
  mdLogCorr = diag(int_number_features),
  ### The correlation matrix of the logged distribution; default is a identity matrix with dimension int_number_features
  fZeroInflate = TRUE,
  ### Turns off Zero inflation if FALSE (default TRUE, zero inflation turned on)
  lSDRel = NA,
  ### If vdSD is not given, SD will be generated by a relationship with Exp parameters using the parameters given here
  lPercentZeroRel = NA,
  ### If vdPercentZero is not given, vdPercentZero will be generated by a relationship with Exp using the parameters given here
  dBetaGrandSD = c_d$BetaGrandSD,
  ### The beta for the relationship between the mu of mus (of feature distributions) and the SD of mus (of feature distributions)
  funcUpdateData = NA,
  ### Function to update the bug matrix before read depth shuffling. If not provided, will not occur.
  fVerbose = FALSE
  ### Turns on logging (typically generates pdfs)
){
  message(paste("func_generate_random_lognormal_matrix START"))
  
  # Preallocating for speed
  mat_bugs = matrix(data=NA,nrow=int_number_features,ncol=int_number_samples)
  mat_bugs_basis = matrix(data=NA,nrow=int_number_features,ncol=int_number_samples)
  
  # Get the initial mu vector for generating features.
  lsInitialDistribution = funcGenerateFeatureParameters(
    int_number_features = int_number_features,
    int_number_samples  = int_number_samples,
    iMinNumberSamples   = iMinNumberSamples,
    iReadDepth          = iReadDepth,
    vdExp               = vdExp,
    vdMu                = vdMu,
    vdSD                = vdSD,
    vdPercentZero       = vdPercentZero,
    lSDRel              = lSDRel,
    lPercentZeroRel     = lPercentZeroRel,
    dBetaGrandSD        = dBetaGrandSD,
    fVerbose            = fVerbose
  )
  
  # Update the Mu, SD and Percent zero bugs and report on distributions
  vdMu = lsInitialDistribution[["mu"]]
  vdSD = lsInitialDistribution[["sd"]]
  vdPercentZero = lsInitialDistribution[["PercentZero"]]
  vdExp = lsInitialDistribution[["exp"]]
  
  # TODO remove
  if(c_f$FreezeSDFeatures)
  { 
    message("Feature SDs are frozen.")
    if(length(vdSD)>0)
    {
      vdMu = rep(iReadDepth/int_number_features, length(vdSD))
      vdSD = rep(1,length(vdSD))
    }
  }
  
  # Number of samples needed to have signal as a constraint
  iNumberSamples = min(int_number_samples, iMinNumberSamples)
  
  # Make features and assign feature samples to samples giving higher counts to lower read depth samples.
  message("func_generate_random_lognormal_matrix: START Making features")
  # Create features
  lFeatureDetails = funcMakeFeature(
    vdMu                = vdMu,
    vdSD                = vdSD,
    vdPercentZero       = vdPercentZero,
    mdLogCorr           = mdLogCorr,
    iNumberSamples     = int_number_samples,
    iMinNumberCounts   = iMinNumberCounts,
    iMinNumberSamples  = iMinNumberSamples,
    vdTruncateThreshold = (c_i$TimesSDIsOutlier*vdSD)+vdMu,
    fZeroInflate       = fZeroInflate,
    fVerbose           = fVerbose
  )
  
  # Need to transpose so that features are rows, columns are samples
  mat_bugs = t(lFeatureDetails$Feature)
  mat_bugs_basis = t(lFeatureDetails$Feature_base)
  
  
  # Shuffle back in removed signal, but only if there are no bug-bug correlations that would be messed up
  if (all(mdLogCorr[upper.tri(mdLogCorr)] == 0)){
    mat_bugs = funcShuffleMatrix(mtrxData=mat_bugs, iTargetReadDepth=iReadDepth)
  }
  
  # Round to counts
  # This round method does not allow value produced lower then the minimal value 
  # This allows control to make sure zeros are produced only if they are not equal to zero.
  # This is appropriate for a zero inflated model but not for a standard lognormal model.
  # So if zero inflation is used then this rounding function is needed otherwise a normal rounding can be performed by not setting the iMinValue.
  mat_bugs = funcRoundMatrix(mtrxData=mat_bugs, fZeroInflated=fZeroInflate)
  mat_basis = funcRoundMatrix(mtrxData=mat_bugs_basis, fZeroInflated=fZeroInflate)
  
  if(c_f$FreezeSDGrandMu ||c_f$FreezeSDFeatures||c_f$PrintLognormalMatrix)
  {
    #    message("bug counts")
    #    message(mat_bugs)
    message("Read Depth")
    message(colSums(mat_bugs))
    message("Average read depth")
    message(mean(colSums(mat_bugs)))
    message("Feature mean")
    message(funcGetRowMetric(mat_bugs,mean))
    message("vdExp")
    message(vdExp)
    message("Sum exp")
    message(sum(vdExp))
    message("Sum exp, should be read depth")
    message(iReadDepth)
    message("Feature mean summary")
    message(summary(funcGetRowMetric(mat_bugs,mean)))
  }
  
  # Truth table for relative abundance log normal data
  mtrxParameters = matrix(data=NA, nrow=6, ncol=1)
  mtrxParameters[1,1] = paste(c_str$SyntheticMicrobiome, c_str$Random, sep='')
  mtrxParameters[2,1] = paste(c_str$NumberOfFeatures, int_number_features)
  mtrxParameters[3,1] = paste(c_str$NumberOfSamples, int_number_samples)
  mtrxParameters[4,1] = paste(c_str$TotalSampleBugOccurrence, iReadDepth)
  mtrxParameters[5,1] = paste(c_str$NumberCounts, iMinNumberCounts)
  mtrxParameters[6,1] = paste(c_str$NumberSamples, iMinNumberSamples)
  
  # Truth table for basis data
  mtrxBasisParameters = matrix(data=NA,nrow=5,ncol=1)
  mtrxBasisParameters[1,1] = paste(c_str$SyntheticMicrobiomeBasis, c_str$Random, sep='')
  mtrxBasisParameters[2,1] = paste(c_str$NumberOfFeatures, int_number_features)
  mtrxBasisParameters[3,1] = paste(c_str$NumberOfSamples, int_number_samples)
  mtrxBasisParameters[4,1] = paste(c_str$NumberCounts, iMinNumberCounts)
  mtrxBasisParameters[5,1] = paste(c_str$NumberSamples, iMinNumberSamples)
  
  message("stop func_generate_random_lognormal_matrix")
  return(list(mat_bugs=mat_bugs, mtrxParameters=mtrxParameters, mat_basis = mat_basis, mtrxBasisParameters = mtrxBasisParameters))
  ### Returns a row major matrix of log-normal data.
}

func_generate_random_lognormal_with_outliers = function(
  ### Generates a random log normal distribution of data as a null matrix
  ### The option of using a matrix passed in as a parameter as the null matrix is provided (mtrxBugs)
  ### A percent of samples are given outliers based on the percent parameter
  int_number_features,
  ### Number of features
  int_number_samples,
  ### Number of samples
  iMinNumberCounts,
  ### Minimum number of counts for a feature to be considered in a sample for the QC filtering
  iMinNumberSamples,
  ### Created bugs must have a minimum number of samples (iMinNumberSamples) that have a minimum number of counts (iMinNumberSamples)
  dMaxPercentOutliers,
  ### The maximum percent of outliers to create in each sample (0 =< dPercent =< 1.0)
  dPercentSamples,
  ### Percent of samples given outliers (0 =< dPercent =< 1.0)
  mtrxBugs,
  ### Precalculated null matrix
  fVerbose = FALSE
){
  message("start func_generate_random_lognormal_with_outliers")
  
  # convert to dataframe
  mat_bugs_dataframe = as.data.frame(mtrxBugs)
  
  # Determine the number of samples
  iNumberSamples = round(int_number_samples*dPercentSamples)
  # Randomly select samples
  liIndices = funcSample( 1:int_number_samples,iNumberSamples,replace=FALSE)
  
  # Detemine the max number of outliers
  iMaxNumberOfOutliers = round(int_number_features*dMaxPercentOutliers)
  
  lviSwapped = vector('list',length=length(liIndices))
  
  if(iMaxNumberOfOutliers>0)
  {
    #for each column
    for(iSample in liIndices)
    {
      # Select the number of outliers for the sample
      iNumberOutliers = funcSample( 1:iMaxNumberOfOutliers,1,replace=FALSE)
      
      # sort the column
      dfSorted = mat_bugs_dataframe[order(mat_bugs_dataframe[,iSample]),]
      
      viAlreadySelected = vector(length=2*iNumberOutliers)
      ivA = 1
      
      iBufferForZeros = 0
      if(c_f$IgnoreZerosInOutliers)
      {
        iBufferForZeros = length(which(dfSorted[[iSample]]==0))
      }
      if((iNumberOutliers+iBufferForZeros)>length(dfSorted[[iSample]]))
      {
        iBufferForZeros = length(dfSorted[[iSample]])-iNumberOutliers
        message(paste("Allowing",iBufferForZeros,"zeros to be selected in outlier swapping so that all",iNumberOutliers,"swaps can be performed in sample",iSample))
      }
      
      for (iOutlier in 1:iNumberOutliers)
      {
        # Get the indices of all Min and max values ties
        iMaxValue = dfSorted[(int_number_features-(iOutlier-1)), iSample]
        iMinValue = dfSorted[iOutlier+iBufferForZeros, iSample]
        viMaxValueIndices = which(mat_bugs_dataframe[[iSample]]==iMaxValue)
        viMinValueIndices = which(mat_bugs_dataframe[[iSample]]==iMinValue)
        iRowMax = setdiff(viMaxValueIndices,viAlreadySelected)
        iRowMin = setdiff(viMinValueIndices,viAlreadySelected)
        
        # Sample from ties
        if(length(iRowMax) > 1)
        {
          iRowMax = funcSample( c(setdiff( viMaxValueIndices, viAlreadySelected ) ), size = 1 )
        }
        if(length(iRowMin) > 1)
        {
          iRowMin = funcSample( c( setdiff( viMinValueIndices, viAlreadySelected ) ), size = 1 )
        }
        
        mtrxBugs[c(iRowMin, iRowMax), iSample] = mtrxBugs[c(iRowMax, iRowMin), iSample]
        viAlreadySelected[c(ivA,ivA+1)] = c( iRowMax, iRowMin )
        ivA = ivA + 2
      }
      # Record which features were swapped with lists as features, basically a transpose
      lviSwapped[[iSample]] = viAlreadySelected
    }
  }
  
  # Count swaps
  iSwapCount = 0
  if(length(lviSwapped)>0)
  {
    for(iIndex in seq_along(lviSwapped))
    {
      iSwapCount = iSwapCount + length(lviSwapped[[iIndex]])
    }
  }
  
  # Truth table for log normal data
  mtrxParameters = vector(length=1e4)
  mtrxParameters[1:5] = c(paste(c_str$SyntheticMicrobiome, c_str$Outlier, sep=''), 
                          paste(c_str$NumberOfFeatures, int_number_features), 
                          paste(c_str$NumberOfSamples, int_number_samples), 
                          paste(c_str$PercentOutliers, dMaxPercentOutliers), 
                          paste(c_str$PercentSampleOutliers, dPercentSamples))
  imP = 6
  
  if(length(lviSwapped)>0)
  {
    for(iIndex in seq_along(lviSwapped))
    {
      for(iItemSwapped in lviSwapped[[iIndex]])
      {
        mtrxParameters[[imP]] = paste(c_str$OutlierParameter,paste(c_str$Feature,c_str$Outlier, iItemSwapped, sep="_"),c_str$SampleParameter,iIndex)
        imP = imP + 1
      }
    }
  }
  
  message("Stop func_generate_random_lognormal_with_outliers")
  # And return
  length( mtrxParameters ) = imP - 1
  return(list(mat_bugs=mtrxBugs, mtrxParameters=mtrxParameters))
}

func_generate_random_lognormal_with_multivariate_spikes = function(
  ### Spike in associations with 1 or more metadata
  int_number_features,
  ### Number of features
  int_number_samples,
  ### Number of samples
  iMinNumberCounts,
  ### Minimum number of counts for a feature to be considered in a sample for the QC filtering
  iMinNumberSamples,
  ### Created bugs must have a minimum number of samples (iMinNumberSamples ) that have a minimum number of counts (iMinNumberSamples)
  percent_spikes,
  ### The percent of features (bugs) to spike
  multiplier,
  ### Used to multiple the metadata before adding to a feature to strengthen the signal of the metadata
  metadata_matrix,
  ### Matrix of metadata (by row) to spike in
  multivariate_parameter,
  ### Number of metadata to spike in
  dMinLevelCountPercent,
  ### Minimum number of samples allowed to be part of the spiked in relationship
  mtrxBugs,
  ### Random log normal matrix
  fZeroInflated,
  ### True indicates it is a zero inflated model.
  lviFrozeMetadataIndices = NULL,
  ### If given, the method must select these specific metadata indicies.
  ### This allow selection of features to be the same when evaluating the multiplier
  liFrozeDataIndicies = NULL,
  ### If given, the method must select these specific data indicies.
  ### This allows the selection of features to be the same when evaluating the multiplier
  lsFrozeLevels = NULL,
  ### If given, the method must select these specific level indicies.
  ### This allows the slection of features to be the same when evaluating the multiplier
  fVerbose = FALSE
){
  message("start func_generate_random_lognormal_with_multivariate_spikes")
  
  # Tracks the bug of interest
  iIndexSpikedFeature = NA
  
  # Initialize frozen levels if need be
  if(is.null(lsFrozeLevels))
  {
    lsFrozeLevels = list()
  }
  
  # looping control
  iLoopingControl = min(c_i$LoopingControlIncrement, int_number_samples*int_number_features)
  
  # Min number of samples in the relationship
  iMinSamples = floor(dMinLevelCountPercent*int_number_samples)
  
  # Creating the truth file which contains true positive spike-ins
  strMult = paste(multiplier)
  strMult = sub(".","_", strMult, fixed=TRUE)
  m_parameter_rec = c(paste(paste(c_str$SyntheticMicrobiome, c_str$Spike,sep=""), "n", multivariate_parameter,"m", strMult, sep='_'))
  m_parameter_rec = c(m_parameter_rec, paste(c_str$NumberOfFeatures, int_number_features))
  m_parameter_rec = c(m_parameter_rec, paste(c_str$NumberOfSamples, int_number_samples))
  m_parameter_rec = c(m_parameter_rec, paste(c_str$PercentSpikes, percent_spikes))
  m_parameter_rec = c(m_parameter_rec, paste(c_str$Multiplier, multiplier))
  m_parameter_rec = c(m_parameter_rec, paste(c_str$MultiplierParameter, multivariate_parameter))
  m_parameter_rec = c(m_parameter_rec, paste(c_str$MinimumSamples, iMinSamples))
  
  # Features to be spiked to select from
  # Optionally sparsity of available features can be controlled for
  viRemainingFeatures = 1:int_number_features
  #TODO new hack code
  #  message("viRemainingFeatures 1")
  #  message(viRemainingFeatures)
  #  if( FALSE )
  #    c_dLowestPercentZeros = -1
  #    c_dGreatestPercentZeros = 0.15
  #    vdPercentZero = funcGetRowMetric( mtrxBugs, function(vdBug) length( vdBug[ vdBug == 0 ] ) / length( vdBug ) )
  #    message("vdPercentZero")
  #    message(vdPercentZero)
  #    viRemainingFeatures = intersect( which( vdPercentZero > c_dLowestPercentZeros ), which( vdPercentZero < c_dGreatestPercentZeros ) )
  #    message("viRemainingFeatures 2")
  #    message(viRemainingFeatures)
  
  
  # Holds the metadata selected for spikin
  lviMetadata = vector(length=1e4)
  ivM = 1
  # Holds the data selected for spikin
  liData = vector(length=1e4)
  ilD = 1
  
  # The number of bugs to spike in
  iSpikeinCount = floor(int_number_features*percent_spikes)
  
  if(iSpikeinCount > 0)
  {
    # For each spiked in bug
    for(iSpikedBug in 1:iSpikeinCount)
    {
      vsCurFrozeLevels = c()
      # Controls breaking the following while loop incase a solution can not be found.
      iSpikeInLoopControl = 1
      # Holds the currently selected metadata indicies
      viSelectedMetadata = NULL
      # Currently selected metadata names
      vstrSpikedMetadata = NA
      # Bug to spike with
      vdCurData = NA
      # Indicator if the spikin passed quality control
      fSpikeInFailed = TRUE
      # Current best failed run
      # List including 
      ## Metadata = matrix of metadata that has been prepped for the spikin (row major),
      ## MetadataNames = vector of strings for the name sof each row of metadata,
      ## SpikinBug = vector of doubles (bug measurements),
      ## Count = vector of non-zeros elements overlapping both the bugs and each metadata,
      ## if not using metadata at a level but using all levels, the overlap is given for each level.
      lxCurrentBestRun = list(Metadata = c(), MetadataNames = c(), SpikinBug = c(), BugIndex = NA, Count = -1)
      
      # Metadata indices that are selected for spikin
      viSelectedMetadata = c()
      
      # Find valid spike-in scenario or best choice
      while(fSpikeInFailed)
      {
        # Get the bug to attempt the association
        # If previously associations have been made
        # liFrozenDataIndices makes the same associations happen here
        # This is so if multiple multipliers are given
        # The different matrices show the differences given increased size of effect
        # Not difference driven by selecting different bugs
        if(!is.null(liFrozeDataIndicies) & length(lviFrozeMetadataIndices)>0)
        {
          iIndexSpikedFeature = liFrozeDataIndicies[[iSpikedBug]]
        } else {
          iIndexSpikedFeature = funcSample( viRemainingFeatures, 1 )
        }
        
        # Select which of the metadatum we will be using to scale
        if(!is.null(lviFrozeMetadataIndices) & length(lviFrozeMetadataIndices)>0)
        {
          viSelectedMetadata = lviFrozeMetadataIndices[[iSpikedBug]]
        } else {
          viSelectedMetadata = funcSample( seq_len(nrow( metadata_matrix )), multivariate_parameter, replace = FALSE )
        }
        
        # Get matrix of metadata and feature
        vdCurMetadata = metadata_matrix[viSelectedMetadata,]
        vdCurData = mtrxBugs[iIndexSpikedFeature,]
        
        # Attempt to spike in a new bug
        
        lsMetadataInfo = funcPrepareMetadata(viSelectedMetadata=viSelectedMetadata, vdCurMetadata=vdCurMetadata, vsFrozeLevels=unlist(lsFrozeLevels[iSpikedBug]))
        vdCurMetadata = lsMetadataInfo[["metadata"]]
        vstrSpikedMetadata = lsMetadataInfo[["names"]]
        vsCurFrozeLevels = lsMetadataInfo[["vsLevels"]]
        
        # Spike in new bug
        vdSpikedBug = funcSpikeNewBug(vdCurMetadata=vdCurMetadata, vdCurData=vdCurData, multiplier=multiplier, fZeroInflated=fZeroInflated)
        
        # Check to see if a failure occured
        lxQCInfo = funcQCSpikin(vdCurMetadata, vdSpikedBug, iMinSamples)
        
        # lxQCInfo has the slots PASS = Boolean, CommonCounts = vector of integers
        fSpikeInFailed = !lxQCInfo[["PASS"]]
        
        # Check to see if the looping must end and no success is given
        # Also if the spikein failed, make sure to update the best failed scenario
        iSpikeInLoopControl = iSpikeInLoopControl + 1
        if(fSpikeInFailed)
        {
          # Keep the best failed scenario so far.
          if(!is.null(lxQCInfo[["CommonCounts"]]) && ( sum(lxQCInfo[["CommonCounts"]]>iMinSamples) > lxCurrentBestRun[["Count"]]))
          {
            lxCurrentBestRun = list(Metadata = vdCurMetadata, MetadataNames = vstrSpikedMetadata, SpikinBug = vdSpikedBug, BugIndex = iIndexSpikedFeature,  Count = sum(lxQCInfo[["CommonCounts"]]>iMinSamples), MetadataIndices = viSelectedMetadata, Levels = vsCurFrozeLevels)
          }
          
          # If we have ran out of iteration, use the best failed scenario and indicate this.
          if(iSpikeInLoopControl > iLoopingControl)
          {
            # Reset the current spikin variables to the best scenario
            vdCurMetadata = lxCurrentBestRun[["Metadata"]]
            vstrSpikedMetadata = lxCurrentBestRun[["MetadataNames"]]
            vdSpikedBug = lxCurrentBestRun[["SpikinBug"]]
            iIndexSpikedFeature = lxCurrentBestRun[["BugIndex"]]
            viSelectedMetadata = lxCurrentBestRun[["MetadataIndices"]]
            if(length(lsFrozeLevels) < iSpikedBug)
            {
              lsFrozeLevels[iSpikedBug] = lxCurrentBestRun[["Levels"]]
            }
            message(paste("While spiking in a relationship between metadata and bug was not able to meet the minimal percentage of spiked-in samples for the relationship. Min sample = ", iMinSamples))
            message(paste("The following spike-in does not pass quality control: Bug=", lxCurrentBestRun$BugIndex," Metadata=",paste(lxCurrentBestRun$MetadataNames,collapse=","),"Multiplier=", multiplier,"Count=", multivariate_parameter))
            # Break while to use the best case scenario
            break
          }
        } else {
          lxCurrentBestRun = list(Metadata = vdCurMetadata, MetadataNames = vstrSpikedMetadata, SpikinBug = vdSpikedBug, BugIndex = iIndexSpikedFeature,  Count = sum(lxQCInfo[["CommonCounts"]]>iMinSamples), MetadataIndices = viSelectedMetadata, Levels = vsCurFrozeLevels)
        }
      }
      
      # If the spike-in was successful, update
      if(!is.na(lxCurrentBestRun[["BugIndex"]]))
      {
        # Update metadata and data indices for the features used for spikins
        lviMetadata[ivM:(ivM-1+length(viSelectedMetadata))] = viSelectedMetadata
        ivM = ivM + length(viSelectedMetadata)
        liData[[ilD]] = iIndexSpikedFeature
        ilD = ilD + 1
        
        # Update spike-in
        mtrxBugs[iIndexSpikedFeature,] = vdSpikedBug
        if(length(lsFrozeLevels) < iSpikedBug)
        {
          lsFrozeLevels[[iSpikedBug]] = vsCurFrozeLevels
        }
        
        # Update the truth table with which feature is spiked with which metadata
        m_parameter_rec = c(m_parameter_rec, paste(c_str$Feature,c_str$Spike,"n", multivariate_parameter,"m", strMult,iIndexSpikedFeature,sep='_'))
        m_parameter_rec = c(m_parameter_rec, vstrSpikedMetadata)
        
        # Remove bug from pool of potential bugs to spike.
        viRemainingFeatures = setdiff(viRemainingFeatures, iIndexSpikedFeature)
      }
    }
  }
  # Floor to counts (spike-ins will be real numbers)
  mtrxBugs = floor(mtrxBugs)
  
  message("stop func_generate_random_lognormal_with_multivariate_spikes")
  
  length(lviMetadata) = ivM - 1
  length(liData) = ilD - 1
  # Return normalized matrix and spike-in list
  return(list(mat_bugs=mtrxBugs, m_parameter_rec=m_parameter_rec, MetadataIndices=lviMetadata, DataIndices=liData, Levels=lsFrozeLevels))
}


# 3 Tests 10/22/2013
funcGetExp = function(
  ### Given the Mu and SD, get the expectation
  ### Will need to 
  # Expects the values not logged
  # Log is base exp(1)
  # Returns a value that target the mean(data) not the logged data.
  dMu,
  dSD
){
  return( exp( log( dMu, exp( 1 ) ) + 0.5 * ( log( dSD, exp( 1 ) )^2 ) ) )
}


# 3 Tests 10/22/2013
funcGetMu = function(
  ### Get the mu given SD and an expectation
  # Expects the Exp to be the mean(data) and SD to be the sd(log(data))
  # so sd(log(rlnorm(10000,3,2))) to calculate the value
  # output from these functions can be directly used in the rlnorm (is the logged mu)
  # Log is base exp(1)
  # To use these values in rlnorm you would log them
  dEx,
  dSD
){
  if( (log(dSD)^2 - log(dEx))>745 ) warning(paste("mu = e^(-x), but x =",round(log(dSD)^2 - log(dEx),2),"will give mu=0."))
  return( exp( -1 * ( ( ( log( dSD, exp( 1 ) )^2 ) /2 ) - log( dEx, exp( 1 ) ) ) ) )
}


# 2 Tests 10/222013
funcGetParamsForReadDepth = function(
  ### From the relationship between the grand mu and grand sd contained in dBetadGradSD
  ### A mu and sd is calculated that should satisfy the given read depth
  dReadDepthPerFeature,
  ### The read depth of interest in the original scale
  dBeta
  ### The beta describing the relationship
){
  # Get the associated SD from the EXP
  dLogSD = funcEstimateGrandSD( dReadDepthPerFeature, dBeta )
  
  # Get the associated mu
  dMu = funcGetMu( dReadDepthPerFeature, exp( dLogSD ) )
  if( dMu == 0 ) warning(paste("dMu is numerically 0.  dReadDepthPerFeature=",dReadDepthPerFeature,"and dBeta =",dBeta,"Try making one smaller."))
  
  return( list( dLogMu = log( dMu ), dLogSD = dLogSD ) )
  ### dLogMu: The logMu which, with the logSD, would give the read depth (exp)
  ### dLogSD: The logMu which, with the logSD, would give the read depth (exp)
}


# 5 Tests 10/22/2013
funcGetRowMetric = function(
  ### Perform function on vector or rows of a matrix
  lxValues,
  ### Values of metadata which may be discrete (in which after the _ is a number) or a numeric
  ### Method to perform on values
  funcMethod
){
  if( is.null( dim( lxValues )[1] ) )
  {
    return( funcMethod( lxValues ) )
  }
  return( apply( lxValues, 1, funcMethod ) )
}


# 3 Tests 10/22/2013
funcGetSD = function(
  ### Get standard deviation for vector or matrix
  lxValues
  ### matrix or vector
){
  if( is.null( dim( lxValues )[ 1 ] ) )
  {
    return( sd( lxValues, na.rm = TRUE ) )
  }
  return( apply( lxValues, 1, sd, na.rm = TRUE ) )
}


# 4 Tests 10/22/2013
funcIsFactorMetadataValid = function(
  ### Check to make sure a level of entrophy is in a discontinuous metadata so it has a minimal level to associate
  vxMetadata,
  ### A metadata that can be factor data
  iMin
  ### Minimum number of instances of a level
){
  vxFactorMetadata = as.factor( vxMetadata )
  lstrLevels = levels( vxFactorMetadata )
  for( strLevel in lstrLevels )
  {
    if( sum( vxFactorMetadata == strLevel ) < iMin )
    {
      return( FALSE )
    }
  }
  return( TRUE )
}


### 3 Test 10/22/2013 (could do more)
### modified by ehs
funcMakeFeature = function(
  ### Create a feature given parameters
  ### Uses a zero inflated model if necessary
  ### Enforces a min number samples with signal if needed.
  vdMu,
  ### Mu of the rnorm distribution associate with the rlnorm draw that will occur, this will be logged
  vdSD,
  ### SD of the rnorm distribution associate with the rlnorm draw that will occur, this will be logged
  vdPercentZero,
  ### Percent of zeros to add to the feature
  iNumberSamples,
  ### Number of measurements for the feature
  iMinNumberCounts,
  ### Minimum number of counts to be considered signal
  iMinNumberSamples,
  ### Minimum number of samples needed to have signal. If this is not fulfilled, signal will be generated and added.
  mdLogCorr = diag(length(vdSD)),
  ### The correlation matrix of the logged distribution; default is a identity matrix with dimension length(vdLogSD)
  vdTruncateThreshold = NA,
  ### Threshold to truncate the underlying distribution
  fZeroInflate = TRUE,
  ### If the feature should be zero inflated
  fVerbose = FALSE
  ### If pdf logging of feaure should occur
){
  # If not zero inflated
  if(!fZeroInflate){vdPercentZero = rep(0, length(vdMu))}
  
  # Check that vdLogMean and vdLogSD are same length
  if (length(vdMu) != length(vdSD)){
    stop("vdMu and vdSD must have equal length")
  }
  # Expectation of the features
  vdExpCal = sapply(seq_along(vdMu), function(i) funcGetExp(vdMu[i],vdSD[i]))
  
  # Generate features
  mdFeature_base = func_zero_inflate(vdLogMean = log(vdMu), vdPercentZeroInflated=vdPercentZero, int_number_samples=iNumberSamples, vdLogSD=log(vdSD), mdLogCorr=mdLogCorr, viThreshold=vdTruncateThreshold)
  
  # Update the distributions to the targeted expectations
  mdFeature = matrix(NA, ncol=ncol(mdFeature_base), nrow=nrow(mdFeature_base))
  for (k in seq_len(ncol(mdFeature))){
    mdFeature[, k] = funcUpdateDistributionToExpectation(vdFeatures=mdFeature_base[, k], dExp = vdExpCal[k] )
  }
  
  # Causes striation, update
  mdFeature = apply(mdFeature, 2, funcForceMinCountsInMinSamples, iMinNumberCounts = iMinNumberCounts, iMinNumberSamples = iMinNumberSamples)
  mdFeature_base = apply(mdFeature_base, 2, funcForceMinCountsInMinSamples, iMinNumberCounts = iMinNumberCounts, iMinNumberSamples = iMinNumberSamples )
  
  # Extra useful measurements, the true and expected means
  vdMean = apply(mdFeature, 2, mean)
  vdMean_base = apply(mdFeature_base, 2, mean)
  
  return(list(Feature = mdFeature, Feature_base = mdFeature_base, Exp=vdMean, ExpCal = vdExpCal, Exp_base = vdMean))
}


# 4 Tests 10/22/2013
funcNormalizeMicrobiome = function(
  ### Normalize the data assuming that the 
  ### Samples are columns and the features (bugs) are rows.
  ### Zero columns are preserved as zero.
  mtrxMicrobiome
  ### Matrix or dataframe microbiome
){
  mtrxRelAb = mtrxMicrobiome / rep( colSums( mtrxMicrobiome ), each = nrow( mtrxMicrobiome ) )
  #unlist( lapply( colSums( mtrxMicrobiome ), function( x ) rep( x, nrow( mtrxMicrobiome ) ) ) )
  viZero = which( apply( mtrxMicrobiome, 2, sum ) == 0 )
  mtrxRelAb[, viZero] = rep( 0, nrow( mtrxRelAb ) )
  return( mtrxRelAb )
}


# 4 Tests 10/22/2013
# Will fail if all are integer if using to detect metadata
# Returns a TRUE (1) for each int value so if
# the return vector's sum is 0, the data is not numeric
funcNumericIsInt = function(
  ### Tests to see if a numeric is an interger
  dValue
  ### The numeric value to check
){
  return( paste( as.integer( dValue ), sep = "" ) == paste( dValue, sep = "" ) )
}


### 16 Test 10/22/2013
funcPrepareMetadata = function(
  ### If the data is continuous just return it and the metadata name
  ### If the data is discontinuous and dummying return one randomly selected level
  ### and give it a name that reflects the level selected.
  viSelectedMetadata,
  ### Indices of metadata to prepare
  vdCurMetadata,
  ### Matrix (row major) of original metadata but only metadata to be spiked
  fDummyData = TRUE,
  vsFrozeLevels = NULL
  ### Holds the name and level of metadata that is discontinuous
  ### Allows one to force the selection of a specific level
){
  # Get the names of the metadata associated with this bug
  # If dummy is occuring this will be written over
  vstrSpikedMetadata = paste(c_str$Metadata, viSelectedMetadata, sep='')
  
  if(is.null(vsFrozeLevels))
  {
    vsFrozeLevels = c()
  }
  vsCurLevels = vector(length=1e4)
  ivC = 1
  
  # Dummy the data if need be
  if(fDummyData)
  {
    if(is.null(nrow(vdCurMetadata)))
    {
      vdCurMetadata = matrix(vdCurMetadata,nrow=1)
    }
    
    # Reset for levels
    vstrSpikedMetadata = vector(length=nrow(vdCurMetadata))
    for(iRowIndex in seq_len(nrow(vdCurMetadata)))
    {
      # Get metadata
      vxMetadata = as.vector(vdCurMetadata[iRowIndex,])
      # If the factor data is not numeric
      if(!sum(sapply(vxMetadata,funcNumericIsInt))==0)
      {
        # Get levels of the metadata
        vsLevels = levels(as.factor(vxMetadata))
        
        # Add factor metadata with levels
        # Select from any level except for the first level
        # This is because the first level is always the reference level
        # in comparisons and will never be significant
        # (the other value compared to it will be).
        strLevel = funcSample( vsLevels[ seq_along( vsLevels )[-1] ], size=1 )
        if(length(vsFrozeLevels)>0)
        {
          strLevel = vsFrozeLevels[iRowIndex]
        } else {
          vsCurLevels[ivC] = strLevel
          ivC = ivC + 1
        }
        # Make the selected level 2 and those levels not selected 1
        vxMetadata[which(vxMetadata != strLevel)]=0
        vxMetadata[which(vxMetadata == strLevel)]=2
        vxMetadata[which(vxMetadata == 0)]=1
        vdCurMetadata[iRowIndex,] = vxMetadata
        vstrSpikedMetadata[iRowIndex] = paste(c_str$Metadata,viSelectedMetadata[iRowIndex],"_",c_str$Level,"_",strLevel,sep="")
        
      } else {
        if(length(vsFrozeLevels)==0)
        {
          vsCurLevels[ivC] = NA
          ivC = ivC + 1
        }
        # Add continuous metadata
        vstrSpikedMetadata[iRowIndex] = paste(c_str$Metadata,viSelectedMetadata[iRowIndex],sep="")
      }
    }
  }
  
  length(vsCurLevels) = ivC - 1
  # Return which level was dummied
  if(!length(vsFrozeLevels)>0)
  {
    vsFrozeLevels = vsCurLevels
  }
  return(list(names=vstrSpikedMetadata, metadata=vdCurMetadata, vsLevels=vsFrozeLevels))
}


### 78 Tests 10/22/2013
funcQCSpikin = function(
  ### Check to make sure the minimal number of spiked-in samples are met
  ### Asks the question is there a minimal number of overlap for metadata features with data features
  ### Works for dummied metadata or not
  ### In level based metadata, checks if there is a minimal # of nonzero bug entries with the level of interest
  ### In factor data not reduced to a level, checks if each level has at least a min of non zero entries >= min # non-zeros / # metadata levels
  vdCurMetadata,
  ### The metadata that was spiked-in, could be a vector or matrices (row major) depending on if the spikin is multivariate or not.
  vdSpikedBug,
  ### The bug feature created with a spiked in relationship
  iMinSpikedSamples,
  ### Minimal number of samples to be spiked in to pass QC
  fDummyFactorData = TRUE
){
  if( is.null( vdSpikedBug ) ){ return( list( PASS = FALSE, CommonCounts = c() ) ) }
  
  # Get which bug samples are 0
  viNonZeroSpikedBugs = which( vdSpikedBug != 0 )
  
  # Get number of row (the matrix is row major)
  # This could be a matrix (if multiple metadata are spiked in)
  # or a vector ( for 1 metadata spikeins )
  dNumberMetadata = nrow( vdCurMetadata )
  if( is.null( dNumberMetadata ) )
  {
    dNumberMetadata = 1
    vdCurMetadata = matrix( vdCurMetadata, nrow = 1 )
  }
  
  # This is to return a pass or fail and the number of nonzero items so
  # even if one fails, the best choice can be made in a series of these calls
  fPass = TRUE
  viCommon = vector(length=1e4)
  ivC = 1
  for( iMetadata in 1:dNumberMetadata )
  {
    vCurrentMetadata = vdCurMetadata[ iMetadata, ]
    
    # Test to see if the data is continuous or not
    if( sum( sapply( vCurrentMetadata, funcNumericIsInt ) ) == 0 )
    {
      iCommon = intersect( viNonZeroSpikedBugs, which( vCurrentMetadata != 0 ) )
      if( length( iCommon ) < iMinSpikedSamples )
      {
        fPass = FALSE
      }
      viCommon[ivC] = length( iCommon )
      ivC = ivC + 1
      
    } else {
      # Check factor data
      # Dummied means there is going to be two values, 1s and 2s 2s being the level of interest.
      # So you are interested in those things not == 1 which is everything but that the level of interest.
      if( fDummyFactorData )
      {
        iCommon = intersect( viNonZeroSpikedBugs, which( vCurrentMetadata != 1 ) )
        if( length( iCommon ) < iMinSpikedSamples )
        {
          fPass = FALSE
        }
        viCommon[ivC] = length( iCommon )
        ivC = ivC + 1
        
      } else {
        # Not dummied
        # Check for each level
        # Store each level and then set diff like normal.
        # At the very end each level will be checked if Non dummying is used.
        vstrLevels = levels( as.factor( vCurrentMetadata ) )
        iMinSamplesPerLevel = ceiling( iMinSpikedSamples / length( vstrLevels ) )
        
        for( strLevel in vstrLevels )
        {
          iCommon = intersect( viNonZeroSpikedBugs, which( vCurrentMetadata == strLevel ) )
          if( length( iCommon ) < iMinSamplesPerLevel )
          {
            fPass = FALSE
          }
          viCommon[ivC] = length( iCommon )
          ivC = ivC + 1
        }
      }
    }
  }
  length(viCommon) = ivC - 1
  # All data / levels passes in this spiked relationship
  return( list( PASS = fPass, CommonCounts = viCommon ) )
}


### 4 Tests 10/22/2013
funcRoundMatrix = function(
  ### Round a sparse matrix. If a value is greater than 0 but less than the minimal number then set to the minimal number and then round.
  ### This keeps the current level of sparsity and does not add more zeros.
  mtrxData,
  ### Matrix of data values to round
  fZeroInflated = FALSE
  ### Indicates if zero inflation is used.
){
  if( fZeroInflated )
  {
    vdSubCount = intersect( which( mtrxData < 1 ), which( mtrxData > 0 ) )
    mtrxData[ vdSubCount ] = 1
  }
  return( round( mtrxData ) )
}


### 8 Tests 10/22/2013
funcSample = function(
  ### Safe sampling on vector data. Sample will sample from 1:N if only one value is in a vector instead of the lenght = 1 vector.
  ### This changes the behavior so only the value is returned or an error is thrown (depending on size and replace and length of x)
  ### if the a vector of N = 1 is given
  x,
  size,
  replace = FALSE,
  prob = NULL
){
  iLength = length(x)
  if( iLength == 0 )
  {
    stop( "funcSample:: Can not sample from length of 0 vector." )
  } else if ( iLength == 1 ) {
    if( size == 1 ){ return( x ) }
    if( replace ){ return( rep( x, size ) ) }
    stop( "funcSample:: Can not create a vector of size > 1 from 1 entry without replacement." )
  } else {
    return( sample( x = x, size = size, replace = replace, prob = prob ) )
  }
}


### 3 Test 10-22-2013
funcShuffleMatrix = function(
  ### Shuffle the martix to make read depth less variable
  mtrxData,
  ### Matrix of values (rows are features, columns are samples)
  iTargetReadDepth
){
  message("Start funcShuffleMatrix")
  
  # Get Read depths
  viReadDepths = colSums( mtrxData )
  dCurDeviance = sum( abs( viReadDepths - iTargetReadDepth ) )
  
  # Measures the increase in the configuration
  dPrevDeviance = dCurDeviance + 1
  
  # Hold the max and min sample and feature info so the last shuffle can be undone
  viMinSample = NA
  viMaxSample = NA
  viMaxFeature = NA
  
  while( dCurDeviance < dPrevDeviance )
  {
    # Get min and max
    viMinSample = which( viReadDepths == min( viReadDepths ) )
    if( length( viMinSample ) > 1 ){ viMinSample  = funcSample( viMinSample, 1 ) }
    viMaxSample = which( viReadDepths == max( viReadDepths ) )
    if( length( viMaxSample ) > 1 ){ viMaxSample  = funcSample( viMaxSample, 1 ) }
    
    # Get feature with max count from max sample
    viSample = mtrxData[,viMaxSample]
    viMaxFeature = which( viSample == max( viSample ) )
    if( length( viMaxFeature ) > 1 ){ viMaxFeature = funcSample( viMaxFeature, 1 ) }
    
    # Shuffle with in feature, switching the max and min samples
    iHold = mtrxData[ viMaxFeature, viMaxSample ]
    mtrxData[ viMaxFeature, viMaxSample ] = mtrxData[ viMaxFeature, viMinSample ]
    mtrxData[ viMaxFeature, viMinSample ] = iHold 
    
    # Update Read depths and deviance
    viReadDepths = colSums( mtrxData )
    dPrevDeviance = dCurDeviance
    dCurDeviance = sum( abs( viReadDepths - iTargetReadDepth ) )
  }
  
  # Undo last move
  if( ! all(is.na( viMinSample ) ) )
  {
    iHold = mtrxData[ viMaxFeature, viMaxSample ]
    mtrxData[ viMaxFeature, viMaxSample ] = mtrxData[ viMaxFeature, viMinSample ]
    mtrxData[ viMaxFeature, viMinSample ] = iHold 
  }
  
  return(mtrxData)
}


### 11 Tests 10/22/2013
funcSpikeNewBug = function(
  ### Combine a bug data and metadata with a certain multiplier to create a known positive association
  vdCurMetadata,
  ### Metadata to spike in
  vdCurData,
  ### Data (bug to use in the spike-in)
  multiplier,
  ### Multiplier to increase the signal of the metadata in the association.
  ### Larger values make metadata more of the signal in the spike-in.
  ### 1 makes an equal parts metadata and data relationship.
  fZeroInflated = TRUE
){
  # Get average and SD of each metadata
  metadata_average = funcGetRowMetric(vdCurMetadata,mean)
  metadata_sigma = funcGetSD(vdCurMetadata)
  
  # Get the average and SD of the feature (ignoring zeros if zero inflated)
  data_average = mean(vdCurData[which(vdCurData>0)])
  data_sigma = sd(vdCurData[which(vdCurData>0)])
  if(is.na(data_sigma)){return(NULL)}
  
  if(!fZeroInflated)
  {
    data_average = mean(vdCurData)
    data_sigma = sd(vdCurData)
  }
  
  metadata_sigma[ metadata_sigma == 0] = 1
  if( data_sigma == 0 ) { data_sigma = 1 }
  
  # Feature as occurence (0,1)
  # This is used so that the features with zero keep those measurements at zero
  # Only used in a zero inflated environment so everything is set to 1 otherwise which cancels it's affect
  liOccurenceFilter = vdCurData
  liOccurenceFilter[liOccurenceFilter>0]=1
  if(!fZeroInflated)
  {
    liOccurenceFilter[seq_along(liOccurenceFilter)] = 1
  }
  
  # Spike in the feature with the metadata (multivariate_parameter == 1 means 1 metadata spike-in for the feature)
  # Make the scaled metadata
  scaled_metadata = (vdCurMetadata - metadata_average)/metadata_sigma
  if(data_sigma!=0){ scaled_metadata = scaled_metadata * data_sigma }
  scaled_metadata = scaled_metadata + data_average
  
  if(!(is.null(nrow(vdCurMetadata))||(nrow(vdCurMetadata)==1)))
  {
    scaled_metadata = colSums(scaled_metadata)
  }
  
  # Spike in the metadata that is now scaled.
  # make the metadata sparse based on the feature ("and" function filter)
  vdSpikedBug = vdCurData + (multiplier*scaled_metadata*liOccurenceFilter)
  
  # Scale the spike in down so the bug count does not increase as a (whole) feature just because adding is happening
  iNumberRows = nrow( vdCurMetadata )
  if( is.null( iNumberRows ) ){ iNumberRows = 1 }
  vdSpikedBug = vdSpikedBug / ( ( iNumberRows * multiplier ) + 1 )
  
  vdSpikedBug[which(vdSpikedBug<0)]=0
  return(vdSpikedBug)
}


### 6 Tests 10-22-2013
### Modified by bor
### Modified by ehs
funcTruncatedRLNorm = function(
  ### Return draws from a random log normal distribution with a truncated tail so outliers are not introduced.
  iNumberMeasurements,
  ### The number of measurements for this distribution
  vdLogMean,
  ### The mean of the logged distribution
  vdLogSD,
  ### The SD of the logged distribution
  mdLogCorr = diag(length(vdLogSD)),
  ### The correlation matrix of the logged distribution; default is a identity matrix with dimension length(vdLogSD)
  viThreshold = NA
  ### The value used to define outliers. 
){
  # Check that vdLogMean and vdLogSD are same length
  if (length(vdLogMean) != length(vdLogSD)){
    stop("vdLogMean and vdLogSD must have equal length")
  }
  # Get truncated normal distribution with rows=samples, cols=features
  mdLogSD = diag(x=vdLogSD, nrow=length(vdLogSD))
  mdLogVar = mdLogSD %*% mdLogCorr %*% mdLogSD
  if (length(viThreshold) == 1 && is.na(viThreshold)){
    viThreshold <- rep(Inf, length(vdLogSD))
  }
  
  mdFeature <- exp(tmvtnorm::rtmvnorm(n=iNumberMeasurements,
                                      mean = vdLogMean,
                                      sigma = mdLogVar,
                                      upper = log(viThreshold),
                                      algorithm="gibbs"))
  
  #  if (ncol(mdFeature)==1) return(as.vector(mdFeature))
  #vdFeature = rlnorm( iNumberMeasurements, dLogMean, dLogSD )
  
  # If a value is given to truncate outliers
  #if( !is.na( iThreshold ) )
  #{
  #  viOutliers = which( vdFeature > iThreshold )
  
  # If there are outliers
  # Change the outliers to the mean of the draws before it is drawn (indices less than it's own)
  #if( length( vioutliers ) )
  #{
  #  for( iindex in vioutliers )
  #  {
  #    # indices of measurements before the outlier
  #    vitomean = 1:( iindex - 1 )
  # if the first value was an outlier then use the means of all the non outlier measurements
  #    if(iindex-1 == 0){ vitomean = setdiff( 1:length( vdfeature ), vioutliers ) }
  # sample two measurments and take the mean, replace the outlier with this mean
  #    vdfeature[ iindex ] = mean( vdfeature[ funcsample( vitomean, 2, replace = true ) ] )
  #  }
  #}
  #}
  #Truncate negatives to zero
  #vdFeature[ vdFeature < 0 ] = 0
  
  return( mdFeature )
  ### vdFeature: The signal (optionally truncated, log normal distribution)
}


### 4 Tests 10/22/2013
funcUpdateDistributionToExpectation = function(
  ### Updates a distribution to the mean while keeping the shape. If any expectation is less than 1, we make them equal to 1.
  vdFeatures,
  ### Distribution of counts
  dExp
  ### The expectation to which to update the distribution
){
  # If the expectation is less than 1 we can get features that are all zero
  # especially after rounding a not zero-inflated matrix
  # So here a control is put in to make sure the expectation is no less than 1
  #if( dExp < 1 ){ dExp = 1 }
  
  # Used to ignore zeros in these calculations
  vdUpdateIndices = which( vdFeatures > 0 )
  # Get the amount of change in signal needed
  dDifference = dExp - mean( vdFeatures )
  
  # If signal is needed to be changed then update or remove signal by count
  if( abs( dDifference ) > 0 )
  {
    # Number of counts to shift
    # Modified by bor
    dCounts = ceiling( abs( dDifference * length( vdFeatures ) ) )
    # Add to the distributio and shift it up
    if(dDifference > 0)
    {
      if( length( vdUpdateIndices ) > 1 )
      {
        vdUpdateIndices = funcSample( vdUpdateIndices, dCounts, replace = TRUE, prob = vdFeatures[ vdUpdateIndices ] / sum( vdFeatures[ vdUpdateIndices ] ) )
      }
      update.freq = table( vdUpdateIndices )
      vdFeatures[as.numeric(names(update.freq))] = as.numeric(update.freq) + vdFeatures[as.numeric(names(update.freq))]
      # for( iIndex in vdUpdateIndices )
      # {
      #   vdFeatures[ iIndex ] = vdFeatures[ iIndex ] + 1
      # }
    } else if( dDifference < 0 )
    {
      # Remove from distribution
      for( iIndex in 1:dCounts )
      {
        viGreaterThan1 = which( vdFeatures > 1 )
        if( length( viGreaterThan1 ) > 1 )
        {
          iUpdateIndex = funcSample( viGreaterThan1, 1, prob = vdFeatures[ viGreaterThan1 ] / sum( vdFeatures[ viGreaterThan1 ] ) )
          vdFeatures[ iUpdateIndex ] = vdFeatures[ iUpdateIndex ] - 1
        }
      }
    }
  }
  return( vdFeatures )
}


# 4 Tests 10-22-2013
funcUpdateDistributionToSum = function(
  ### Update the distribution to the sum (read depth) requested.
  ### If the current read depth is not as large as needed sample the difference using
  ### a multinomial model based on the current taxa means (values) (as proportions)
  ### If the current read depth is larger than needed
  ### Then remove values in the same manner as described above for adding values.
  ### Except that a count of 1 can not be removed. This would throw off the sparsity.
  vdDistribution,
  ### A vector of numeric data that will be forced to sum to iTargetSum (+- less than 1)
  iTargetSum
  ### The target sum to update the distribution to
){
  # Update vdExp to target
  dCurReadDepthDifference = iTargetSum - sum( vdDistribution )
  
  # Multinomial probability
  vdProbabilities = vdDistribution / sum( vdDistribution )
  
  # Need more counts for the sum
  if(dCurReadDepthDifference > 0)
  {
    # Use a multinomial model to select indices
    viSamplingIndices = funcSample( x = seq_along( vdDistribution ), size = abs( dCurReadDepthDifference ), replace = TRUE, prob = vdProbabilities )
    # Update indices
    sample.update = table( viSamplingIndices )
    vdDistribution[as.numeric(names(sample.update))] = vdDistribution[as.numeric(names(sample.update))] + as.numeric( sample.update )
    # for(iIndex in viSamplingIndices)
    # {
    #   vdDistribution[iIndex] = vdDistribution[iIndex]+1
    # }
    # Need less counts for sum
  } else if( dCurReadDepthDifference < 0 ) {
    for( iIndex in seq_len(round( abs( dCurReadDepthDifference ) )) )
    {
      # Find those entries that are greater than one (do not want to make zeros here)
      viGreaterThanOne = which( vdDistribution > 1 )
      if(length( viGreaterThanOne ) > 0)
      {
        # If there is only one entry that is not 1 then select it otherwise sample.
        iGreaterThanOne = viGreaterThanOne
        if( length( iGreaterThanOne ) > 1 ){ iGreaterThanOne = funcSample( x = iGreaterThanOne, size = 1, prob = vdProbabilities[viGreaterThanOne] ) }
        vdDistribution[ iGreaterThanOne ] = vdDistribution[ iGreaterThanOne ] - 1
      }
    }
  }
  return( vdDistribution )
}


# 9 Tests 10/22/2013
### Modified by ehs
func_zero_inflate = function(
  ### Create a zero inflated log normal distribution with a specified mean and percentage of zeros.
  ### If you want to get the original values of mu and sd used in rlnorm, use mean(log(func_zero_inflate()))
  ### and sd(log(func_zero_inflate()))
  vdLogMean,
  ### Mean of the distribution (logged)
  vdPercentZeroInflated,
  ### Percentage of return which is zero
  int_number_samples,
  ### The number of samples to create
  vdLogSD,
  ### The sd of the distribution (logged)
  mdLogCorr = diag(length(vdLogSD)),
  ### The correlation matrix of the logged distribution; default is a identity matrix with dimension length(vdLogSD)
  viThreshold = NA
  ### The threshold for outliers
){
  # Get feature given distribution parameters; returns matrix with rows=samples, cols=features
  mdFeature = funcTruncatedRLNorm( int_number_samples, vdLogMean, vdLogSD, mdLogCorr = mdLogCorr, viThreshold = viThreshold )
  
  # Zero inlate
  # modified by bor
  miZeroLocations = as.logical( sapply( vdPercentZeroInflated, rbinom, n=int_number_samples, size=1 ) )
  mdFeature[miZeroLocations] = 0
  
  if( any(colSums( mdFeature ) == 0 ) ){
    zero_cols <- which( colSums( mdFeature ) == 0 )
    one_rows  <- sample( nrow(mdFeature), length(zero_cols), replace = TRUE )
    mdFeature[cbind(one_rows[zero_cols],zero_cols)] = 1 
    # for ( k in seq_along( zero_cols ) ){
    #     mdFeature[one_rows[k], zero_cols[k]] = 1
    # }
    #    vdFeature[sample(1:length(vdFeature),1)] = 1
  }
  
  #  if (ncol(mdFeature)==1) return(as.vector(mdFeature))
  #viZeroLocations = funcSample( 1:int_number_samples, floor( int_number_samples * dPercentZeroInflated ), replace = FALSE )
  #if( length( viZeroLocations ) ){ vdFeature[ viZeroLocations ] = 0 }
  
  # Return zero-inflated truncated lognormal feature
  return( mdFeature )
}
