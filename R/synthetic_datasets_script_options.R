#!/usr/bin/env Rscript
#######################################################################################
# This file is provided under the Creative Commons Attribution 3.0 license.
#
# You are free to share, copy, distribute, transmit, or adapt this work
# PROVIDED THAT you attribute the work to the authors listed below.
# For more information, please see the following web page:
# http://creativecommons.org/licenses/by/3.0/
#
# This file is a component of the MaAsLin (Multivariate Associations Using Linear Models), 
# authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Timothy Tickle, ttickle@hsph.harvard.edu).
#######################################################################################

library(optparse)

# Defaults for options (these will also be used as defaults to function call)
option_default = list(
    strNormalizedFileName = "SyntheticMicrobiome.pcl",
    strCountFileName = "SyntheticMicrobiome-Counts.pcl",
    parameter_filename = "SyntheticMicrobiomeParameterFile.txt",
    bugs_to_spike = 0,
    calibrate = NA,
    datasetCount = 1,
    read_depth = 8030,
    number_features = 300,
    bugBugCorr = "0.5",
    spikeCount = "1",
    lefse_file = NULL,
    percent_spiked = 0.03,
    minLevelPercent = 0.1,
    max_domain_bugs = 2,
    number_samples = 50,
    max_percent_outliers = 0.05,
    number_metadata = 5,
    spikeStrength = "1.0",
    seed = NA,
    percent_outlier_spikins = 0.05,
    minOccurence = 0,
    verbose = TRUE,
    minSample = 0,
    scalePercentZeros = 1,
    association_type = "linear",
    spikeFile = NULL,
    noZeroInflate = FALSE,
    noRunMetadata = FALSE,
    runBugBug = FALSE)

option_list = list(
      make_option(
          c("-b","--bugs_to_spike"),
          type="integer",
          default= option_default[['bugs_to_spike']],
          help="Number of bugs to correlate with others.  A non-negative integer value is expected."
          ),
      make_option(
          c("-c","--calibrate"),
          type="character",
          default=option_default[['calibrate']],
          help="Calibration file for generating the random log normal data. TSV file (column = feature)"
          ),
      make_option(
          c("-d", "--datasetCount"),
          type="integer",
          default = option_default[['datasetCount']],
          help="The number of bug-bug spiked datasets to generate.  A positive integer value is expected."
          ),
      make_option(
          c("-e","--read_depth"),
          type="integer",
          default=option_default[['read_depth']],
          help="Simulated read depth for counts. A positive integer value is expected."
          ),
      make_option(
          c("-f","--number_features"),
          type="integer",
          default=option_default[['number_features']],
          help="The number of features per sample to create. A positive integer value is expected."
          ),
      make_option(
          c("-g","--bugBugCorr"),
          type="character",
          default = option_default[['bugBugCorr']],
          help = paste(
              "A vector of string separated values for the correlation ",
              "values of the pairwise bug-bug associations. This ",
              "is the correlation of the log-counts. ",
              "Values are comma-separated; for example: ",
              "0.7,0.5. Default is %default",
              sep=""
              )
          ),
      make_option(
          c("-i","--spikeCount"),
          type = "character",
          default = option_default[['spikeCount']],
          help = paste(
              "Counts of spiked metadata used in the spike-in dataset",
              "These values should be comma delimited values, in the order of the spikeStrength values (if given)",
              "Can be one value, in this case the value will be repeated to pair with the spikeCount values (if multiple are present)",
              "Example 1,2,3",
              sep = ". "
              )
          ),
      make_option(
          c("-j","--lefse_file"),
          type="character",
          default=option_default[['lefse_file']],
          help="Folder containing lefSe inputs."
          ),
      make_option(
          c("-k","--percent_spiked"),
          type="double",
          default=option_default[['percent_spiked']],
          help="The percent of features spiked-in. A real number between 0 and 1 is expected."
          ),
      make_option(
          c("-l","--minLevelPercent"),
          type="double",
          default=option_default[['minLevelPercent']],
          help=paste(
              "Minimum percent of measurements out of the total a level can have in a discontinuous metadata (Rounded up to the nearest count)",
              "A real number between 0 and 1 is expected.",
              sep = ". "
              )
          ),
      make_option(
          c("-n","--number_samples"),
          type="integer",
          default=option_default[['number_samples']],
          help="The number of samples to generate. A positive integer greater than 0 is expected."
          ),
      make_option(
          c("-o","--max_percent_outliers"),
          type="double",
          default=option_default[['max_percent_outliers']],
          help="The maximum percent of outliers to spike into a sample. A real number between 0 and 1 is expected."
          ),
      make_option(
          c("-p","--number_metadata"),
          type="integer",
          default=option_default[['number_metadata']],
          dest='number_metadata',
          help=paste(
              "Indicates how many metadata are created",
              "number_metadata*2 = number continuous metadata, number_metadata = number binary metadata, number_metadata = number quaternary metadata",
              "A positive integer greater than 0 is expected.",
              sep = ". "
              )
          ),
      make_option(
          c("-r","--spikeStrength"),
          type = "character",
          default = option_default[['spikeStrength']],
          help = paste(
              "Strength of the metadata association with the spiked-in feature",
              "These values should be comma delimited and in the order of the spikeCount values (if given)",
              "Can be one value, in this case the value wil be repeated to pair with the spikeStrength values (if multiple are present)",
              "Example 0.2,0.3,0.4.",
              sep = ". "
              )
          ),
      make_option(
          c("-s","--seed"),
          type="integer",
          default=option_default[['seed']],
          help=paste(
              "A seed to freeze the random generation of counts/relative abundance",
              "If left as default (NA), generation is random",
              "If seeded, data generation will be random within a run but identical if ran again under the same settings",
              "An integer is expected.",
              sep = ". "
              )
          ),
      make_option(
          c("-t","--percent_outlier_spikins"),
          type="double",
          default=option_default[['percent_outlier_spikins']],
          help="The percent of samples to spike in outliers. A real number between 0 to 1 is expected."
          ),
      make_option(
          c("-u","--minOccurence"),
          type="integer",
          default=option_default[['minOccurence']],
          help=paste(
              "Minimum counts a bug can have for the ocurrence quality control filter used when creating bugs",
              "( Filtering minimum number of counts in a minimum number of samples)",
              "A positive integer is expected.",
              sep = ". "
              )
          ),
      make_option(
          c("-v","--verbose"),
          action="store_false",
          default = option_default[['verbose']],
          help=paste(
              "If True logging and plotting is made by the underlying methodology",
              "This is a flag, it is either included or not included in the commandline, no value needed.",
              sep = ". "
              )
          ),
      make_option(
          c("-w","--minSample"),
          type="integer",
          default=option_default[['minSample']],
          help=paste(
              "Minimum samples a bug can be in for the ocurrence quality control filter used when creating bugs",
              "( Filtering minimum number of counts in a minimum number of samples)",
              "A positive integer is expected.",
              sep = ". "
              )
          ),
      make_option(
          c("-x","--scalePercentZeros"),
          type="double",
          default=option_default[['scalePercentZeros']],
          help=paste(
              "A scale used to multiply the percent zeros of all features across the sample after it is derived from the relatiohships with it and the feature abundance or calibration file",
              "Requires a number greater than 0",
              "A number greater than 1 increases sparsity, a number less than 1 decreases sparsity",
              "O removes sparsity, 1 (default) does not change the value and the value.",
              sep = ". "
              )
          ),
      make_option(
          c("-y","--association_type"),
          type="character",
          default = option_default[['association_type']],
          help=paste(
              "The type of association to generate",
              "Options are 'linear' or 'rounded_linear'.",
              sep = ". "
              )
          ),
      make_option(
          c("-z","--noZeroInflate"),
          action="store_true",
          default = option_default[['noZeroInflate']],
          help=paste(
              "If given, zero inflation is not used when generating a feature",
              "This is a flag, it is either included or not included in the commandline, no value needed.",
              sep = ". "
              )
          ),
      make_option(
          c("--noRunMetadata"),
          action="store_true",
          default = option_default[['noRunMetadata']],
          help=paste(
              "If given, no metadata files are generated",
              "This is a flag, it is either included or not included in the commandline, no value needed.",
              sep = ". "
              )
          ),
      make_option(
          c("--runBugBug"),
          action="store_true",
          default = option_default[['runBugBug']],
          help=paste(
              "If given, bug-bug interaction files are generated in addition to any metadata files",
              "This is a flag, it is either included or not included in the commandline, no value needed.",
              sep = ". "
              )
          ),
      make_option(
          c("--spikeFile"),
          type="character",
          default = option_default[["spikeFile"]],
          help=paste0("If specified, a file that delineates all the ",
              "pairwise correlations and correlated feature pairs. The ",
              "file must have columns `Domain`, `Range`, and ",
              "`Correlations`. This option overrides `bugs_to_spike` and ",
              "`bugbugCorr`.")
          )
    )
