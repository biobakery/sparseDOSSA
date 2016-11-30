# Integer constants
c_i <- list(
  CountTypesOfMetadata = 4,
  ### There are 4 types of metadata supported in metadata generation
  TimesSDIsOutlier = 4,
  ### The number of deviations until a value is truncated as an outlier
  LoopingControlIncrement = 1000
  ### Max number of loops in looping functions
)

# float constants
c_d <- list(
  RunifMin = .01,
  RunifMax = .99,
  ### Continuous metadata are draw from uniform distributions, these are the bounds
  MinBinary = .3,
  MaxBinary = .7,
  ### Min and max probabilities for binary metadata levels
  ALittleMoreThanZero = 0.00001,
  SDBeta = 0.1251,
  SDIntercept = 1.121,
  ### The estimate for the relationship between SD and exp
  BetaZero = -0.6197338,
  #Beta2Zero = -0.0111924
  InterceptZero = 2.3536094,
  ### The estimate for the relationship between exp and zero percent
  ### modified by bor to logistic regression results on the IBD study
  BetaGrandSD = 0.04982219
  ### The estimate for the relationship between the mu of mus (of feature distributions) and the SD of mus (of feature distributions)
)

# Labeling Constants
c_str <- list(
  Level = "Level",
  Feature = "Feature",
  Metadata = "Metadata",
  Random = "Lognormal",
  Outlier = "Outlier",
  Spike = "spike",
  Null = "null",
  BugBugAssociations = "BugToBugAssociations",
  CorrDomainBugs = "Number of bugs each correlated bug is correlated with:",
  CorrDomainBugsIdx = "Indices of the bugs each correlated bug is correlated with:",
  CorrRangeBugsIdx = "Indices of bugs correlated with others:",
  DistributionParameters = "DistributionParameters",
  ExpVector = "Expected value vector (of lognormal):",
  MaxCorrDomainBugs = "Maximum number of bugs with which one bug is correlated:",
  MinimumSamples = "Minimum Spiked-in Samples:",
  MuVector = "Mu vector (of normal):",
  NoiseScaling  = "Scaling parameter for variance of noise:",
  NumberDatasets = "Number of datasets generated:",
  NumberOfAssociations = "Number of associations (bugs correlated with others):",
  LogCorrValues = "Specified correlation values of the log-counts:",
  DirOfAssociations    = "Direction of associations:",
  PercentZeroVector = "Percent zeros vector:",
  SDVector = "SD vector (of normal):",
  SyntheticMicrobiome = "SyntheticMicrobiome",
  SyntheticMicrobiomeBasis = "SyntheticMicrobiomeBasis",
  NumberOfFeatures = 'Number of features:',
  NumberOfSamples = 'Number of samples:',
  PercentSpikes = 'Percent spikes:',
  Multiplier = 'Multiplier:',
  MultiplierParameter = 'Multivariate Parameter:',
  TotalSampleBugOccurrence = "Total Reads per Sample:",
  NumberCounts = "Minimum Number of Counts:",
  NumberSamples = "in Minimum Number of Samples:",
  PercentOutliers = "Max Percent Outliers in a Sample:",
  PercentSampleOutliers = "Percent Samples with Outliers:",
  OutlierParameter = "Outlier Swap:",
  SampleParameter = "Sample:",
  Continuous = "Continuous",
  Factor = "Factor Levels",
  MetadataDetails = "Metadata: Details"
)

# Boolean constants
c_f <- list(
  IgnoreZerosInOutliers = TRUE,
  ### Will not allow zeros to be used as the min value in swapping unless they are needed to fulfill the number of
  ### swaps (if there are a whole bunch of zeros, some zeros may be needed or no swapping can be performed).
  FreezeSDFeatures = FALSE,
  FreezeSDGrandMu = FALSE,
  PrintLognormalMatrix = FALSE
)
