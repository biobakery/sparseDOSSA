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

# Plot constants
c_strDefaultMarkerColor = "black"
c_strOutlierColor = "red"
c_strFeatureOutlier = "cyan"

# Control flow constants
iLoopingControlIncrement = 1000
### Max number of looping in looping functions
c_iCountTypesOfMetadata = 4
### There are 4 types of metadata supported in metadata generation
c_fIgnoreZerosInOutliers = TRUE
### Will not allow zeros to be used as the min value in swapping unless they are needed to fulfill the number of
### swaps (if there are a whole bunch of zeros, some zeros may be needed or no swapping can be performed).
c_dRunifMin = .01
c_dRunifMax = .99
### Continuous metadata are draw from uniform distributions, these are the bounds
c_dMinBinary = .3
c_dMaxBinary = .7
### Min and max probabilities for binary metadata levels
c_iTimesSDIsOutlier = 4
### The number of deviations until a value is truncated as an outlier
c_dALittleMoreThanZero = 0.00001

# Labeling Constants
c_strLevel = "Level"
c_strFeature = "Feature"
c_strMetadata = "Metadata"
c_strRandom = "Lognormal"
c_strOutlier = "Outlier"
c_strSpike = "spike"
c_strNull = "null"
c_strBugBugAssociations = "BugToBugAssociations"
c_strMinimumSamples = "Minimum Spiked-in Samples:"
c_strSyntheticMicrobiome = "SyntheticMicrobiome"
c_strSyntheticMicrobiomeBasis = "SyntheticMicrobiomeBasis"
c_strNumberOfFeatures = 'Number of features:'
c_strNumberOfSamples = 'Number of samples:'
c_strPercentSpikes = 'Percent spikes:'
c_strMultiplier = 'Multiplier:'
c_strMultiplierParameter = 'Multivariate Parameter:'
c_strTotalSampleBugOccurrence = "Total Reads per Sample:"
c_strNumberCounts = "Minimum Number of Counts:"
c_strNumberSamples = "in Minimum Number of Samples:"
c_strPercentOutliers = "Max Percent Outliers in a Sample:"
c_strPercentSampleOutliers = "Percent Samples with Outliers:"
c_strOutlierParameter = "Outlier Swap:"
c_strSampleParameter = "Sample:"
c_strContinuous = "Continuous"
c_strFactor = "Factor Levels"
c_strMetadataDetails = "Metadata: Details"

### For Bug-bug spikin matrix (in alphabetical order)
c_strCorrDomainBugs = "Number of bugs each correlated bug is correlated with:"
c_strCorrDomainBugsIdx ="Indices of the bugs each correlated bug is correlated with:"
c_strCorrRangeBugsIdx = "Indices of bugs correlated with others:"
c_strMaxCorrDomainBugs = "Maximum number of bugs with which one bug is correlated:"
c_strNumberOfAssociations = "Number of associations (bugs correlated with others):"
c_strNoiseScaling  = "Scaling parameter for variance of noise:"

# Temporary control flags
c_dfFreezeSDFeatures = FALSE
c_dfFreezeSDGrandMu = FALSE
c_fPrintLognormalMatrix = FALSE

# Define the relationship between different feature properties
# These variables are associated with settings for
# Calculating the SD and percent Zero based on the mu or expectation.
# The constants for these variables have been estimated from
# a real (IBD) data set. These constants are used unless a
# calibration file is given which estimates the values in the
# same manner as the constants were estimated.
c_dSDBeta = 0.1251
c_dSDIntercept = 1.121
### The estimate for the relationship between SD and exp
c_dBetaZero = -0.6197338
#c_dBeta2Zero = -0.0111924
c_dInterceptZero = 2.3536094
### The estimate for the relationship between exp and zero percent
### modified by bor to logistic regression results on the IBD study
c_dBetaGrandSD = 0.04982219
### The estimate for the relationship between the mu of mus (of feature distributions) and the SD of mus (of feature distributions)
