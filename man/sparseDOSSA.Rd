\name{sparseDOSSA}
\alias{sparseDOSSA}
\title{
Sparse Data Observations for Simulating Synthetic Abundance}
\description{
Sparse Data Observations for Simulating Synthetic Abundance}
\usage{
sparseDOSSA( strNormalizedFileName = "SyntheticMicrobiome.pcl",
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
             runBugBug =  FALSE )
}
\arguments{
  \item{runBugBug}{
  If given, bug-bug interaction files are generated in addition to any metadata files. This is a flag, it is either included or not included in the command line, no value needed.}
}

\value{
A list contains the names of the output files.
}

\author{
Boyu Ren<bor158@mail.harvard.edu>, Emma Schwager<eschwager@hsph.harvard.edu>, Timothy Tickle<ttickle@hsph.harvard.edu>, Curtis Huttenhower <chuttenh@hsph.harvard.edu>
}
\examples{
sparseDOSSA(strNormalizedFileName = "SyntheticMicrobiome.pcl",
	strCountFileName = "SyntheticMicrobiome-Counts.pcl",
	parameter_filename = "SyntheticMicrobiomeParameterFile.txt",
	bugs_to_spike = 0,
	calibrate = NA,
	datasetCount = 1,
	read_depth = 8030,
	number_features = 300,
	spikeCount = "1",
	lefse_file = NA,
	percent_spiked = 0.03,
	minLevelPercent =  0.1,
	number_samples = 50, 
	max_percent_outliers = 0.05,
	number_metadata = 5,
	spikeStrength =  "1.0",
	seed =  1,
	percent_outlier_spikins = 0.05,
	minOccurence =  0,
	verbose =  TRUE,
	minSample =  0,
	association_type =  "linear",
	noZeroInflate =  FALSE,
	noRunMetadata = FALSE,
	runBugBug =  FALSE)
}
 
