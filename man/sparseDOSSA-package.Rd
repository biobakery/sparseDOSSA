\name{sparseDOSSA-package}
\alias{sparseDOSSA-package}
\docType{package}
\title{
Sparse Data Observations for Simulating Synthetic Abundance
}
\description{
Sparse Data Observations for Simulating Synthetic Abundance
}
\details{
\tabular{ll}{
Package: \tab sparseDOSSA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-27\cr
License: \tab MIT\cr
}
Sparse Data Observations for Simulating Synthetic Abundance can be used by invoking the \cr
command ./synthetic_datasets_script.R  and overriding the input parms
}
\author{
Timothy Tickle<ttickle@hsph.harvard.edu>,\cr
Curtis Huttenhower <chuttenh@hsph.harvard.edu>\cr
Maintainer:  Boyu Ren<bor158@mail.harvard.edu>, Emma Schwager <emma.schwager@gmail.com>, George Weingart<george.weingart@gmail.com>  
}
\examples{
sparseDOSSA (
strNormalizedFileName = "SyntheticMicrobiome.pcl",
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
runBugBug =  FALSE
)
}
