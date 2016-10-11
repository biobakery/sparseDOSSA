library(testthat)
library(sparseDOSSA)



setwd("..") 
 

expected_sparsedossa_results <- read.csv("tests/expected_SyntheticMicrobiome.pcl",  row.names=1,header=TRUE,sep="\t")

sparseDOSSA (
	bugs_to_spike = 0,
	calibrate = NA,
	datasetCount = 1,
	read_depth = 8030,
	number_features = 300,
	bugBugCorr =  "0.5",
	spikeCount = "1",
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

sparsedossa_results <- read.csv("SyntheticMicrobiome.pcl", row.names=1, header=TRUE,sep="\t")

expect_that(colnames(expected_sparsedossa_results),equals(colnames(sparsedossa_results)))
expect_that(rownames(expected_sparsedossa_results),equals(rownames(sparsedossa_results)))
