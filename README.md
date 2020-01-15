[![Travis-CI Build Status](https://travis-ci.com/tk382/HIPPO.svg?branch=master)](https://travis-ci.org/tk382/HIPPO)
[![CRAN status](https://www.r-pkg.org/badges/version/DynamicCorrelation)](https://cran.r-project.org/package=DynamicCorrelation)

# HIPPO <img src="https://github.com/tk382/HIPPO/blob/master/hippo_image.png" width="60">


Single cell UMI analysis tool that focuses on zero-inflation to detect biological heterogeneity

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

HIPPO works on the SingleCellExperiment object. You can download the library like the following.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
```

### Installing

HIPPO is under review for Bioconductor and CRAN submission. You can download the developer version as below. Please allow up to 5 minutes to completely compile the vignette.

```
devtools::install_github("tk382/HIPPO", build_vignettes = TRUE)
```

## Example Analysis

You can see a full analysis example through the vignette. 

```
browseVignettes("HIPPO")
```

Here is a brief tutorial to get you started.

### Prepare the data

Many UMI data sets are already SingleCellExperiment objects. We have an example data set uploaded on github "sce_Zhengmix4eq.rds" around 6Mb. 

```
sce = readRDS("sce_Zhengmix4eq.rds")
```

We also have an example data set that is only count matrix: "zhengmix4eq_counts.rds". In that case, you can use the count matrix to create a SingleCellExperiment object.

```
X = readRDS("zhengmix4eq_counts.rds")
sce = SingleCellExperiment(assays=list(counts = X))
```

### Run hippo

HIPPO assumes that the count matrix is placed in sce@assays@data$counts. Some objects that we found online have the count matrix in sce@assays$data$counts. In this case, HIPPO will throw an error because it cannot found a count matrix. In this case, you have to create another SingleCellExperiment object to assign the count matrix in the correct slot.

Next, you can run hippo function to do the pre-processing that simutlaneously conducts feature selection and hierarchcial clustering. There are three arguments that help you decide the stopping criterion of clustering procedure.

1. K is the maximum number of clusters that you want. HIPPO will return the clustering results for all k = 2, 3, ..., K, so you can overestimate the number of potential clusters. The default is 10, but users are highly recommended to adjust this.

2. z_threshold is the feature selection criterion. For each round of hierarchical clustering, hippo will find outlier genes where the z-value of significance is greater than the threshold. For example, if you would like to select genes with p-values less than 0.05, z_threshold would be 1.96. The default threshold is 2, but users can use their discretion to change this value.

3. outlier_proportion is the number of outlier genes to allow. The default is 0.01 (1\%) which means the clustering procedure will automatically stop if there are less than 1\% of genes remain as important features. With the example data set, the default choice has empirically worked well.

```
sce = hippo(sce, K=4, z_threshold = 2, outlier_proportion = 0.01)
```

### Dimension Reduction and Visualization

## Authors

* [Tae Kim](https://github.com/tk382)

## Acknowledgments

* [Mengjie Chen](http://www.mengjiechen.com) provided guidance in methodology development.
* [Yong Peng](https://github.com/bigdataage) contributed in packaging the code to meet the Bioconductor requirements.
* The hippo icon comes from [here](https://www.needpix.com/photo/178308/hippo-head-cartoon-cute-grey-zoo-wildlife)
