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

## 



## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* [Tae Kim](https://github.com/tk382)

## Acknowledgments

* [Mengjie Chen](http://www.mengjiechen.com) provided guidance in methodology development.
* [Yong Peng](https://github.com/bigdataage) contributed in packaging the code to meet the Bioconductor requirements.
* The hippo icon comes from [here](https://www.needpix.com/photo/178308/hippo-head-cartoon-cute-grey-zoo-wildlife)
