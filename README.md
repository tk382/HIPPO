[![Travis-CI Build Status](https://travis-ci.com/tk382/HIPPO.svg?branch=master)](https://travis-ci.org/tk382/HIPPO)
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/DynamicCorrelation)](https://cran.r-project.org/package=DynamicCorrelation)-->

# HIPPO <img src="https://github.com/tk382/HIPPO/blob/master/readme/hippo_image.png" width="40">


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

### Plot diagnostics

This plot shows the zero inflation compared to the expected Poisson line. If most genes don't align with the black line, it shows that there is cell heterogeneity driving the zero inflation. 

```
hippo_diagnostic_plot(sce)
```
<img src="https://github.com/tk382/HIPPO/blob/master/readme/diagnostic_plot.png" width="350">

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

We offer two dimension reduction methods: umap and tsne.

```
sce = dimension_reduction(sce, method = "umap")
sce = dimension_reduction(sce, method = "tsne")
```

And we offer two separate visualization functions.

```
hippo_umap_plot(sce)
```
<img src="https://github.com/tk382/HIPPO/blob/master/readme/umap.png" width="800">
 
```
hippo_tsne_plot(sce)
```
<img src="https://github.com/tk382/HIPPO/blob/master/readme/tsne.png" width="800">

### Diagnostic plots after HIPPO

We also recommend users to check the diagnostic plots of zero inflation for each round of hippo as below. The zero proportion of each gene in each clustered group aligns better and better after each round of hierarchical clustering. 

```
zero_inflation_plot(sce)
```
<img src="https://github.com/tk382/HIPPO/blob/master/readme/zero_inflation.png" width="800">

### Differential Expression with HIPPO

We also offer a differential expression analysis tool. 

*diffexp* has arguments called *switch_to_hgnc* and *ref*. These aim to provide the users an option to change the gene names from ENSG IDs to HGNC symbols for ease of understanding. Many SingleCellExperiment objects have such data embedded in *rowData(sce)*. Users can create a data frame with *ensg* and *hgnc* columns for the genes, and HIPPO will automatically switch the row names of the count matrix from ENSG IDs to HGNC symbols. The default is set to FALSE, assuming that the row names are already HGNC symbols.

*top.n* argument lets users choose how many top genes to show in the box plot. The default is 5.

```
ref = data.frame(hgnc = rowData(sce)$symbol, ensg = rowData(sce)$id)
sce = diffexp(sce, top.n = 5, switch_to_hgnc = TRUE, ref = ref)
```
<img src="https://github.com/tk382/HIPPO/blob/master/readme/diffexp.png" width="600">

The labels of boxplots are not quite straightfoward, as we look at different cell groups at each round of HIPPO. This must be interepreted simultaneously with the hierarchical clustering plot such as t-SNE plot above. 

First, in the first round, K moves from 1 to 2, and the red group is separated. This group is Monocytes in this particular data set. The box plot always shows the "separated" group in the green box, and hence the group 1. The red boxes represent the remaining cells, so groups 2, 3, and 4 combined.

In the second round, K moves from 2 to 3, and as shown in the t-SNE plot, groups 3 and 4 are separated from the group 2 because they're assigned a new color. Group 2 in this data set is B cells. (Note that this is different from group 2 is separated. Group 2 remains green.) Therefore, in the second round of box plots, the green boxes represent groups 3 and 4, and the red box represents group 2. In this round, the first group of cells have been removed from the samples. 

In the last round, K moves from 3 to 4, and group 4 has been assigned a new color of violet. Hence, the green boxes represent group 4, which is Regulatory T cells, while the red boxes represent the remaining cells of group 3: Naive T cells.



Each round of differential expression test results are also saved in the list of data frames.

```
sce@int_metadata$hippo$diffexp$result_table[[1]]
```

Above code will show the list of genes in the order of significance that differentiates the first group from the rest. The second element of the list will show the list of genes in the order of significance that differentiates the third and fourth group from the second group.

## Authors

* [Tae Kim](https://github.com/tk382)

## Acknowledgments

* [Mengjie Chen](http://www.mengjiechen.com) provided guidance in methodology development.
* [Yong Peng](https://github.com/bigdataage) contributed in packaging the code to meet the Bioconductor requirements.
* The hippo icon is from [here](https://www.needpix.com/photo/178308/hippo-head-cartoon-cute-grey-zoo-wildlife)
