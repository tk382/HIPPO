#' Expected zero proportion under Poisson
#'
#' @param lambda numeric vector of means of Poisson
#' @return numeric vector of expected proportion of zeros for each lambda
#' @examples
#' pois_prob_zero(3)
#' @export
pois_prob_zero = function(lambda){
  exp(-lambda)
}

#' Expected zero proportion under Negative Binomial
#'
#' @param lambda numeric vector of means of negative binomial
#' @param theta numeric vector of the dispersion parameter for negative binomial, 0 if poisson
#' @return numeric vector of expected zero proportion under Negative Binomial
#' @examples
#' nb_prob_zero(3, 1.1)
#' @export
nb_prob_zero = function(lambda, theta){
  if(theta==0){
    return(exp(-lambda))
  }else{
    return((1/(lambda * theta + 1))^(1/theta))
  }
}
#' Expected zero proportion under Negative Binomial
#'
#' @param lambda gene mean
#' @param theta dispersion parameter, 0 if zero-inflated poisson
#' @param pi zero inflation, 0 if negative binomial
#' @return Expected zero proportion under Zero-Inflated Negative Binomial
#' @examples
#' zinb_prob_zero(3, 1.1, 0.1)
#' @export
zinb_prob_zero = function(lambda, theta, pi){
  if(theta==0){
    return((1-pi) * exp(-lambda) + pi)
  }else{
    return(pi + (1-pi) * (1/(lambda * theta + 1))^(1/theta))
  }
}

#' Preprocess UMI data without cell label so that each row contains information about each gene
#'
#' @param sce SingleCellExperiment object with counts data
#' @return data frame with one row for each gene.
#' @examples
#' library(SingleCellExperiment)
#' library(HIPPO)
#' X = matrix(rpois(1000, 10), nrow = 100) # create random count matrix from poisson(10)
#' rownames(X) = paste0('gene',1:100)
#' colnames(X) = paste0('cell',1:10)
#' sce = SingleCellExperiment(assays = list(counts = X)) #create SingleCellExperiment object
#' df = preprocess_heterogeneous(sce) #get gene information
#' @export
preprocess_heterogeneous = function(sce){
  if(class(sce)=="SingleCellExperiment"){
    X = sce@assays@data$counts
  }else{
    stop("input must be a SingleCellExperiment object")
  }
  gene_mean = rowMeans(X)
  zero_proportion = rowMeans(X==0)
  where = which(gene_mean > 0)
  gene_var = NA
  gene_var[where] = matrixStats::rowVars(X[where,])

  df = data.frame(gene = rownames(X),
                  gene_mean = rowMeans(X),
                  zero_proportion = rowMeans(X==0),
                  gene_var = gene_var)
  df$samplesize = ncol(X)
  return(df)
}

#' Preprocess UMI data with inferred or known labels
#'
#' @param sce SingleCellExperiment object with counts data
#' @param label a numeric or character vector of inferred or known label
#' @param normalize boolean whether to normalize each cell to have the same sequencing depth. Default as FALSE
#' @return data frame with one row for each gene.
#' @examples
#' library(SingleCellExperiment)
#' library(HIPPO)
#' X = matrix(rpois(2000, 10), nrow = 100) # create random count matrix from poisson(10)
#' rownames(X) = paste0('gene',1:100)
#' colnames(X) = paste0('cell',1:20)
#' sce = SingleCellExperiment(assays = list(counts = X)) #create SingleCellExperiment object
#' label = sample(1:2, size=20, replace = TRUE) # create fake cell type label
#' label = as.factor(label)
#' df = preprocess_homogeneous(sce, label = label) #get gene information
#' @export
preprocess_homogeneous = function(sce, label, normalize = FALSE){
  if(class(sce)=="SingleCellExperiment"){
    X = sce@assays@data$counts
  }else{
    stop("input must be a SingleCellExperiment object")
  }
  sf = median(colSums(X))
  if(normalize){
    X = apply(X, 2, function(x) x/sum(x)*sf)
  }
  labelnames = as.character(unique(label))
  zero_proportion = matrix(NA, nrow(X), length(labelnames))
  gene_mean = matrix(NA, nrow(X), length(labelnames))
  positive_mean = matrix(NA, nrow(X), length(labelnames))
  gene_var = matrix(NA, nrow(X), length(labelnames))
  samplesize = table(label)
  for (i in 1:length(labelnames)){
    ind = which(label==labelnames[i])
    zero_proportion[,i] = rowMeans(X[,ind]==0)
    gene_mean[,i] = rowMeans(X[, ind])
    where = gene_mean[,i] != 0
    gene_var[where,i] = matrixStats::rowVars(X[where, ind])
  }
  colnames(zero_proportion) =
    colnames(gene_mean) =
    colnames(gene_var) = labelnames
  gene_mean = as.data.frame(gene_mean)
  zero_proportion = as.data.frame(zero_proportion)
  gene_var = as.data.frame(gene_var)
  gene_mean$id = zero_proportion$id  = gene_var$id = rownames(X)
  mgm = reshape2::melt(gene_mean, id = "id")
  mdor = reshape2::melt(zero_proportion, id = "id")
  mgv = reshape2::melt(gene_var, id = "id")
  df = data.frame(gene = rownames(X),
                  gene_mean = mgm$value,
                  gene_var = mgv$value,
                  zero_proportion = mdor$value,
                  celltype = mgm$variable)
  df$samplesize = NA
  for (i in names(samplesize)){
    df[df$celltype == i, "samplesize"] = samplesize[i]
  }
  rownames(df) = c()
  return(df)
}

#' Conduct feature selection by computing test statistics for each gene
#'
#' @param df pre-processed data frame
#' @return data frame with added columns with test results
#' @examples
#' library(SingleCellExperiment)
#' library(HIPPO)
#' X = matrix(rpois(1000, 10), nrow = 100) # create random count matrix from poisson(10)
#' rownames(X) = paste0('gene',1:100)
#' colnames(X) = paste0('cell',1:10)
#' sce = SingleCellExperiment(assays = list(counts = X)) #create SingleCellExperiment object
#' df = preprocess_heterogeneous(sce) #get gene information
#' df = compute_test_statistic(df)
#' @export
compute_test_statistic = function(df){
  ind = which(df$gene_mean==0)
  if(length(ind)>0){
    df = df[-ind,]
  }
  ind = grep("^MT-", df$gene)
  if(length(ind) > 0){
    df = df[-grep("^MT-", df$gene), ]
  }
  df = df %>% dplyr::mutate(expected_pi= pmin(2*.data$samplesize/(2*.data$samplesize-1) * exp(-.data$gene_mean), 1 - 1e-10)) %>%
    dplyr::mutate(se = sqrt(.data$expected_pi * (1-.data$expected_pi) / (.data$samplesize-1.25))) %>%
    dplyr::mutate(minus_logp = -pnorm(.data$zero_proportion, .data$expected_pi, .data$se, log.p = TRUE, lower.tail=FALSE)) %>%
    dplyr::mutate(minus_logp = pmin(.data$minus_logp, 500)) %>%
    dplyr::mutate(zvalue = -qnorm(exp(-.data$minus_logp)))
  df$gene = as.character(df$gene)
  return(df)
}

#' Clustering of one-step
#' @param subX a matrix that is gene by cell that needs clustering
#' @param z_threshold z-value threshold for feature selection
#' @return a list of clustering result with relevant features, PCA reslut, and kmeans result
one_level_clustering = function(subX, z_threshold){
  subdf = preprocess_heterogeneous(subX)
  subdf = compute_test_statistic(subdf)
  features = subdf$gene[subdf$zvalue>z_threshold]
  if(length(features)<10){
    return(list(features = NA, pcs = NA, km = NA))
  }
  pcs = irlba::irlba(log(subX[features, ]+1), 10)$v
  unscaledpc = irlba::prcomp_irlba(log(t(subX[features,])+1), n = 10, scale.=FALSE, center=FALSE)$x[,1:10]
  km = kmeans(pcs, 2, nstart = 10, iter.max = 50)
  return(list(features = features, pcs = pcs, km = km, unscaled_pcs = unscaledpc, subdf = subdf))
}

#' HIPPO's hierarchical clustering
#'
#' @param sce SingleCellExperiment object
#' @param K number of clusters to ultimately get
#' @param z_threshold numeric > 0 as a z-value threshold for selecting the features
#' @param outlier_proportion numeric between 0 and 1, a cut-off so that when the proportion of important features reach this number, the clustering terminates
#' @return a list of clustering result for each level of k=1, 2, ... K.
#' @examples
#' library(SingleCellExperiment)
#' library(HIPPO)
#' X = matrix(rpois(50000, 3), nrow = 1000) # create random count matrix from poisson(10)
#' X[X%in%c(1,2)] = 0
#' rownames(X) = paste0('gene',1:1000)
#' colnames(X) = paste0('cell',1:50)
#' sce = SingleCellExperiment(assays = list(counts = X)) #create SingleCellExperiment object
#' sce = hippo(sce, K = 3)
#' @export
hippo = function(sce, K=10, z_threshold = 3, outlier_proportion = 0.01){
  if(class(sce)=="SingleCellExperiment"){
    X = sce@assays@data$counts
    ref = SingleCellExperiment::rowData(sce)
  }else if (class(sce)=="matrix"){
    sce = SingleCellExperiment::SingleCellExperiment(assays = list(counts = sce))
    X = sce@assays@data$counts
  }else{
    stop("input must be either matrix or SingleCellExperiment object")
  }
  if(outlier_proportion >= 1 | outlier_proportion <= 0){
    stop("Outlier_proportion must be a number between 0 and 1. Default is 5%")
  }
  outlier_number = nrow(X) * outlier_proportion
  labelmatrix = matrix(NA, ncol(X), K); labelmatrix[,1] = 1
  eachlevel = list();   subX = X
  subXind = 1:ncol(X); withinss = rep(0, K)
  oldk = 1
  features = list();   featuredata = list()
  for (k in 2:K){
    print(paste0("K = ", k, ".."))
    thisk = one_level_clustering(subX, z_threshold)
    if(length(thisk$features) < outlier_number){
      print("not enough important features left; terminating the procedure")
      labelmatrix = labelmatrix[,1:(k-1)]
      break
    }
    labelmatrix[,k] = labelmatrix[,k-1]
    labelmatrix[subXind[thisk$km$cluster==2], k] = k
    withinss[oldk] = sum(apply(thisk$unscaled_pcs[thisk$km$cluster==1, ], 1, var)^2)
    withinss[k] = sum(apply(thisk$unscaled_pcs[thisk$km$cluster==2, ], 1, var)^2)
    oldk = which.max(withinss[1:k])
    subX = X[thisk$features, which(labelmatrix[,k]==oldk)]
    subXind = which(labelmatrix[,k]==oldk)
    features[[k-1]] = thisk$features
    featuredata[[k-1]] = thisk$subdf
  }
  sce@int_metadata$hippo = list(X = X, features = features, labelmatrix = labelmatrix)
  return(sce)
}

#' visualize each round of hippo through zero proportion plot
#' @param sce SingleCellExperiment object with hippo element in it
#' @return a ggplot object that shows the zero proportions for each round
#' @examples
#' library(SingleCellExperiment)
#' library(HIPPO)
#' X = matrix(rpois(50000, 3), nrow = 1000) # create random count matrix from poisson(10)
#' X[X%in%c(1,2)] = 0
#' rownames(X) = paste0('gene',1:1000)
#' colnames(X) = paste0('cell',1:50)
#' sce = SingleCellExperiment(assays = list(counts = X)) #create SingleCellExperiment object
#' sce = hippo(sce, K = 3)
#' zero_proportion_plot(sce)
#' @export
zero_proportion_plot = function(sce){
  hippo_object = sce@int_metadata$hippo
  plist = list()
  dflist = list()
  K = ncol(hippo_object$labelmatrix)
  df = preprocess_heterogeneous(hippo_object$X)
  df$celltype = "combined"
  df$selected_feature = FALSE
  df$K = 1
  dflist[[1]] = df
  for (i in 2:K){
    df = preprocess_homogeneous(sce, label = hippo_object$labelmatrix[,i])
    df$selected_feature = df$gene %in% hippo_object$features[[i-1]]
    df$K = i
    dflist[[i]] = df
  }
  df = do.call(rbind, dflist)
  df = df[sample(nrow(df)), ]
  zero_plot = ggplot2::ggplot(df, ggplot2::aes(x = .data$gene_mean, y = .data$zero_proportion, col = .data$celltype)) +
    ggplot2::geom_point(size = 0.4, alpha = 0.5) +
    ggplot2::facet_wrap(~.data$K) +
    ggplot2::geom_line(aes(x = .data$gene_mean, y = exp(-.data$gene_mean)), col = 'black') +
    ggplot2::xlim(c(0,10))+
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme_bw() +
    ggplot2::ylab("zero proportion") +
    ggplot2::xlab("gene mean") +
    ggplot2::theme(legend.title = ggplot2::element_blank())
}

#' compute t-SNE or umap of each round of HIPPO
#' @param sce SingleCellExperiment object with hippo object in it.
#' @param method a string that determines the method for dimension reduction: either "umap" or 'tsne
#' @param perplexity numeric perplexity parameter for Rtsne function
#' @return a data frame of dimension reduction result for each k in 1, ..., K
#' @examples
#' library(SingleCellExperiment)
#' library(HIPPO)
#' X = matrix(rpois(50000, 3), nrow = 1000) # create random count matrix from poisson(10)
#' X[X%in%c(1,2)] = 0
#' rownames(X) = paste0('gene',1:1000)
#' colnames(X) = paste0('cell',1:50)
#' sce = SingleCellExperiment(assays = list(counts = X)) #create SingleCellExperiment object
#' sce = hippo(sce, K = 3)
#' sce = dimension_reduction(sce, method = "tsne", perplexity = 2)
#' @export
dimension_reduction = function(sce, method = c("umap", "tsne"), perplexity = 30){
  hippo_object = sce@int_metadata$hippo
  dflist = list()
  K = ncol(hippo_object$labelmatrix)
  sce@int_metadata$hippo_object$umap = NA
  sce@int_metadata$hippo_object$tsne = NA
  if (method=="umap"){
    um = umap::umap(log(t(hippo_object$X[hippo_object$features[[1]], ])+1))
    um = as.data.frame(um$layout)
    umdf = data.frame()
    for (i in 2:K){
      df = preprocess_homogeneous(sce, label = hippo_object$labelmatrix[,i])
      df$selected_feature = df$gene %in% hippo_object$features[[i-1]]
      df$K = i
      dflist[[i]] = df
      umdf = rbind(umdf, data.frame(umap1 = um$V1,
                                    umap2 = um$V2,
                                    K = i,
                                    label = hippo_object$labelmatrix[,i]))
    }
    umdf$label = as.factor(umdf$label)
    sce@int_metadata$hippo$umap = umdf
    return(sce)
  }else if(method == "tsne"){
    tsne = Rtsne::Rtsne(log(t(hippo_object$X[hippo_object$features[[1]], ])+1), perplexity = perplexity)
    tsnedf = data.frame()
    for (i in 2:K){
      df = preprocess_homogeneous(sce, label = hippo_object$labelmatrix[,i])
      df$selected_feature = df$gene %in% hippo_object$features[[i-1]]
      df$K = i
      dflist[[i]] = df
      tsnedf = rbind(tsnedf, data.frame(tsne1 = tsne$Y[,1],
                                        tsne2 = tsne$Y[,2],
                                        K=i,
                                        label = hippo_object$labelmatrix[,i]))
    }
    tsnedf$label = as.factor(tsnedf$label)
    sce@int_metadata$hippo$tsne = tsnedf
    return(sce)
  }
}


#' visualize each round of hippo through UMAP
#'
#' @param sce SingleCellExperiment object with hippo and UMAP result in it
#' @return ggplot object for umap in each round
#' @examples
#' library(SingleCellExperiment)
#' library(HIPPO)
#' X = matrix(rpois(50000, 3), nrow = 1000) # create random count matrix from poisson(10)
#' X[X%in%c(1,2)] = 0
#' rownames(X) = paste0('gene',1:1000)
#' colnames(X) = paste0('cell',1:50)
#' sce = SingleCellExperiment(assays = list(counts = X)) #create SingleCellExperiment object
#' sce = hippo(sce, K = 3)
#' sce = dimension_reduction(sce, method="umap")
#' hippo_umap_plot(sce)
#' @export
hippo_umap_plot = function(sce){
  umdf = sce@int_metadata$hippo$umap
  if(!is.na(umdf[1])){
    umap_plot = ggplot2::ggplot(umdf, ggplot2::aes(x = .data$umap1, y= .data$umap2, col = .data$label)) +
      ggplot2::facet_wrap(~.data$K) +
      ggplot2::geom_point(size = 0.4, alpha = 0.5) +
      ggplot2::theme_bw() +
      ggplot2::ylab("umap2") +
      ggplot2::xlab("umap1")+
      ggplot2::theme(legend.title = ggplot2::element_blank())
  }else{
    stop("use dimension_reduction to compute umap first")
  }
}

#' visualize each round of hippo through t-SNE
#' @param sce SincleCellExperiment object with hippo and t-SNE result in it
#' @return ggplot object for t-SNE in each round
#' @examples
#' library(SingleCellExperiment)
#' library(HIPPO)
#' X = matrix(rpois(50000, 3), nrow = 1000) # create random count matrix from poisson(10)
#' X[X%in%c(1,2)] = 0
#' rownames(X) = paste0('gene',1:1000)
#' colnames(X) = paste0('cell',1:50)
#' sce = SingleCellExperiment(assays = list(counts = X)) #create SingleCellExperiment object
#' sce = hippo(sce, K = 3)
#' sce = dimension_reduction(sce, method = "tsne", perplexity = 2)
#' hippo_tsne_plot(sce)
#' @export
hippo_tsne_plot = function(sce){
  tsnedf = sce@int_metadata$hippo$tsne
  if(!is.na(tsnedf[1])){
    tsne_plot = ggplot2::ggplot(tsnedf, ggplot2::aes(x = .data$tsne1, y=.data$tsne2, col=.data$label)) +
      ggplot2::facet_wrap(~.data$K) +
      ggplot2::geom_point(size=0.4, alpha = 0.5) +
      ggplot2::theme_bw() +
      ggplot2::ylab("tsne2") +
      ggplot2::xlab("tsne1")+
      ggplot2::theme(legend.title = ggplot2::element_blank())
  }else{
    stop("use dimension_reduction to compute tsne first")
  }
  return(tsne_plot)
}

#' HIPPO's differential expression
#'
#' @param sce SingleCellExperiment object with hippo
#' @param top.n number of markers to return
#' @param switch_to_hgnc if the current gene names are ensemble ids, and would like to switch to hgnc
#' @param ref a data frame with columns "hgnc" and "ensg" to match each other, only required when switch_to_hgnc is set to TRUE
#' @return list of differential expression result
#' @examples
#' library(SingleCellExperiment)
#' library(HIPPO)
#' X = matrix(rpois(50000, 3), nrow = 1000) # create random count matrix from poisson(10)
#' X[X%in%c(1,2)] = 0
#' rownames(X) = paste0('gene',1:1000)
#' colnames(X) = paste0('cell',1:50)
#' sce = SingleCellExperiment(assays = list(counts = X)) #create SingleCellExperiment object
#' sce = hippo(sce, K = 3)
#' sce = diffexp(sce)
#' @export
diffexp = function(sce, top.n = 5, switch_to_hgnc=FALSE, ref = NA){
  if(switch_to_hgnc & is.na(ref)){
    stop("A reference must be provided in order to match ENSG ids to HGNC symbols")
  }
  hippo_object = sce@int_metadata$hippo
  K = ncol(hippo_object$labelmatrix)
  featureind = list()
  cellind = list()
  plist = list()
  featureind[[1]] = 1:nrow(hippo_object$X)
  cellind[[1]] = 1:ncol(hippo_object$X)
  labelmatrix = hippo_object$labelmatrix
  result = list()
  count = hippo_object$X
  for (k in 2:K){
    features = hippo_object$features[[k-1]]

    cellind = which(labelmatrix[,k-1] == labelmatrix[which(labelmatrix[,k-1] != labelmatrix[,k])[1], k-1])
    types = unique(hippo_object$labelmatrix[cellind, k])
    cellgroup1 = which(hippo_object$labelmatrix[,k] == types[1])
    cellgroup2 = which(hippo_object$labelmatrix[,k] == types[2])
    rowdata = data.frame(genes = features)
    rowdata$meandiff = rowMeans(count[features,cellgroup1]) - rowMeans(count[features,cellgroup2])
    rowdata$sd = sqrt(rowMeans(count[features,cellgroup1])/length(cellgroup1) +
                        rowMeans(count[features,cellgroup2])/length(cellgroup2))
    rowdata$z = rowdata$meandiff/rowdata$sd
    rowdata = rowdata[order(rowdata$z, decreasing=TRUE), ]
    rowdata$genes = as.character(rowdata$genes)
    newcount = t(log(cbind(count[rowdata$genes[1:top.n], cellgroup1],
                           count[rowdata$genes[1:top.n], cellgroup2])+1))
    topgenes = rowdata$genes[1:top.n]
    if(switch_to_hgnc){
      features_hgnc = ref$hgnc[match(topgenes, ref$ensg)]
    }
    newcount = as.data.frame(newcount)
    colnames(newcount) = rowdata$genes[1:top.n]
    if(switch_to_hgnc){
      colnames(newcount) = features_hgnc
      rowdata$hgnc = ref$hgnc[match(rowdata$genes, ref$ensg)]
    }
    newcount$celltype = c(rep(types[1], length(cellgroup1)), rep(types[2], length(cellgroup2)))
    newcount = reshape2::melt(newcount, id="celltype")
    newcount$celltype = as.factor(newcount$celltype)
    colnames(newcount) = c("celltype", "gene", "logcount")
    g = ggplot(newcount, ggplot2::aes(x = .data$gene, y = .data$logcount, col = .data$celltype)) +
      geom_boxplot(outlier.size = 0.2) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      ggtitle(paste("Round", k-1)) +
      xlab("") +
      theme(legend.position='none') +
      ylab("log count")

    plist[[k-1]] = g
    result[[k-1]] = rowdata
  }
  bpl = gridExtra::grid.arrange(grobs = plist, ncol = 2)
  sce@int_metadata$hippo$diffexp$result_table = result
  sce@int_metadata$hippo$diffexp$plot = bpl
  return(sce)
}
