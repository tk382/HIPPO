RowVar <- function(x) {
  Matrix::rowSums((x - Matrix::rowMeans(x))^2)/(ncol(x) - 1)
}

compute_test_statistic = function(df) {
  ind = which(df$gene_mean == 0)
  if (length(ind)) {
    df = df[-ind, ]
  }
  ind = grep("^MT-", df$gene)
  if (length(ind)) {
    df = df[-grep("^MT-", df$gene), ]
  }
  df = df %>%
    dplyr::mutate(const1 = 2*.data$samplesize/(2*.data$samplesize-1)) %>%
    dplyr::mutate(expected_pi1 = .data$const1*exp(-.data$gene_mean)) %>%
    dplyr::mutate(expected_pi = pmin(.data$expected_pi1, 1 - 1e-10)) %>%
    dplyr::mutate(const2 = (1-.data$expected_pi)/(.data$samplesize-1.25)) %>%
    dplyr::mutate(se = sqrt(.data$expected_pi*.data$const2)) %>%
    dplyr::mutate(propdiff = .data$zero_proportion - .data$expected_pi) %>%
    dplyr::mutate(zvalue = .data$propdiff/.data$se) %>%
    dplyr::mutate(zvalue = pmin(.data$zvalue, 500)) %>%
    dplyr::mutate(minus_logp = -pnorm(.data$zero_proportion,
                                      .data$expected_pi,
                                      .data$se,
                                      log.p = TRUE,
                                      lower.tail = FALSE)) %>%
    dplyr::mutate(minus_logp = pmin(.data$minus_logp, 500))
  df$gene = as.character(df$gene)
  return(df)
}


one_level_clustering = function(subX, z_threshold) {
  subdf = preprocess_heterogeneous(subX)
  subdf = compute_test_statistic(subdf)
  features = subdf[subdf$zvalue > z_threshold, ]
  nullfeatures = data.frame(matrix(ncol = 11, nrow = 0))
  colnames(nullfeatures) = c("gene", "gene_mean", "zero_proportion",
                             "gene_var", "samplesize", "expected_pi", "se",
                             "minus_logp","zvalue", "subsetK", "K")
  if (nrow(features) < 10) {
    return(list(features = nullfeatures, pcs = NA, km = NA))
  }
  if (nrow(features) < 10) {
    return(list(features = nullfeatures, pcs = NA, km = NA,
                unscaled_pcs = NA,subdf = NA))
  }
  pcs = tryCatch(expr = {
    irlba::irlba(log(subX[features$gene, ] + 1), min(9, nrow(features) -
                                                       1, ncol(subX) - 1))$v
  }, error = function(e) {
    NA
  }, warning = function(w) {
    NA
  })
  if (is.na(pcs[1])) {
    return(list(features = nullfeatures, pcs = NA, km = NA,
                unscaled_pcs = NA,
                subdf = NA))
  } else {
    unscaledpc = irlba::prcomp_irlba(log(Matrix::t((subX[features$gene,])) + 1),
                                     n = min(9, nrow(features) - 1,
                                             ncol(subX) - 1),
                                     scale. = FALSE, center = FALSE)$x
    km = kmeans(pcs, 2, nstart = 10, iter.max = 50)
  }
  return(list(features = features,
              pcs = pcs,
              km = km,
              unscaled_pcs = unscaledpc,
              subdf = subdf))
}


#' Expected zero proportion under Poisson
#'
#' @param lambda numeric vector of means of Poisson
#' @return numeric vector of expected proportion of zeros for each lambda
#' @examples
#' pois_prob_zero(3)
#' @export
pois_prob_zero = function(lambda) {
  exp(-lambda)
}

#' Expected zero proportion under Negative Binomial
#'
#' @param lambda numeric vector of means of negative binomial
#' @param theta numeric vector of the dispersion parameter
#' for negative binomial, 0 if poisson
#' @return numeric vector of expected zero proportion under
#' Negative Binomial
#' @examples
#' nb_prob_zero(3, 1.1)
#' @export
nb_prob_zero = function(lambda, theta) {
  if (theta == 0) {
    return(exp(-lambda))
  } else {
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
zinb_prob_zero = function(lambda, theta, pi) {
  return(pi + (1-pi) * nb_prob_zero(lambda, theta))
}

#' Preprocess UMI data without cell label so that each row contains
#' information about each gene
#'
#' @param X a matrix object with counts data
#' @return data frame with one row for each gene.
#' @examples
#' data(toydata)
#' df = preprocess_heterogeneous(get_data_from_sce(toydata))
#' @export
preprocess_heterogeneous = function(X) {
  gene_mean = Matrix::rowMeans(X)
  zero_proportion = Matrix::rowMeans(X == 0)
  where = which(gene_mean > 0)
  gene_var = rep(NA, nrow(X))
  gene_var[where] = RowVar(X[where, ])
  df = data.frame(gene = rownames(X), gene_mean = Matrix::rowMeans(X),
                  zero_proportion = Matrix::rowMeans(X == 0),
                  gene_var = gene_var)
  df$samplesize = ncol(X)
  df = compute_test_statistic(df)
  return(df)
}

#' Preprocess UMI data with inferred or known labels
#'
#' @param sce SingleCellExperiment object with counts data
#' @param label a numeric or character vector of inferred or known label
#' @return data frame with one row for each gene.
#' @examples
#' data(toydata)
#' labels = SingleCellExperiment::colData(toydata)$phenoid
#' df = preprocess_homogeneous(toydata, label = labels)
#' @export
preprocess_homogeneous = function(sce, label) {
  if (is(sce, "SingleCellExperiment")) {
    X = sce@assays@data$counts
  } else if (is(sce, "matrix")){
    sce = SingleCellExperiment::SingleCellExperiment(assays=list(counts = sce))
    X = sce@assays@data$counts
  }else{
    stop("input must be a SingleCellExperiment object")
  }
  labelnames = as.character(unique(label))
  zero_proportion = matrix(NA, nrow(X), length(labelnames))
  gene_mean = matrix(NA, nrow(X), length(labelnames))
  positive_mean = matrix(NA, nrow(X), length(labelnames))
  gene_var = matrix(NA, nrow(X), length(labelnames))
  samplesize = table(label)
  for (i in seq(length(labelnames))) {
    ind = which(label == labelnames[i])
    if(length(ind) >= 2){
      zero_proportion[, i] = Matrix::rowMeans(X[, ind] == 0)
      gene_mean[, i] = Matrix::rowMeans(X[, ind])
      where = gene_mean[, i] != 0
      gene_var[where, i] = RowVar(X[where, ind])
    }else{
      zero_proportion[, i] = mean(X[, ind] == 0)
      gene_mean[, i] = mean(X[, ind])
      where = gene_mean[, i] != 0
      gene_var[where, i] = var(X[where, ind])
    }

  }
  colnames(zero_proportion) = colnames(gene_mean) =
    colnames(gene_var) = labelnames
  gene_mean = as.data.frame(gene_mean)
  zero_proportion = as.data.frame(zero_proportion)
  gene_var = as.data.frame(gene_var)
  gene_mean$id = zero_proportion$id = gene_var$id = rownames(X)
  mgm = reshape2::melt(gene_mean, id = "id")
  mdor = reshape2::melt(zero_proportion, id = "id")
  mgv = reshape2::melt(gene_var, id = "id")
  df = data.frame(gene = rownames(X), gene_mean = mgm$value,
                  gene_var = mgv$value, zero_proportion = mdor$value,
                  celltype = mgm$variable)
  df$samplesize = NA
  for (i in names(samplesize)) {
    df[df$celltype == i, "samplesize"] = samplesize[i]
  }
  rownames(df) = NULL
  return(df)
}

#' Conduct feature selection by computing test statistics for each gene
#'
#' @param sce SingleCellExperiment object with count matrix
#' @param show_outliers boolean to indicate whether to circle the outliers
#' with given zvalue_thresh
#' @param zvalue_thresh a numeric v for defining outliers
#' @examples
#' data(toydata)
#' hippo_diagnostic_plot(toydata, show_outliers=TRUE, zvalue_thresh = 2)
#' @return a diagnostic plot that shows genes with zero inflation
#' @export
hippo_diagnostic_plot = function(sce,
                                 show_outliers = FALSE,
                                 zvalue_thresh = 10) {
  df = preprocess_heterogeneous(sce@assays@data$counts)
  df = compute_test_statistic(df)
  subset = df[which(df$zvalue > zvalue_thresh), ]
  g = ggplot2::ggplot(df, ggplot2::aes(x = .data$gene_mean,
                                       y = .data$zero_proportion)) +
    ggplot2::geom_point(size = 0.4, alpha = 0.5, na.rm=TRUE) +
    ggplot2::geom_line(ggplot2::aes(x = .data$gene_mean,
                                    y = exp(-.data$gene_mean)),
                       col = "black",
                       na.rm=TRUE) +
    ggplot2::xlim(c(0, 10)) +
    ggplot2::theme_bw() +
    ggplot2::ylab("zero proportion") +
    ggplot2::xlab("gene mean")
  if (show_outliers) {
    g = g + ggplot2::geom_point(data = subset,
                                ggplot2::aes(x = .data$gene_mean,
                                             y = .data$zero_proportion),
                                shape = 21, col = "red",
                                na.rm=TRUE)
  }
  gridExtra::grid.arrange(g, nrow = 1, ncol = 1)
}

#' Access hippo object from SingleCellExperiment object.
#'
#' @param sce SingleCellExperiment object
#' @examples
#' data(toydata)
#' set.seed(20200321)
#' toydata = hippo(toydata,K = 10,z_threshold = 1,outlier_proportion = 0.01)
#' hippo_object = get_hippo(toydata)
#' @return hippo object embedded in SingleCellExperiment object
#' @export
get_hippo = function(sce) {
  if ("hippo" %in% names(sce@int_metadata)) {
    return(sce@int_metadata$hippo)
  } else {
    stop("hippo object does not exist")
  }
}

#' Access data from SCE object
#' @param sce SingleCellExperiment object
#' @examples
#' data(toydata)
#' X = get_data_from_sce(toydata)
#' @return count matrix
#' @export
get_data_from_sce = function(sce){
  return(sce@assays@data$counts)
}


#' HIPPO's hierarchical clustering
#'
#' @param sce SingleCellExperiment object
#' @param K number of clusters to ultimately get
#' @param z_threshold numeric > 0 as a z-value threshold
#' for selecting the features
#' @param outlier_proportion numeric between 0 and 1, a cut-off
#' so that when the proportion of important features reach this
#' number, the clustering terminates
#' @param verbose if set to TRUE, it shows progress of the algorithm
#' @examples
#' data(toydata)
#' toydata = hippo(toydata,K = 10,z_threshold = 1,outlier_proportion = 0.01)
#' @return a list of clustering result for each level of k=1, 2, ... K.
#' @export
hippo = function(sce, K = 20,
                 z_threshold = 2,
                 outlier_proportion = 0.001,
                 verbose = TRUE) {
  if (is(sce, "SingleCellExperiment")) {
    X = sce@assays@data$counts
  } else if (is(sce, "matrix")) {
    sce = SingleCellExperiment::SingleCellExperiment(assays=list(counts = sce))
    X = sce@assays@data$counts
  } else {
    stop("input must be either matrix or SingleCellExperiment object")
  }
  if (outlier_proportion > 1 | outlier_proportion < 0) {
    stop("Outlier_proportion must be a number between 0 and 1.
         Default is 5%")
  }
  param = list(z_threshold = z_threshold,
               outlier_proportion = outlier_proportion,
               maxK = K)
  outlier_number = nrow(X) * outlier_proportion
  labelmatrix = matrix(NA, ncol(X), K)
  labelmatrix[, 1] = 1
  eachlevel = list()
  subX = X
  subXind = seq(ncol(X))
  withinss = rep(0, K)
  oldk = 1
  features = list()
  featuredata = list()
  for (k in 2:K) {
    thisk = one_level_clustering(subX, z_threshold)
    if (is.na(thisk$features$gene[1])) {
      if(verbose){
        message("not enough important features left; terminate the procedure")
      }
      labelmatrix = labelmatrix[, seq((k - 1))]
      break
    }
    if (nrow(thisk$features) < outlier_number) {
      if(verbose){
        message("not enough important features; terminate the procedure")
      }
      labelmatrix = labelmatrix[, seq((k - 1))]
      break
    }
    # if(min(table(thisk$km$cluster)) <= 1){
    #   if(verbose){
    #     message("only one cell in a cluster: terminate procedure")
    #   }
    #   labelmatrix = labelmatrix[, seq((k - 1))]
    #   break
    # }
    if (verbose) {message(paste0("K = ", k, ".."))}
    labelmatrix[, k] = labelmatrix[, k - 1]
    labelmatrix[subXind[thisk$km$cluster == 2], k] = k
    oneind = thisk$km$cluster == 1
    twoind = thisk$km$cluster == 2
    if(sum(oneind) >= 2){
      withinss[oldk] = sum(apply(thisk$unscaled_pcs[oneind, ],1, var)^2)
    }else{
      withinss[oldk] = var(as.numeric(thisk$unscaled_pcs[oneind,]))^2
    }
    if(sum(twoind) >= 2){
      withinss[k] = sum(apply(thisk$unscaled_pcs[twoind, ], 1, var)^2)
    }else{
      withinss[k] = var(as.numeric(thisk$unscaled_pcs[twoind,]))^2
    }


    ind = which(table(thisk$km$cluster) <= 5)
    if (length(ind) >= 1){
      valid_indices = seq(k-1)[-ind]
      oldk = which(withinss == max(withinss[valid_indices]))
    }else{
      oldk = which.max(withinss[seq(k-1)])
    }
    if (sum(labelmatrix[, k] == oldk) < 2) {
      if(verbose){
        message("too few cells in one cluster; terminating the procedure")
      }
      labelmatrix = labelmatrix[, seq(k)]
      break
    }
    subX = X[thisk$features$gene, which(labelmatrix[, k] == oldk)]
    subXind = which(labelmatrix[, k] == oldk)
    thisk$features$subsetK = oldk
    thisk$features$K = k
    features[[k - 1]] = thisk$features
  }
  sce@int_metadata$hippo = list(X = X,
                                features = features,labelmatrix = labelmatrix,
                                z_threshold = z_threshold, param = param,
                                outlier_proportion = outlier_proportion)
  return(sce)
  }


#' visualize each round of hippo through zero proportion plot
#' @param sce SingleCellExperiment object with hippo element in it
#' @param switch_to_hgnc boolean argument to indicate whether to change the gene
#'  names from ENSG IDs to HGNC symbols
#' @param ref a data frame with hgnc column and ensg column
#' @param k select rounds of clustering that you would like to see result.
#' Default is 1 to K
#' @param plottitle Title of your plot output
#' @param top.n number of top genes to show the name
#' @param pointsize size of the ggplot point
#' @param pointalpha transparency level of the ggplot point
#' @param textsize text size of the resulting plot
#' @examples
#' data(toydata)
#' set.seed(20200321)
#' toydata = hippo(toydata,K = 10,z_threshold = 1,outlier_proportion = 0.01)
#' data(ensg_hgnc)
#' zero_proportion_plot(toydata, switch_to_hgnc = TRUE, ref = ensg_hgnc)
#' @return a ggplot object that shows the zero proportions for each round
#' @export
zero_proportion_plot = function(sce,
                                switch_to_hgnc = FALSE,
                                ref = NA,
                                k = NA,
                                plottitle = "",
                                top.n = 5,
                                pointsize = 0.5,
                                pointalpha = 0.5,
                                textsize = 3) {
  df = do.call(rbind, sce@int_metadata$hippo$features)
  topz = df %>% dplyr::group_by(K) %>% dplyr::arrange(desc(zvalue)) %>%
    dplyr::slice(seq_len(5))
  featurelength = as.numeric(table(df$K))
  df$featurecount = featurelength[df$K - 1]
  if (is.na(k[1])) {
    k = 2:ncol(sce@int_metadata$hippo$labelmatrix)
  } else {
    df = df[df$K %in% k, ]
    topz = topz[topz$K %in% k, ]
  }
  topz$hgnc = topz$gene
  if (switch_to_hgnc) {
    topz$hgnc = as.character(ref$hgnc[match(topz$gene, ref$ensg)])
    topz$hgnc[is.na(topz$hgnc)] = topz$gene[is.na(topz$hgnc)]
  }
  g = ggplot2::ggplot(df,ggplot2::aes(x = .data$gene_mean,
                                      y = .data$zero_proportion)) +
    ggplot2::geom_point(size = pointsize, alpha = pointalpha) +
    ggplot2::facet_wrap(~.data$K, ncol = 4) +
    ggplot2::geom_line(ggplot2::aes(x = .data$gene_mean,
                                    y = exp(-.data$gene_mean)),
                       col = "black") +
    ggplot2::xlim(c(0,10)) +
    ggrepel::geom_label_repel(data = topz,
                              ggplot2::aes(label = .data$hgnc),
                              size = textsize, col = "black") +
    ggplot2::geom_text(ggplot2::aes(label =
                                      paste0(featurecount,"genes"),
                                    x = 8, y = 0.8),
                       check_overlap = TRUE,col = "red",size = textsize) +
    ggplot2::theme(legend.position = "none") + ggplot2::theme_bw() +
    ggplot2::ylab("Zero Proportion of Selected Features") +
    ggplot2::xlab("Gene Mean") +
    ggplot2::guides(colour =
                      ggplot2::guide_legend(override.aes =
                                              list(size = 5,alpha = 1),
                                            shape = 19)) +
    ggplot2::theme(axis.text.x =
                     ggplot2::element_text(angle = 45,hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "none", strip.placement = "inside") +
    ggplot2::ggtitle(plottitle) +
    ggplot2::scale_color_manual(values = c("black", "red"))
  gridExtra::grid.arrange(g, nrow = 1, ncol = 1)
}



#' compute t-SNE or umap of each round of HIPPO
#' @param sce SingleCellExperiment object with hippo object in it.
#' @param method a string that determines the method for dimension
#' reduction: either 'umap' or 'tsne
#' @param perplexity numeric perplexity parameter for Rtsne function
#' @param featurelevel the round of clustering that you will extract
#' features to reduce the dimension
#' @return a data frame of dimension reduction result for each
#' k in 1, ..., K
#' @examples
#' data(toydata)
#' set.seed(20200321)
#' set.seed(20200321)
#' toydata = hippo(toydata,K = 10,z_threshold = 1,outlier_proportion = 0.01)
#' toydata = hippo_dimension_reduction(toydata, method="tsne")
#' hippo_tsne_plot(toydata)
#' @export
hippo_dimension_reduction = function(sce, method = c("umap", "tsne"),
                                     perplexity = 30,
                                     featurelevel = 1) {
  hippo_object = sce@int_metadata$hippo
  dflist = list()
  K = ncol(hippo_object$labelmatrix)
  if (method == "umap"){
    dimred = umap::umap(log(t(hippo_object$X[hippo_object$features[[1]]$gene,
                                             ]) + 1))$layout
  }else{
    dimred = tsne = Rtsne::Rtsne(log(t(hippo_object$X[hippo_object$features[[1]]$gene,
                                                      ]) + 1), perplexity = perplexity,
                                 check_duplicates = FALSE)$Y
  }
  dimred = as.data.frame(dimred)
  dimreddf = data.frame()
  for (i in 2:K){
    df = preprocess_homogeneous(sce, label = hippo_object$labelmatrix[,i])
    df$selected_feature = df$gene %in% hippo_object$features[[i -1]]
    df$K = i
    dflist[[i]] = df
    dimreddf = rbind(dimreddf,
                     data.frame(dim1 = dimred$V1, dim2 = dimred$V2,
                                K = i,
                                label = hippo_object$labelmatrix[, i]))
  }
  if (method == "umap"){
    sce@int_metadata$hippo$umap = NA
    colnames(dimreddf) = c("umap1", "umap2", "K", "label")
    dimreddf$label = as.factor(dimreddf$label)
    sce@int_metadata$hippo$umap = dimreddf
  }else{
    sce@int_metadata$hippo$tsne = NA
    colnames(dimreddf) = c("tsne1", "tsne2", "K", "label")
    dimreddf$label = as.factor(dimreddf$label)
    sce@int_metadata$hippo$tsne = dimreddf
  }
  return(sce)
}

#' visualize each round of hippo through UMAP
#'
#' @param sce SingleCellExperiment object with hippo and
#' UMAP result in it
#' @param k number of rounds of clustering that you'd like to see result.
#' Default is 1 to K
#' @param pointsize size of the point for the plot (default 0.5)
#' @param pointalpha transparency level of points for the plot (default 0.5)
#' @param plottitle title of the resulting plot
#' @return ggplot object for umap in each round
#' @examples
#' data(toydata)
#' set.seed(20200321)
#' toydata = hippo(toydata,K = 10,z_threshold = 1,outlier_proportion = 0.01)
#' toydata = hippo_dimension_reduction(toydata, method="umap")
#' hippo_umap_plot(toydata)
#' @export
hippo_umap_plot = function(sce,
                           k = NA,
                           pointsize = 0.5,
                           pointalpha = 0.5,
                           plottitle = "") {
  if (is.na(k[1])) {
    k = seq(1, ncol(sce@int_metadata$hippo$labelmatrix))
  }
  umdf = sce@int_metadata$hippo$umap
  umdf = umdf %>% dplyr::filter(K %in% k)
  if (length(umdf)) {
    g = ggplot2::ggplot(umdf,
                        ggplot2::aes(x = .data$umap1,y = .data$umap2,
                                     col = .data$label)) +
      ggplot2::facet_wrap(~.data$K, ncol = 4) +
      ggplot2::geom_point(size = pointsize,
                          alpha = pointalpha) +
      ggplot2::theme_bw() +
      ggplot2::ylab("umap2") + ggplot2::xlab("umap1") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0,hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black"),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     legend.position = "none", strip.placement = "inside") +
      ggplot2::guides(colour =
                        ggplot2::guide_legend(override.aes =
                                                list(size = 5,alpha = 1))) +
      ggplot2::ggtitle(plottitle)
    gridExtra::grid.arrange(g, nrow = 1, ncol = 1)
  } else {
    stop("use dimension_reduction to compute umap first")
  }
}

#' visualize each round of hippo through t-SNE
#' @param sce SincleCellExperiment object with hippo and t-SNE
#' result in it
#' @param k number of rounds of clustering that you'd like to see result.
#' Default is 1 to K
#' @param pointsize size of the point for the plot (default 0.5)
#' @param pointalpha transparency level of points for the plot (default 0.5)
#' @param plottitle title for the ggplot output
#' @return ggplot object for t-SNE in each round
#' @examples
#' data(toydata)
#' set.seed(20200321)
#' toydata = hippo(toydata,K = 10,z_threshold = 1,outlier_proportion = 0.01)
#' toydata = hippo_dimension_reduction(toydata, method="tsne")
#' hippo_tsne_plot(toydata)
#' @export
hippo_tsne_plot = function(sce,
                           k = NA,
                           pointsize = 0.5,
                           pointalpha = 0.5,
                           plottitle = "") {
  if (is.na(k[1])) {
    k = seq(1, ncol(sce@int_metadata$hippo$labelmatrix))
  }
  tsnedf = sce@int_metadata$hippo$tsne
  tsnedf = tsnedf %>% dplyr::filter(K %in% k)
  if (length(tsnedf)) {
    g = ggplot2::ggplot(tsnedf,
                        ggplot2::aes(x = .data$tsne1, y = .data$tsne2,
                                     col = .data$label)) +
      ggplot2::facet_wrap(~.data$K, ncol = 4) +
      ggplot2::geom_point(size = pointsize, alpha = pointalpha) +
      ggplot2::theme_bw() +
      ggplot2::ylab("tsne2") + ggplot2::xlab("tsne1") +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::guides(colour =
                        ggplot2::guide_legend(override.aes = list(size = 5,
                                                                  alpha = 1))) +
      ggplot2::xlab("TSNE1") + ggplot2::ylab("TSNE2") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0,hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black"),
                     legend.title = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     legend.position = "none", strip.placement = "inside") +
      ggplot2::ggtitle(plottitle)
    gridExtra::grid.arrange(g, nrow = 1, ncol = 1)
  } else {
    stop("use dimension_reduction to compute tsne first")
  }
}

#' visualize each round of hippo through t-SNE
#' @param sce SincleCellExperiment object with hippo and t-SNE result in it
#' @param k number of rounds of clustering that you'd like to see result.
#' Default is 1 to K
#' @param pointsize size of the point for the plot (default 0.5)
#' @param pointalpha transparency level of points for the plot (default 0.5)
#' @param plottitle title for the ggplot
#' @return ggplot for pca in each round
#' @examples
#' data(toydata)
#' set.seed(20200321)
#' toydata = hippo(toydata, K = 10,z_threshold = 1)
#' hippo_pca_plot(toydata, k = 2:3)
#' @export
hippo_pca_plot = function(sce,
                          k = NA,
                          pointsize = 0.5,
                          pointalpha = 0.5,
                          plottitle = "") {
  if (is.na(k[1])) {
    k = seq(2, ncol(sce@int_metadata$hippo$labelmatrix))
  }
  hippo_object = sce@int_metadata$hippo
  counts = get_data_from_sce(sce)
  pc = irlba::irlba(log(counts[hippo_object$features[[1]]$gene,
                               ] + 1), v = 2)$v
  pcadf = data.frame()
  for (kk in k) {
    pcadf = rbind(pcadf, data.frame(PC1 = pc[, 1], PC2 = pc[, 2],K = kk,
                                    label =
                                      sce@int_metadata$hippo$labelmatrix[, kk]))
  }
  pcadf$label = as.factor(pcadf$label)
  pcadf$K = as.factor(pcadf$K)
  g = ggplot2::ggplot(pcadf,
                      ggplot2::aes(x = .data$PC1,
                                   y = .data$PC2,
                                   col = .data$label)) +
    ggplot2::facet_wrap(~.data$K, ncol = 4) +
    ggplot2::geom_point(size = pointsize, alpha = pointalpha) +
    ggplot2::theme_bw() +
    ggplot2::ylab("PC2") + ggplot2::xlab("PC1") +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::guides(colour =
                      ggplot2::guide_legend(override.aes =
                                              list(size = 5,alpha = 1))) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"),
                   legend.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "none", strip.placement = "inside") +
    ggplot2::ggtitle(plottitle)
  gridExtra::grid.arrange(g, nrow = 1, ncol = 1)
}


#' HIPPO's differential expression
#'
#' @param sce SingleCellExperiment object with hippo
#' @param top.n number of markers to return
#' @param switch_to_hgnc if the current gene names are ensemble ids, and would
#' like to switch to hgnc
#' @param ref a data frame with columns 'hgnc' and 'ensg' to match each other,
#' only required when switch_to_hgnc is set to TRUE
#' @param k number of rounds of clustering that you'd like to see result.
#' Default is 1 to K
#' @param plottitle title of the resulting plot
#' @return list of differential expression result
#' @examples
#' data(toydata)
#' set.seed(20200321)
#' toydata = hippo(toydata,K = 10,z_threshold = 1,outlier_proportion = 0.01)
#' result = hippo_diffexp(toydata)
#' @export
hippo_diffexp = function(sce,
                         top.n = 5,
                         switch_to_hgnc = FALSE,
                         ref = NA,
                         k = NA,
                         plottitle = "") {
  if (switch_to_hgnc & length(ref) < 2) {
    stop("A reference must be provided in order to match
         ENSG ids to HGNC symbols")
  }
  hippo_object = sce@int_metadata$hippo
  if (is.na(k[1])) {k = seq(2,ncol(hippo_object$labelmatrix))}
  param = hippo_object$param
  featureind = cellind = result = list()
  featureind[[1]] = seq(nrow(hippo_object$X))
  cellind[[1]] = seq(ncol(hippo_object$X))
  labelmatrix = hippo_object$labelmatrix
  count = hippo_object$X
  finalnewcount = data.frame()
  ind = 1
  for (kk in k) {
    features = hippo_object$features[[ind]]
    cellind = which(labelmatrix[, kk - 1] ==
                      labelmatrix[which(labelmatrix[,kk - 1] !=
                                          labelmatrix[, kk])[1], kk - 1])
    types = unique(hippo_object$labelmatrix[cellind, kk])
    cellgroup1 = which(hippo_object$labelmatrix[, kk] == types[1])
    cellgroup2 = which(hippo_object$labelmatrix[, kk] == types[2])
    rowdata = diffexp_subfunction(count,features, cellgroup1, cellgroup2)
    topgenes = rowdata$genes[seq(top.n)]
    tmpx = cbind(count[rowdata$genes[seq(top.n)],cellgroup1],
                 count[rowdata$genes[seq(top.n)], cellgroup2])
    newcount = as.data.frame(Matrix::t(log(tmpx+1))[,topgenes])
    if (switch_to_hgnc) {
      colnames(newcount) = ref$hgnc[match(colnames(newcount), ref$ensg)]
    }
    newcount$celltype = c(rep(types[1], length(cellgroup1)),
                          rep(types[2],length(cellgroup2)))
    newcount = reshape2::melt(newcount, id = "celltype")
    newcount$celltype = as.factor(newcount$celltype)
    colnames(newcount) = c("celltype", "gene", "logcount")
    newcount$round = paste0("K = ", kk)
    finalnewcount = rbind(finalnewcount, newcount)
    result[[ind]] = rowdata
    ind = ind + 1
  }
  sce@int_metadata$hippo$diffexp$result_table = result
  finalnewcount$round = factor(finalnewcount$round,
                               levels = paste0("K = ", k))
  g = ggplot2::ggplot(finalnewcount,
                      ggplot2::aes(x = .data$gene,
                                   y = exp(.data$logcount) - 1,
                                   col = .data$celltype)) +
    ggplot2::facet_wrap(~round,scales = "free", ncol = 4) +
    ggplot2::geom_boxplot(outlier.size = 0.2) +
    ggplot2::theme_classic()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,hjust=1),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"),
                   strip.background = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::ylab("UMI count") + ggplot2::xlab("") +
    ggplot2::scale_y_continuous(trans = "log1p",breaks = c(0, 10, 100, 1000)) +
    ggplot2::ggtitle(plottitle)
  gridExtra::grid.arrange(g, nrow = 1, ncol = 1)
  sce@int_metadata$hippo$diffexp$plot = g
  return(sce)
  }

diffexp_subfunction = function(count, features, cellgroup1, cellgroup2){
  rowdata = data.frame(genes = features$gene)
  tmpcount1 = count[features$gene, cellgroup1]
  tmpcount2 = count[features$gene, cellgroup2]
  if(length(cellgroup1) == 1){
    tmpmean1 = mean(tmpcount1)
  }else{
    tmpmean1 = Matrix::rowMeans(tmpcount1)
  }
  if(length(cellgroup2) == 1){
    tmpmean2 = mean(tmpcount2)
  }else{
    tmpmean2 = Matrix::rowMeans(tmpcount2)
  }
  rowdata$meandiff = tmpmean1 - tmpmean2
  rowdata$sd = sqrt(tmpmean1/length(cellgroup1) + tmpmean2/length(cellgroup2))
  rowdata$z = rowdata$meandiff/rowdata$sd
  rowdata = rowdata[order(rowdata$z, decreasing = TRUE), ]
  rowdata$genes = as.character(rowdata$genes)
  return(rowdata)
}

#' HIPPO's feature heatmap
#'
#' @param sce SingleCellExperiment object with hippo
#' @param top.n number of markers to return
#' @param switch_to_hgnc if the current gene names are ensemble ids, and would
#' like to switch to hgnc
#' @param ref a data frame with columns 'hgnc' and 'ensg' to match each other,
#' only required when switch_to_hgnc is set to TRUE
#' @param kk integer for the round of clustering that you'd like to see result.
#' Default is 2
#' @param plottitle title for the plot
#' @return list of differential expression result
#' @examples
#' data(toydata)
#' set.seed(20200321)
#' toydata = hippo(toydata,K = 10,z_threshold = 1,outlier_proportion = 0.01)
#' hippo_feature_heatmap(toydata)
#' @export
hippo_feature_heatmap = function(sce,
                                 switch_to_hgnc = FALSE,
                                 ref = NA,
                                 top.n = 50,
                                 kk = 2,
                                 plottitle = "") {
  if (switch_to_hgnc & length(ref) < 2) {
    stop("A reference must be provided to match ENSG ids to HGNC symbols")
  }
  hippo_object = sce@int_metadata$hippo
  labelmatrix = as.data.frame(hippo_object$labelmatrix)
  labelmatrix$barcode = colnames(hippo_object$X)
  bigX = data.frame()
  feat = hippo_object$features[[kk - 1]]
  feat = feat %>% dplyr::arrange(desc(zvalue))
  feat = feat[seq(top.n), ]
  lab = sce@int_metadata$hippo$labelmatrix[, kk]
  tmp = log(hippo_object$X[feat$gene, ] + 1)
  tmp = as.data.frame(tmp)
  tmp$gene = feat$gene
  tmp$hgnc = tmp$gene
  if (switch_to_hgnc) {
    tmp$hgnc = as.character(ref$hgnc[match(tmp$gene, ref$ensg)])
    tmp$hgnc[is.na(tmp$hgnc)] = tmp$gene[is.na(tmp$hgnc)]
  }
  X = reshape2::melt(tmp, id = c("hgnc", "gene"))
  X$label = labelmatrix[match(X$variable, labelmatrix$barcode), kk]
  X$K = kk - 1
  X$value = as.numeric(X$value)
  g = ggplot2::ggplot(X,
                      ggplot2::aes(x = .data$variable, y = .data$hgnc,
                                   fill = .data$value)) +
    ggplot2::facet_grid(~.data$label, scales = "free_x") +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(high = "darkred",low = "white") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   strip.placement = "inside",
                   legend.position = "none",
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::ggtitle(plottitle)
  gridExtra::grid.arrange(g, nrow = 1, ncol = 1)
}


ensg_to_hgnc = function(ensg) {
  data(ensg_hgnc)
  maps = ensg_hgnc
  maps2 = data.frame(ensg = ensg, hgnc = maps$hgnc[match(ensg, maps$ensembl)])
  maps2$ensg = as.character(maps2$ensg)
  maps2$hgnc = as.character(maps2$hgnc)
  ind_na = which(is.na(maps2$hgnc))
  ind_blank = which(maps2$hgnc == "")
  hgnc = maps2$hgnc
  hgnc[c(ind_na, ind_blank)] = maps2$ensg[c(ind_na, ind_blank)]
  return(hgnc)
}



#' Return hippo_diffexp object
#'
#' @param sce SingleCellExperiment object with hippo
#' @param k integer round of result of interest
#' @return data frame of differential expression test
#' @examples
#' data(toydata)
#' set.seed(20200321)
#' toydata = hippo(toydata,K = 10,z_threshold = 1,outlier_proportion = 0.01)
#' toydata = hippo_diffexp(toydata)
#' result1 = get_hippo_diffexp(toydata)
#' @export
get_hippo_diffexp = function(sce, k=1){
  hippo_object = get_hippo(sce)
  return(hippo_object$diffexp$result_table[[k]])
}
