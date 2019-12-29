#' Expected zero proportion under Poisson
#'
#' @param lambda gene mean
#' @return expected proportion of zero under poisson
pois_prob_zero = function(lambda){
  exp(-lambda)
}
#' Expected zero proportion under Negative Binomial
#'
#' @param lambda gene mean
#' @param theta dispersion parameter, 0 if poisson
#' @param theta expected proportion of zero under poisson
#' @return Expected zero proportion under Negative Binomial
nb_prob_zero = function(lambda, theta){
  if(theta==0){
    return(exp(-lambda))
  }else{
    return((1/(mu * theta + 1))^(1/theta))
  }
}
#' Expected zero proportion under Negative Binomial
#'
#' @param lambda gene mean
#' @param theta dispersion parameter, 0 if zero-inflated poisson
#' @param pi zero inflation, 0 if negative binomial
#' @return Expected zero proportion under Zero-Inflated Negative Binomial
zinb_prob_zero = function(lambda, theta, pi){
  if(theta==0){
    return((1-pi) * exp(-lambda) + pi)
  }else{
    return(pi + (1-pi) * (1/(lambda * theta + 1))^(1/theta))
  }
}

#' Preprocess UMI data so that each row contains information about each gene
#'
#' @param X gene by cell matrix
#' @return data frame with one row for each gene.
preprocess_heterogeneous = function(X){
  X = t(X)
  Y = X
  X = X>0 * 1
  numcounts = colSums(X)
  ind = which(numcounts == 0)
  if(length(ind) > 0){
    X = X[, -ind]
    Y = Y[, -ind]
  }
  rm(numcounts)
  det_rate = colMeans(X)
  gene_mean = colMeans(Y)
  gene_var = matrixStats::colVars(Y)

  zero_proportion = 1-det_rate
  rm(Y); gc()

  gene_mean = as.data.frame(gene_mean)
  det_rate = as.data.frame(det_rate)
  zero_proportion = as.data.frame(zero_proportion)
  gene_var = as.data.frame(gene_var)

  df = data.frame(gene = colnames(X),
                  det_rate = det_rate,
                  gene_mean = gene_mean,
                  gene_var = gene_var,
                  zero_proportion = zero_proportion)
  df$samplesize = nrow(X)
  return(df)

}

#' Preprocess UMI data with inferred or known labels
#'
#' @param X gene by cell matrix
#' @param label inferred or known label in factor
#' @param normalize normalize each cell to have the same sequencing depth. Default as FALSE
#' @return data frame with one row for each gene.
preprocess_homogeneous = function(X, label, normalize = FALSE){
  X = t(X)
  sf = median(rowSums(X))
  if(normalize){
    X = X/rowSums(X) * sf
  }
  Y = X
  X = X>0 * 1
  numcounts = colSums(X)
  ind = which(numcounts == 0)
  if(length(ind) > 0){
    X = X[, -ind]
    Y = Y[, -ind]
  }
  rm(numcounts)
  labelnames = as.character(unique(label))
  det_rate = matrix(NA, ncol(X), length(labelnames))
  gene_mean = matrix(NA, ncol(X), length(labelnames))
  positive_mean = matrix(NA, ncol(X), length(labelnames))
  gene_var = matrix(NA, ncol(X), length(labelnames))
  samplesize = table(label)

  for (i in 1:length(labelnames)){
    ind = which(label==labelnames[i])
    det_rate[,i] = colMeans(X[ind, ])
    gene_mean[,i] = colMeans(Y[ind, ])
    gene_var[,i] = matrixStats::colVars(Y[ind,])
  }
  colnames(det_rate) =
    colnames(gene_mean) =
    colnames(gene_var) = labelnames

  zero_proportion = 1-det_rate

  gene_mean = as.data.frame(gene_mean)
  det_rate = as.data.frame(det_rate)
  zero_proportion = as.data.frame(zero_proportion)
  gene_var = as.data.frame(gene_var)

  gene_mean$id = det_rate$id = zero_proportion$id  = gene_var$id = colnames(X)

  mgm = melt(gene_mean, id = "id")
  mdr = melt(det_rate, id = "id")
  mdor = melt(zero_proportion, id = "id")
  mgv = melt(gene_var, id = "id")

  df = data.frame(gene = colnames(X),
                  det_rate = mdr$value,
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

visualize_hippo = function(hippo_object){
  plist = list()
  dflist = list()
  K = ncol(hippo_object$labelmatrix)
  df = preprocess_heterogeneous(hippo_object$X)
  df$celltype = "combined"
  df$selected_feature = FALSE
  df$K = 1
  dflist[[1]] = df
  um = umap::umap(log(t(hippo_object$X[hippo_object$features[[1]], ])+1))
  um = as.data.frame(um$layout)
  tsne = Rtsne::Rtsne(log(t(hippo_object$X[hippo_object$features[[1]], ])+1))
  umdf = data.frame()
  tsnedf = data.frame()
  for (i in 2:K){
    df = preprocess_homogeneous(hippo_object$X, label = hippo_object$labelmatrix[,i])
    df$selected_feature = df$gene %in% hippo_object$features[[i-1]]
    df$K = i
    dflist[[i]] = df
    umdf = rbind(umdf, data.frame(umap1 = um$V1,
                                    umap2 = um$V2,
                                    K = i,
                                    label = hippo_object$labelmatrix[,i]))
    tsnedf = rbind(tsnedf, data.frame(tsne1 = tsne$Y[,1],
                                      tsne2 = tsne$Y[,2],
                                      K=i,
                                      label = hippo_object$labelmatrix[,i]))


  }
  df = do.call(rbind, dflist)
  df = df[sample(nrow(df)), ]
  zero_plot = ggplot(df, aes(x = gene_mean, y = zero_proportion, col = celltype)) +
    geom_point(size = 0.8, alpha = 0.4) +
    facet_wrap(~K) +
    geom_line(aes(x = gene_mean, y = exp(-gene_mean)), col = 'black') +
    xlim(c(0,10))+
    theme(legend.position = "none") +
    theme_bw()

  umap_plot = ggplot(umdf, aes(x = umap1, y= umap2, col = factor(label))) +
    facet_wrap(~K) +
    geom_point(size = 0.8, alpha = 0.4)

  tsne_plot = ggplot(tsnedf, aes(x = tsne1, y=tsne2, col=factor(label))) +
    facet_wrap(~K) +
    geom_point(size=0.8, alpha = 0.4)
  return(list(zero_plot = zero_plot,
              umap_plot = umap_plot,
              tsne_plot = tsne_plot))
}

#' Likelihood ratio test for dispersion parameter = 0
#'
#' @param y a vector of each gene across cells
#' @return p-value for the significance for non-zero dispersion parameter
pois_vs_nb = function(y){
  require(MASS)
  pois = sum(dpois(y, mean(y), log = TRUE))
  if(mean(y) > var(y)){
    return(1)
  }else{
    nb = fitdistrplus::fitdist(y, "nbinom",
                               start = list(size = 0.1, prob = 0.1),
                               lower = c(0, 0), upper = c(10^10, 1))
    chi = 2 * (logLik(nb) - pois)
    p = pchisq(chi, df = 1, lower.tail = FALSE)
  }
  return(p)
}

#' Change ENSG id's to HGNC symbols
#'
#' @param ensg a vector of ENSG names
#' @return a dataframe of both ENSG id's and HGNC symbol
ensg_to_hgnc = function(ensg){
  maps = read.table("~/Work/SC/data/Annotations/hgnc_ensembl.txt", header=TRUE, stringsAsFactors = FALSE)
  maps2 = data.frame(ensg = ensg,
                     hgnc = maps$hgnc[match(ensg, maps$ensembl)])
  maps2$ensg = as.character(maps2$ensg)
  maps2$hgnc = as.character(maps2$hgnc)
  ind_na = which(is.na(maps2$hgnc))
  ind_blank = which(maps2$hgnc=="")
  hgnc = maps2$hgnc
  hgnc[c(ind_na, ind_blank)] = maps2$ensg[c(ind_na, ind_blank)]
  return(hgnc)
}


#' Conduct feature selection
#'
#' @param df pre-processed data frame
#' @return data frame with added columns with test results
compute_test_statistic = function(df){
  ind = which(df$gene_mean==0)
  if(length(ind)>0){
    df = df[-ind,]
  }
  ind = grep("^MT-", df$gene)
  if(length(ind) > 0){
    df = df[-grep("^MT-", df$gene), ]
  }
  require(dplyr)
  df = df %>% mutate(expected_pi= 2*samplesize/(2*samplesize-1) * exp(-gene_mean)) %>%
      mutate(se = sqrt(expected_pi * (1-expected_pi) / (samplesize-1.25))) %>%
      mutate(minus_logp = -pnorm(zero_proportion, expected_pi, se, log.p = TRUE, lower.tail=FALSE)) %>%
      mutate(minus_logp = pmin(minus_logp, 500)) %>%
      mutate(zvalue = -qnorm(exp(-minus_logp)))
  df$gene = as.character(df$gene)
  return(df)
}




#' Clustering of one-step
#'
#' @param subX gene by cell matrix that needs clustering
#' @param z_threshold z-value threshold for feature selection
#' @return a list of clustering result. km object contains the k-means result
one_level_clustering = function(subX, z_threshold){
  subdf = preprocess_heterogeneous(subX)
  subdf = compute_test_statistic(subdf)
  features = subdf$gene[subdf$zvalue>z_threshold]
  if(length(features)<10){
    return(list(features = NA, pcs = NA, km = NA))
  }
  pcs = irlba::irlba(log(subX[features, ]+1), 10)$v
  unscaledpc = prcomp(log(t(subX[features,])+1), scale.=FALSE, center=FALSE)$x[,1:10]
  km = kmeans(pcs, 2, nstart = 500)
  return(list(features = features, pcs = pcs, km = km, unscaled_pcs = unscaledpc, subdf = subdf))
}

#' HIPPO's hierarchical clustering
#'
#' @param X gene by cell matrix
#' @param K number of clusters to ultimately get
#' @return a list of clustering result for each level of k=1, 2, ... K.
hippo = function(X, K=10, z_threshold = 20){
  labelmatrix = matrix(NA, ncol(X), K)
  labelmatrix[,1] = 1
  eachlevel = list()
  subX = X
  subXind = 1:ncol(X)
  withinss = rep(0, K)
  oldk = 1
  features = list()
  featuredata = list()
  for (k in 2:K){
    thisk = one_level_clustering(subX, z_threshold)
    if(is.na(thisk$features[1])){
      print("ran out of important features")
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
  return(list(X = X, features = features, labelmatrix = labelmatrix))
}


#' HIPPO's differential expression
#'
#' @param sce single cell experiment object
#' @param clust hierarchical_clustering object
#' @param top.n number of markers to return
#' @param ref a data frame with columns "hgnc" and "ensg" to match each other
#' @return list of differential expression result
diffexp = function(hippo_object, top.n = 10, switch_to_hgnc=FALSE, ref = NA){
  if(switch_to_hgnc & is.na(ref)){
    stop("A reference must be provided in order to match ENSG ids to HGNC symbols")
  }
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
    }
    newcount$celltype = c(rep(types[1], length(cellgroup1)), rep(types[2], length(cellgroup2)))
    newcount = reshape::melt(newcount, id="celltype")
    newcount$celltype = as.factor(newcount$celltype)
    colnames(newcount) = c("celltype", "gene", "logcount")
    g = ggplot(newcount, aes(x = gene, y = logcount, col = celltype)) +
      geom_boxplot(outlier.size = 0.2) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      ggtitle(paste("Round", k-1)) +
      xlab("") +
      theme(legend.position='bottom')
    plist[[k-1]] = g
    result[[k-1]] = rowdata
  }
  bpl = gridExtra::grid.arrange(grobs = plist, ncol = 2)

  return(list(result_table = result,
              result_boxplot = bpl))
}


#' HIPPO's differential expression
#'
#' @param sce single cell experiment object
#' @param clust hierarchical_clustering object
#' @param top.n number of markers to return
#' @return list of differential expression result
hierarchical_dimred = function(sce, clust, top.n = 10){
  K = ncol(clust$labelmatrix)
  featureind = list()
  cellind = list()
  featureind[[1]] = 1:nrow(rowData(sce))
  cellind[[1]] = 1:nrow(colData(sce))
  labelmatrix = clust$labelmatrix
  rowd = rowData(sce)
  result = list()
  count = sce@assays$data$counts
  for (k in 2:K){
    features = clust$intermediate_data[[k-1]]$features
    cellind = which(labelmatrix[,k-1] == labelmatrix[which(labelmatrix[,k-1] != labelmatrix[,k])[1], k-1])
    types = unique(clust$labelmatrix[cellind, k])
    cellgroup1 = which(clust$labelmatrix[,k] == types[1])
    cellgroup2 = which(clust$labelmatrix[,k] == types[2])
    newset = count[features, c(cellgroup1, cellgroup2)]
    tsne = Rtsne(newset)

    rowdata = data.frame(gene = features)
    rowdata$meandiff = rowMeans(count[features,cellgroup1]) - rowMeans(count[features,cellgroup2])
    rowdata$sd = sqrt(rowMeans(count[features,cellgroup1])/length(cellgroup1) +
                        rowMeans(count[features,cellgroup2])/length(cellgroup2))
    rowdata$z = rowdata$meandiff/rowdata$sd
    rowdata = rowdata[order(rowdata$z, decreasing=TRUE), ]
    newcount = t(log(cbind(count[rowdata$id[1:top.n], cellgroup1],
                           count[rowdata$id[1:top.n], cellgroup2])+1))
    newcount = as.data.frame(newcount)
    colnames(newcount) = rowdata$gene[1:top.n]
    newcount$celltype = c(rep(types[1], length(cellgroup1)), rep(types[2], length(cellgroup2)))
    newcount = reshape::melt(newcount, id="celltype")
    newcount$celltype = as.factor(newcount$celltype)
    colnames(newcount) = c("celltype", "gene", "logcount")
    g = ggplot(newcount, aes(x = gene, y = logcount, col = celltype)) +
      geom_boxplot(outlier.size = 0.2) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      ggtitle(paste("Round", k-1)) +
      xlab("") +
      theme(legend.position='bottom')
    result[[k-1]] = list(features = rowdata,
                         boxplot = g)
  }
  return(result)
}
