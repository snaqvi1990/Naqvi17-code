# given matrix of gene-pct values, data frame of genes and classes, performs nrep gene samples, and calculates aggregate, median, or mean pct
# for each sampling. Returns information on either quantiles or standard deviations of the resulting resampled distributions
samplePcts <- function(gene_pcts, gene_classes,nrep,agg_mean,quantile_stdev) {
  gene_classes = na.omit(gene_classes)
  class_estimates = matrix(nrow = length(unique(gene_classes[,2])),ncol = 4,dimnames = list(as.character(unique(gene_classes[,2])),c("Estimate","Lo","Hi","Wilcox_vs_below")))
  class_vecs = matrix(nrow = length(unique(gene_classes[,2])),ncol = nrep)
  #remove all NA genes/rows
  if (sum(apply(is.na(gene_pcts[,-1]),1,sum)==ncol(gene_pcts[,-1])) > 0) {
    gene_pcts = gene_pcts[-which(apply(is.na(gene_pcts[,-1]),1,sum)==ncol(gene_pcts[,-1])),]
  }
  for (j in 1:nrep){
    for (i in 1:length(unique(gene_classes[,2]))) {
      class = as.character(unique(gene_classes[,2]))[i]
      genes_inclass = gene_classes[which((gene_classes[,2] == class)),1]
      genes_inclass_subsamp = sample(x=genes_inclass,size=length(genes_inclass),replace = TRUE)
      genes_inclass_subsamp_pcts = as.matrix(gene_pcts[which(gene_pcts[,1] %in% genes_inclass_subsamp),-1])
      if (agg_mean == 'mean_pct') {
        class_vecs[i,j] = mean(as.numeric(c(genes_inclass_subsamp_pcts)),na.rm = TRUE)
      }
      if (agg_mean == 'med_pct') {
        class_vecs[i,j] = median(as.numeric(c(genes_inclass_subsamp_pcts)),na.rm = TRUE)
      }
      if (agg_mean == 'agg_pct') {
        vec = c(genes_inclass_subsamp_pcts)[which(!is.na(c(genes_inclass_subsamp_pcts)))]
        #print(vec)
        class_vecs[i,j] = 1 - prod(1-vec)
      }
      if (!(agg_mean %in% c("mean_pct","med_pct","agg_pct"))) {stop()}
    }
  }
  if (quantile_stdev == "quantile") {
    #print(head(class_vecs))
    class_estimates[,"Estimate"] = apply(class_vecs,1,median,na.rm=TRUE)
    class_estimates[,"Lo"] = apply(class_vecs,1,quantile,probs=0.05)
    class_estimates[,"Hi"] = apply(class_vecs,1,quantile,probs=0.95)
  }
  if (quantile_stdev == "stdev") {
    class_estimates[,"Estimate"] = apply(class_vecs,1,mean)
    class_estimates[,"Lo"] = class_estimates[,"Estimate"] - apply(class_vecs,1,sd)
    class_estimates[,"Hi"] = class_estimates[,"Estimate"] + apply(class_vecs,1,sd)
  }
  return(class_estimates)
}