# reads in gene count and abundance files, along with sample info, and uses limma/voom to calculate chicken male/female ratios
# (normalized to ancestral human and anolis male/female ratios) in heart, liver, kidney, and brain

library(edgeR)
library(limma)

counts = read.table("human.chicken.anolis.genecounts.txt",header=TRUE,stringsAsFactors = FALSE)
abundance = read.table("human.chicken.anolis.genetpm.txt",header=TRUE,stringsAsFactors = FALSE)
combined_table = read.table("human.chicken.anolis.sampleinfo.txt",header = TRUE,stringsAsFactors = FALSE)
table_heart = subset(combined_table,source_name_s=="Heart")
table_brain = subset(combined_table,source_name_s=="Brain")
table_liver = subset(combined_table,source_name_s=="Liver")
table_kidney = subset(combined_table,source_name_s=="Kidney")

heart_tmm = DGEList(counts[which(apply(abundance[,table_heart$Run_s] > 1,1,sum) >= 12),table_heart$Run_s])
heart_tmm = calcNormFactors(heart_tmm)
v = voom(heart_tmm,design=model.matrix(~Sex_s*Species,table_heart),plot=TRUE)
heart_fit = eBayes(lmFit(v,design=model.matrix(~Sex_s*Species,table_heart)))
heart_fit = eBayes(contrasts.fit(heart_fit,as.matrix(c(0,0,0,0,1,-1))))
heart_results = topTable(heart_fit,coef=1,number = Inf)
colnames(heart_results)[1] = "logFC_chickenheart_anolishuman"

brain_tmm = DGEList(counts[which(apply(abundance[,table_brain$Run_s] > 1,1,sum) >= 11),table_brain$Run_s])
brain_tmm = calcNormFactors(brain_tmm)
v = voom(brain_tmm,design=model.matrix(~Sex_s*Species,table_brain),plot=TRUE)
brain_fit = eBayes(lmFit(v,design=model.matrix(~Sex_s*Species,table_brain)))
brain_fit = eBayes(contrasts.fit(brain_fit,as.matrix(c(0,0,0,0,1,-1))))
brain_results = topTable(brain_fit,coef=1,number = Inf)
colnames(brain_results)[1] = "logFC_chickenbrain_anolishuman"

kidney_tmm = DGEList(counts[which(apply(abundance[,table_kidney$Run_s] > 1,1,sum) >= 12),table_kidney$Run_s])
kidney_tmm = calcNormFactors(kidney_tmm)
v = voom(kidney_tmm,design=model.matrix(~Sex_s*Species,table_kidney),plot=TRUE)
kidney_fit = eBayes(lmFit(v,design=model.matrix(~Sex_s*Species,table_kidney)))
kidney_fit = eBayes(contrasts.fit(kidney_fit,as.matrix(c(0,0,0,0,1,-1))))
kidney_results = topTable(kidney_fit,coef=1,number = Inf)
colnames(kidney_results)[1] = "logFC_chickenkidney_anolishuman"

liver_tmm = DGEList(counts[which(apply(abundance[,table_liver$Run_s] > 1,1,sum) >= 12),table_liver$Run_s])
liver_tmm = calcNormFactors(liver_tmm)
v = voom(liver_tmm,design=model.matrix(~Sex_s*Species,table_liver),plot=TRUE)
liver_fit = eBayes(lmFit(v,design=model.matrix(~Sex_s*Species,table_liver)))
liver_fit = eBayes(contrasts.fit(liver_fit,as.matrix(c(0,0,0,0,1,-1))))
liver_results = topTable(liver_fit,coef=1,number = Inf)
colnames(liver_results)[1] = "logFC_chickenliver_anolishuman"


