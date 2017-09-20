# reads in Supplemental Table 3 or 4 (in txt format), and fits multinomial or binomial logistic regression models to X- or Z-linked genes, with or without 
# Pct score as a predictor
library(nnet)
ancX_forlogit = read.table("PATH_TO_SUPPTABLE_S3",sep="\t",fill=TRUE,header = TRUE,stringsAsFactors = FALSE)
ancX_forlogit$evo_class = NA
ancX_forlogit[which(ancX_forlogit$gene %in% subset(ancX_forlogit,Human.Y==TRUE)$gene),"evo_class"] = "XY"
ancX_forlogit[which(ancX_forlogit$gene %in% subset(ancX_forlogit,(Human.Y==FALSE)&(XCI.Status =="I"))$gene),"evo_class"] = "I"
ancX_forlogit[which(ancX_forlogit$gene %in% subset(ancX_forlogit,(Human.Y==FALSE)&(XCI.Status  %in% c("E","VE")))$gene),"evo_class"] = "E"
mn_x_nomir = multinom(evo_class~Haploinsufficiency.probability+Human.expression.breadth+Human-mouse.dN/dS,data=na.omit(ancX_forlogit[,c("evo_class","Haploinsufficiency.probability","Human.expression.breadth","Mean.Pct","Human-mouse.dN/dS")]))
mn_x_pctmir = multinom(evo_class~Haploinsufficiency.probability+Human.expression.breadth+Human-mouse.dN/dS+Mean.Pct,data=na.omit(ancX_forlogit[,c("evo_class","Haploinsufficiency.probability","Human.expression.breadth","Mean.Pct","Human-mouse.dN/dS")]))

ancZ_forlogit = read.table("PATH_TO_SUPPTABLE_S3",sep="\t",fill=TRUE,header = TRUE,stringsAsFactors = FALSE)
ancZ_forlogit$evo_class = NA
ancZ_forlogit[which(ancZ_forlogit$gene %in% subset(ancZ_forlogit,chicken.W=="Y")$gene),"evo_class"] = "ZW"
ancZ_forlogit[which(ancZ_forlogit$gene %in% subset(ancZ_forlogit,chicken.W=="N")$gene),"evo_class"] = "AncZ"
mn_z_nomir = multinom(evo_class~Haploinsufficiency.probability+Human.expression.breadth+Chicken-flycatcher.dN/dS,data=na.omit(ancZ_forlogit[,c("evo_class","Haploinsufficiency.probability","Human.expression.breadth","Mean.Pct","Chicken-flycatcher.dN/dS")]))
mn_z_pctmir = multinom(evo_class~Haploinsufficiency.probability+Human.expression.breadth+Chicken-flycatcher.dN/dS+Mean.Pct,data=na.omit(ancZ_forlogit[,c("evo_class","Haploinsufficiency.probability","Human.expression.breadth","Mean.Pct","Chicken-flycatcher.dN/dS")]))
