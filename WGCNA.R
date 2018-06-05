
## This is an r script for WGCNA analysis.
## Author: Du Yang
## Date: 2018-6-3

library(WGCNA)
library(CancerSubtypes)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

protDataPath<-'/mnt/dellfs/projects/gc_proteogenomics/MS/Assay/Result/2018_5_30/'

data<-read.table(paste0(protDataPath,'afterBE/Tumor.meanRep.annotated.0.5_Imputed.normalized.BEremoved.data.mean.txt')
	,header=T,sep="\t")

row.names(data)<-data[,1]
data<-data[,-1]

data.removed<-read.table("/mnt/dellfs/projects/gc_proteogenomics/MS/Assay/metafile/removed_sample.txt",header=T,sep="\t")
row.names(data.removed)<-data.removed[,1]
d<-data[,!(colnames(data) %in% row.names(data.removed[data.removed$Integrative_clustering==1,]))]
d<-d[,!(colnames(d) %in% row.names(data.removed[data.removed$Adjust==2,]))]

d=FSbyVar(d, cut.type="topk",value=2504)  ## top 80% by variance. 

datExpr<-t(d)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut=0.9,verbose = 5)

par(mfrow = c(1,2));
cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

net = blockwiseModules(
				 datExpr,
				 power = sft$powerEstimate,
				 TOMType = "unsigned", minModuleSize = 30,
				 reassignThreshold = 0, mergeCutHeight = 0.15,
				 numericLabels = TRUE, pamRespectsDendro = FALSE,
				 saveTOMs = FALSE,
				 verbose = 3
 )
 table(net$colors) 

 sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
pdf("Top2504genesByVar.pdf")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                  "Module colors",
                  dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05)
dev.off()

