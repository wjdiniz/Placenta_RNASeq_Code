---
title: "DEA_Placenta_Forage_concentrate_Project_SDSU"
author: "Wellison Diniz"
date: "01/26/2026"
---

Code used to analyse the data for Manuscript: 
---
Disentangling Dietary Substrate from Energy Intake Level: 
Effects of Iso-Energetic Limit-Feeding on Placental Gene 
Expression, Lactation, and Performance in Beef Cows 
By Brander et al. 2026
---

#libraries 

library(edgeR)
library(limma)
library(DESeq2)
library(factoextra)
library(pheatmap)
library(plyr)
library(rafalib)
library(biomaRt)
library(janitor)

#Datasets


#Quality Control**

#Set directory

#setwd("~/path/")


#Expression Set


mRNAExp <- read.table("~/Read_counts_placenta.txt", header = T, row.names = 1)

dim(mRNAExp)

mRNAExp[1:5,1:5]



#ColData object

#Create a metadata file with the info for the samples. 


ColData <- read.table("~/ColData.txt", header=T, row.names = 1, sep = '\t')
dim(ColData)
ColData[1:length(ColData)]


#Match mRNA with group Data 


#We examine the count matrix(mRNAExp) and column from ColData data to see if they are consistent in terms of sample order.


all(rownames(ColData) %in% colnames(mRNAExp))


#Are they on the same order?

all(rownames(ColData) == colnames(mRNAExp))



#if false, let’s order it


order.sample.cols<- match(rownames(ColData), colnames(mRNAExp))

mRNAExp<-mRNAExp[, order.sample.cols]

#Are they on the same order?

all(rownames(ColData) == colnames(mRNAExp))


#Distribution of the animals per treatment

table(ColData$Treatment)
tabyl(ColData, Calf_Sex, Treatment)


#Filtering to remove low counts

design <- model.matrix(~ColData$Treatment)
keep <- filterByExpr(mRNAExp, design)
table(keep)
mRNAExp.fil<-mRNAExp[keep, ]
dim(mRNAExp.fil)


#PCA

pca.counts <- prcomp(t(cpm(mRNAExp.fil, log=TRUE)),scale=TRUE)

pdf(file = "PCA_cpm.pdf", width = 10, height = 8)
par(mfrow = c(1, 1))
par(mar = c(8,6,10,6))
par (cex=2)
par(las=2)

fviz_pca_ind(pca.counts, 
             col.var= "contrib", #color by the contributions to the PC
             gradient.colors= c("#00AFBB", "E7B800","FC4EO7"),
             repel=FALSE,habillage =ColData$Treatment, title="PCA by Treatment", addEllipses = F)

fviz_pca_ind(pca.counts, 
             col.var= "contrib", #color by the contributions to the PC
             gradient.colors= c("#00AFBB", "E7B800","FC4EO7"),
             repel=FALSE,habillage =ColData$Calf_Sex, title="PCA by Calf Sex", addEllipses = F)
dev.off()


#Converting treatment to factors

ColData$Treatment<- as.factor(ColData$Treatment)


#DEA - Calf_sex in the model 

#Create object to work with DESeq2


dds<- DESeqDataSetFromMatrix(countData = mRNAExp.fil, colData = ColData, design = ~Calf_Sex+Treatment)
dds


#Run DESeq function 

dds <- DESeq(dds)

#Annotation

#Retrieving  gene annotations from ensembl biomaRt

ensembl = useEnsembl(biomart="ensembl", dataset="btaurus_gene_ensembl")
Annot<-getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype'), 
             filters = 'ensembl_gene_id', values =rownames(dds), mart = ensembl)
dim(Annot)

#Run contrasts

Res.HFvsHC <- results(dds, contrast=c("Treatment","HF" ,"HC"), alpha = 0.05, pAdjustMethod="BH", tidy = F)

summary(Res.HFvsHC)

Res.HFvsHC


Res.HFvsHC <- results(dds, contrast=c("Treatment","HF" ,"HC"), alpha = 0.05, pAdjustMethod="BH", tidy = T)

#Annotation

idx <- match(Res.HFvsHC$row, Annot$ensembl_gene_id)

Res.HFvsHC$Gene.Name<-Annot$external_gene_name[ idx ]

Res.HFvsHC$Gene.Biotype<-Annot$gene_biotype[ idx ]

dim(Res.HFvsHC)


#Write down the results

write.csv(Res.HFvsHC,"DEGs_Res.HFvsHC_Placenta.csv")

#Ranking genes for GSEA

#FROM RESULTS DESEQ2 OBJECT 

GSEA_transcriptome<-Res.HFvsHC[,c(1, 3, 6)]

#Ranking equation

rankings <- sign(Res.HFvsHC$log2FoldChange)*(-log10(Res.HFvsHC$pvalue))
names(rankings) <- Res.HFvsHC$row
head(rankings)

rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

#Adding the info
rankings<-as.data.frame(rankings)

write.csv(rankings,"~/RanK_GSEA.csv")


sessionInfo()


