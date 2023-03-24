#!/usr/bin/env Rscript

options(warn=-1)

library(withr)
library(httr)
library(curl)
library(devtools)
library(RUVSeq)

args = commandArgs(TRUE)
matrix <- args[1]
treatments <- args[2]
controls <- args[3]

countData <- read.table(matrix,sep="\t",header=TRUE,row.names=1,check.names=F)

# Find only the genes with the correct format (exclude mRNA or NoGeneId) and convert to gene###### format)
genes <- rownames(countData)
filter <- grepl("XB.+g.+",genes)
countData <- countData[filter,]
genes <- genes[filter]
genes <- sub(".*[^0-9](\\d+)$", "\\1", genes)
genes <- paste0("gene",genes)
rownames(countData) <- genes
countData <- as.matrix(countData)

controls <- strsplit(controls,",")[[1]]
treatments <- strsplit(treatments,",")[[1]]
countData <- countData[,c(controls,treatments)]
x <- as.factor(c(rep("Control",length(controls)),rep("Treatment",length(treatments))))
phenoData <- data.frame(x, row.names=colnames(countData))
set <- newSeqExpressionSet(as.matrix(countData), phenoData = phenoData)
set <- betweenLaneNormalization(set, which="upper")

filename <- paste0(sort(controls)[1],"_vs_",sort(treatments)[1],".txt")

if (length(controls) == 1 && length(treatments) == 1) {
  design <- model.matrix(~x, data=pData(set))
  y <- DGEList(counts=counts(set), group=x)
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, method="deviance", robust=TRUE, subset=NULL)
  y <- estimateGLMTagwiseDisp(y)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef=2)
  results <- topTags(lrt, n=dim(y)[1]+1, adjust.method="BH", sort.by="logFC")
  results <- results$table
  Gene <-  rownames(results)
  results <- cbind(Gene,results)
  write.table(results,file = filename,sep = "\t",quote=FALSE,row.names=F)
} else {
  design <- model.matrix(~x, data=pData(set))
  y <- DGEList(counts=counts(set), group=x)
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef=2)
  results <- topTags(lrt, n=dim(y)[1]+1, adjust.method="BH", sort.by="logFC")
  results <- results$table
  Gene <-  rownames(results)
  results <- cbind(Gene,results)
  write.table(results,file = filename,sep = "\t",quote=FALSE,row.names=F)
}