#!/usr/bin/env Rscript

library(tidyr)
library(tools)

args <- commandArgs(TRUE)

df <- read.table(args[1],sep="\t",header=T)
# Replace all pipes with underscores
df <- data.frame(lapply(df, function(x) {gsub("\\|", "_", x)}))
# Filter out unsupported assay types
df <- df[df$ASSAYTYPE %in% c("RNA-SEQ","ATAC","ChIP-TF","ChIP-Epigenetic"),]
# Correct the library if necessary
ends <- df$LIBRARYLAYOUT
ends[ends == "SINGLE"] <- "SINGLE-END"
ends[ends == "PAIRED"] <- "PAIRED-END"
df$LIBRARYLAYOUT <- ends
# Create a copy of SRR_TREATMENT to split and use for downloading SRRs
df$SRR_COPY <- df$SRR_TREATMENT
# Modify controls to exclude everything but RNA-SEQ, also explicity set those without controls to NA
controls <- df$SRR_CONTROL
controls[controls == ""] <- "NA"
controls[df$ASSAYTYPE != "RNA-SEQ"] <- "NA"
df$SRR_CONTROL <- controls
# There are rare instances where a GSE contains multiple species. GSE_SPECIES is used to track this
df$GSE_SPECIES <- paste0(df$GSE,"_",df$SPECIES)
# MERGE_BIO is used to combine biological replicates, AKA the samples listed in SRR_TREATMENT, for bigwigs/peaks
df$MERGE_BIO <- sapply(strsplit(df$SRR_COPY,","), function(x) length(x))
# Rows are split so each corresponds to a single sample
df_split <- df %>% separate_rows(SRR_COPY,sep=",")
df_split$SAMPLE <- df_split$SRR_COPY
# Shortening the names to the first technical replicate
df_split$SAMPLE <- sapply(strsplit(df_split$SAMPLE,"_"), function(x) x[1])
srr <- strsplit(df_split$SRR_TREATMENT,",")
srr_final <- c()
for (x in srr) {
  step1 <- strsplit(x,"_")
  step2 <- sapply(step1,function(x) x[1])
  step3 <- paste(step2,collapse=",")
  srr_final <- c(srr_final,step3)
}
df_split$SRR_TREATMENT <- srr_final
srr <- strsplit(df_split$SRR_CONTROL,",")
srr_final <- c()
for (x in srr) {
  step1 <- strsplit(x,"_")
  step2 <- sapply(step1,function(x) x[1])
  step3 <- paste(step2,collapse=",")
  srr_final <- c(srr_final,step3)
}
df_split$SRR_CONTROL <- srr_final
# MERGE_GSE is used to combine different rows when creating the TPM/counts matrices. Only applies to RNA-seq
gses <- df_split$GSE_SPECIES
gses[df_split$ASSAYTYPE != "RNA-SEQ"] = "GSE"
merge_gse <- sapply(gses, function(x) length(gses[gses == x]))
merge_gse[df_split$ASSAYTYPE != "RNA-SEQ"] = 0
df_split$MERGE_GSE <- merge_gse
# Rows are split further if they contain technical replicates
df_split <- df_split %>% separate_rows(SRR_COPY,sep="_")
# MERGE_TECH is used to merge FASTQs of the technical replicates
samples <- df_split$SAMPLE
df_split$MERGE_TECH <- sapply(samples, function(x) length(samples[samples == x]))
# Table is reordered and renamed, then saved
df_split$FILENAME <- sapply(df_split$BIGWIG, function(x) file_path_sans_ext(x))
df_split <- df_split[,c("SAMPLE","MERGE_TECH","MERGE_BIO","MERGE_GSE","SRR_COPY","GSE","GSE_SPECIES","SPECIES","ASSAYTYPE","LIBRARYLAYOUT","BIGWIG","FILENAME","SRR_TREATMENT","SRR_CONTROL")]
colnames(df_split) <- c("SAMPLE","MERGE_TECH","MERGE_BIO","MERGE_GSE","SRR","GSE","GSE_SPECIES","SPECIES","ASSAY","LIBRARY","BIGWIG","FILENAME","TREATMENT","CONTROL")
write.table(df_split,file="srr_table.txt",sep="\t",quote=F, row.names=F)