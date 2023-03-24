#!/usr/bin/env Rscript

library(tximport)

args <- commandArgs(TRUE)

# Nextflow sends the isoforms as individual arguments, this selects all of them
isoforms <- head(args,-2)
# The final two arguments are the work directory and transcript/gene map
remainder <- tail(args,2)
work_dir <- remainder[1]
map <- read.table(remainder[2],sep="\t",header=T)
isoforms <- paste0(work_dir,"/",isoforms)
# Use tximport to automatically create isoform matrices from files
names(isoforms) <- sub(".isoforms.results","",basename(isoforms))
txi <- tximport(isoforms, type = "rsem", txIn = TRUE, txOut = TRUE)

# Extract counts, reorder, use map to find gene name, combine with transcript name
counts <- as.data.frame(txi$counts)
if (ncol(counts) > 1) {
    srr <- colnames(counts)
    srr_ordered <- sort(sapply(srr, function(x) substr(x,4,nchar(x))))
    counts <- counts[,names(srr_ordered)]
}
counts[counts == ""] <- NA
transcripts <- rownames(counts)
genes <- map$GeneSym[match(transcripts,map$TranscriptId)]
genes_copy <- genes
genes_copy[is.na(genes_copy)] <- "NoGeneId"
"Gene<Transcript>" <-  paste0(genes_copy,"<",transcripts,">")
counts <- cbind(`Gene<Transcript>`,counts)
write.table(counts, "Isoforms_Counts_Matrix.txt", sep = '\t', row.names = F, col.names = T, quote = F)

# Remove the "NoGeneId" rows from the counts matrix, sum expression of all transcripts for each gene
if (ncol(counts) > 2) {
    gene_counts <- counts[!is.na(genes),c(2:ncol(counts))]
} else {
    gene_counts <- counts[!is.na(genes),]
}
genes <- genes[!is.na(genes)]
gene_counts <- rowsum(gene_counts, genes, F)
Gene <- rownames(gene_counts)
gene_counts <- cbind(Gene,gene_counts)
write.table(gene_counts, "Genes_Counts_Matrix.txt", sep = '\t', row.names = F, col.names = T, quote = F)

# Extract tpm, reorder, use map to find gene name, combine with transcript name
tpm <- as.data.frame(txi$abundance)
if (ncol(tpm) > 1) {
    srr <- colnames(tpm)
    srr_ordered <- sort(sapply(srr, function(x) substr(x,4,nchar(x))))
    tpm <- tpm[,names(srr_ordered)]
}
tpm[tpm == ""] <- NA
transcripts <- rownames(tpm)
genes <- map$GeneSym[match(transcripts,map$TranscriptId)]
genes_copy <- genes
genes_copy[is.na(genes_copy)] <- "NoGeneId"
"Gene<Transcript>" <-  paste0(genes_copy,"<",transcripts,">")
tpm <- cbind(`Gene<Transcript>`,tpm)
write.table(tpm, "Isoforms_TPM_Matrix.txt", sep = '\t', row.names = F, col.names = T, quote = F)

# Remove the "NoGeneId" rows from the tpm matrix, sum expression of all transcripts for each gene
if (ncol(tpm) > 2) {
    gene_tpm <- tpm[!is.na(genes),c(2:ncol(tpm))]
} else {
    gene_tpm <- tpm[!is.na(genes),]
}
genes <- genes[!is.na(genes)]
gene_tpm <- rowsum(gene_tpm, genes, F)
Gene <- rownames(gene_tpm)
gene_tpm <- cbind(Gene,gene_tpm)
write.table(gene_tpm, "Genes_TPM_Matrix.txt", sep = '\t', row.names = F, col.names = T, quote = F)