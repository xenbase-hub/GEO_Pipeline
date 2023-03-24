#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {processRunFile; createReadMe} from './modules/xenbase'
include {downloadSRR; mergeFASTQ} from './modules/fastq'
include {alignRNA; processRNA; createMatrices; performDE} from './modules/rna'
include {alignDNA; processDNASingle; processDNAMerged; callPeaks} from './modules/dna'
include {cleanWorkFiles as cleanRNA1; cleanWorkFiles as cleanDNA1; cleanWorkFiles as cleanRNA2; cleanWorkFiles as cleanDNA2} from './modules/cleanup'

params.runFile = null
params.outDir = null

results = file("${params.outDir}").mkdir()

workflow {
	// Reads in run file, reformats it into a table, specific columns are selected for each channel
	// "sample_data" has a group key to handle technical replicates, rest are needed for alignment
	// "metadata" is used for combining biological replicates, performing DE, and saving files to the correct folders

	Channel.fromPath(params.runFile, checkIfExists:true) | set{runfile}
	runfile | processRunFile | splitCsv(header: true, sep: '\t') \
	| multiMap {row -> 
		sample_data: tuple(groupKey(row.SAMPLE, row.MERGE_TECH as int), row.MERGE_TECH as int, row.SRR, row.SPECIES, row.ASSAY, row.LIBRARY)
		metadata: tuple(row.SAMPLE, row.GSE, row.GSE_SPECIES, row.MERGE_GSE as int, row.BIGWIG, row.MERGE_BIO as int, row.FILENAME)
		deg: tuple(row.GSE_SPECIES, row.TREATMENT, row.CONTROL)
		readme: tuple(row.GSE, row.SPECIES, row.ASSAY)
	} | set{runs}

	runs.deg | unique | filter{GSE_SPECIES, TREATMENT, CONTROL -> CONTROL != "NA"} | set{deg}

	runfile | combine(runs.readme.unique()) | createReadMe

	// FASTQs for each sample are downloaded, and those that are technical replicates are grouped together

	downloadSRR(runs.sample_data) 
	downloadSRR.out.raw_data | groupTuple \
	| branch {SAMPLE, MERGE_TECH, SRR, SPECIES, ASSAY, LIBRARY, FASTQ ->
		to_merge: MERGE_TECH[0] > 1
			return tuple(SAMPLE,SPECIES[0],ASSAY[0],LIBRARY[0],FASTQ.flatten())
		not_merged: true
			return tuple(SAMPLE,SPECIES[0],ASSAY[0],LIBRARY[0],FASTQ.flatten())
	} | set{fastqs}

	// Technical replicates are merged and recombined with the other FASTQs. FASTQs are then divided by assay

	mergeFASTQ(fastqs.to_merge)
	fastqs.not_merged | mix(mergeFASTQ.out.merged) \
	| branch { SAMPLE, SPECIES, ASSAY, LIBRARY, FASTQ ->
		RNA: ASSAY == "RNA-SEQ"
		DNA: true
	} | set{assays}

	// Run alignment for RNA, combine with metadata, split into two channels
	// processRNA will merge biological replicates if necessary and generate bigwig files
	// createMatrices will combine samples with the same GSE to create matrices from isoform results
	// performDE takes the matrices to calculate DEGs

	alignRNA(assays.RNA) 
	downloadSRR.out.raw_data | groupTuple | mix(mergeFASTQ.out.merged) | join(alignRNA.out.alignment) | flatten | filter{it =~ /.fastq$/} | cleanRNA1
	alignRNA.out.alignment | join(runs.metadata) | multiMap{SAMPLE, SPECIES, ASSAY, LIBRARY, ALIGNMENT, INDEX, ISOFORMS, GSE, GSE_SPECIES, MERGE_GSE, BIGWIG, MERGE_BIO, FILENAME -> 
		bigwigs: tuple(groupKey(BIGWIG, MERGE_BIO), MERGE_BIO, ALIGNMENT, INDEX, GSE, SPECIES)
		matrices: tuple(groupKey(GSE_SPECIES, MERGE_GSE), GSE, SPECIES, ISOFORMS)
		removal: tuple(BIGWIG, ALIGNMENT)
	} | set{dataRNA}
	dataRNA.bigwigs | groupTuple | processRNA
	dataRNA.removal | join(processRNA.out.complete) | flatten | filter{it =~ /.bam$/} | cleanRNA2
	dataRNA.matrices | groupTuple | createMatrices
	createMatrices.out.counts | combine(deg, by: 0) | performDE

	// Run alignment for DNA, combine with metadata, group and check if there are biological replicates
	// Two processDNAs are needed to generate bigwig files (unlike RNA, bam files are needed downstream)
	// callPeaks takes the bams to generate peak files

	alignDNA(assays.DNA) 
	downloadSRR.out.raw_data | groupTuple | mix(mergeFASTQ.out.merged)| join(alignDNA.out.alignment) | flatten | filter{it =~ /.fastq$/} | cleanDNA1
	alignDNA.out.alignment | join(runs.metadata) | multiMap{SAMPLE, SPECIES, ASSAY, LIBRARY, ALIGNMENT, INDEX, GSE, GSE_SPECIES, MERGE_GSE, BIGWIG, MERGE_BIO, FILENAME -> 
		process: tuple(groupKey(BIGWIG, MERGE_BIO), MERGE_BIO, ALIGNMENT, INDEX, GSE, SPECIES, ASSAY, LIBRARY, FILENAME)
		removal: tuple(BIGWIG, ALIGNMENT)
	} | set{dataDNA}
	dataDNA.process | groupTuple | branch { BIGWIG, MERGE_BIO, ALIGNMENT, INDEX, GSE, SPECIES, ASSAY, LIBRARY, FILENAME ->
		to_merge: MERGE_BIO[0] > 1
		not_merged: true
	} | set{mergeDNA}
	processDNASingle(mergeDNA.not_merged)
	processDNAMerged(mergeDNA.to_merge)
	processDNASingle.out.peaks | mix(processDNAMerged.out.peaks) | callPeaks
	dataDNA.removal | mix(processDNAMerged.out.peaks) | join(callPeaks.out.complete) | flatten | filter{it =~ /.bam$/} | cleanDNA2
}
