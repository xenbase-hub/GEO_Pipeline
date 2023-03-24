#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.bowtie2_path = "/usr/local/bowtie2/2.1.0/bin"
params.rsem_path = "/usr/local/rsem/1.3.0/bin/rsem-calculate-expression"
params.rsem_laevis = "$baseDir/data/RSEM_LAEVIS_v10/Laevis"
params.rsem_trop = "$baseDir/data/RSEM_TROP_v10/Trop"
params.map_laevis = "$baseDir/data/RSEM_LAEVIS_v10/Laevis-v10_Map.txt"
params.map_trop = "$baseDir/data/RSEM_TROP_v10/Trop-v10_Map.txt"
params.gsize_laevis = 1700000000
params.gsize_trop = 3100000000
params.gsize_laevis_2 = "1.7e9"
params.gsize_trop_2 = "3.1e9"

process alignRNA {

    // Aligns FASTQ files using RSEM and Bowtie2, followed by sorting and indexing bam
    // Outputs are bams and isoform results

    tag "$SAMPLE"

    module "deeptools:samtools/1.8.0"

    input:
    tuple val(SAMPLE), val(SPECIES), val(ASSAY), val(LIBRARY), path(FASTQ)

    output:
    tuple val(SAMPLE), val(SPECIES), val(ASSAY), val(LIBRARY), path("*.bam"), path("*.bai"), path("*.isoforms.results"), emit: alignment

    script:
    if(LIBRARY == "PAIRED-END" && SPECIES == "xenopus-trop-v10")
        """
        ${params.rsem_path} --output-genome-bam --paired-end --num-threads 4 --phred33-quals --bowtie2 --bowtie2-path ${params.bowtie2_path} ${FASTQ[0]} ${FASTQ[1]} ${params.rsem_trop} ${SAMPLE}
        rm ${SAMPLE}.genes.results
        rm ${SAMPLE}.transcript.bam
        nreads=`samtools view -@ 4 -c -F 4 ${SAMPLE}.genome.bam`
        frac=`gawk -v total="\$nreads" 'BEGIN {frac=500000000/total; if (frac > 1) {print 0.999} else {print frac}}'`
        samtools view -@ 4 -bs \$frac ${SAMPLE}.genome.bam -o ${SAMPLE}.genome.downsample.bam
        rm ${SAMPLE}.genome.bam
        samtools sort -@ 4 -m 12G ${SAMPLE}.genome.downsample.bam -o ${SAMPLE}.genome.sorted.bam
        rm ${SAMPLE}.genome.downsample.bam
        samtools index -@ 4 ${SAMPLE}.genome.sorted.bam
        """
    else if(LIBRARY == "PAIRED-END" && SPECIES == "xenopus-laevis-v10")
        """
        ${params.rsem_path} --output-genome-bam --paired-end --num-threads 4 --phred33-quals --bowtie2 --bowtie2-path ${params.bowtie2_path} ${FASTQ[0]} ${FASTQ[1]} ${params.rsem_laevis} ${SAMPLE}
        rm ${SAMPLE}.genes.results
        rm ${SAMPLE}.transcript.bam
        nreads=`samtools view -@ 4 -c -F 4 ${SAMPLE}.genome.bam`
        frac=`gawk -v total="\$nreads" 'BEGIN {frac=500000000/total; if (frac > 1) {print 0.999} else {print frac}}'`
        samtools view -@ 4 -bs \$frac ${SAMPLE}.genome.bam -o ${SAMPLE}.genome.downsample.bam
        rm ${SAMPLE}.genome.bam
        samtools sort -@ 4 -m 12G ${SAMPLE}.genome.downsample.bam -o ${SAMPLE}.genome.sorted.bam
        rm ${SAMPLE}.genome.downsample.bam
        samtools index -@ 4 ${SAMPLE}.genome.sorted.bam
        """
    else if(LIBRARY == "SINGLE-END" && SPECIES == "xenopus-trop-v10")
        """
        ${params.rsem_path} --output-genome-bam --num-threads 4 --phred33-quals --bowtie2 --bowtie2-path ${params.bowtie2_path} ${FASTQ[0]} ${params.rsem_trop} ${SAMPLE}
        rm ${SAMPLE}.genes.results
        rm ${SAMPLE}.transcript.bam
        nreads=`samtools view -@ 4 -c -F 4 ${SAMPLE}.genome.bam`
        frac=`gawk -v total="\$nreads" 'BEGIN {frac=500000000/total; if (frac > 1) {print 0.999} else {print frac}}'`
        samtools view -@ 4 -bs \$frac ${SAMPLE}.genome.bam -o ${SAMPLE}.genome.downsample.bam
        rm ${SAMPLE}.genome.bam
        samtools sort -@ 4 -m 12G ${SAMPLE}.genome.downsample.bam -o ${SAMPLE}.genome.sorted.bam
        rm ${SAMPLE}.genome.downsample.bam
        samtools index -@ 4 ${SAMPLE}.genome.sorted.bam
        """
    else if(LIBRARY == "SINGLE-END" && SPECIES == "xenopus-laevis-v10")
        """
        ${params.rsem_path} --output-genome-bam --num-threads 4 --phred33-quals --bowtie2 --bowtie2-path ${params.bowtie2_path} ${FASTQ[0]} ${params.rsem_laevis} ${SAMPLE}
        rm ${SAMPLE}.genes.results
        rm ${SAMPLE}.transcript.bam
        nreads=`samtools view -@ 4 -c -F 4 ${SAMPLE}.genome.bam`
        frac=`gawk -v total="\$nreads" 'BEGIN {frac=500000000/total; if (frac > 1) {print 0.999} else {print frac}}'`
        samtools view -@ 4 -bs \$frac ${SAMPLE}.genome.bam -o ${SAMPLE}.genome.downsample.bam
        rm ${SAMPLE}.genome.bam
        samtools sort -@ 4 -m 12G ${SAMPLE}.genome.downsample.bam -o ${SAMPLE}.genome.sorted.bam
        rm ${SAMPLE}.genome.downsample.bam
        samtools index -@ 4 ${SAMPLE}.genome.sorted.bam
        """
}

process processRNA {

    // Completes processing of bams from alignRNA
    // Biological replicates are merged and reprocessed
    // Outputs are bams

    publishDir "${params.outDir}", mode: 'copy', pattern: '*.bw', saveAs: {filename ->
        if(SPECIES[0] == "xenopus-trop-v10") "${GSE[0]}/XENTR_10.0/RNA-Seq/BigWigs/$filename"
        else if(SPECIES[0] == "xenopus-laevis-v10") "${GSE[0]}/XENLA_10.1/RNA-Seq/BigWigs/$filename"
    }
    
    module "deeptools:samtools/1.8.0"

    input:
    tuple val(BIGWIG), val(MERGE_BIO), path(ALIGNMENT), path(INDEX), val(GSE), val(SPECIES)

    output:
    tuple val(BIGWIG), val(GSE), emit: complete
    path("*.bw")

    script:

    if(SPECIES[0] == "xenopus-trop-v10" && MERGE_BIO[0] > 1)
        """
        samtools merge -@ 4 combined.bam ${ALIGNMENT.join(" ")}
        nreads=`samtools view -@ 4 -c -F 4 combined.bam`
        frac=`gawk -v total="\$nreads" 'BEGIN {frac=500000000/total; if (frac > 1) {print 0.999} else {print frac}}'`
        samtools view -@ 4 -bs \$frac combined.bam -o combined.downsample.bam
        rm combined.bam
        samtools sort -@ 4 -m 25G combined.downsample.bam -o combined.sorted.bam
        rm combined.downsample.bam
        samtools index -@ 4 combined.sorted.bam
        bamCoverage --bam combined.sorted.bam -o ${BIGWIG} --effectiveGenomeSize ${params.gsize_trop} --normalizeUsing BPM --binSize 20 --smoothLength 60 --numberOfProcessors 4 --verbose
        rm combined.sorted.bam
        rm combined.sorted.bam.bai
        """
    else if(SPECIES[0] == "xenopus-trop-v10" && MERGE_BIO[0] == 1)
        """
        bamCoverage --bam ${ALIGNMENT[0]} -o ${BIGWIG} --effectiveGenomeSize ${params.gsize_trop} --normalizeUsing BPM --binSize 20 --smoothLength 60 --numberOfProcessors 4 --verbose
        """
    else if(SPECIES[0] == "xenopus-laevis-v10" && MERGE_BIO[0] > 1)
        """
        samtools merge -@ 4 combined.bam ${ALIGNMENT.join(" ")}
        nreads=`samtools view -@ 4 -c -F 4 combined.bam`
        frac=`gawk -v total="\$nreads" 'BEGIN {frac=500000000/total; if (frac > 1) {print 0.999} else {print frac}}'`
        samtools view -@ 4 -bs \$frac combined.bam -o combined.downsample.bam
        rm combined.bam
        samtools sort -@ 4 -m 25G combined.downsample.bam -o combined.sorted.bam
        rm combined.downsample.bam
        samtools index -@ 4 combined.sorted.bam
        bamCoverage --bam combined.sorted.bam -o ${BIGWIG} --effectiveGenomeSize ${params.gsize_laevis} --normalizeUsing BPM --binSize 20 --smoothLength 60 --numberOfProcessors 4 --verbose
        rm combined.sorted.bam
        rm combined.sorted.bam.bai
        """
    else if(SPECIES[0] == "xenopus-laevis-v10" && MERGE_BIO[0] == 1)
        """
        bamCoverage --bam ${ALIGNMENT[0]} -o ${BIGWIG} --effectiveGenomeSize ${params.gsize_laevis} --normalizeUsing BPM --binSize 20 --smoothLength 60 --numberOfProcessors 4 --verbose
        """
}

process createMatrices {

    // Create matrices from isoforms.results files from RSEM
    // Output are the counts/TPM matrices for genes/isoforms

    publishDir "${params.outDir}", mode: 'copy', pattern: '*.txt', saveAs: {filename ->
        if(SPECIES[0] == "xenopus-trop-v10") "${GSE[0]}/XENTR_10.0/RNA-Seq/ExpressionFiles/$filename"
        else if(SPECIES[0] == "xenopus-laevis-v10") "${GSE[0]}/XENLA_10.1/RNA-Seq/ExpressionFiles/$filename"
    }

    module 'R/4.1.1'

    input:
    tuple val(GSE_SPECIES), val(GSE), val(SPECIES), path(ISOFORMS)

    output:
    tuple val(GSE_SPECIES), val(GSE), val(SPECIES), path("Genes_Counts_Matrix.txt"), emit: counts
    path("*.txt"), emit: matrices

    script:
    if(SPECIES[0] == "xenopus-trop-v10")
        """
        Rscript ${workflow.projectDir}/bin/generate_tpm_counts.R ${ISOFORMS} \$PWD ${params.map_trop}
        """
    else if(SPECIES[0] == "xenopus-laevis-v10")
        """
        Rscript ${workflow.projectDir}/bin/generate_tpm_counts.R ${ISOFORMS} \$PWD ${params.map_laevis}
        """
}

process performDE {

    // Perform DE analysis using RUVSeq and UQ normalization

    publishDir "${params.outDir}", mode: 'copy', pattern: '*.txt', saveAs: {filename ->
        if(SPECIES[0] == "xenopus-trop-v10") "${GSE[0]}/XENTR_10.0/RNA-Seq/DE_Analysis/$filename"
        else if(SPECIES[0] == "xenopus-laevis-v10") "${GSE[0]}/XENLA_10.1/RNA-Seq/DE_Analysis/$filename"
    }

    module 'R/4.1.1'

    input:
    tuple val(GSE_SPECIES), val(GSE), val(SPECIES), path(MATRIX), val(TREATMENT), val(CONTROL)

    output:
    path("*.txt")

    script:
    """
    Rscript ${workflow.projectDir}/bin/perform_de.R ${MATRIX} ${TREATMENT} ${CONTROL}
    """

}