#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.bowtie_laevis = "$baseDir/data/Bowtie_indexes_LAEVIS_v10/Laevis"
params.bowtie_trop = "$baseDir/data/Bowtie_indexes_TROP_v10/Trop"
params.gsize_laevis = 1700000000
params.gsize_trop = 3100000000
params.gsize_laevis_2 = "1.7e9"
params.gsize_trop_2 = "3.1e9"

process alignDNA {

    // Aligns FASTQ files using Bowtie2, followed by marking duplicates, sorting and indexing bam
    // Paired end also uses fixmate
    // Outputs are bams

    tag "$SAMPLE"
    
    module "bowtie2/2.1.0:samtools/1.8.0"

    input:
    tuple val(SAMPLE), val(SPECIES), val(ASSAY), val(LIBRARY), path(FASTQ)

    output:
    tuple val(SAMPLE), val(SPECIES), val(ASSAY), val(LIBRARY), path("*.bam"), path("*.bai"), emit: alignment

    script:
    if(LIBRARY == "PAIRED-END" && SPECIES == "xenopus-trop-v10")
        """
        bowtie2 -x ${params.bowtie_trop} -1 ${FASTQ[0]} -2 ${FASTQ[1]} -S ${SAMPLE}.sam --threads 4
        samtools view -@ 4 -b -o ${SAMPLE}.bam ${SAMPLE}.sam
        rm ${SAMPLE}.sam
        samtools fixmate -@ 4 -m ${SAMPLE}.bam ${SAMPLE}.temp.bam
        rm ${SAMPLE}.bam
        nreads=`samtools view -@ 4 -c -F 4 ${SAMPLE}.temp.bam`
        frac=`gawk -v total="\$nreads" 'BEGIN {frac=500000000/total; if (frac > 1) {print 0.999} else {print frac}}'`
        samtools view -@ 4 -bs \$frac ${SAMPLE}.temp.bam -o ${SAMPLE}.downsample.bam
        rm ${SAMPLE}.temp.bam
        samtools sort -@ 4 -m 12G ${SAMPLE}.downsample.bam -o ${SAMPLE}.sorted.bam
        rm ${SAMPLE}.downsample.bam
        samtools markdup -@ 4 -r ${SAMPLE}.sorted.bam ${SAMPLE}.final.bam
        rm ${SAMPLE}.sorted.bam
        samtools sort -@ 4 -m 12G ${SAMPLE}.final.bam -o ${SAMPLE}.final.sorted.bam
        rm ${SAMPLE}.final.bam
        samtools index -@ 4 ${SAMPLE}.final.sorted.bam
        """
    else if(LIBRARY == "PAIRED-END" && SPECIES == "xenopus-laevis-v10")
        """
        bowtie2 -x ${params.bowtie_laevis} -1 ${FASTQ[0]} -2 ${FASTQ[1]} -S ${SAMPLE}.sam --threads 4
        samtools view -@ 4 -b -o ${SAMPLE}.bam ${SAMPLE}.sam
        rm ${SAMPLE}.sam
        samtools fixmate -@ 4 -m ${SAMPLE}.bam ${SAMPLE}.temp.bam
        rm ${SAMPLE}.bam
        nreads=`samtools view -@ 4 -c -F 4 ${SAMPLE}.temp.bam`
        frac=`gawk -v total="\$nreads" 'BEGIN {frac=500000000/total; if (frac > 1) {print 0.999} else {print frac}}'`
        samtools view -@ 4 -bs \$frac ${SAMPLE}.temp.bam -o ${SAMPLE}.downsample.bam
        rm ${SAMPLE}.temp.bam
        samtools sort -@ 4 -m 12G ${SAMPLE}.downsample.bam -o ${SAMPLE}.sorted.bam
        rm ${SAMPLE}.downsample.bam
        samtools markdup -@ 4 -r ${SAMPLE}.sorted.bam ${SAMPLE}.final.bam
        rm ${SAMPLE}.sorted.bam
        samtools sort -@ 4 -m 12G ${SAMPLE}.final.bam -o ${SAMPLE}.final.sorted.bam
        rm ${SAMPLE}.final.bam
        samtools index -@ 4 ${SAMPLE}.final.sorted.bam
        """
    else if(LIBRARY == "SINGLE-END" && SPECIES == "xenopus-trop-v10")
        """
        bowtie2 -x ${params.bowtie_trop} -q -U ${FASTQ[0]} -S ${SAMPLE}.sam --threads 4
        samtools view -@ 4 -b -o ${SAMPLE}.bam ${SAMPLE}.sam
        rm ${SAMPLE}.sam
        nreads=`samtools view -@ 4 -c -F 4 ${SAMPLE}.bam`
        frac=`gawk -v total="\$nreads" 'BEGIN {frac=500000000/total; if (frac > 1) {print 0.999} else {print frac}}'`
        samtools view -@ 4 -bs \$frac ${SAMPLE}.bam -o ${SAMPLE}.downsample.bam
        rm ${SAMPLE}.bam
        samtools sort -@ 4 -m 12G ${SAMPLE}.downsample.bam -o ${SAMPLE}.sorted.bam
        rm ${SAMPLE}.downsample.bam
        samtools markdup -@ 4 -r ${SAMPLE}.sorted.bam ${SAMPLE}.final.bam
        rm ${SAMPLE}.sorted.bam
        samtools sort -@ 4 -m 12G ${SAMPLE}.final.bam -o ${SAMPLE}.final.sorted.bam
        rm ${SAMPLE}.final.bam
        samtools index -@ 4 ${SAMPLE}.final.sorted.bam
        """
    else if(LIBRARY == "SINGLE-END" && SPECIES == "xenopus-laevis-v10")
        """
        bowtie2 -x ${params.bowtie_laevis} -q -U ${FASTQ[0]} -S ${SAMPLE}.sam --threads 4
        samtools view -@ 4 -b -o ${SAMPLE}.bam ${SAMPLE}.sam
        rm ${SAMPLE}.sam
        nreads=`samtools view -@ 4 -c -F 4 ${SAMPLE}.bam`
        frac=`gawk -v total="\$nreads" 'BEGIN {frac=500000000/total; if (frac > 1) {print 0.999} else {print frac}}'`
        samtools view -@ 4 -bs \$frac ${SAMPLE}.bam -o ${SAMPLE}.downsample.bam
        rm ${SAMPLE}.bam
        samtools sort -@ 4 -m 12G ${SAMPLE}.downsample.bam -o ${SAMPLE}.sorted.bam
        rm ${SAMPLE}.downsample.bam
        samtools markdup -@ 4 -r ${SAMPLE}.sorted.bam ${SAMPLE}.final.bam
        rm ${SAMPLE}.sorted.bam
        samtools sort -@ 4 -m 12G ${SAMPLE}.final.bam -o ${SAMPLE}.final.sorted.bam
        rm ${SAMPLE}.final.bam
        samtools index -@ 4 ${SAMPLE}.final.sorted.bam
        """
}

process processDNASingle {

    // Create bigwigs from bam files
    // Unlike RNA, could not process DNA in a single process because bams are needed downstream

    publishDir "${params.outDir}", mode: 'copy', pattern: '*.bw', saveAs: {filename ->
        if(SPECIES[0] == "xenopus-trop-v10" && ASSAY[0] == "ChIP-Epigenetic") "${GSE[0]}/XENTR_10.0/ChIP-Seq/BigWigs/$filename"
        else if(SPECIES[0] == "xenopus-laevis-v10" && ASSAY[0] == "ChIP-Epigenetic") "${GSE[0]}/XENLA_10.1/ChIP-Seq/BigWigs/$filename"
        else if(SPECIES[0] == "xenopus-trop-v10" && ASSAY[0] == "ChIP-TF") "${GSE[0]}/XENTR_10.0/ChIP-Seq/BigWigs/$filename"
        else if(SPECIES[0] == "xenopus-laevis-v10" && ASSAY[0] == "ChIP-TF") "${GSE[0]}/XENLA_10.1/ChIP-Seq/BigWigs/$filename"
        else if(SPECIES[0] == "xenopus-trop-v10" && ASSAY[0] == "ATAC") "${GSE[0]}/XENTR_10.0/ATAC-Seq/BigWigs/$filename"
        else if(SPECIES[0] == "xenopus-laevis-v10" && ASSAY[0] == "ATAC") "${GSE[0]}/XENLA_10.1/ATAC-Seq/BigWigs/$filename"
    }
    
    module "deeptools:samtools/1.8.0"

    input:
    tuple val(BIGWIG), val(MERGE_BIO), path(ALIGNMENT), path(INDEX), val(GSE), val(SPECIES), val(ASSAY), val(LIBRARY), val(FILENAME)

    output:
    tuple val(BIGWIG), val(GSE), val(SPECIES), val(ASSAY), val(LIBRARY), val(FILENAME), path(ALIGNMENT), emit: peaks
    path("*.bw"), emit: bigwigs

    script:

    if(SPECIES[0] == "xenopus-trop-v10")
        """
        bamCoverage --bam ${ALIGNMENT[0]} -o ${BIGWIG} --binSize 20 --normalizeUsing RPGC --effectiveGenomeSize ${params.gsize_trop}  --extendReads 200 --smoothLength 60 --numberOfProcessors 4 --verbose
        """
    else if(SPECIES[0] == "xenopus-laevis-v10")
        """
        bamCoverage --bam ${ALIGNMENT[0]} -o ${BIGWIG} --binSize 20 --normalizeUsing RPGC --effectiveGenomeSize ${params.gsize_laevis}  --extendReads 200 --smoothLength 60 --numberOfProcessors 4 --verbose
        """
}

process processDNAMerged {

    // Merge biological replicates and create bigwig files

    publishDir "${params.outDir}", mode: 'copy', pattern: '*.bw', saveAs: {filename ->
        if(SPECIES[0] == "xenopus-trop-v10" && ASSAY[0] == "ChIP-Epigenetic") "${GSE[0]}/XENTR_10.0/ChIP-Seq/BigWigs/$filename"
        else if(SPECIES[0] == "xenopus-laevis-v10" && ASSAY[0] == "ChIP-Epigenetic") "${GSE[0]}/XENLA_10.1/ChIP-Seq/BigWigs/$filename"
        else if(SPECIES[0] == "xenopus-trop-v10" && ASSAY[0] == "ChIP-TF") "${GSE[0]}/XENTR_10.0/ChIP-Seq/BigWigs/$filename"
        else if(SPECIES[0] == "xenopus-laevis-v10" && ASSAY[0] == "ChIP-TF") "${GSE[0]}/XENLA_10.1/ChIP-Seq/BigWigs/$filename"
        else if(SPECIES[0] == "xenopus-trop-v10" && ASSAY[0] == "ATAC") "${GSE[0]}/XENTR_10.0/ATAC-Seq/BigWigs/$filename"
        else if(SPECIES[0] == "xenopus-laevis-v10" && ASSAY[0] == "ATAC") "${GSE[0]}/XENLA_10.1/ATAC-Seq/BigWigs/$filename"
    }
    
    module "deeptools:samtools/1.8.0"

    input:
    tuple val(BIGWIG), val(MERGE_BIO), path(ALIGNMENT), path(INDEX), val(GSE), val(SPECIES), val(ASSAY), val(LIBRARY), val(FILENAME)

    output:
    tuple val(BIGWIG), val(GSE), val(SPECIES), val(ASSAY), val(LIBRARY), val(FILENAME), path("*.bam"), emit: peaks
    path("*.bw"), emit: bigwigs

    script:

    if(SPECIES[0] == "xenopus-trop-v10")
        """
        samtools merge -@ 4 combined.bam ${ALIGNMENT.join(" ")}
        nreads=`samtools view -@ 4 -c -F 4 combined.bam`
        frac=`gawk -v total="\$nreads" 'BEGIN {frac=500000000/total; if (frac > 1) {print 0.999} else {print frac}}'`
        samtools view -@ 4 -bs \$frac combined.bam -o combined.downsample.bam
        rm combined.bam
        samtools sort -@ 4 -m 25G combined.downsample.bam -o combined.sorted.bam
        rm combined.downsample.bam
        samtools index -@ 4 combined.sorted.bam
        samtools markdup -@ 4 -r combined.sorted.bam combined.final.bam
        rm combined.sorted.bam
        samtools sort -@ 4 -m 25G combined.final.bam -o combined.final.sorted.bam
        rm combined.final.bam
        samtools index -@ 4 combined.final.sorted.bam
        bamCoverage --bam combined.final.sorted.bam -o ${BIGWIG} --binSize 20 --normalizeUsing RPGC --effectiveGenomeSize ${params.gsize_trop}  --extendReads 200 --smoothLength 60 --numberOfProcessors 4 --verbose
        """
    else if(SPECIES[0] == "xenopus-laevis-v10")
        """
        samtools merge -@ 4 combined.bam ${ALIGNMENT.join(" ")}
        nreads=`samtools view -@ 4 -c -F 4 combined.bam`
        frac=`gawk -v total="\$nreads" 'BEGIN {frac=500000000/total; if (frac > 1) {print 0.999} else {print frac}}'`
        samtools view -@ 4 -bs \$frac combined.bam -o combined.downsample.bam
        rm combined.bam
        samtools sort -@ 4 -m 25G combined.downsample.bam -o combined.sorted.bam
        rm combined.downsample.bam
        samtools index -@ 4 combined.sorted.bam
        samtools markdup -@ 4 -r combined.sorted.bam combined.final.bam
        rm combined.sorted.bam
        samtools sort -@ 4 -m 25G combined.final.bam -o combined.final.sorted.bam
        rm combined.final.bam
        samtools index -@ 4 combined.final.sorted.bam
        bamCoverage --bam combined.final.sorted.bam -o ${BIGWIG} --binSize 20 --normalizeUsing RPGC --effectiveGenomeSize ${params.gsize_laevis}  --extendReads 200 --smoothLength 60 --numberOfProcessors 4 --verbose
        """
}

process callPeaks {

    // Produce bed peak files using MACS2 

    publishDir "${params.outDir}", mode: 'copy', pattern: '*Peaks.bed', saveAs: {filename ->
        if(SPECIES[0] == "xenopus-trop-v10" && ASSAY[0] == "ChIP-Epigenetic") "${GSE[0]}/XENTR_10.0/ChIP-Seq/Called_Peaks/$filename"
        else if(SPECIES[0] == "xenopus-laevis-v10" && ASSAY[0] == "ChIP-Epigenetic") "${GSE[0]}/XENLA_10.1/ChIP-Seq/Called_Peaks/$filename"
        else if(SPECIES[0] == "xenopus-trop-v10" && ASSAY[0] == "ChIP-TF") "${GSE[0]}/XENTR_10.0/ChIP-Seq/Called_Peaks/$filename"
        else if(SPECIES[0] == "xenopus-laevis-v10" && ASSAY[0] == "ChIP-TF") "${GSE[0]}/XENLA_10.1/ChIP-Seq/Called_Peaks/$filename"
        else if(SPECIES[0] == "xenopus-trop-v10" && ASSAY[0] == "ATAC") "${GSE[0]}/XENTR_10.0/ATAC-Seq/Called_Peaks/$filename"
        else if(SPECIES[0] == "xenopus-laevis-v10" && ASSAY[0] == "ATAC") "${GSE[0]}/XENLA_10.1/ATAC-Seq/Called_Peaks/$filename"
    }

    module 'MACS/2.1.0'

    input:
    tuple val(BIGWIG), val(GSE), val(SPECIES), val(ASSAY), val(LIBRARY), val(FILENAME), path(ALIGNMENT)

    output:
    tuple val(BIGWIG), val(GSE), emit: complete
    path("*.bed")

    script:
    if(ASSAY[0] == "ChIP-Epigenetic" && LIBRARY[0] == "SINGLE-END" && SPECIES[0] == "xenopus-trop-v10")
        """
        macs2 callpeak -t ${ALIGNMENT[0]} --buffer-size 3000 -f BAM -g ${params.gsize_trop_2} -n ${FILENAME[0]} --broad --broad-cutoff 0.1
        mv ${FILENAME[0]}_peaks.broadPeak ${FILENAME[0]}_broadPeaks.bed
        """
    else if(ASSAY[0] == "ChIP-Epigenetic" && LIBRARY[0] == "SINGLE-END" && SPECIES[0] == "xenopus-laevis-v10")
        """
        macs2 callpeak -t ${ALIGNMENT[0]} --buffer-size 3000 -f BAM -g ${params.gsize_laevis_2} -n ${FILENAME[0]} --broad --broad-cutoff 0.1
        mv ${FILENAME[0]}_peaks.broadPeak ${FILENAME[0]}_broadPeaks.bed
        """
    else if(ASSAY[0] == "ChIP-Epigenetic" && LIBRARY[0] == "PAIRED-END" && SPECIES[0] == "xenopus-trop-v10")
        """
        macs2 callpeak -t ${ALIGNMENT[0]} --buffer-size 3000 -f BAMPE -g ${params.gsize_trop_2} -n ${FILENAME[0]} --broad --broad-cutoff 0.1
        mv ${FILENAME[0]}_peaks.broadPeak ${FILENAME[0]}_broadPeaks.bed
        """
    else if(ASSAY[0] == "ChIP-Epigenetic" && LIBRARY[0] == "PAIRED-END" && SPECIES[0] == "xenopus-laevis-v10")
        """
        macs2 callpeak -t ${ALIGNMENT[0]} --buffer-size 3000 -f BAMPE -g ${params.gsize_laevis_2} -n ${FILENAME[0]} --broad --broad-cutoff 0.1
        mv ${FILENAME[0]}_peaks.broadPeak ${FILENAME[0]}_broadPeaks.bed
        """
    else if(ASSAY[0] == "ChIP-TF" && LIBRARY[0] == "SINGLE-END" && SPECIES[0] == "xenopus-trop-v10")
        """
        macs2 callpeak -t ${ALIGNMENT[0]} --buffer-size 3000 -f BAM -g ${params.gsize_trop_2} -n ${FILENAME[0]} --call-summits
        mv ${FILENAME[0]}_peaks.narrowPeak ${FILENAME[0]}_narrowPeaks.bed
        """
    else if(ASSAY[0] == "ChIP-TF" && LIBRARY[0] == "SINGLE-END" && SPECIES[0] == "xenopus-laevis-v10")
        """
        macs2 callpeak -t ${ALIGNMENT[0]} --buffer-size 3000 -f BAM -g ${params.gsize_laevis_2} -n ${FILENAME[0]} --call-summits
        mv ${FILENAME[0]}_peaks.narrowPeak ${FILENAME[0]}_narrowPeaks.bed
        """
    else if(ASSAY[0] == "ChIP-TF" && LIBRARY[0] == "PAIRED-END" && SPECIES[0] == "xenopus-trop-v10")
        """
        macs2 callpeak -t ${ALIGNMENT[0]} --buffer-size 3000 -f BAMPE -g ${params.gsize_trop_2} -n ${FILENAME[0]} --call-summits
        mv ${FILENAME[0]}_peaks.narrowPeak ${FILENAME[0]}_narrowPeaks.bed
        """
    else if(ASSAY[0] == "ChIP-TF" && LIBRARY[0] == "PAIRED-END" && SPECIES[0] == "xenopus-laevis-v10")
        """
        macs2 callpeak -t ${ALIGNMENT[0]} --buffer-size 3000 -f BAMPE -g ${params.gsize_laevis_2} -n ${FILENAME[0]} --call-summits
        mv ${FILENAME[0]}_peaks.narrowPeak ${FILENAME[0]}_narrowPeaks.bed
        """
    else if(ASSAY[0] == "ATAC" && LIBRARY[0] == "SINGLE-END" && SPECIES[0] == "xenopus-trop-v10")
        """
        macs2 callpeak -t ${ALIGNMENT[0]} --buffer-size 3000 -f BAM -g ${params.gsize_trop_2} -n ${FILENAME[0]} --nomodel --shift 37 --extsize 73 --call-summits
        mv ${FILENAME[0]}_peaks.narrowPeak ${FILENAME[0]}_narrowPeaks.bed
        """
    else if(ASSAY[0] == "ATAC" && LIBRARY[0] == "SINGLE-END" && SPECIES[0] == "xenopus-laevis-v10")
        """
        macs2 callpeak -t ${ALIGNMENT[0]} --buffer-size 3000 -f BAM -g ${params.gsize_laevis_2} -n ${FILENAME[0]} --nomodel --shift 37 --extsize 73 --call-summits
        mv ${FILENAME[0]}_peaks.narrowPeak ${FILENAME[0]}_narrowPeaks.bed
        """
    else if(ASSAY[0] == "ATAC" && LIBRARY[0] == "PAIRED-END" && SPECIES[0] == "xenopus-trop-v10")
        """
        macs2 callpeak -t ${ALIGNMENT[0]} --buffer-size 3000 -f BAMPE -g ${params.gsize_trop_2} -n ${FILENAME[0]} --call-summits
        mv ${FILENAME[0]}_peaks.narrowPeak ${FILENAME[0]}_narrowPeaks.bed
        """
    else if(ASSAY[0] == "ATAC" && LIBRARY[0] == "PAIRED-END" && SPECIES[0] == "xenopus-laevis-v10")
        """
        macs2 callpeak -t ${ALIGNMENT[0]} --buffer-size 3000 -f BAMPE -g ${params.gsize_laevis_2} -n ${FILENAME[0]} --call-summits
        mv ${FILENAME[0]}_peaks.narrowPeak ${FILENAME[0]}_narrowPeaks.bed
        """
}