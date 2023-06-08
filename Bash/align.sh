#!bin/bash

process alignRNA {

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
