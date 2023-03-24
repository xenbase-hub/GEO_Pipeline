#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process downloadSRR {

    // Takes in SRR values, returns FASTQs
    // prefetch + fasterq-dump is the preferred method

    tag "$SRR"

    module "sratoolkit/3.0.0"

    input:
    tuple val(KEY), val(MERGE_TECH), val(SRR), val(SPECIES), val(ASSAY), val(LIBRARY)

    output:
    tuple val(KEY), val(MERGE_TECH), val(SRR), val(SPECIES), val(ASSAY), val(LIBRARY), path("*.fastq"), emit: raw_data

    script:
    """
    echo ${SRR}
    prefetch --max-size 1000G ${SRR}
    fasterq-dump ${SRR}
    rm -r ${SRR}
    """
}

process mergeFASTQ {

    // Concatenates FASTQ files

    input:
    tuple val(SAMPLE), val(SPECIES), val(ASSAY), val(LIBRARY), path(FASTQ)
        
    output:
    tuple val(SAMPLE), val(SPECIES), val(ASSAY), val(LIBRARY), path("*.fastq"), emit: merged

    script:

    if (LIBRARY == "PAIRED-END")
        """
        cat ${FASTQ.findAll{it.name.endsWith("_1.fastq")}.join(" ")} >> merged_1.fastq
        cat ${FASTQ.findAll{it.name.endsWith("_2.fastq")}.join(" ")} >> merged_2.fastq
        """
    else if (LIBRARY == "SINGLE-END")
        """
        cat ${FASTQ.join(" ")} >> merged.fastq
        """
}