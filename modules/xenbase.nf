#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process processRunFile {

    module 'R/4.1.1'

    input:
    path SRR_TABLE

    output:
    path "srr_table.txt"

    script:
    """
    Rscript ${workflow.projectDir}/bin/create_srr_table.R ${SRR_TABLE}
    """
}

process createReadMe {

    publishDir "${params.outDir}", mode: 'copy', pattern: '*.txt', saveAs: {filename ->
        if(SPECIES == "xenopus-trop-v10" && ASSAY == "RNA-SEQ") "${GSE}/XENTR_10.0/RNA-Seq/$filename"
        else if(SPECIES == "xenopus-laevis-v10" && ASSAY == "RNA-SEQ") "${GSE}/XENLA_10.1/RNA-Seq/$filename"
        else if(SPECIES == "xenopus-trop-v10" && ASSAY == "ChIP-Epigenetic") "${GSE}/XENTR_10.0/ChIP-Seq/$filename"
        else if(SPECIES == "xenopus-laevis-v10" && ASSAY == "ChIP-Epigenetic") "${GSE}/XENLA_10.1/ChIP-Seq/$filename"
        else if(SPECIES == "xenopus-trop-v10" && ASSAY == "ChIP-TF") "${GSE}/XENTR_10.0/ChIP-Seq/$filename"
        else if(SPECIES == "xenopus-laevis-v10" && ASSAY == "ChIP-TF") "${GSE}/XENLA_10.1/ChIP-Seq/$filename"
        else if(SPECIES == "xenopus-trop-v10" && ASSAY == "ATAC") "${GSE}/XENTR_10.0/ATAC-Seq/$filename"
        else if(SPECIES == "xenopus-laevis-v10" && ASSAY == "ATAC") "${GSE}/XENLA_10.1/ATAC-Seq/$filename"
    }

    module 'python3'

    input:
    tuple path(SRR_TABLE), val(GSE), val(SPECIES), val(ASSAY) 

    output:
    path("*.txt")

    script:
    """
    python ${workflow.projectDir}/bin/create_readme.py ${SRR_TABLE} ${GSE} ${SPECIES} ${ASSAY}
    """
}