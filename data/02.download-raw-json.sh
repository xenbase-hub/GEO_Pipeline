#!/bin/bash
#PROJECT_ID="PRJNA759592"
PROJECT_ID=$1

OUTFILE=$PROJECT_ID".raw.json"
curl -o $OUTFILE "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$PROJECT_ID&result=read_run&fields=study_accession,sample_accession,experiment_accession,submission_accession,tax_id,scientific_name,instrument_model,library_layout,library_strategy,read_count,study_alias,experiment_alias,fastq_ftp,sample_alias,sample_title,study_title&format=json&download=true&limit=0"
