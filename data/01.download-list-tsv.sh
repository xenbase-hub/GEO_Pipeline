#!/bin/bash
DATE=$(date +"%Y-%m-%d")

XENLA_LIST="XENLA_study_list."$DATE".tsv"
XENTR_LIST="XENTR_study_list."$DATE".tsv"

curl -o $XENLA_LIST -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=read_study&query=tax_eq(8355)&fields=tax_id%2Cstudy_accession%2Cstudy_alias%2Csubmission_accession%2Clibrary_strategy%2Crun_accession%2Cfastq_ftp%2Cstudy_title&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search"
curl -o $XENTR_LIST -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=read_study&query=tax_eq(8364)&fields=tax_id%2Cstudy_accession%2Cstudy_alias%2Csubmission_accession%2Clibrary_strategy%2Crun_accession%2Cfastq_ftp%2Cstudy_title&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search"
