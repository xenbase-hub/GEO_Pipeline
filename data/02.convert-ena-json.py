#!/usr/bin/env python3
import sys
import json

# print(json.dumps({'4': 5, '6': 7}, sort_keys=True, indent=4))

filename_raw_json = sys.argv[1]

f_raw = open(filename_raw_json, 'r')

out_json = dict()
out_json['study_accession'] = 'NA'
out_json['study_alias'] = 'NA'
out_json['tax_id'] = 'NA'
out_json['scientific_name'] = 'NA'
out_json['run_list'] = dict()

raw_data = json.load(f_raw)
for tmp in raw_data:
    out_json['study_accession'] = tmp['study_accession']
    out_json['study_alias'] = tmp['study_alias']
    out_json['study_title'] = tmp['study_title']
    out_json['tax_id'] = tmp['tax_id']
    out_json['scientific_name'] = tmp['tax_id']

    tmp_run_acc = tmp['run_accession']
    out_json['run_list'][tmp_run_acc] = dict()
    out_json['run_list'][tmp_run_acc]['sample_accession'] = tmp['sample_accession']
    out_json['run_list'][tmp_run_acc]['experiment_accession'] = tmp['experiment_accession']
    out_json['run_list'][tmp_run_acc]['library_layout'] = tmp['library_layout']
    out_json['run_list'][tmp_run_acc]['library_strategy'] = tmp['library_strategy']
    out_json['run_list'][tmp_run_acc]['read_count'] = tmp['read_count']
    out_json['run_list'][tmp_run_acc]['sample_alias'] = tmp['sample_alias']
    out_json['run_list'][tmp_run_acc]['sample_title'] = tmp['sample_title']


print(json.dumps(out_json, sort_keys=True, indent=4))
