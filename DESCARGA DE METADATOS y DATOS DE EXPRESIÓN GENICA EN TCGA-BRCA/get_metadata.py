#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 15:52:44 2023

@author: laura
"""

import requests
import json

fields = [
    "cases.submitter_id",
    "cases.sample_ids",
    "cases.case_id",
    "cases.primary_site",
    "cases.project.name",
    "cases.samples.sample_type",
    "cases.samples.sample_id", 
    "diagnoses.vital_status", 
    "file_id", 
    "file_name"
 ]

fields = ",".join(fields)

files_endpt = "https://api.gdc.cancer.gov/files"

filters = {
    "op": "and",
    "content":[
        {
        "op": "in",
        "content":{
            "field": "cases.project.project_id",
            "value": ["TCGA-BRCA"]
            }
        },
        {
        "op": "in",
        "content":{
            "field":"files.data_type",
            "value":"Gene Expression Quantification"
            }
        }
    ]
}

# With a GET request, the filters parameter needs to be converted
# from a dictionary to JSON-formatted string

params = {
    "filters": json.dumps(filters),
    "fields": fields,
    "format": "TSV" ,
    "size": 50000
    }


# The parameters are passed to 'json' rather than 'params' in this case
response = requests.get(files_endpt, params = params)

file_name = "TCGA-BRCA_metadata.txt"

with open(file_name, "wb") as output_file:
    output_file.write(response.content)
