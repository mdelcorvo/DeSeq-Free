# Copyright 2024
# Marcello Del Corvo
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms.

##### Modules #####
include: "rules/common.smk"
include: "rules/qc.smk"
include: "rules/align-fastq.smk"
include: "rules/cna.smk"
if config['variant_calling']:
    include: "rules/variants_analysis.smk"


rule all:
   input:
       f'{final}/report/Processed_data_qc.html',
       f'{final}/qc/qc.pdf',
       f'{final}/cna/plasma-tumor-shared.cna.txt',
       expand(f'{final}/variant_calling/{{sample_calling}}/{{sample_calling}}.plasma-tumor.summary.txt',sample_calling=s2['sample']) if config['variant_calling'] else []