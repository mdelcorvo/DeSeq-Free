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
include: "rules/variants_analysis.smk"


rule all:
   input:
       f'{final}/cna/plasma-tumor-shared.cna.txt',
       expand(f'{derived}/qc/{{sample_type}}.sequencing-qc.txt',sample_type=s1["sample_type"]),
       expand(f'{derived}/qc/coverage/{{sample_type}}.mosdepth.global.dist.txt',sample_type=s1["sample_type"])
