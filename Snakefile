# Copyright 2024
# Marcello Del Corvo
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms.

##### Modules #####
include: "rules/common.smk"
include: "rules/align-fastq.smk"
include: "rules/align-fastq_after_merge.smk"


rule all:
   input:
       expand(f'{derived}/recal/{{sample_type}}.bam',sample_type=s1["sample_type"]),
       expand(f'{derived}/recal/{{sample_type}}.notCombined.final.bam',sample_type=s1[s1["type"]=='plasma']['sample_type']),
       expand(f'{derived}/recal/{{sample_type}}.extendedFrags.final.bam',sample_type=s1[s1["type"]=='plasma']['sample_type'])