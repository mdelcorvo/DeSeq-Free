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
#include: "rules/align-fastq_after_merge.smk"


rule all:
   input:
       expand(f'{derived}/cna/{{sample_type}}.cna.seg',sample_type=s1[s1["type"]=='plasma']['sample_type']),
       expand(f'{derived}/cna/{{sample_type}}.cna.seg',sample_type=s1[s1["type"]== 'tumor']['sample_type']),
       expand(f'{derived}/qc/{{sample_type}}.sequencing-qc.txt',sample_type=s1["sample_type"]),
       expand(f'{derived}/qc/coverage/{{sample_type}}.mosdepth.global.dist.txt',sample_type=s1["sample_type"]),
       expand(f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.Varscan2.filtered.vcf.gz',sample_calling=s2["sample"]),
       expand(f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.Varscan2.filtered.vcf.gz',sample_calling=s2["sample"]),
       expand(f'{derived}/variant_calling/LoFreq/{{sample_calling}}/{{sample_calling}}-plasma.LoFreq.filtered.vcf.gz',sample_calling=s2["sample"]),
       expand(f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.LoFreq.filtered.vcf.gz',sample_calling=s2["sample"])