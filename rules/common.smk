configfile: "config.yaml"

import os
import pandas as pd
import numpy as np
import re
import glob

from snakemake.utils import validate
from snakemake.utils import min_version

def check_fastq(path):
    if os.path.isfile(path):
       return 1
    else:
       return len(glob.glob1(path,"*.fastq*"))

def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

def get_fastq(wildcards):
    """Get fastq files of given sample."""
    fastqs = samples.loc[(wildcards.id), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def is_single_end(id):
    """Return True if sample is single end."""
    return pd.isnull(samples.loc[(id), "fq2"])

def is_single_end1(sample_type):
    """Return True if sample is single end."""
    return pd.isnull(s1.loc[(sample_type), "fq2"])

def is_plasma_sample(sample_type):
    """Return True if sample is single end."""
    return s1.loc[(sample_type),'type']=='plasma'

###########
#For unpack
###########

def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        r1 = f'{derived}/pre_processing/temp/{{id}}.R1.fastq.gz'
        r2 = f'{derived}/pre_processing/temp/{{id}}.R2.fastq.gz'
        return {"r1": r1, "r2": r2}
    r1 = f'{derived}/pre_processing/temp/{{id}}.fastq.gz'
    return r1

#############
# For bwa mem
#############

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample."""
    if not is_single_end1(**wildcards):
        # paired-end sample
        return expand(f'{derived}/pre_processing/final/{{sample_type}}-paired.R{{group}}.fastq.gz',
                      group=[1, 2], **wildcards)
    # single end sample
    return f'{derived}/pre_processing/final/{{sample_type}}-single.fastq.gz'.format(**wildcards)

def get_trimmed_reads_notCombined(wildcards):

    if is_plasma_sample(**wildcards):
        # paired-end sample
        return expand(f'{derived}/pre_processing/final/{{sample_type}}.notCombined_{{group}}.fastq.gz',
                      group=[1, 2], **wildcards)

def get_trimmed_reads_extendedFrags(wildcards):

    if is_plasma_sample(**wildcards):
        # paired-end sample
        return f'{derived}/pre_processing/final/{{sample_type}}.extendedFrags.fastq.gz'.format(**wildcards)

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample_type}\tSM:{sample_type}\tPL:{sample_type}'".format(
        sample_type=wildcards.sample_type,
        platform='ILLUMINA')



samples = pd.read_excel(config['meta'], dtype=str)

samples = samples[~samples['fq1'].apply(check_fastq) != 1]
samples = samples[~samples['fq2'].apply(check_fastq) != 1]

samples['type']=samples['type'].astype(str).astype(int)
samples['type'] = np.where(samples['type']==0, 'control',
 np.where(samples['type']==1, 'plasma', 'tumor'))
samples['fastq']= np.where(samples['fq2'].isnull(), 'single', 'paired')

samples["id"] = samples['sample'] +"-"+ samples["type"] +"-"+ samples["lane"]+"-"+ samples["fastq"]
samples['sample_type'] = samples['sample'] +"-"+ samples["type"]
samples = samples.set_index(["id"], drop=False)

s1=samples.drop_duplicates('sample_type')
s1=s1.set_index(["sample_type"], drop=False)
#plasma=s1.loc[s1["type"]=='plasma']
#plasma['plasma']= plasma['sample_type']
#plasma=plasma.set_index(["plasma"], drop=False)
#del plasma["sample_type"]

####### HELPER VARIABLES #######
# Directories
outputdir = getpath(config["output"])
derived = outputdir + config["derived_data"]
final = outputdir + config["final_data"]
logs = outputdir + config["logs"]
tmpdir = config["tmpdir"]
benchmarks = outputdir + config["benchmarks"]

genome_data = config["genome_data"]

