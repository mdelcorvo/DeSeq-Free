configfile: "config.yaml"

import os
import pandas as pd
import numpy as np
import re
import glob
import sys
from pathlib import Path

from snakemake.utils import validate
from snakemake.utils import min_version

def detect_file_type(file_path):
    if not os.path.exists(file_path):
        return 'File does not exist'
    _, extension = os.path.splitext(file_path)
    if extension.lower() == '.csv':
        try:
            samples = pd.read_csv(file_path, dtype=str)
            return samples
        except pd.errors.ParserError:
            return 'Invalid CSV format'
        except Exception:
            return 'Error reading CSV file'
    elif extension.lower() == '.xlsx':
        try:
            samples = pd.read_excel(file_path, dtype=str,engine="openpyxl")
            return samples
        except Exception:
            return 'Error reading Excel (XLSX) file'
    else:
        return 'Unknown format'

def make_absolute(path):
    """Convert a relative path to an absolute path."""
    if not os.path.isabs(path):
        return os.path.abspath(path)
    return path

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

def is_single_end(id):
    """Return True if sample is single end."""
    return pd.isnull(samples.loc[(id), "fq2"])

def is_single_end1(sample_type):
    """Return True if sample is single end."""
    return pd.isnull(s1.loc[(sample_type), "fq2"])

def is_plasma_sample(sample_type):
    """Return True if sample is single end."""
    return s1.loc[(sample_type),'type']=='plasma'

def is_tumor_sample(sample_type):
    """Return True if sample is single end."""
    return s1.loc[(sample_type),'type']=='tumor'

def is_control_sample(sample_type):
    """Return True if sample is single end."""
    return s1.loc[(sample_type),'type']=='control'

###########
#For unpack
###########

def get_fastq(wildcards):
    """Get fastq files of given sample."""
    fastqs = samples.loc[(wildcards.id), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def fastq(wildcards):
    if not is_single_end1(**wildcards):
        fq1= f'{derived}/pre_processing/final/{{sample_type}}-paired.R1.fastq.gz'
        fq2 = f'{derived}/pre_processing/final/{{sample_type}}-paired.R2.fastq.gz'
        return {"r1": fq1, "r2": fq2}
    fq1 = f'{derived}/pre_processing/final/{{sample_type}}.fastq.gz'
    return fq1

def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        r1 = f'{derived}/pre_processing/temp/{{id}}.R1.fastq.gz'
        r2 = f'{derived}/pre_processing/temp/{{id}}.R2.fastq.gz'
        return {"r1": r1, "r2": r2}
    r1 = f'{derived}/pre_processing/temp/{{id}}.fastq.gz'
    return r1

def get_pileup(wildcards):
    plasma = f'{derived}/variant_calling/Varscan2/pileup/{{sample_calling}}-plasma.pileup'
    tumor = f'{derived}/variant_calling/Varscan2/pileup/{{sample_calling}}-tumor.pileup'
    control = f'{derived}/variant_calling/Varscan2/pileup/{{sample_calling}}-control.pileup'
    return {"plasma": plasma, "tumor": tumor, "control": control}

def get_bam(wildcards):
    plasma = f'{derived}/recal/{{sample_calling}}-plasma.bam'
    tumor = f'{derived}/recal/{{sample_calling}}-tumor.bam'
    control = f'{derived}/recal/{{sample_calling}}-control.bam'
    return {"plasma": plasma, "tumor": tumor, "control": control}

######################
# For align-fastq rule
######################

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample."""
    if not is_single_end1(**wildcards):
        # paired-end sample
        return expand(f'{derived}/pre_processing/final/{{sample_type}}-paired.R{{group}}.fastq.gz',
                      group=[1, 2], **wildcards)
    # single end sample
    return f'{derived}/pre_processing/final/{{sample_type}}-single.fastq.gz'.format(**wildcards)

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample_type}\tSM:{sample_type}\tPL:{sample_type}'".format(
        sample_type=wildcards.sample_type,
        platform='ILLUMINA')

#################
# For CNA calling
#################

def get_IchorCNA(wildcards):
    """Denote sample name and platform in read group."""
    return (r"{gcwig} {mapwig} "
            r"{centromere} "
            r"{minMapScore}").format(
        gcwig=gcwig,
        mapwig=mapwig,
        centromere=centromere,
        minMapScore=minMapScore)

def get_bwa2_index(wildcards):
    """Denote sample name and platform in read group."""
    return (r"{genome} {genome}.amb {genome}.ann {genome}.bwt.2bit.64 {genome}.pac ").format(
        genome=genome)

def get_plasma_wig(wildcards):
    plasma = f'{derived}/cna/wig/{{sample_calling}}-plasma.wig'
    return {"plasma": plasma}

def get_tumor_wig(wildcards):
    tumor = f'{derived}/cna/wig/{{sample_calling}}-tumor.wig'
    return {"tumor": tumor}

samples = detect_file_type(make_absolute(config['input']))
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

s2=s1.drop_duplicates('sample')['sample'].to_frame().set_index(["sample"], drop=False)
s2.index.names = ['sample_calling']

####### HELPER VARIABLES #######
# Directories
outputdir = getpath(config["output"])
derived = outputdir + config["derived_data"]
final = outputdir + config["final_data"]
logs = outputdir + config["logs"]
benchmarks = outputdir + config["benchmarks"]
genome = config["genome"]
genome_fai=genome+'.fai'
genome_dict = os.path.splitext(genome)[0] + ".dict"
exon_transcripts_hg38 = os.path.join(config["annotation"], "canonical_exon_transcripts_hg38.bed")

binSize = config["binSize"]
qual= config["qual"]
chrs= config["chrs"]

gcwig=config["ichorCNA_gcWig"]
mapwig=config["ichorCNA_mapWig"]
centromere=config["ichorCNA_centromere"]
minMapScore=config["ichorCNA_minMapScore"]
