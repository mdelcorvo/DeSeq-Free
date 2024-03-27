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


def get_IchorCNA(wildcards):
    """Denote sample name and platform in read group."""
    return (r"--id {sample_type} --gcWig {gcwig} --mapWig {mapwig} --maxCN {maxCN} "
            r"--includeHOMD {includeHOMD} --chrs '{chrs_ichorCNA}' --chrTrain '{chrTrain}' "
            r"--genomeStyle {genomeStyle} --estimateNormal {estimateNormal} "
            r"--estimatePloidy {estimatePloidy} --estimateScPrevalence {estimateClonality} "
            r"--scStates '{scStates}' --centromere {centromere} --genomeBuild hg38 "
            r"--txnE {txnE} --txnStrength {txnStrength} "
            r" --minMapScore  {minMapScore} --fracReadsInChrYForMale  {fracReadsChrYMale} "
            r"--maxFracGenomeSubclone {maxFracGenomeSubclone} --maxFracCNASubclone {maxFracCNASubclone} "
            r" --plotFileType {plotFileType} --plotYLim '{plotYlim}' ").format(
        sample_type=wildcards.sample_type,
        gcwig=gcwig,mapwig=mapwig,
        maxCN=maxCN,includeHOMD=includeHOMD,
        chrs_ichorCNA=chrs_ichorCNA,chrTrain=chrTrain,
        genomeStyle=genomeStyle,
        estimatePloidy=estimatePloidy,estimateNormal=estimateNormal,
        estimateClonality=estimateClonality,scStates=scStates,
        centromere=centromere,txnE=txnE,txnStrength=txnStrength,
        minMapScore=minMapScore,fracReadsChrYMale=fracReadsChrYMale,
        maxFracGenomeSubclone=maxFracGenomeSubclone,
        maxFracCNASubclone=maxFracCNASubclone,
        plotFileType=plotFileType,plotYlim=plotYlim)

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

s2=s1.drop_duplicates('sample')['sample'].to_frame().set_index(["sample"], drop=False)
s2.index.names = ['sample_calling']

####### HELPER VARIABLES #######
# Directories
outputdir = getpath(config["output"])
derived = outputdir + config["derived_data"]
final = outputdir + config["final_data"]
logs = outputdir + config["logs"]
tmpdir = config["tmpdir"]
benchmarks = outputdir + config["benchmarks"]
genome_data = config["genome_data"]

binSize = config["binSize"]
qual= config["qual"]
chrs= config["chrs"]

gcwig=config["ichorCNA_gcWig"]
mapwig=config["ichorCNA_mapWig"]
maxCN=config["ichorCNA_maxCN"]
includeHOMD=config["ichorCNA_includeHOMD"]
chrs_ichorCNA=config["ichorCNA_chrs"]
chrTrain=config["ichorCNA_chrTrain"]
genomeStyle=config["ichorCNA_genomeStyle"]
estimatePloidy=config["ichorCNA_estimatePloidy"]
estimateNormal=config["ichorCNA_estimateNormal"]
estimateClonality=config["ichorCNA_estimateClonality"]
scStates=config["ichorCNA_scStates"]
centromere=config["ichorCNA_centromere"]
txnE=config["ichorCNA_txnE"]
txnStrength=config["ichorCNA_txnStrength"]
minMapScore=config["ichorCNA_minMapScore"]
fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"]
maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"]
maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"]
plotFileType=config["ichorCNA_plotFileType"]
plotYlim=config["ichorCNA_plotYlim"]