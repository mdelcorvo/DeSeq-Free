# DeSeq-Free <img src="img/DeSeq-Free_logo.png" width="200" align="right" />

**DeSeq-Free** (Whole genome **De**ep **Seq**uencing analysis of Cell **Free** tumor DNA ) is a [Snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html), 
aimed to analyze WGS of circulating cell-free DNA (cfDNA) in the plasma of cancer patients in a reproducible, automated, and partially contained manner. 
It is implemented such that alternative or similar analysis can be added or removed. 

## Contents

- [Using the DeSeq-Free workflow](#using-the-snhichip-workflow)
- [Documentation](#documentation)
- [Input files](#input-files)
- [Output files](#output-files)
- [Test Dataset](#test-dataset)
 

## Using the DeSeq-Free workflow

We assume that you already have conda installed, otherwise you can easily install it:

To install conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

In order to ease the use of DeSeq-Free, we provide a yml file for conda with all required tools, including Snakemake. 
```
To use DeSeq-Free:

git clone https://github.com/mdelcorvo/DeSeq-Free.git
cd DeSeq-Free && conda env create -f envs/workflow.yaml
conda activate DeSeq-Free_workflow

#edit config and meta file
snakemake -j 8
```


## Documentation

The pipeline leverages several tools to QC DeSeq-Free library, create statistics/interactive report and calculate/annotate interaction matrices at different bin size: 
[bwa mem](https://bio-bwa.sourceforge.net/bwa.shtml), [pairtools](https://pairtools.readthedocs.io/en/latest/index.html), [juicer](https://github.com/aidenlab/juicer), [cooler](https://cooler.readthedocs.io/en/latest/index.html),
[pairix](https://github.com/4dn-dcic/pairix), [Macs2](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) and [FitHiChIP](https://ay-lab.github.io/FitHiChIP/html/index.html).

You will need to specify the location of the `reference genome` (hg38) in fasta/fa format with bwa index.
Use the parameter `genome_data` in the config file to add it.

## Input files

Users are required to provide a metadata file for running the **DeSeq-Free** workflow:

- **metadata file**  – a _tab-delimited_  text  file  listing  the  name  of  the  samples,  the  sequencing  technology  and  the paths to raw paired FASTQ files

| sample        | platform      |fq1     |fq2    |
| ------------- |:-------------:| :-----:|:-----:|
| Sample1       | ILLUMINA      | data/S1_1.fastq.gz |data/S1_2.fastq.gz |
| Sample2       | ILLUMINA      | data/S2_1.fastq.gz |data/S2_2.fastq.gz |
| Sample3       | ILLUMINA      | data/S3_1.fastq.gz |data/S3_2.fastq.gz |

- **configuration file**

The configuration file (`config.yaml`) contains all the paths to input, output and reference files and additional parameters to customize the pipeline and the performed tests. All of these need to be carefully specified in accordance with the specific experiment.

**Important**: ALL relative paths will be interpreted relative to the directory where the Snakefile is located. Alternatively, you can use absolute paths.

- **reference in a fasta file format**, e.g. hg38 with bwa index

## Output files

- Somatic variant analysis
- Variant allele frequency
- Annotation of somatic variants
- Somatic signatures
- Analysis of somatic CNAs
- Fragment size analysis


```
