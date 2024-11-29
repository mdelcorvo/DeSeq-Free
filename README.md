# DeSeq-Free <img src="img/DeSeq-Free_logo.png" width="200" align="right" />

**DeSeq-Free** (Whole genome **De**ep **Seq**uencing analysis of Cell **Free** tumor DNA ) is a [Snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html), 
aimed to analyze WGS of circulating cell-free DNA (cfDNA) in the plasma of cancer patients together with their matched germline and tumour samples in a reproducible, automated, and partially contained manner. 
It is implemented such that alternative or similar analysis can be added or removed. 


## Contents

- [Using the DeSeq-Free workflow](#using-the-deseq-free-workflow)
- [Documentation](#documentation)
- [Input files](#input-files)
- [Output files](#output-files)
- [Test Dataset](#test-dataset)
 

## Using the DeSeq-Free workflow

We assume that you already have conda installed, otherwise you can easily install it:

To install conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

* Input:
  
  metafile  (can be .xlsx or .csv) with raw fastq.gz data that looks as follows:
  ```
  sample, lane, fq1, fq2, type
  
  Sample1, lane1, S1_L001_R1_001.fastq.gz, S1_L001_R2_001.fastq.gz, 0
  Sample1, lane2, S1_L002_R1_001.fastq.gz, S1_L002_R2_001.fastq.gz, 0
  ```
  Each row represents a single-end fastq file. Rows with the same sample identifier are considered technical replicates and will be automatically merged. ``` type ``` refers to sample type (0=
  buffy coat, 1= plasma, 2=tumor).

  - Reference genome<br />
   
    Before starting, a user need to download reference genome. 

    Download from [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/), [Ensembl](https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/), or any other autorities
    ```
    wget https://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
    ```
    
    - Index reference genome for bwa-mem2<br />
  
      Prepare indexed genome for bwa-mem2 to boost mapping.  Refer to the [bwa-mem2 instruction](https://github.com/bwa-mem2/bwa-mem2).<br />
      
      Example code:
      ```
      ./bwa-mem2 index <in.fasta>
      Where 
      <in.fasta> is the path to reference sequence fasta file and 
      ```
      
* Code:

  ```
  git clone https://github.com/mdelcorvo/DeSeq-Free.git
  cd DeSeq-Free && conda env create -f envs/workflow.yaml
  conda activate DeSeq-Free_workflow

  snakemake --use-conda \
  --config \
  input=inputfile.xlsx \
  output=output_directory \
  genome=genome.fasta
  ```

## Output files

- Somatic variant analysis
- Variant allele frequency
- Annotation of somatic variants
- Somatic signatures
- Analysis of somatic CNAs
- Fragment size analysis


```
