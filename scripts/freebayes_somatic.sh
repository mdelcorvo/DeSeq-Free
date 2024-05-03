#!/bin/bash
if [[ $# -eq 0 ]] ; then
    echo """
********************************************************
Welcome to freebayes somatic calling with matched normal
********************************************************

Copyright 2024 IEO - Data Science Unit
Marcello Del Corvo, PhD

Description:
Freebayes parallel somatic calling (inspired by speedseq https://github.com/hall-lab/speedseq/)
Default parameters are set to maximize somatic variant calling

Usage: freebayes_somatic.sh --tumor <tumor bam file> --normal <normal bam file> --ref <reference genome> [OPTIONS]

* You must have an active Conda environment to automatically install the required software *

Required argument:

  -t, --tumor  tumor bam file
  -n, --normal normal bam file
  -r, --ref   Reference fasta file

Options:

  -t  --tumor   tumor bam file
  -n  --normal  normal bam file
  -f  --fasta-reference Reference fasta file
  -d  --tmp Folder path for temporary data (default: ./tmp)
  -v  --vcf Output VCF-format results to FILE. (default: tumor-name.vcf.gz)
  -P  --pon panels of normals
  -m  --min-mapping-quality Exclude alignments from analysis if they have a mapping quality less than Q.  default: 1
  -q  --min-base-quality Exclude alleles from analysis if their supporting base quality is less than Q.  default: 0
  -F  --min-alternate-fraction   Require at least this fraction of observations supporting an alternate allele within a single individual in order to evaluate the position.  default: 0.02
  -C  --min-alternate-count Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position.  default: 2
  -E  --min-repeat-entropy To detect interrupted repeats, build across sequence until it has entropy > N bits per bp. Set to 0 to turn off. (default: 1)
  -p  --ploidy Sets the default ploidy for the analysis to N.  default: 2
  -e  --bed   Bed file with specifics position for statistics (default: NONE)
  -c  --threads   Number of cores (default: 1)


"""
    exit 0
fi

while [[ "$#" -gt 0 ]]
do case $1 in
    -t|--tumor) tumor="$2"
    shift;;
    -n|--normal) normal="$2"
    shift;;
    -f|--fasta-reference) ref="$2"
    shift;;
    -d|--tmp) tmp="$2"
    shift;;
    -v|--vcf) vcf="$2"
    shift;;
    -P|--pon) pon="$2"
    shift;;
    -m|--min-mapping-quality) mmq="$2"
    shift;;
    -q|--min-base-quality) mbq="$2"
    shift;;
    -F|--min-alternate-fraction) maf="$2"
    shift;;
    -C|--min-alternate-count) mac="$2"
    shift;;
    -E|--min-repeat-entropy) mre="$2"
    shift;;
    -p|--ploidy) ploidy="$2"
    shift;;
    -e|--bed) bed="$2"
    shift;;
    -c|--threads) threads="$2"
    shift;;
    *) echo "Unknown parameter passed: $1"
    exit 1;;
esac
shift
done

if [ -z "$tumor" ]; then echo "tumor bam is required!"; exit; fi
if [ -z "$normal" ]; then echo "normal bam is required!"; exit; fi
if [[ $tumor != *.bam ]]; then echo "the tumor input file don't seem to be in the right format"; exit; fi
if [[ $normal != *.bam ]]; then echo "the tumor input file don't seem to be in the right format"; exit; fi
if [ -z "$pon" ]; then echo "panels of normals is required!"; exit; fi
if [ -z "$ref" ]; then echo "need to provide a reference genome"; exit; fi

if [ -z "$tmp" ]; then mkdir ./tmp;tmp=tmp; fi
if [ -z "$vcf" ]; then tumor_name="${tumor%.*}" vcf=${tumor_name}.vcf.gz; fi

if [ -z "$mmq" ]; then mmq=1; fi
if [ -z "$mbq" ]; then mbq=0; fi
if [ -z "$maf" ]; then maf=0.02; fi
if [ -z "$mac" ]; then mac=2; fi
if [ -z "$mre" ]; then mre=1; fi
if [ -z "$ploidy" ]; then ploidy=2; fi

if [ -z "$threads" ]; then threads=1; fi

##########
#functions
##########

# deletes the temp directory
function cleanup {
  rm -rf "$tmp"
  echo "Deleted temp working directory $tmp"
}

function somatic_filter() {
    awk -v MINQUAL="$1" -v SSC_THRES="$2" -v ONLY_SOMATIC="$3" 'BEGIN {NORMAL=10; TUMOR=11; GL_IDX=0;}
    {
        if ($0~"^#") { print ; next; }
        if (! GL_IDX) {
            split($9,fmt,":")
            for (i=1;i<=length(fmt);++i) { if (fmt[i]=="GL") GL_IDX=i }
        }
        split($NORMAL,N,":");
        split(N[GL_IDX],NGL,",");
        split($TUMOR,T,":");
        split(T[GL_IDX],TGL,",");
        LOD_NORM=NGL[1]-NGL[2];
        LOD_TUMOR_HET=TGL[2]-TGL[1];
        LOD_TUMOR_HOM=TGL[3]-TGL[1];
        if (LOD_TUMOR_HET > LOD_TUMOR_HOM) { LOD_TUMOR=LOD_TUMOR_HET }
        else { LOD_TUMOR=LOD_TUMOR_HOM }
        DQUAL=LOD_TUMOR+LOD_NORM;
        if (DQUAL>=SSC_THRES && $NORMAL~"^0/0") {
            $7="PASS"
            $8="SSC="DQUAL";"$8
            print
        }
        else if (!ONLY_SOMATIC && $6>=MINQUAL && $10~"^0/0" && ! match($11,"^0/0")) {
            $8="SSC="DQUAL";"$8
            print
        }
    }' OFS="\t"
}

export -f somatic_filter

#####################################

#####################################

tumor_name="`basename $tumor`"
normal_name="`basename $normal`"
tname="${tumor_name%.*}"
nname="${normal_name%.*}"

if [ -z "$bed" ]; then

#Phase 1: raw somatic variant calling

  printf "%s\n" $(seq -f "chr%g" 22) chrX chrY \
  | parallel -j ${threads} \
  "freebayes \
  -f $ref \
  --region {} \
  --allele-balance-priors-off \
  --pooled-discrete --pooled-continuous \
  --genotype-qualities \
  --min-repeat-entropy ${mre} \
  --min-alternate-fraction ${maf} \
  --min-alternate-count ${mac} \
  --ploidy ${ploidy} \
  ${tumor} \
  ${normal} \
  | bgzip -c > ${tmp}/{}.raw.${tname}.vcf.gz"

  printf "%s\n" $(seq -f "chr%g" 22) chrX chrY \
  | parallel -j ${threads} \
  "tabix -p vcf ${tmp}/{}.raw.${tname}.vcf.gz"


#Phase 2: Somatic filtering

  printf "%s\n" $(seq -f "chr%g" 22) chrX chrY \
  | parallel -j ${threads} \
  "vcfallelicprimitives \
  ${tmp}/{}.raw.${tname}.vcf.gz \
  -t DECOMPOSED --keep-geno \
  | vcfsamplediff -s VT ${nname} ${tname} - \
  | bcftools filter \
  --soft-filter 'Reject' -e '(FORMAT/DP < 25) || (SAF < 0 && SAR < 0) || (RPR < 0 && RPL < 0) || (AF[0] < 0.05 && AF[0] > 0.5)' -m '+' \
  | vcffilter -f 'VT = somatic' \
  | SnpSift annotate ${pon} | awk '( \$3 < 0.001 || \$3 == \"ID\" )' \
  | somatic_filter 1e-4 18 0 \
  | bgzip -c > ${tmp}/{}.filt.${tname}.vcf.gz"

  printf "%s\n" $(seq -f "chr%g" 22) chrX chrY \
  | parallel -j ${threads} \
  "tabix -p vcf ${tmp}/{}.filt.${tname}.vcf.gz"

else

#Phase 1: raw somatic variant calling

  if [[ $bed != *.bed ]]; then echo "the bed file don't seem to be in the right format"; exit; fi
  mkdir -p ${tmp}/bed/${tname}
  Rscript --vanilla scripts/split.bed.R ${bed} ${tmp}/bed/${tname} ${threads}
  cd ${tmp}/bed/${tname}
  ls | sed -n 's/\.bed$//p' \
  | parallel -j ${threads} \
  "freebayes \
  -f $ref \
  --targets {} \
  --allele-balance-priors-off \
  --pooled-discrete --pooled-continuous \
  --genotype-qualities \
  --min-repeat-entropy ${mre} \
  --min-alternate-fraction ${maf} \
  --min-alternate-count ${mac} \
  --ploidy ${ploidy} \
  ${tumor} \
  ${normal} \
  | bgzip -c > ${tmp}/{}.raw.${tname}.vcf.gz"

  ls | sed -n 's/\.bed$//p' \
  | parallel -j ${threads} \
  "tabix -p vcf ${tmp}/{}.raw.${tname}.vcf.gz"

#Phase 2: Somatic filtering

  ls | sed -n 's/\.bed$//p' \
  | parallel -j ${threads} \
  "vcfallelicprimitives \
  ${tmp}/{}.raw.${tname}.vcf.gz \
  -t DECOMPOSED --keep-geno \
  | vcfsamplediff -s VT ${nname} ${tname} - \
  | bcftools filter \
  --soft-filter 'Reject' -e '(FORMAT/DP < 25) || (SAF < 0 && SAR < 0) || (RPR < 0 && RPL < 0) || (AF[0] < 0.05 && AF[0] > 0.5)' -m '+' \
  | vcffilter -f 'VT = somatic' \
  | SnpSift annotate ${pon} | awk '( \$3 < 0.01 || \$3 == \"ID\" )' \
  | somatic_filter 1e-4 18 0 \
  | bgzip -c > ${tmp}/{}.filt.${tname}.vcf.gz"

  ls | sed -n 's/\.bed$//p' \
  | parallel -j ${threads} \
  "tabix -p vcf ${tmp}/{}.filt.${tname}.vcf.gz"

fi

vcfcat  ${tmp}/*.filt.${tname}.vcf.gz \
| vcffixup - | vcfuniq  | bgzip -c > ${vcf}

tabix -p vcf ${vcf}

# register the cleanup function to be called on the EXIT signal
trap cleanup EXIT