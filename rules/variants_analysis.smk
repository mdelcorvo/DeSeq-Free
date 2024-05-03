###### Rules for somatic variant analysis ######
#     Somatic variant analysis with VarScan2 / LoFreq
#     # Varscan2
#     1. Generating mpileup.
#     2. VarScan2
#     3. Filtering somatic calls from Varscan2
#     4. Extraction of shared variants (i.e. variants shared between plasma and tumour)
#     5. Extraction of unique variants (i.e. variants not shared between plasma and tumour)
##
#     # LoFreq
#     1. Calling somatic variants: lofreq somatic
#
#     Freebayes
#     1. Calling somatic variants with Freebayes
##########################################################################################

##########
# Varscan2
##########


rule mpileup:
    input:
        bam=f'{derived}/recal/{{sample_type}}.bam',
        ref=f'{genome_data}/GRCh38.d1.vd1.fa',
    output:
        pileup=f'{derived}/variant_calling/Varscan2/pileup/{{sample_type}}.pileup'
    log: f'{logs}/{{sample_type}}/01_Varscan2.log'
    benchmark: f'{benchmarks}/{{sample_type}}_Varscan2.txt'
    threads: 1
    group: "variant_analysis"
    priority: 40
    shell:
        """
        samtools mpileup \
        -q 2 -f {input.ref} \
        {input.bam} > {output.pileup}

        """

rule VarScan2:
    input:
        unpack(get_pileup)
    params:
        varscan = "tools/VarScan.jar",
        plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma',
        tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor',
    output:
        snp_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.snp',
        snp_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.snp',
        indel_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.indel',
        indel_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.indel'
    threads: 1
    priority: 40
    shell:
        """
        java -Xmx16g -jar {params.varscan} somatic  {input.control}  {input.plasma} {params.plasma} --min-var-freq 0.01
        
        java -Xmx16g -jar {params.varscan} somatic  {input.control}  {input.tumor} {params.tumor} --min-var-freq 0.01

        """


rule VarScan2_process:
    input:
        snp_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.snp',
        snp_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.snp',
        indel_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.indel',
        indel_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.indel'
    params:
        varscan = "tools/VarScan.jar"
    output:
        snp_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.snp.Somatic.hc',
        snp_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.snp.Somatic.hc',
        indel_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.indel.Somatic.hc',
        indel_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.indel.Somatic.hc',
        snp_plasma_list=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.snp.Somatic.hc.list',
        snp_tumor_list=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.snp.Somatic.hc.list',
        indel_plasma_list=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.indel.Somatic.hc.list',
        indel_tumor_list=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.indel.Somatic.hc.list'
    threads: 1
    priority: 40
    shell:
        """
        java -Xmx16g -jar {params.varscan} processSomatic {input.snp_plasma} --min-tumor-freq 0.01 --max-normal-freq 0.00
        java -Xmx16g -jar {params.varscan} processSomatic {input.snp_tumor} --min-tumor-freq 0.01 --max-normal-freq 0.00
        java -Xmx16g -jar {params.varscan} processSomatic {input.indel_plasma} --min-tumor-freq 0.01 --max-normal-freq 0.00
        java -Xmx16g -jar {params.varscan} processSomatic {input.indel_tumor} --min-tumor-freq 0.01 --max-normal-freq 0.00
        
        awk 'BEGIN{{OFS="\t"}} NR!=1 {{print $1, $2, $2}}' {output.snp_plasma} > {output.snp_plasma_list}
        awk 'BEGIN{{OFS="\t"}} NR!=1 {{print $1, $2, $2}}' {output.snp_tumor} > {output.snp_tumor_list}
        awk 'BEGIN{{OFS="\t"}} NR!=1 {{print $1, $2, $2}}' {output.indel_plasma} > {output.indel_plasma_list}
        awk 'BEGIN{{OFS="\t"}} NR!=1 {{print $1, $2, $2}}' {output.indel_tumor} > {output.indel_tumor_list}
        
        """


rule bam_readcount:
    input:
        snp_plasma_list=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.snp.Somatic.hc.list',
        snp_tumor_list=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.snp.Somatic.hc.list',
        indel_plasma_list=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.indel.Somatic.hc.list',
        indel_tumor_list=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.indel.Somatic.hc.list',
        bam_plasma=f'{derived}/recal/{{sample_calling}}-plasma.bam',
        bam_tumor=f'{derived}/recal/{{sample_calling}}-tumor.bam',
        ref=f'{genome_data}/GRCh38.d1.vd1.fa'
    params:
        varscan="tools/VarScan.jar"
    output:
        snp_rc_plasma=f'{derived}/variant_calling/Varscan2/bam_readcount/{{sample_calling}}-plasma.somatic_hc.snp.readcount',
        snp_rc_tumor=f'{derived}/variant_calling/Varscan2/bam_readcount/{{sample_calling}}-tumor.somatic_hc.snp.readcount',
        indel_rc_plasma=f'{derived}/variant_calling/Varscan2/bam_readcount/{{sample_calling}}-plasma.somatic_hc.indel.readcount',
        indel_rc_tumor=f'{derived}/variant_calling/Varscan2/bam_readcount/{{sample_calling}}-tumor.somatic_hc.indel.readcount'
    conda:
        "../envs/bam_read-count.yaml"
    threads: 1
    priority: 40
    shell:
        """
        bam-readcount -q 2 -f {input.ref} \
        -l {input.snp_plasma_list} {input.bam_plasma} \
        > {output.snp_rc_plasma}

        bam-readcount -q 2 -f {input.ref} \
        -l {input.snp_tumor_list} {input.bam_tumor} \
        > {output.snp_rc_tumor}
        
        bam-readcount -q 2 -f {input.ref} \
        -l {input.indel_plasma_list} {input.bam_plasma} \
        > {output.indel_rc_plasma}

        bam-readcount -q 2 -f {input.ref} \
        -l {input.indel_tumor_list} {input.bam_tumor} \
        > {output.indel_rc_tumor}
        """

rule VarScan2_filtering:
    input:
        snp_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.snp.Somatic.hc',
        snp_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.snp.Somatic.hc',
        indel_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.indel.Somatic.hc',
        indel_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.indel.Somatic.hc',
        snp_rc_plasma=f'{derived}/variant_calling/Varscan2/bam_readcount/{{sample_calling}}-plasma.somatic_hc.snp.readcount',
        snp_rc_tumor=f'{derived}/variant_calling/Varscan2/bam_readcount/{{sample_calling}}-tumor.somatic_hc.snp.readcount',
        indel_rc_plasma=f'{derived}/variant_calling/Varscan2/bam_readcount/{{sample_calling}}-plasma.somatic_hc.indel.readcount',
        indel_rc_tumor=f'{derived}/variant_calling/Varscan2/bam_readcount/{{sample_calling}}-tumor.somatic_hc.indel.readcount'
    params:
        varscan="tools/VarScan.jar"
    output:
        snp_plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.snp.fpfiltered',
        snp_tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.snp.fpfiltered',
        indel_plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.indel.fpfiltered',
        indel_tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.indel.fpfiltered',
        snp_plasma_edit=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.snp.fpfiltered.edit',
        snp_tumor_edit=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.snp.fpfiltered.edit'
    threads: 1
    priority: 40
    shell:
        """
        java -jar {params.varscan} fpfilter {input.snp_plasma} {input.snp_rc_plasma} --output-file {output.snp_plasma} --min-var-freq 0.01 --min-strandedness 0.0
        java -jar {params.varscan} fpfilter {input.snp_tumor} {input.snp_rc_tumor} --output-file {output.snp_tumor} --min-var-freq 0.01 --min-strandedness 0.0
        
        java -jar {params.varscan} fpfilter {input.indel_plasma} {input.indel_rc_plasma} --output-file {output.indel_plasma} --min-var-freq 0.01 --min-strandedness 0.0
        java -jar {params.varscan} fpfilter {input.indel_tumor} {input.indel_rc_tumor} --output-file {output.indel_tumor} --min-var-freq 0.01 --min-strandedness 0.0
        
        awk '{{if($10 >= 5 && $6 == 0 && $5 + $6 >= 10 && $9 + $10 >= 10) {{print $0}} }}' {output.snp_plasma} > {output.snp_plasma_edit}
        awk '{{if($10 >= 5 && $6 == 0 && $5 + $6 >= 10 && $9 + $10 >= 10) {{print $0}} }}' {output.snp_tumor} > {output.snp_tumor_edit}
        
        """

rule VarScan2_shared_variants:
    input:
        snp_plasma_edit=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.snp.fpfiltered.edit',
        snp_tumor_edit=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.snp.fpfiltered.edit'
    output:
        snp_shared=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma-tumor.shared.snp',
        snp_unique=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma-tumor.unique.snp'
    threads: 1
    priority: 40
    shell:
        """
        awk 'FNR==NR{{a[$1,$2]=$0;next}}{{if(b=a[$1,$2]){{print b}} }}' \
        {input.snp_plasma_edit} {input.snp_tumor_edit} \
        > {output.snp_shared}
        
        awk 'FNR==NR{{a[$1,$2]++}}FNR!=NR && !a[$1,$2]{{print}}' \
        {input.snp_tumor_edit} {input.snp_plasma_edit} \
        > {output.snp_unique}
        
        """

rule Varscan2_snp2vcf:
    input:
        snp_shared = f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma-tumor.shared.snp',
        snp_unique = f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma-tumor.unique.snp',
        snp_plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.snp.fpfiltered',
        snp_tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.snp.fpfiltered',
        indel_plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.indel.fpfiltered',
        indel_tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.indel.fpfiltered',
    params:
        vs_format_converter="scripts/vs_format_converter.py"
    output:
        snp_plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.snp.fpfiltered.vcf.gz',
        snp_tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.snp.fpfiltered.vcf.gz',
        indel_plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.indel.fpfiltered.vcf.gz',
        indel_tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.indel.fpfiltered.vcf.gz',
        plasma=f'{final}/variant_calling/plasma/{{sample_calling}}-plasma.Varscan2.final.vcf.gz',
        tumor=f'{final}/variant_calling/tumor/{{sample_calling}}-tumor.Varscan2.final.vcf.gz'
    threads: 1
    priority: 40
    shell:
        """
        python {params.vs_format_converter} {input.snp_plasma} | bgzip -c > {output.snp_plasma}
        
        python {params.vs_format_converter} {input.snp_tumor}  | bgzip -c > {output.snp_tumor}
        
        python {params.vs_format_converter} {input.indel_plasma} | bgzip -c > {output.indel_plasma}
        
        python {params.vs_format_converter} {input.indel_tumor}  | bgzip -c > {output.indel_tumor}
        
        bcftools concat -O z -o {output.plasma} {output.snp_plasma} {output.indel_plasma}
        
        bcftools concat -O z -o {output.tumor} {output.snp_tumor} {output.indel_tumor}
        
        """

########
# LoFreq
########

rule LoFreq:
    input:
        unpack(get_bam),
        ref=f'{genome_data}/GRCh38.d1.vd1.fa',
        known=config['dbsnp']
    params:
        plasma=f'{derived}/variant_calling/LoFreq/{{sample_calling}}/{{sample_calling}}-plasma.lofreq_',
        tumor=f'{derived}/variant_calling/LoFreq/{{sample_calling}}/{{sample_calling}}-tumor.lofreq_',
    output:
        snp_plasma=f'{derived}/variant_calling/LoFreq/{{sample_calling}}/{{sample_calling}}-plasma.lofreq_somatic_final_minus-dbsnp.snvs.vcf.gz',
        snp_tumor=f'{derived}/variant_calling/LoFreq/{{sample_calling}}/{{sample_calling}}-tumor.lofreq_somatic_final_minus-dbsnp.snvs.vcf.gz',
        indel_plasma=f'{derived}/variant_calling/LoFreq/{{sample_calling}}/{{sample_calling}}-plasma.lofreq_somatic_final_minus-dbsnp.indels.vcf.gz',
        indel_tumor=f'{derived}/variant_calling/LoFreq/{{sample_calling}}/{{sample_calling}}-tumor.lofreq_somatic_final_minus-dbsnp.indels.vcf.gz',
        plasma=f'{final}/variant_calling/plasma/{{sample_calling}}-plasma.LoFreq.final.vcf.gz',
        tumor=f'{final}/variant_calling/tumor/{{sample_calling}}-tumor.LoFreq.final.vcf.gz'
    conda:
        "../envs/lofreq.yaml"
    threads: config["ncores"]
    priority: 40
    shell:
        """
        lofreq somatic -n {input.control} -t {input.plasma} -f {input.ref} \
        --threads {threads} -o {params.plasma} -d {input.known}
        
        lofreq somatic -n {input.control} -t {input.tumor} -f {input.ref} \
        --threads {threads} -o {params.tumor} -d {input.known}
        
        bcftools concat -O z -o {output.plasma} {output.snp_plasma} {output.indel_plasma}
        
        bcftools concat -O z -o {output.tumor} {output.snp_tumor} {output.indel_tumor}
        
        """

###########
# Freebayes
###########

rule freebayes:
    input:
        unpack(get_bam),
        ref=f'{genome_data}/GRCh38.d1.vd1.fa',
        known=config['dbsnp']
    params:
        freebayes_somatic = "scripts/freebayes_somatic.sh",
        tmp_plasma=f'{derived}/variant_calling/Freebayes/{{sample_calling}}/tmp_plasma',
        tmp_tumor=f'{derived}/variant_calling/Freebayes/{{sample_calling}}/tmp_tumor',
    output:
        plasma=f'{derived}/variant_calling/Freebayes/{{sample_calling}}/{{sample_calling}}-plasma.Freebayes.filtered.vcf.gz',
        tumor=f'{derived}/variant_calling/Freebayes/{{sample_calling}}/{{sample_calling}}-tumor.Freebayes.filtered.vcf.gz',
        plasma_pass=f'{final}/variant_calling/plasma/{{sample_calling}}-plasma.Freebayes.final.vcf.gz',
        tumor_pass=f'{final}/variant_calling/tumor/{{sample_calling}}-tumor.Freebayes.final.vcf.gz'
    conda:
        "../envs/fb_somatic.yaml"
    threads: config["ncores"]
    priority: 40
    shell:
        """
        mkdir {params.tmp_plasma}
        mkdir {params.tmp_tumor}
        
        {params.freebayes_somatic}  \
        -t {input.tumor} -n {input.control} \
        -f {input.ref} -d {params.tmp_tumor} -P {input.known} \
        -c {threads} -v {output.tumor}
        
        {params.freebayes_somatic}  \
        -t {input.plasma} -n {input.control} \
        -f {input.ref} -d {params.tmp_plasma} -P {input.known} \
        -c {threads} -v {output.plasma}
        
        bcftools view -f PASS {output.tumor} | bgzip -c > {output.tumor_pass}
        
        bcftools view -f PASS {output.plasma} | bgzip -c > {output.plasma_pass}
        """
