###### Rules for somatic variant analysis ######
#     Somatic variant analysis with VarScan2 / LoFreq
#     # Varscan2
#     1. Generating mpileup.
#     2. VarScan2
#     3. Filtering somatic calls from Varscan2
#     4. Extraction of shared variants (i.e. variants shared between plasma and tumour)
#     5. Extraction of unique variants (i.e. variants not shared between plasma and tumour)
#
#     # LoFreq
##########################################################################################

rule mpileup:
    input:
        bam=f'{derived}/recal/{{sample_type}}.bam',
        ref=f'{genome_data}/GRCh38.d1.vd1.fa',
    output:
        pileup=f'{derived}/variant_calling/{{sample_type}}.pileup'
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
        plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma',
        tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor',
    output:
        snp_plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma.snp',
        snp_tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor.snp'
    threads: 1
    priority: 40
    shell:
        """
        java -Xmx16g -jar {params.varscan} somatic  {input.control}  {input.plasma} {params.plasma} --min-var-freq 0.01
        
        java -Xmx16g -jar {params.varscan} somatic  {input.control}  {input.tumor} {params.tumor} --min-var-freq 0.01
        
        """


rule VarScan2_process:
    input:
        snp_plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma.snp',
        snp_tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor.snp'
    params:
        varscan = "tools/VarScan.jar"
    output:
        snp_plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma.snp.Somatic.hc',
        snp_tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor.snp.Somatic.hc',
        snp_plasma_list=f'{derived}/variant_calling/{{sample_calling}}-plasma.snp.Somatic.hc.list',
        snp_tumor_list=f'{derived}/variant_calling/{{sample_calling}}-tumor.snp.Somatic.hc.list'
    threads: 1
    priority: 40
    shell:
        """
        
        java -Xmx16g -jar {params.varscan} processSomatic {input.snp_plasma} --min-tumor-freq 0.01 --max-normal-freq 0.00
        
        java -Xmx16g -jar {params.varscan} processSomatic {input.snp_tumor} --min-tumor-freq 0.01 --max-normal-freq 0.00
        
        awk 'BEGIN{{OFS="\t"}} NR!=1 {{print $1, $2, $2}}' {output.snp_plasma} > {output.snp_plasma_list}
        
        awk 'BEGIN{{OFS="\t"}} NR!=1 {{print $1, $2, $2}}' {output.snp_tumor} > {output.snp_tumor_list}
     
        """


rule bam_readcount:
    input:
        snp_plasma_list=f'{derived}/variant_calling/{{sample_calling}}-plasma.snp.Somatic.hc.list',
        snp_tumor_list=f'{derived}/variant_calling/{{sample_calling}}-tumor.snp.Somatic.hc.list',
        bam_plasma=f'{derived}/recal/{{sample_calling}}-plasma.bam',
        bam_tumor=f'{derived}/recal/{{sample_calling}}-tumor.bam',
        ref=f'{genome_data}/GRCh38.d1.vd1.fa'
    params:
        varscan="tools/VarScan.jar"
    conda:
        "../envs/bam_read-count.yaml"
    output:
        rc_plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma.somatic_hc.readcount',
        rc_tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor.somatic_hc.readcount'
    threads: 1
    priority: 40
    shell:
        """

        bam-readcount -q 2 -f {input.ref} \
        -l {input.snp_plasma_list} {input.bam_plasma} \
        > {output.rc_plasma}

        bam-readcount -q 2 -f {input.ref} \
        -l {input.snp_tumor_list} {input.bam_tumor} \
        > {output.rc_tumor}
        """

rule VarScan2_filtering:
    input:
        snp_plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma.snp.Somatic.hc',
        snp_tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor.snp.Somatic.hc',
        rc_plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma.somatic_hc.readcount',
        rc_tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor.somatic_hc.readcount'
    params:
        varscan="tools/VarScan.jar",
        plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma',
        tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor',
    output:
        plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma.fpfilter',
        tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor.fpfilter',
        snp_plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma.fpfilter.snp',
        snp_tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor.fpfilter.snp'
    threads: 1
    priority: 40
    shell:
        """

        java -jar VarScan.v2.4.4.jar fpfilter {input.snp_plasma} {input.rc_plasma} --output-file {params.plasma} --min-var-freq 0.01 --min-strandedness 0.0

        java -jar VarScan.v2.4.4.jar fpfilter {input.snp_tumor} {input.rc_tumor} --output-file {params.tumor} --min-var-freq 0.01 --min-strandedness 0.0
        
        awk '{{if($10 >= 5 && $6 == 0 && $5 + $6 >= 10 && $9 + $10 >= 10) {{print $0}} }}' {output.plasma} > {output.snp_plasma}

        awk '{{if($10 >= 5 && $6 == 0 && $5 + $6 >= 10 && $9 + $10 >= 10) {{print $0}} }}' {output.tumor} > {output.snp_tumor}
        """

rule shared_variants:
    input:
        snp_plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma.fpfilter.snp',
        snp_tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor.fpfilter.snp',
    output:
        shared=f'{derived}/variant_calling/{{sample_calling}}-plasma-tumor.shared.snp',
        unique=f'{derived}/variant_calling/{{sample_calling}}-plasma-tumor.unique.snp'
    threads: 1
    priority: 40
    shell:
        """

        awk 'FNR==NR{{a[$1,$2]=$0;next}}{{if(b=a[$1,$2]){{print b}} }}' \
        {input.snp_plasma} {input.snp_tumor} \
        > {output.shared}
        
        awk 'FNR==NR{{a[$1,$2]++}}FNR!=NR && !a[$1,$2]{{print}}' \
        {input.snp_tumor} {input.snp_plasma} \
        > {output.unique}
        """

rule lofreq:
    input:
        unpack(get_bam),
        ref=f'{genome_data}/GRCh38.d1.vd1.fa',
        known=config['dbsnp']
    params:
        plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma.lofreq_',
        tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor.lofreq_',
    output:
        snp_plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma.lofreq_somatic_final_minus-dbsnp.snvs.vcf.gz',
        snp_tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor.lofreq_somatic_final_minus-dbsnp.snvs.vcf.gz'
    threads: config["ncores"]
    priority: 40
    shell:
        """
        lofreq somatic -n {input.control} -t {input.plasma} -f {input.ref} \
        --threads {threads} -o {params.plasma} -d {input.known}
        
        lofreq somatic -n {input.control} -t {input.tumor} -f {input.ref} \
        --threads {threads} -o {params.tumor} -d {input.known}
        """