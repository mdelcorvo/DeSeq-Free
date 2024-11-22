###### Rules for somatic variant analysis ######
#
#     * Calling somatic variants with Varscan2
#     1. Generating mpileup.
#     2. VarScan2
#     3. Filtering somatic calls from Varscan2
#     4. Extraction of shared variants (i.e. variants shared between plasma and tumour)
#     5. Extraction of unique variants (i.e. variants not shared between plasma and tumour)
#
##########################################################################################

##########
# Varscan2
##########

rule mpileup:
    input:
        bam=f'{derived}/recal/{{sample_type}}.bam',
        ref=genome,
    output:
        pileup=f'{derived}/variant_calling/Varscan2/pileup/{{sample_type}}.pileup'
    conda: "../envs/qc/samtools.yaml"
    log: f'{logs}/variant_calling/{{sample_type}}/mpileup.log'
    benchmark: f'{benchmarks}/{{sample_type}}_mpileup.txt'
    threads: 1
    group: "variant_analysis"
    priority: 40
    shell:
        """
        samtools mpileup \
        -q 2 -f {input.ref} \
        {input.bam} > {output.pileup} 2> {log}

        """

rule VarScan2:
    input:
        unpack(get_pileup)
    params:
        plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma',
        tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor',
    output:
        snp_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.snp',
        snp_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.snp',
        indel_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.indel',
        indel_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.indel'
    conda: "../envs/variant_calling/varscan.yaml"
    log: f'{logs}/variant_calling/{{sample_calling}}/VarScan2.log'
    benchmark: f'{benchmarks}/{{sample_calling}}_VarScan2.txt'
    threads: 1
    priority: 40
    shell:
        """
        varscan somatic  {input.control}  {input.plasma} {params.plasma} --min-var-freq 0.01
        
        varscan somatic  {input.control}  {input.tumor} {params.tumor} --min-var-freq 0.01

        """

rule VarScan2_process:
    input:
        snp_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.snp',
        snp_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.snp',
        indel_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.indel',
        indel_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.indel'
    output:
        snp_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.snp.Somatic.hc',
        snp_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.snp.Somatic.hc',
        indel_plasma=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.indel.Somatic.hc',
        indel_tumor=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.indel.Somatic.hc',
        snp_plasma_list=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.snp.Somatic.hc.list',
        snp_tumor_list=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.snp.Somatic.hc.list',
        indel_plasma_list=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-plasma.indel.Somatic.hc.list',
        indel_tumor_list=f'{derived}/variant_calling/Varscan2/raw/{{sample_calling}}/{{sample_calling}}-tumor.indel.Somatic.hc.list'
    conda: "../envs/variant_calling/varscan.yaml"
    log: f'{logs}/variant_calling/{{sample_calling}}/VarScan2_process.log'
    benchmark: f'{benchmarks}/{{sample_calling}}_VarScan2_process.txt'
    threads: 1
    priority: 40
    shell:
        """
        varscan processSomatic {input.snp_plasma} --min-tumor-freq 0.01 --max-normal-freq 0.00
        varscan processSomatic {input.snp_tumor} --min-tumor-freq 0.01 --max-normal-freq 0.00
        varscan processSomatic {input.indel_plasma} --min-tumor-freq 0.01 --max-normal-freq 0.00
        varscan processSomatic {input.indel_tumor} --min-tumor-freq 0.01 --max-normal-freq 0.00
        
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
        ref=genome
    output:
        snp_rc_plasma=f'{derived}/variant_calling/Varscan2/bam_readcount/{{sample_calling}}-plasma.somatic_hc.snp.readcount',
        snp_rc_tumor=f'{derived}/variant_calling/Varscan2/bam_readcount/{{sample_calling}}-tumor.somatic_hc.snp.readcount',
        indel_rc_plasma=f'{derived}/variant_calling/Varscan2/bam_readcount/{{sample_calling}}-plasma.somatic_hc.indel.readcount',
        indel_rc_tumor=f'{derived}/variant_calling/Varscan2/bam_readcount/{{sample_calling}}-tumor.somatic_hc.indel.readcount'
    conda: "../envs/variant_calling/bam_read-count.yaml"
    log: f'{logs}/variant_calling/{{sample_calling}}/bam_readcount.log'
    benchmark: f'{benchmarks}/{{sample_calling}}_bam_readcount.txt'
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
    output:
        snp_plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.snp.fpfiltered',
        snp_tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.snp.fpfiltered',
        indel_plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.indel.fpfiltered',
        indel_tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.indel.fpfiltered',
        snp_plasma_edit=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.snp.fpfiltered.edit',
        snp_tumor_edit=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.snp.fpfiltered.edit'
    conda: "../envs/variant_calling/varscan.yaml"
    log: f'{logs}/variant_calling/{{sample_calling}}/VarScan2_filtering.log'
    benchmark: f'{benchmarks}/{{sample_calling}}_VarScan2_filtering.txt'
    threads: 1
    priority: 40
    shell:
        """
        varscan fpfilter {input.snp_plasma} {input.snp_rc_plasma} --output-file {output.snp_plasma} --min-var-freq 0.01 --min-strandedness 0.0
        varscan fpfilter {input.snp_tumor} {input.snp_rc_tumor} --output-file {output.snp_tumor} --min-var-freq 0.01 --min-strandedness 0.0
        
        varscan fpfilter {input.indel_plasma} {input.indel_rc_plasma} --output-file {output.indel_plasma} --min-var-freq 0.01 --min-strandedness 0.0
        varscan fpfilter {input.indel_tumor} {input.indel_rc_tumor} --output-file {output.indel_tumor} --min-var-freq 0.01 --min-strandedness 0.0
        
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
    log: f'{logs}/variant_calling/{{sample_calling}}/VarScan2_shared_variants.log'
    benchmark: f'{benchmarks}/{{sample_calling}}_VarScan2_shared_variants.txt'
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
        vs_format_converter="scripts/variant_calling/vs_format_converter.py"
    output:
        snp_plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.snp.fpfiltered.vcf.gz',
        snp_tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.snp.fpfiltered.vcf.gz',
        indel_plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.indel.fpfiltered.vcf.gz',
        indel_tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.indel.fpfiltered.vcf.gz',
        plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.Varscan2.final.vcf.gz',
        tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.Varscan2.final.vcf.gz'
    conda: "../envs/variant_calling/tabix.yaml"
    log: f'{logs}/variant_calling/{{sample_calling}}/Varscan2_snp2vcf.log'
    benchmark: f'{benchmarks}/{{sample_calling}}_Varscan2_snp2vcf.txt'
    threads: 1
    priority: 40
    shell:
        """
        python {params.vs_format_converter} {input.snp_plasma} | bgzip -c > {output.snp_plasma}
        
        python {params.vs_format_converter} {input.snp_tumor}  | bgzip -c > {output.snp_tumor}
        
        python {params.vs_format_converter} {input.indel_plasma} | bgzip -c > {output.indel_plasma}
        
        python {params.vs_format_converter} {input.indel_tumor}  | bgzip -c > {output.indel_tumor}
        
        tabix -p vcf {output.snp_plasma}
        tabix -p vcf {output.indel_plasma}
        
        tabix -p vcf {output.snp_tumor}
        tabix -p vcf {output.indel_tumor}
        
        vcf-concat {output.snp_plasma} {output.indel_plasma} | bgzip -c > {output.plasma}
        
        vcf-concat {output.snp_tumor} {output.indel_tumor} | bgzip -c > {output.tumor}
        
        """

rule comparison_Varscan2:
    input:
        plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.Varscan2.final.vcf.gz',
        tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.Varscan2.final.vcf.gz'
    output:
        tabulate=f'{derived}/variant_calling/Varscan2/{{sample_calling}}/{{sample_calling}}.Varscan2.snp.table.txt',
        summary=f'{derived}/variant_calling/Varscan2/{{sample_calling}}/{{sample_calling}}.Varscan2.snp.summary.txt',
    conda: "../envs/variant_calling/vcftoolz.yaml"
    log: f'{logs}/variant_calling/{{sample_calling}}/comparison_Varscan2.log'
    benchmark: f'{benchmarks}/{{sample_calling}}_comparison_Varscan2.txt'
    threads: 1
    group: "variant_analysis"
    shell:
        """
        vcftoolz compare \
        {input.plasma} \
        {input.tumor} \
        --tabulate {output.tabulate} > {output.summary}

        """

rule gather_calling_results:
    input:
        plasma=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-plasma.Varscan2.final.vcf.gz',
        tumor=f'{derived}/variant_calling/Varscan2/filtered/{{sample_calling}}/{{sample_calling}}-tumor.Varscan2.final.vcf.gz',
        tabulate=f'{derived}/variant_calling/Varscan2/{{sample_calling}}/{{sample_calling}}.Varscan2.snp.table.txt',
        summary=f'{derived}/variant_calling/Varscan2/{{sample_calling}}/{{sample_calling}}.Varscan2.snp.summary.txt'
    output:
        plasma=f'{final}/variant_calling/{{sample_calling}}/{{sample_calling}}-plasma.vcf.gz',
        tumor=f'{final}/variant_calling/{{sample_calling}}/{{sample_calling}}-tumor.vcf.gz',
        tabulate=f'{final}/variant_calling/{{sample_calling}}/{{sample_calling}}.plasma-tumor.snp.table.txt',
        summary=f'{final}/variant_calling/{{sample_calling}}/{{sample_calling}}.plasma-tumor.summary.txt',
    shell:
        """
        cp {input.plasma} {output.plasma}
        cp {input.tumor} {output.tumor}
        cp {input.tabulate} {output.tabulate}
        cp {input.summary} {output.summary}
        """