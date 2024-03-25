###### Rules for somatic variant analysis ######
#     Somatic variant analysis with VarScan2
#     1. Generating mpileup.
#     2. VarScan2
#     3. Filtering somatic calls from Varscan2
#     4. Extraction of shared variants (i.e. variants shared between plasma and tumour)
#     5. Extraction of unique variants (i.e. variants not shared between plasma and tumour)
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

rule VarScan2_plasma:
    input:
        plasma=f'{derived}/variant_calling/{{sample_calling}}-plasma.pileup',
        control=f'{derived}/variant_calling/{{sample_calling}}-control.pileup'
    params:
        varscan = "tools/VarScan.jar"
    output:
        snp=f'{derived}/variant_calling/{{sample_calling}}-plasma.snp'
    threads: 1
    priority: 40
    shell:
        """
        java -Xmx16g -jar {params.varscan} somatic  {input.control}  {input.plasma} {output.snp} --min-var-freq 0.01
        """

rule VarScan2_tumor:
    input:
        tumor=f'{derived}/variant_calling/{{sample_calling}}-tumor.pileup',
        control=f'{derived}/variant_calling/{{sample_calling}}-control.pileup'
    params:
        varscan = "tools/VarScan.jar"
    output:
        snp=f'{derived}/variant_calling/{{sample_calling}}-tumor.snp'
    threads: 1
    priority: 40
    shell:
        """
        java -Xmx16g -jar {params.varscan} somatic  {input.control}  {input.tumor} {output.snp} --min-var-freq 0.01
        """