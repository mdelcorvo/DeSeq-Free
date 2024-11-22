###### Rules for somatic CNAs analysis ######
#     There are 2 main steps in the CNAs workflow:
#     1. Generating read count coverage information using readCounter from the HMMcopy suite.
#     2. Copy number analysis and prediction of tumor fraction using ichorCNA R package
#
## 1. Readcounts for IchorCNA analysis are generated with HMMcopy readcounter option (https://github.com/shahcompbio/hmmcopy_utils).
#     parameters recommended: read counts generation: 1 Mb window, reads with mapping quality greater or equal to 20
## 2. CNAs in both tumour and plasma samples are assessed using IchorCNA (https://github.com/broadinstitute/ichorCNA)
#     For plasma samples, tumour content is expected to be low,
#     therefore estimated tumour fractions of 5, 1, 0.5 and 0.1% are suggested
#     and ploidy set to diploid.
#     For tumour samples estimated tumour fractions of 50, 60, 70, 80, 90% are suggested
#     and ploidy set to 2 and 3
##########################################################################################

rule readCounter:
    input:
        bam=f'{derived}/recal/{{sample_type}}.bam',
        bai=f'{derived}/recal/{{sample_type}}.bam.bai'
    params:
        window= binSize,
        qual= qual,
        chrs=chrs
    output:
        wig=f'{derived}/cna/wig/{{sample_type}}.wig'
    conda: "../envs/cna/hmmcopy.yaml"
    log: f'{logs}/cna/{{sample_type}}/readCounter.log'
    benchmark: f'{benchmarks}/{{sample_type}}_readCounter.txt'
    threads: 1
    group: "cna_analysis"
    priority: 40
    shell:
        """
        readCounter --window {params.window} --quality {params.qual} \
        --chromosome {params.chrs} \
        {input.bam} > {output.wig} \
        2> {log}

        """

rule normal_list:
    input:
        expand(f'{derived}/cna/wig/{{sample_type}}.wig',sample_type=s1[s1["type"]=='control']['sample_type'])
    output:
        f'{derived}/cna/normal_samples.txt'
    threads: 1
    priority: 40
    shell:
        """
        ls {input} > {output}
        """

rule IchorCNA_normal_panel:
    input:
        f'{derived}/cna/normal_samples.txt'
    params:
        ichor="scripts/cna/normal_panels.R",
        gcwig=config["ichorCNA_gcWig"],
        mapwig=config["ichorCNA_mapWig"],
        centromere=config["ichorCNA_centromere"],
        normal_sample=f'{derived}/cna/normal_samples'
    output:
        f'{derived}/cna/normal_samples_median.rds'
    conda: "../envs/cna/ichorCNA.yaml"
    log: f'{logs}/cna/IchorCNA_normal_panel_creation.log'
    benchmark: f'{benchmarks}/IchorCNA_normal_panel_creation.txt'
    threads: config['ncores']
    group: "cna_analysis"
    priority: 40
    shell:
        """
        Rscript --vanilla {params.ichor} \
        {input} \
        {params.gcwig} {params.mapwig} \
        {params.centromere} \
        {params.normal_sample} \
        {log}
        """

rule IchorCNA_plasma:
    input:
        unpack(get_plasma_wig),
        normalpanel=f'{derived}/cna/normal_samples_median.rds'
    params:
        ichor="scripts/cna/IchorCNA.R",
        ichor_par=get_IchorCNA,
        outDir=f'{derived}/cna/seg/',
        id="{sample_calling}-plasma"
    output:
        cna=f'{derived}/cna/seg/{{sample_calling}}-plasma.cna.seg'
    conda: "../envs/cna/ichorCNA.yaml"
    log: f'{logs}/cna/{{sample_calling}}-plasma/IchorCNA_run_matched_plasma.log'
    benchmark: f'{benchmarks}/{{sample_calling}}-plasma/02_IchorCNA_run_matched_plasma.txt'
    threads: config['ncores']
    priority: 40
    shell:
        """
        Rscript --vanilla {params.ichor} \
        {params.id} \
        {input.plasma} \
        \"c(2)\" \
        \"c(0.95, 0.99, 0.995, 0.999)\" \
        3 \
        plasma \
        \"c()\" \
        {input.normalpanel} \
        {params.outDir} \
        {threads} \
        {params.ichor_par} \
        2> {log} 
        """

rule IchorCNA_tumor:
    input:
        unpack(get_tumor_wig),
        normalpanel=f'{derived}/cna/normal_samples_median.rds'
    params:
        ichor="scripts/cna/IchorCNA.R",
        ichor_par=get_IchorCNA,
        outDir=f'{derived}/cna/seg/',
        id="{sample_calling}-tumor"
    output:
        cna=f'{derived}/cna/seg/{{sample_calling}}-tumor.cna.seg'
    conda: "../envs/cna/ichorCNA.yaml"
    log: f'{logs}/cna/{{sample_calling}}-tumor/IchorCNA_run_matched_tumor.log'
    benchmark: f'{benchmarks}/{{sample_calling}}-tumor/03_IchorCNA_run_matched_tumor.txt'
    threads: 1
    priority: 40
    shell:
        """
        Rscript --vanilla {params.ichor} \
        {params.id} \
        {input.tumor} \
        \"c(2,3)\"  \
        \"c(0.5,0.6,0.7,0.8,0.9)\" \
        5 \
        tumor \
        \"c(1,3)\" \
        {input.normalpanel} \
        {params.outDir} \
        {threads} \
        {params.ichor_par} \
        2> {log} 
        """

rule cna_plasma_tumor_comparison:
    input:
        plasma_cna=f'{derived}/cna/seg/{{sample_calling}}-plasma.cna.seg',
        tumor_cna=f'{derived}/cna/seg/{{sample_calling}}-tumor.cna.seg'
    params:
        compare="scripts/cna/IchorCNA_result_comp.R",
        annot=exon_transcripts_hg38
    output:
        shared_cna=f'{derived}/cna/results/{{sample_calling}}-plasma-tumor-shared.cna.txt'
    conda: "../envs/cna/ichorCNA.yaml"
    log: f'{logs}/cna/{{sample_calling}}/IchorCNA.plasma_tumor_comparison.log'
    benchmark: f'{benchmarks}/{{sample_calling}}-plasma-tumor/04_IchorCNA.plasma_tumor_comparison.txt'
    threads: 1
    priority: 40
    shell:
            """
            Rscript --vanilla {params.compare} \
            {input.plasma_cna} \
            {input.tumor_cna} \
            {params.annot} \
            {output.shared_cna} \
            {log} 

            """

rule list_results:
    input:
        expand(f'{derived}/cna/results/{{sample_calling}}-plasma-tumor-shared.cna.txt', sample_calling=s2["sample"])
    output:
        f'{derived}/cna/output_samples.txt'
    threads: 1
    priority: 40
    shell:
            """
            ls {input} > {output}
            """

rule cna_merge_results:
    input:
        f'{derived}/cna/output_samples.txt'
    output:
        f'{final}/cna/plasma-tumor-shared.cna.txt'
    params:
        merge="scripts/cna/merge_results.R"
    conda: "../envs/cna/ichorCNA.yaml"
    threads: 1
    priority: 40
    shell:
            """
            Rscript --vanilla {params.merge} \
            {input} \
            {output} \
            """