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
        rc = "tools/readCounter",
        window= binSize,
        qual= qual,
        chrs=chrs
    output:
        wig=f'{derived}/cna/wig/{{sample_type}}.wig'
    log: f'{logs}/{{sample_type}}/01_readCounter.log'
    benchmark: f'{benchmarks}/{{sample_type}}_readCounter.txt'
    threads: 1
    group: "cna_analysis"
    priority: 40
    shell:
        """
        {params.rc} --window {params.window} --quality {params.qual} \
        --chromosome {params.chrs} \
        {input.bam} > {output.wig} \
        2> {log}

        """

if config['ichorCNA_Nornal_Panel']:

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
            ichor="scripts/createPanelOfNormals.R",
            gcwig=config["ichorCNA_gcWig"],
            mapwig=config["ichorCNA_mapWig"],
            centromere=config["ichorCNA_centromere"],
            normal_sample=f'{derived}/cna/normal_samples'
        output:
            f'{derived}/cna/normal_samples_median.rds'
        log: f'{logs}/{{sample_type}}/02_IchorCNA_normal_panel_creation.log'
        benchmark: f'{benchmarks}/{{sample_type}}/02_IchorCNA_normal_panel_creation.txt'
        threads: 1
        group: "cna_analysis"
        priority: 40
        shell:
            """
            Rscript --vanilla {params.ichor} \
            --filelist {input} \
            --gcWig {params.gcwig} --mapWig {params.mapwig} \
            --centromere {params.centromere} \
            --outfile {params.normal_sample} > {log}
    
            """

    rule IchorCNA_with_normals:
        input:
            wig=f'{derived}/cna/wig/{{sample_type}}.wig',
            normalpanel=f'{derived}/cna/normal_samples_median.rds'
        params:
            ichor="scripts/runIchorCNA.R",
            ichor_par=get_IchorCNA,
            outDir=f'{derived}/cna/seg/',
            id="{sample_type}"
        output:
            cna=f'{derived}/cna/seg/{{sample_type}}.cna.seg'
        log: f'{logs}/{{sample_type}}/03_IchorCNA_run_with_normal_panel.log'
        benchmark: f'{benchmarks}/{{sample_type}}/03_IchorCNA_run_with_normal_panel.txt'
        threads: 1
        priority: 40
        run:
            if 'plasma' in wildcards.sample_type:
                cmd = "Rscript --vanilla {params.ichor} \
                --id {params.id} \
                --WIG {input.wig} \
                --ploidy \"c(2)\" \
                --normal \"c(0.95, 0.99, 0.995, 0.999)\" \
                --maxCN 3 \
                --estimateScPrevalence FALSE --scStates \"c()\" \
                --normalPanel {input.normalpanel} \
                --outDir {params.outDir} \
               {params.ichor_par} > {log} "

            elif 'tumor' in wildcards.sample_type:
                cmd = "Rscript --vanilla {params.ichor} \
                --id {params.id} \
                --WIG {input.wig} \
                --ploidy \"c(2,3)\"  \
                --normal \"c(0.5,0.6,0.7,0.8,0.9)\" \
                --maxCN 5 \
                --estimateScPrevalence TRUE --scStates \"c(1,3)\" \
                --normalPanel {input.normalpanel} \
                --outDir {params.outDir} \
                {params.ichor_par} > {log} "

            shell(cmd)

else:

    rule IchorCNA_plasma:
        input:
            unpack(get_wig)
        params:
            ichor="scripts/runIchorCNA.R",
            ichor_par=get_IchorCNA,
            outDir=f'{derived}/cna/seg/',
            id="{sample_calling}-plasma"
        output:
            cna=f'{derived}/cna/seg/{{sample_calling}}-plasma.cna.seg',
        log: f'{logs}/{{sample_calling}}-plasma/02_IchorCNA_run_matched_plasma.log'
        benchmark: f'{benchmarks}/{{sample_calling}}-plasma/02_IchorCNA_run_matched_plasma.txt'
        threads: 1
        priority: 40
        shell:
            """
            Rscript --vanilla {params.ichor} \
            --id {params.id} \
            --WIG {input.plasma} \
            --NORMWIG {input.control} \
            --ploidy \"c(2)\" \
            --normal \"c(0.95, 0.99, 0.995, 0.999)\" \
            --maxCN 4 \
            --estimateScPrevalence FALSE --scStates \"c()\" \
            --outDir {params.outDir} \
            {params.ichor_par} > {log}

            """

    rule IchorCNA_tumor:
        input:
            unpack(get_wig)
        params:
            ichor="scripts/runIchorCNA.R",
            ichor_par=get_IchorCNA,
            outDir=f'{derived}/cna/seg/',
            id="{sample_calling}-tumor"
        output:
            cna=f'{derived}/cna/seg/{{sample_calling}}-tumor.cna.seg'
        log: f'{logs}/{{sample_calling}}-tumor/03_IchorCNA_run_matched_tumor.log'
        benchmark: f'{benchmarks}/{{sample_calling}}-tumor/03_IchorCNA_run_matched_tumor.txt'
        threads: 1
        priority: 40
        shell:
            """
            Rscript --vanilla {params.ichor} \
            --id {params.id} \
            --WIG {input.tumor} \
            --NORMWIG {input.control} \
            --ploidy \"c(2,3)\"  \
            --normal \"c(0.5,0.6,0.7,0.8,0.9)\" \
            --maxCN 5 \
            --estimateScPrevalence TRUE --scStates \"c(1,3)\" \
            --outDir {params.outDir} \
            {params.ichor_par} > {log} 

            """

rule plasma_tumor_comparison:
    input:
        plasma_cna=f'{derived}/cna/seg/{{sample_calling}}-plasma.cna.seg',
        tumor_cna=f'{derived}/cna/seg/{{sample_calling}}-tumor.cna.seg'
    params:
        compare="scripts/IchorCNA_result_comp.R",
        annot=exon_transcripts_hg38
    output:
        shared_cna=f'{final}/cna/{{sample_calling}}-plasma-tumor-shared.cna.txt'
    log: f'{logs}/{{sample_calling}}/04_IchorCNA.plasma_tumor_comparison.log'
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