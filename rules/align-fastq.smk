if not os.path.exists(genome_dict):
    rule dict:
        input:
            genome=genome
        output:
            dict=genome_dict
        log: f'{logs}/genome/genome_dict.log'
        benchmark: f'{benchmarks}/genome/genome_dict.txt'
        conda: "../envs/qc/picard.yaml"
        group: "preprocess_align"
        threads: 1
        shell:
            """
            picard CreateSequenceDictionary \
            -R {input.genome} \
            -O {output.dict}
            """

rule trimmomatic_se:
    input:
        unpack(get_fastq)
    output:
        temp(f'{derived}/pre_processing/temp/{{id}}.fastq.gz')
    log: f'{logs}/pre_processing/{{id}}/trimmomatic_se.log'
    benchmark: f'{benchmarks}/{{id}}_trimmomatic_se.txt'
    params: **config["params"]["trimmomatic"]["se"]
    threads: config["ncores"]
    group: "preprocess_align"
    priority: 40
    resources:
        mem_mb=4096
    wrapper:
        "v5.1.0/bio/trimmomatic/se"

rule trimmomatic_pe:
    input:
        unpack(get_fastq)
    output:
        r1=f'{derived}/pre_processing/temp/{{id}}.R1.fastq.gz',
        r2=f'{derived}/pre_processing/temp/{{id}}.R2.fastq.gz',
        r1_unpaired=temp(f'{derived}/pre_processing/unpaired/{{id}}.R1.fastq.gz'),
        r2_unpaired=temp(f'{derived}/pre_processing/unpaired/{{id}}.R2.fastq.gz'),
    log: f'{logs}/pre_processing/{{id}}/trimmomatic_pe.log'
    benchmark: f'{benchmarks}/{{id}}_trimmomatic_pe.txt'
    params: **config["params"]["trimmomatic"]["pe"],
    threads: config["ncores"]
    group: "preprocess_align"
    priority: 40
    resources:
        mem_mb=4096
    wrapper:
        "v5.1.0/bio/trimmomatic/pe"

rule fastx_trimmer_se:
    input:
        unpack(get_trimmed)
    output:
        f'{derived}/pre_processing/{{id}}.fastq.gz'
    conda: "../envs/align-fastq/fastx.yaml"
    log: f'{logs}/pre_processing/{{id}}/fastx_trimmer.se.log'
    benchmark: f'{benchmarks}/{{id}}_fastx_trimmer.txt'
    threads: 1
    group: "preprocess_align"
    priority: 40
    shell:
        """
        fastx_trimmer -z -l 145 -i {input} -o {output}

        """

rule fastx_trimmer_pe:
    input:
        unpack(get_trimmed)
    output:
        r1=f'{derived}/pre_processing/{{id}}.R1.fastq.gz',
        r2=f'{derived}/pre_processing/{{id}}.R2.fastq.gz',
    conda: "../envs/align-fastq/fastx.yaml"
    log: f'{logs}/pre_processing/{{id}}/fastx_trimmer.pe.log'
    benchmark: f'{benchmarks}/{{id}}_fastx_trimmer.txt'
    threads: 1
    group: "preprocess_align"
    priority: 40
    shell:
        """
        zcat {input.r1} | fastx_trimmer -z -l 145 -o {output.r1}
        zcat {input.r2} | fastx_trimmer -z -l 145 -o {output.r2}

        """

rule match_list_se:
    input:
        expand(f'{derived}/pre_processing/{{id}}.fastq.gz',id=samples["id"])
    output:
        f'{derived}/pre_processing/samples.txt'
    threads: 1
    priority: 40
    shell:
        """
        ls {input} > {output}
        """

rule match_list_pe:
    input:
        r1=expand(f'{derived}/pre_processing/{{id}}.R1.fastq.gz',id=samples["id"]),
        r2=expand(f'{derived}/pre_processing/{{id}}.R2.fastq.gz',id=samples["id"])
    output:
        r1=f'{derived}/pre_processing/samples_R1.txt',
        r2=f'{derived}/pre_processing/samples_R2.txt'
    threads: 1
    priority: 40
    shell:
        """
        ls {input.r1} > {output.r1}
        ls {input.r2} > {output.r2}
        """

rule cat_se:
        input:
            f'{derived}/pre_processing/samples.txt'
        output:
            f'{derived}/pre_processing/final/{{sample_type}}-single.fastq.gz'
        params:
            sample=f'{derived}/pre_processing/{{sample_type}}'
        threads: 1
        priority: 40
        shell:
            """
            cat {params.sample}-*-single.fastq.gz  > {output}   
            """

rule cat_pe:
        input:
            r1=f'{derived}/pre_processing/samples_R1.txt',
            r2=f'{derived}/pre_processing/samples_R2.txt'
        output:
            r1= f'{derived}/pre_processing/final/{{sample_type}}-paired.R1.fastq.gz',
            r2= f'{derived}/pre_processing/final/{{sample_type}}-paired.R2.fastq.gz'
        params:
            sample=f'{derived}/pre_processing/{{sample_type}}'
        threads: 1
        priority: 40
        shell:
            """
            cat {params.sample}-*-paired.R1.fastq.gz  > {output.r1}
            cat {params.sample}-*-paired.R2.fastq.gz  > {output.r2}   
            """

rule bwa_mem:
        input:
            reads=get_trimmed_reads,
            genome=genome
        output:
            bam=temp(f'{derived}/alignments/{{sample_type}}.bam')
        conda: "../envs/align-fastq/bwa.yaml"
        log: f'{logs}/alignments/{{sample_type}}/bwa_mem_alignment.log'
        benchmark: f'{benchmarks}/{{sample_type}}_bwa_mem_alignment.txt'
        params:
            extra=get_read_group,
            sort_extra=""
        threads: config['ncores']
        priority: 40
        shell:
            """
            bwa-mem2 mem -t {threads} \
            {params.extra} \
            {input.genome} \
            {input.reads} \
            | samblaster  | sambamba view -S -f bam /dev/stdin -t {threads} \
            | sambamba sort /dev/stdin -t {threads} \
            -o {output.bam} \
            2> {log}   
            """

if config['variant_calling']:
    rule realignertargetcreator:
        input:
            bam=f'{derived}/alignments/{{sample_type}}.bam',
            ref=genome,
            fai=genome_fai,
            dict=genome_dict,
            known=config['dbsnp'],
            known_idx=config['dbsnp_tbi']
        output:
            intervals=f'{derived}/alignments/{{sample_type}}.intervals',
        log: f'{logs}/recal/{{sample_type}}/GATK_realignertargetcreator.log'
        benchmark: f'{benchmarks}/{{sample_type}}_04_GATK_realignertargetcreator.txt'
        threads: config['ncores']
        priority: 40
        resources:
            mem_mb=4096,
        wrapper:
            "v5.1.0/bio/gatk3/realignertargetcreator"

    rule indelrealigner:
        input:
            bam=f'{derived}/alignments/{{sample_type}}.bam',
            ref=genome,
            fai=genome_fai,
            dict=genome_dict,
            known=config['dbsnp'],
            known_idx=config['dbsnp_tbi'],
            target_intervals=f'{derived}/alignments/{{sample_type}}.intervals'
        output:
            bam=temp(f'{derived}/alignments/{{sample_type}}.realigned.bam'),
            bai=temp(f'{derived}/alignments/{{sample_type}}.realigned.bai'),
        log: f'{logs}/recal/{{sample_type}}/GATK_indelrealigner.log'
        benchmark: f'{benchmarks}/{{sample_type}}_05_GATK_indelrealigner.txt'
        params:
            extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
        threads: config['ncores']
        priority: 40
        resources:
            mem_mb=4096,
        wrapper:
            "v5.1.0/bio/gatk3/indelrealigner"

    rule gatk_baserecalibrator:
        input:
            bam=f'{derived}/alignments/{{sample_type}}.realigned.bam',
            ref=genome,
            dict=genome_dict,
            known=config['dbsnp']
        output:
            recal_table=f'{derived}/alignments/{{sample_type}}.recal_data_table',
        log: f'{logs}/recal/{{sample_type}}/GATK_baserecalibrator.log'
        benchmark: f'{benchmarks}/{{sample_type}}_06_GATK_baserecalibrator.txt'
        params:
            extra="",  # optional
            java_opts="",  # optional
        priority: 40
        resources:
            mem_mb=4096,
        wrapper:
            "v5.1.0/bio/gatk/baserecalibrator"

    rule gatk_applybqsr:
        input:
            bam=f'{derived}/alignments/{{sample_type}}.realigned.bam',
            ref=genome,
            dict=genome_dict,
            recal_table=f'{derived}/alignments/{{sample_type}}.recal_data_table'
        output:
            bam=f'{derived}/recal/{{sample_type}}.bam'
        log: f'{logs}/recal/{{sample_type}}/GATK_APPLYBQSR.log'
        benchmark: f'{benchmarks}/{{sample_type}}_07_GATK_APPLYBQSR.txt'
        params:
            extra="",  # optional
            java_opts="",  # optional
            embed_ref=True,  # embed the reference in cram output
        priority: 40
        resources:
            mem_mb=4096,
        wrapper:
            "v5.1.0/bio/gatk/applybqsr"
else:
    rule no_bsqr:
        input:
            bam=f'{derived}/alignments/{{sample_type}}.bam'
        output:
            bam=f'{derived}/recal/{{sample_type}}.bam'
        threads: 1
        priority: 40
        shell:
            """
            cp {input.bam} {output.bam}
            """

rule samtools_index:
    input:
        bam=f'{derived}/recal/{{sample_type}}.bam'
    output:
        bai=f'{derived}/recal/{{sample_type}}.bam.bai'
    log: f'{logs}/recal/{{sample_type}}/samtools_index.log'
    benchmark: f'{benchmarks}/{{sample_type}}.samtools_index.txt'
    conda: "../envs/qc/samtools.yaml"
    threads: 1
    priority: 40
    shell:
        """

        samtools index {input.bam}

        """
