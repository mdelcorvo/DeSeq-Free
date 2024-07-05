rule fastp_report:
    input:
        unpack(fastq)
    output:
        html = f'{derived}/qc/{{sample_type}}.fastp.html',
        json = f'{derived}/qc/{{sample_type}}.fastp.json'
    log: f'{logs}/{{sample_type}}/01_fastp-QC.log'
    benchmark: f'{benchmarks}/{{sample_type}}.fastp.benchmark.txt'
    group: "preprocess_align"
    threads:   4
    priority:  40
    shell:
        """
        fastp --thread {threads} -i {input.r1} -I {input.r2} \
        -h {output.html} -j {output.json} 2> {log}
        """

rule qc_sequencing_fastq:
    input:
        rep_json = f'{derived}/qc/{{sample_type}}.fastp.json'
    output:
        report = f'{derived}/qc/{{sample_type}}.sequencing-qc.txt'
    params:
        encoding       = "Illumina Basecalls",
        sample_name    = "{sample_type}",
        read_filenames = fastq,
        file_type      = "FASTQ"
    log:      f'{logs}/{{sample_type}}/02_fastq-QC.log',
    group:   "preprocess_align"
    threads:  1
    priority: 40
    script:
        "../scripts/extract-fastq-qc.py"

        # file_type = "Conventional base calls",
        # encoding = "Sanger/Illumina 1.9",


rule mosdepth:
    input:
        bam=f'{derived}/recal/{{sample_type}}.bam',
        bai=f'{derived}/recal/{{sample_type}}.bam.bai'
    output:
        dist=f'{derived}/qc/coverage/{{sample_type}}.mosdepth.global.dist.txt'
    params:
        sample=f'{derived}/qc/coverage/{{sample_type}}'
    log:
        f'{logs}/{{sample_type}}/03_mosdepth.log',
    threads: config['ncores']
    shell:
        """
         mosdepth \
         --threads {threads} \
         --fast-mode \
         {params.sample} \
         {input.bam} \
         2> {log}
        """