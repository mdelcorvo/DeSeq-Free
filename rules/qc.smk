rule fastp_report:
    input:
        unpack(fastq)
    output:
        html = f'{derived}/qc/fastp/{{sample_type}}.fastp.html',
        json = f'{derived}/qc/fastp/{{sample_type}}.fastp.json'
    log: f'{logs}/qc/{{sample_type}}/fastp-QC.log'
    benchmark: f'{benchmarks}/{{sample_type}}.fastp.benchmark.txt'
    group: "preprocess_align"
    conda:  "../envs/qc/fastp.yaml"
    threads:   4
    priority:  40
    shell:
        """
        fastp --thread {threads} -i {input.r1} -I {input.r2} \
        -h {output.html} -j {output.json} 2> {log}
        """

rule gather_fastp_stats:
    input:
        rep_json = f'{derived}/qc/fastp/{{sample_type}}.fastp.json'
    output:
        stats_json = f'{derived}/qc/fastp/{{sample_type}}.sequencing-qc.txt'
    params:
        encoding       = "Illumina Basecalls",
        sample_name    = "{sample_type}",
        read_filenames = fastq,
        file_type      = "FASTQ"
    log:      f'{logs}/qc/{{sample_type}}/fastq-QC.stats.log',
    group:   "preprocess_align"
    threads:  1
    priority: 40
    script:
        "../scripts/qc/extract-fastq-qc.py"

        # file_type = "Conventional base calls",
        # encoding = "Sanger/Illumina 1.9",

rule samtools_stats:
    input:
        bam=f'{derived}/recal/{{sample_type}}.bam',
    output:
        stats=f'{derived}/qc/samtools/{{sample_type}}.stats',
        edit_stats=f'{derived}/qc/samtools/{{sample_type}}.edit.stats'
    log: f'{logs}/qc/{{sample_type}}/samtools_stats.log'
    benchmark: f'{benchmarks}/{{sample_type}}.samtools_stats.txt'
    conda: "../envs/qc/samtools.yaml"
    threads: max(1,int(config["ncores"] // 3))
    shell:
        """
        samtools stats \
        -@ {threads} \
        {input.bam} \
        > {output.stats} 2> {log}
        grep '^SN' {output.stats} | sed 's/ \{{2,\}}/\t/g' | sed 's/ /./g' | awk '{{ $1=""; $4=""; print $0 }}' > {output.edit_stats} 2>> {log}
        """

rule samtools_flagstat:
    input:
        bam=f'{derived}/recal/{{sample_type}}.bam'
    output:
        flagstats=f'{derived}/qc/samtools/{{sample_type}}.flagstats',
        edit_flagstats=f'{derived}/qc/samtools/{{sample_type}}.edit.flagstats'
    log: f'{logs}/qc/{{sample_type}}/samtools_flagstat.log'
    benchmark: f'{benchmarks}/{{sample_type}}.samtools_flagstat.txt'
    conda: "../envs/qc/samtools.yaml"
    threads: max(1,int(config["ncores"] // 3))
    shell:
        """
        samtools flagstat \
        -@ {threads} \
        {input.bam} \
        > {output.flagstats} 2> {log}

        awk '{{ print $1 }}' {output.flagstats} > {output.edit_flagstats}
        """

rule mosdepth:
    input:
        bam=f'{derived}/recal/{{sample_type}}.bam',
        bai=f'{derived}/recal/{{sample_type}}.bam.bai'
    output:
        dist=f'{derived}/qc/mosdepth/{{sample_type}}.mosdepth.global.dist.txt'
    params:
        sample=f'{derived}/qc/mosdepth/{{sample_type}}',
    conda: "../envs/qc/mosdepth.yaml"
    log: f'{logs}/qc/{{sample_type}}/mosdepth.log'
    benchmark: f'{benchmarks}/{{sample_type}}.mosdepth.txt'
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

rule mosdepth_plot:
    input:
        dist=expand(f'{derived}/qc/mosdepth/{{sample_type}}.mosdepth.global.dist.txt',sample_type=s1['sample_type'])
    output:
        dist_single=f'{derived}/qc/mosdepth/mosdepth.dist.tiff',
        dist_grouped=f'{derived}/qc/mosdepth/mosdepth.dist.grouped.tiff',
    params:
        plot_script="scripts/qc/coverage_plot.R",
        sample_dir=f'{derived}/qc/mosdepth/'
    conda: "../envs/qc/cfdnapro.yaml"
    log: f'{logs}/qc/mosdepth_plot.log'
    benchmark: f'{benchmarks}/mosdepth_plot.txt'
    threads: 1
    shell:
        """
        Rscript --vanilla {params.plot_script} \
        {params.sample_dir} \
        {output.dist_single} \
        {output.dist_grouped} \
        {log}
        """

rule collect_insert_size_metrics:
    input:
        unpack(get_bam)
    output:
        metrics_plasma=f'{derived}/qc/picard/plasma/{{sample_calling}}.plasma.insert_size_metrics.txt',
        metrics_tumor=f'{derived}/qc/picard/tumor/{{sample_calling}}.tumor.insert_size_metrics.txt',
        metrics_control=f'{derived}/qc/picard/control/{{sample_calling}}.control.insert_size_metrics.txt',
        histo_plasma=f'{derived}/qc/picard/plasma/{{sample_calling}}.plasma.insert_size_histogram.pdf',
        histo_tumor=f'{derived}/qc/picard/tumor/{{sample_calling}}.tumor.insert_size_histogram.pdf',
        histo_control=f'{derived}/qc/picard/control/{{sample_calling}}.control.insert_size_histogram.pdf'
    conda: "../envs/qc/picard.yaml"
    log: f'{logs}/qc/{{sample_calling}}/picard.log'
    benchmark: f'{benchmarks}/{{sample_calling}}.picard.txt'
    threads: 1
    priority: 40
    shell:
        """
        picard CollectInsertSizeMetrics \
         -I {input.plasma} \
         -O {output.metrics_plasma} \
         -H {output.histo_plasma} \
         2> {log}

        picard CollectInsertSizeMetrics \
         -I {input.tumor} \
         -O {output.metrics_tumor} \
         -H {output.histo_tumor} \
         2>> {log}
         
         picard CollectInsertSizeMetrics \
         -I {input.control} \
         -O {output.metrics_control} \
         -H {output.histo_control} \
         2>> {log}

        """

rule cfDNAPro :
    input:
        metrics=expand(f'{derived}/qc/picard/{{experiment}}/{{sample_calling}}.{{experiment}}.insert_size_metrics.txt',sample_calling=s2['sample'],experiment=["plasma", "tumor", "control"]),
    output:
        fragment_size_dist=f'{derived}/qc/cfDNAPro/fragment_size_dist.tiff',
        median_size=f'{derived}/qc/cfDNAPro/median_size_dist.tiff'
    params:
        cfDNAPro="scripts/qc/cfDNAPro.R",
        sample_dir=f'{derived}/qc/picard/'
    conda: "../envs/qc/cfdnapro.yaml"
    log: f'{logs}/qc/cfDNAPro.log'
    benchmark: f'{benchmarks}/cfDNAPro.txt'
    threads: 1
    shell:
        """
        Rscript --vanilla {params.cfDNAPro} \
        {params.sample_dir} \
        {output.fragment_size_dist} \
        {output.median_size} \
        {log}
        """

rule multiqc:
    input:
        fragment_size_dist = f'{derived}/qc/cfDNAPro/fragment_size_dist.tiff',
        stats_json = expand(f'{derived}/qc/fastp/{{sample_type}}.sequencing-qc.txt',sample_type=s1["sample_type"]),
        stats=expand(f'{derived}/qc/samtools/{{sample_type}}.stats',sample_type=s1["sample_type"]),
        flagstats=expand(f'{derived}/qc/samtools/{{sample_type}}.flagstats',sample_type=s1["sample_type"]),
        fastp=expand(f'{derived}/qc/fastp/{{sample_type}}.fastp.json',sample_type=s1["sample_type"]),
        mosdepth=expand(f'{derived}/qc/mosdepth/{{sample_type}}.mosdepth.global.dist.txt',sample_type=s1["sample_type"])
    output:
        report(f'{final}/report/Processed_data_qc.html',caption="multiqc.rst",category="Quality control")
    params:
        outdir=f'{final}/report/',
        sample_proc='Processed_data_qc'
    log: f'{logs}/report/Multiqc_processed.log'
    benchmark: f'{benchmarks}/report/Multiqc_processed.txt'
    conda: "../envs/qc/multiqc.yaml"
    threads: 1
    shell:
        """
        multiqc --force \
        {input.stats} \
        {input.flagstats} \
        {input.fastp} \
        {input.mosdepth} \
        --filename {params.sample_proc} \
        -o {params.outdir} 2> {log}
        """

rule gather_qc:
    input:
        fragment_size_dist=f'{derived}/qc/cfDNAPro/fragment_size_dist.tiff',
        median_size=f'{derived}/qc/cfDNAPro/median_size_dist.tiff',
        dist_single=f'{derived}/qc/mosdepth/mosdepth.dist.tiff',
        dist_grouped=f'{derived}/qc/mosdepth/mosdepth.dist.grouped.tiff',
    output:
        plot_qc=f'{final}/qc/qc.pdf'
    params:
        sample_dir=f'{derived}/qc'
    log: f'{logs}/qc/QC_report.log'
    benchmark: f'{benchmarks}/qc/QC_report.txt'
    conda: "../envs/qc/qc_report.yaml"
    shell:
        """
        cd {params.sample_dir}
        magick -density 150 */*.tiff {output.plot_qc}

        """