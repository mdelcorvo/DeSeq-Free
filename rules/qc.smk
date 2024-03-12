rule fastp_report:
    input:
        unpack(fastq)
    output:
        html = f'{derived}/qc_reports/fastp/{{sample}}/{{sample}}.fastp.html',
        json = f'{derived}/qc_reports/fastp/{{sample}}/{{sample}}.fastp.json'
    log: f'{logs}/{{sample}}/01_fastp-QC.log'
    benchmark: f'{benchmarks}/{{sample}}.fastp.benchmark.txt'
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
        rep_json = f'{derived}/qc_reports/fastp/{{sample}}/{{sample}}.fastp.json'
    output:
        report = f'{derived}/qc_reports/sequencing/{{sample}}/{{sample}}.sequencing-qc.txt'
    params:
        encoding       = "Illumina Basecalls",
        sample_name    = "{sample}",
        read_filenames = fastq,
        file_type      = "FASTQ"
    log:      f'{logs}/{{sample}}/02_fastq-QC.log',
    group:   "preprocess_align"
    threads:  1
    priority: 40
    script:
        "../scripts/extract-fastq-qc.py"

        # file_type = "Conventional base calls",
        # encoding = "Sanger/Illumina 1.9",

rule paired_tools_report:
    input:
        stats = f'{derived}/qc_reports/pairtools/{{sample}}/{{sample}}.stats'
    output:
        stats = f'{derived}/qc_reports/pairtools/{{sample}}/{{sample}}.paired-tools.stats.txt'
    log: f'{logs}/{{sample}}/03_paired_tools-QC.log',
    group:   "preprocess_align"
    threads:  1
    priority: 40
    script:
        "../scripts/get_qc.py"

if peaks:        
	rule enrichment_report:
		input:
			mapped = f'{derived}/alignments/pairtools/{{sample}}/{{sample}}.mapped.PT.bam',
			ref = f'{genome_data}/GRCh38.genome',
			peak = get_peak_calls,
			script = "scripts/enrichment_stats.sh" 
		output:
			stats = f'{derived}/qc_reports/enrichment_stats/{{sample}}/{{sample}}.CTCF_hichip_qc_metrics.txt'
		params:
			out_dir = f'{derived}/qc_reports/enrichment_stats/{{sample}}',
			sample =  '{sample}'
		log: f'{logs}/{{sample}}/04_enrichment-QC.log',
		group:   "preprocess_align"
		threads:  1
		priority: 40
		shell:
			"""
			{input.script} \
			-g {input.ref} \
			-b {input.mapped} \
			-p {input.peak} \
			-t {threads} \
			-x CTCF \
			-o {params.out_dir} \
			-s {params.sample}
			"""      

	rule enrichment_plot:
		input:
			mapped = f'{derived}/alignments/pairtools/{{sample}}/{{sample}}.mapped.PT.bam',
			peak = get_peak_calls,
			stats = f'{derived}/qc_reports/enrichment_stats/{{sample}}/{{sample}}.CTCF_hichip_qc_metrics.txt',
			script = "scripts/plot_chip_enrichment_bed.py"
		output:
			enrichment = f'{derived}/qc_reports/enrichment_stats/{{sample}}/{{sample}}.enrichment.png'
		log: f'{logs}/{{sample}}/05_enrichment-plot.log'
		group:   "preprocess_align"
		threads:  1
		priority: 40
		shell:
			"""
			
			{input.script} \
			-bam {input.mapped} \
			-peaks {input.peak} \
			-output {output.enrichment}
			
			"""      

rule mosdepth:
    input:
        bam=f'{derived}/alignments/pairtools/{{sample}}/{{sample}}.mapped.PT.bam'
    output:     
        dist = f'{derived}/qc_reports/mosdepth/{{sample}}/{{sample}}.mosdepth.region.dist.txt' if peaks else f'{derived}/qc_reports/mosdepth/{{sample}}/{{sample}}.mosdepth.global.dist.txt',
        perbase_cov = f'{derived}/qc_reports/mosdepth/{{sample}}/{{sample}}.per-base.bed.gz',
        summary = f'{derived}/qc_reports/mosdepth/{{sample}}/{{sample}}.mosdepth.summary.txt'       
    log: f'{logs}/{{sample}}/06_mosdepth.log'
    benchmark: f'{benchmarks}/{{sample}}.mosdepth.benchmark.txt'
    params:
        extra=get_bed,
        mosdepth_prefix = f'{derived}/qc_reports/mosdepth/{{sample}}/{{sample}}'
    threads: config["ncores"]
    shell: """	
		mosdepth -t {threads} \
        -m {params.extra} \
        {params.mosdepth_prefix}  \
        {input.bam}
		"""   
        
        
rule multiqc_report:
    input:
        json = expand(f'{derived}/qc_reports/fastp/{{sample}}/{{sample}}.fastp.json',sample=samples["sample"]),
        stats = expand(f'{derived}/qc_reports/pairtools/{{sample}}/{{sample}}.stats',sample=samples["sample"]),
        mosdepth = expand(f'{derived}/qc_reports/mosdepth/{{sample}}/{{sample}}.mosdepth.region.dist.txt',sample=samples["sample"]) if peaks else expand(f'{derived}/qc_reports/mosdepth/{{sample}}/{{sample}}.mosdepth.global.dist.txt',sample=samples["sample"])
    output:
        report(f'{derived}/results/multiqc_report.html', caption="multiqc.rst", category="Quality control")
    params:
        outdir= f'{derived}/results'  
    log: f'{logs}/MultiQC.log'
    group: "preprocess_align"
    threads:   2
    priority:  40
    shell:
        """
        multiqc {input.mosdepth} {input.stats} {input.json} \
        -o {params.outdir} 2> {log}
        """

rule create_summary_xls:
    input:
        report = f'{derived}/qc_reports/sequencing/{{sample}}/{{sample}}.sequencing-qc.txt',
        stats = f'{derived}/qc_reports/pairtools/{{sample}}/{{sample}}.stats',
        sum_stats = f'{derived}/qc_reports/pairtools/{{sample}}/{{sample}}.paired-tools.stats.txt',
        mosdepth = f'{derived}/qc_reports/mosdepth/{{sample}}/{{sample}}.mosdepth.summary.txt'
    output:
        summary_report = f'{derived}/results/{{sample}}.summary.xlsx'
    log:        f"{logs}/{{sample}}/21_convert_summary_xls.log"
    benchmark:  f"{benchmarks}/{{sample}}.convert_summary_xls.benchmark.txt"
    threads: 1
    group: "blood_typing"
    script:
        "../scripts/txt2xlsx.py"

rule gather_reports:
    input:
        rep = f'{derived}/results/{{sample}}.summary.xlsx'
    output:
        rep = f'{final}/{{sample}}.summary.xlsx'
    shell:
        """
        cp {input.rep} {output.rep}
        """

rule gather_html:
    input:
        html = f'{derived}/results/multiqc_report.html'
    output:
        html = f'{final}/Quality_report.html'
    shell:
        """
        cp {input.html} {output.html}
        """
