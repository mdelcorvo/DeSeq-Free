rule samtools_macs2:
    input:
        mapped = f'{derived}/alignments/pairtools/{{sample}}/{{sample}}.mapped.PT.bam'
    output:
        bed = f'{derived}/alignments/pairtools/{{sample}}/{{sample}}.mapped.PT.bed'
    log: f'{logs}/{{sample}}/16_samtools_view_pre_macs2.log'
    benchmark: f'{benchmarks}/{{sample}}.samtools_view_pre_macs2.benchmark.txt'
    group:   "preprocess_align"
    threads:  1
    priority: 40
    shell:
        """
        samtools view -h -F 0x904 {input.mapped} \
        | bedtools bamtobed -i stdin \
        > {output.bed}
        
        """ 
        
rule macs2:
    input:
        bed = f'{derived}/alignments/pairtools/{{sample}}/{{sample}}.mapped.PT.bed'
    output:
        peak = f'{derived}/loop_calling/macs2/{{sample}}/{{sample}}_peaks.narrowPeak',
        xls = f'{derived}/loop_calling/macs2/{{sample}}/{{sample}}_peaks.xls'
    params:
        sample_name = f'{derived}/loop_calling/macs2/{{sample}}/{{sample}}'    
    log: f'{logs}/{{sample}}/17_macs2.log'
    benchmark: f'{benchmarks}/{{sample}}.macs2.benchmark.txt'
    group:   "preprocess_align"
    threads:  1
    priority: 40
    shell:
        """       
        macs2 callpeak \
        --treatment {input.bed} \
        -n {params.sample_name} 
        """ 
          
rule fithichip_prep_input:
    input:
        mapped = f'{derived}/alignments/pairtools/{{sample}}/{{sample}}.mapped.pairs'
    output:
        hicpro = f'{derived}/alignments/pairtools/{{sample}}/{{sample}}_hicpro_mapped.pairs.gz'
    log: f'{logs}/{{sample}}/18_fithichip_prep_input.log'
    benchmark: f'{benchmarks}/{{sample}}.fithichip_prep_input.benchmark.txt'
    group:   "preprocess_align"
    threads:  1
    priority: 40
    shell:
        """       
        grep -v '#' {input.mapped} \
        | awk -F"\t" '{{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}}' \
        | gzip -c > {output.hicpro}
        """ 

rule fithichip_prep_config:
    input:
        hicpro = f'{derived}/alignments/pairtools/{{sample}}/{{sample}}_hicpro_mapped.pairs.gz',
        peak = f'{derived}/loop_calling/macs2/{{sample}}/{{sample}}_peaks.narrowPeak',
        config = "tools/FitHiChIP_config.txt"
    output:
        config = f'{derived}/loop_calling/fithichip/{{sample}}/{{sample}}_FitHiChIP_config.txt'
    params:
        validpairs = get_validpairs,
        peakfile = get_peakfile,
        outdir = get_outdir,
        inttype = get_inttype,
        binsize = get_binsize,
        usep2pbackgrnd = get_usep2pbackgrnd,
        biastype = get_biastype,
        qvalue = get_qvalue,
        chrsizefile = get_chrsizefile,
        prefix = get_prefix
    log: f'{logs}/{{sample}}/19_fithichip_prep_config.log'
    benchmark: f'{benchmarks}/{{sample}}.fithichip_prep_config.benchmark.txt'
    group:   "preprocess_align"
    threads:  1
    priority: 40
    shell:
        """       
        cp {input.config} {output.config}
        sed -i "11i {params.validpairs}" {output.config}
        sed -i "31i {params.peakfile}" {output.config}
        sed -i "34i {params.outdir}" {output.config}
        sed -i "37i {params.inttype}" {output.config}
        sed -i "40i {params.binsize}" {output.config}
        sed -i "53i {params.usep2pbackgrnd}" {output.config}
        sed -i "57i {params.biastype}" {output.config}
        sed -i "64i {params.qvalue}" {output.config}
        sed -i "67i {params.chrsizefile}" {output.config}
        sed -i "70i {params.prefix}" {output.config}
        """ 
        
rule fithichip_run:
    input:
        config = f'{derived}/loop_calling/fithichip/{{sample}}/{{sample}}_FitHiChIP_config.txt',
        tool = "tools/FitHiChIP_HiCPro.sh"
    output:
        res = f'{derived}/loop_calling/fithichip/{{sample}}/Summary_results_FitHiChIP.html'
    log: f'{logs}/{{sample}}/20_fithichip_run.log'
    benchmark: f'{benchmarks}/{{sample}}.fithichip_run.benchmark.txt'
    group:   "preprocess_align"
    threads:  1
    priority: 40
    shell:
        """       
        {input.tool} -C {input.config}
        """ 
