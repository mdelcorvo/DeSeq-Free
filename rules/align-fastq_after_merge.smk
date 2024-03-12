rule bwa_mem_plasma_notCombined:
    input:
        reads=get_trimmed_reads_notCombined,
        idx=multiext(f'{genome_data}/GRCh38.d1.vd1.fa',".amb",".ann",".bwt",".pac",".sa"),
        fai=f'{genome_data}/GRCh38.d1.vd1.fa.fai'
    output:
        bam=temp(f'{derived}/alignments/{{sample_type}}.notCombined.bam')
    log: f'{logs}/{{sample_type}}/03_bwa_mem_alignment.notCombined.log'
    benchmark: f'{benchmarks}/{{sample_type}}_bwa_mem_alignment.notCombined.txt'
    params:
        extra=get_read_group,
        sort_extra=""
    threads: config['ncores']
    priority: 50
    wrapper:
        "v3.4.1/bio/bwa/mem-samblaster"

rule bwa_mem_plasma_extendedFrags:
    input:
        reads=get_trimmed_reads_extendedFrags,
        idx=multiext(f'{genome_data}/GRCh38.d1.vd1.fa',".amb",".ann",".bwt",".pac",".sa"),
        fai=f'{genome_data}/GRCh38.d1.vd1.fa.fai'
    output:
        bam=temp(f'{derived}/alignments/{{sample_type}}.extendedFrags.bam')
    log: f'{logs}/{{sample_type}}/03_bwa_mem_alignment.extendedFrags.log'
    benchmark: f'{benchmarks}/{{sample_type}}_bwa_mem_alignment.extendedFrags.txt'
    params:
        extra=get_read_group,
        samblaster_extra='--ignoreUnmated',
        sort_extra=""
    threads: config['ncores']
    priority: 50
    wrapper:
        "v3.4.1/bio/bwa/mem-samblaster"

rule realignertargetcreator_plasma_notCombined:
    input:
        bam=f'{derived}/alignments/{{sample_type}}.notCombined.bam',
        ref = f'{genome_data}/GRCh38.d1.vd1.fa',
        fai=f'{genome_data}/GRCh38.d1.vd1.fa.fai',
        dict=f'{genome_data}/GRCh38.d1.vd1.dict',
        known=config['dbsnp'],
        known_idx=config['dbsnp_tbi']
    output:
        intervals=f'{derived}/alignments/{{sample_type}}.notCombined.GATK3.intervals',
    log: f'{logs}/{{sample_type}}/04_GATK_realignertargetcreator.notCombined.log'
    benchmark: f'{benchmarks}/{{sample_type}}_04_GATK_realignertargetcreator.notCombined.txt'
    threads: config['ncores']
    priority: 50
    resources:
        mem_mb=4096,
    wrapper:
        "v3.4.1/bio/gatk3/realignertargetcreator"

rule realignertargetcreator_plasma_extendedFrags:
    input:
        bam=f'{derived}/alignments/{{sample_type}}.extendedFrags.bam',
        ref = f'{genome_data}/GRCh38.d1.vd1.fa',
        fai=f'{genome_data}/GRCh38.d1.vd1.fa.fai',
        dict=f'{genome_data}/GRCh38.d1.vd1.dict',
        known=config['dbsnp'],
        known_idx=config['dbsnp_tbi']
    output:
        intervals=f'{derived}/alignments/{{sample_type}}.extendedFrags.GATK3.intervals',
    log: f'{logs}/{{sample_type}}/04_GATK_realignertargetcreator.extendedFrags.log'
    benchmark: f'{benchmarks}/{{sample_type}}_04_GATK_realignertargetcreator.extendedFrags.txt'
    threads: config['ncores']
    priority: 50
    resources:
        mem_mb=4096,
    wrapper:
        "v3.4.1/bio/gatk3/realignertargetcreator"

rule indelrealigner_plasma_notCombined:
    input:
        bam=f'{derived}/alignments/{{sample_type}}.notCombined.bam',
        ref=f'{genome_data}/GRCh38.d1.vd1.fa',
        fai=f'{genome_data}/GRCh38.d1.vd1.fa.fai',
        dict=f'{genome_data}/GRCh38.d1.vd1.dict',
        known=config['dbsnp'],
        known_idx=config['dbsnp_tbi'],
        target_intervals=f'{derived}/alignments/{{sample_type}}.notCombined.GATK3.intervals'
    output:
        bam=temp(f'{derived}/alignments/{{sample_type}}.notCombined.rea.bam'),
        bai=temp(f'{derived}/alignments/{{sample_type}}.notCombined.rea.bai'),
    log: f'{logs}/{{sample_type}}/05_GATK_indelrealigner.notCombined.log'
    benchmark: f'{benchmarks}/{{sample_type}}_05_GATK_indelrealigner.notCombined.txt'
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    threads: config['ncores']
    resources:
        mem_mb=4096,
    wrapper:
        "v3.4.1/bio/gatk3/indelrealigner"

rule indelrealigner_plasma_extendedFrags:
    input:
        bam=f'{derived}/alignments/{{sample_type}}.extendedFrags.bam',
        ref=f'{genome_data}/GRCh38.d1.vd1.fa',
        fai=f'{genome_data}/GRCh38.d1.vd1.fa.fai',
        dict=f'{genome_data}/GRCh38.d1.vd1.dict',
        known=config['dbsnp'],
        known_idx=config['dbsnp_tbi'],
        target_intervals=f'{derived}/alignments/{{sample_type}}.extendedFrags.GATK3.intervals'
    output:
        bam=temp(f'{derived}/alignments/{{sample_type}}.extendedFrags.rea.bam'),
        bai=temp(f'{derived}/alignments/{{sample_type}}.extendedFrags.rea.bai'),
    log: f'{logs}/{{sample_type}}/05_GATK_indelrealigner.extendedFrags.log'
    benchmark: f'{benchmarks}/{{sample_type}}_05_GATK_indelrealigner.extendedFrags.txt'
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    threads: config['ncores']
    resources:
        mem_mb=4096,
    wrapper:
        "v3.4.1/bio/gatk3/indelrealigner"


rule gatk_baserecalibrator_plasma_notCombined:
    input:
        bam=f'{derived}/alignments/{{sample_type}}.notCombined.rea.bam',
        ref=f'{genome_data}/GRCh38.d1.vd1.fa',
        dict=f'{genome_data}/GRCh38.d1.vd1.dict',
        known=config['dbsnp']
    output:
        recal_table=f'{derived}/alignments/{{sample_type}}.notCombined.rt',
    log: f'{logs}/{{sample_type}}/06_GATK_baserecalibrator.notCombined.log'
    benchmark: f'{benchmarks}/{{sample_type}}_06_GATK_baserecalibrator.notCombined.txt'
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=4096,
    wrapper:
        "v3.4.1/bio/gatk/baserecalibrator"

rule gatk_baserecalibrator_plasma_extendedFrags:
    input:
        bam=f'{derived}/alignments/{{sample_type}}.extendedFrags.rea.bam',
        ref=f'{genome_data}/GRCh38.d1.vd1.fa',
        dict=f'{genome_data}/GRCh38.d1.vd1.dict',
        known=config['dbsnp']
    output:
        recal_table=f'{derived}/alignments/{{sample_type}}.extendedFrags.rt',
    log: f'{logs}/{{sample_type}}/06_GATK_baserecalibrator.extendedFrags.log'
    benchmark: f'{benchmarks}/{{sample_type}}_06_GATK_baserecalibrator.extendedFrags.txt'
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=4096,
    wrapper:
        "v3.4.1/bio/gatk/baserecalibrator"

rule gatk_applybqsr_plasma_notCombined:
    input:
        bam=f'{derived}/alignments/{{sample_type}}.notCombined.rea.bam',
        ref=f'{genome_data}/GRCh38.d1.vd1.fa',
        dict=f'{genome_data}/GRCh38.d1.vd1.dict',
        recal_table=f'{derived}/alignments/{{sample_type}}.notCombined.rt'
    output:
        bam=f'{derived}/recal/{{sample_type}}.notCombined.final.bam'
    log: f'{logs}/{{sample_type}}/07_GATK_APPLYBQSR.notCombined.log'
    benchmark: f'{benchmarks}/{{sample_type}}_07_GATK_APPLYBQSR.notCombined.txt'
    params:
        extra="",  # optional
        java_opts="",  # optional
        embed_ref=True,  # embed the reference in cram output
    resources:
        mem_mb=4096,
    wrapper:
        "v3.4.1/bio/gatk/applybqsr"

rule gatk_applybqsr_plasma_extendedFrags:
    input:
        bam=f'{derived}/alignments/{{sample_type}}.extendedFrags.rea.bam',
        ref=f'{genome_data}/GRCh38.d1.vd1.fa',
        dict=f'{genome_data}/GRCh38.d1.vd1.dict',
        recal_table=f'{derived}/alignments/{{sample_type}}.extendedFrags.rt'
    output:
        bam=f'{derived}/recal/{{sample_type}}.extendedFrags.final.bam'
    log: f'{logs}/{{sample_type}}/07_GATK_APPLYBQSR.extendedFrags.log'
    benchmark: f'{benchmarks}/{{sample_type}}_07_GATK_APPLYBQSR.extendedFrags.txt'
    params:
        extra="",  # optional
        java_opts="",  # optional
        embed_ref=True,  # embed the reference in cram output
    resources:
        mem_mb=4096,
    wrapper:
        "v3.4.1/bio/gatk/applybqsr"