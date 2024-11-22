rule filter_and_index_FrEIA:
    input:
        bam=f'{derived}/recal/{{sample_type}}.bam'
    output:
        bam= f'{derived}/fragment_analysis/{{sample_type}}.bam',
        bai= f'{derived}/fragment_analysis/{{sample_type}}.bam.bai'
    threads: config['ncores']
    params:
        qual = config["FiltQual"]
    conda: "../envs/samtools.yaml"
    log: f'{logs}/fragment_analysis/{{sample_type}}/filter_and_index_FrEIA.log'
    benchmark: f'{benchmarks}/{{sample_type}}_filter_and_index_FrEIA.log'
    shell:
        """
        samtools view {input.bam} \
        -h `# include header` \
        -q {params.qual} `#> only reads with high quality` \
        -F 4 `# not unmapped` \
        -F 256 `# no secondary alignment; thus primary alignment` \
        -F 1024 `# no PCR duplicate` \
        -F 2048 `# no supplementary (=chimeric) alignment` \
        -@ {threads} \
        -b | tee {output.bam} | \
        \
        samtools index - `# stdin` \
        -@ {threads} \
        {output.bai}
        """

rule FragmEndSeq:
    input:
        bam= f'{derived}/freia/{{sample_type}}.bam'
    output:
        pq= report(f'{derived}/freia/{{sample_type}}.pq',
               category="Fragment end retrieval")
    threads: config["ThreadNr"]
    params:
        freia_script= "scripts/1_extract_fragment_ends.py",
        nBase = config["nBase"],
        contigs = chrs
    log: f'{logs}/fragment_analysis/{{sample_type}}/FragmEndSeq.log'
    benchmark: f'{benchmarks}/{{sample_type}}_FragmEndSeq.log'
    shell:
        """
        python3 {params.freia_script} \
        -b {input.bam} \
        -o {output.pq} \
        -c {params.contigs} \
        -t {threads} \
        -n {params.nBase} \
        2> {log}
        """

# Calculating fragment end sequence diversity and creating outputs,
# containing all the samples.
rule compare_groups:
    input:
        expand(config["OutPath"] + "/" + ProjDirName +
               "/4_FrEIA/3_Abundances/{group}/{lvl}/{prefix}__{sample}.pq",
               zip,
               group=np.repeat(SampleGroups, (len(Prefixes) * len(Levels))),
               lvl=np.repeat(Levels, len(Prefixes)).tolist() * len(SampleGroups),
               prefix=Prefixes * len(Levels) * len(SampleGroups),
               sample=np.repeat(Samplesheet["sample_name"],
                                (len(Prefixes) * len(Levels))).tolist() * len(SampleGroups))
    output:
        config["OutPath"] + "/" + ProjDirName +
        "/4_FrEIA/4_Compare/Data/Dat_GM__sample.csv",
        config["OutPath"] + "/" + ProjDirName +
        "/4_FrEIA/4_Compare/Data/Dat_GT__MDS_sample.csv",
        config["OutPath"] + "/" + ProjDirName +
        "/4_FrEIA/4_Compare/Data/Dat_GT__sample.csv"
    threads: config["ThreadNr"]
    conda: "../../envs/FrEIA_env.yaml"
    params:
        freia_script = "scripts/3_compare_groups.py",
        outPath = config["OutPath"] + "/" + ProjDirName + "/",
        sampTable = config["Samplesheet"],
        subsResults = config["SubsetResults"],
        regroup = config["Regroup"]
    conda: "../../envs/FrEIA_env.yaml"
    log:
        (config["OutPath"] + "/logs/" + ProjDirName +
         "/4_FrEIA/4_Compare/FrEIA_compare.log")
    shell:
        """
        python3 ../../scripts/1_FrEIA/3_compare_groups.py \
        -i {params.outPath} \
        -st {params.sampTable} \
        -t {threads} \
        -sN {params.subsResults} \
        -rgr {params.regroup}
        """

# Calculating base and motif fractions per sample.
rule data_transformation:
    input:
        config["OutPath"] + "/" + ProjDirName +
        "/4_FrEIA/1_extract_fragment_ends/{sample}.pq"
    output:
        config["OutPath"] + "/" + ProjDirName +
        "/4_FrEIA/3_Abundances/{group}/{lvl}/M__{sample}.pq",
        config["OutPath"] + "/" + ProjDirName +
        "/4_FrEIA/3_Abundances/{group}/{lvl}/T__{sample}.pq",
    params:
        outPath = config["OutPath"] + "/" + ProjDirName + "/",
        sampTable = config["Samplesheet"],
        fraSizeMin = config["FragmSizeMin"],
        fraSizeMax = config["FragmSizeMax"],
        subSamp = config["SubSampleRate"],
        bsSampNr = config["BsSampNr"]
    conda: "../../envs/FrEIA_env.yaml"
    log:
        (config["OutPath"] + "/logs/" + ProjDirName +
         "/4_FrEIA/3_Abundances/{group}/{lvl}/M__{sample}.log"),
        (config["OutPath"] + "/logs/" + ProjDirName +
         "/4_FrEIA/3_Abundances/{group}/{lvl}/T__{sample}.log")
    shell:
        """
        python3 ../../scripts/1_FrEIA/2_data_transformation.py \
        -i {input} \
        -o {params.outPath} \
        -st {params.sampTable} \
        -fsmin {params.fraSizeMin} \
        -fsmax {params.fraSizeMax} \
        -subs {params.subSamp} \
        --bootstrap_sample {params.bsSampNr}
        """

# Compute the FrEIA score.
rule FrEIA_score:
    input:
        config["OutPath"] + "/" + ProjDirName +
        "/4_FrEIA/4_Compare/Data/Dat_GM__sample.csv",
        config["OutPath"] + "/" + ProjDirName +
        "/4_FrEIA/4_Compare/Data/Dat_GT__MDS_sample.csv",
        config["OutPath"] + "/" + ProjDirName +
        "/4_FrEIA/4_Compare/Data/Dat_GT__sample.csv"
    output:
        config["OutPath"] + "/" + ProjDirName +
        "/4_FrEIA/5_FrEIA_score/" + config["ProjName"] + "_logFC_p.csv",
        config["OutPath"] + "/" + ProjDirName +
        "/4_FrEIA/5_FrEIA_score/" + config["ProjName"] + "_FrEIA_score.csv"
    params:
        inPath = config["OutPath"] + "/" + ProjDirName +
                 "/4_FrEIA/4_Compare/Data/",
        metaPath = config["MetaPath"],
        ctrName = config["CtrName"],
        casName= config["CasName"],
        panelMed = config["PanelMed"],
        tncCtr = config["TncCtr"],
        tncCas= config["TncCas"],
        outPath = config["OutPath"] + "/" + ProjDirName +
                  "/4_FrEIA/5_FrEIA_score/" + config["ProjName"],
        batch = config["BatchIDs"],
        threads = config["ThreadNr"]
    conda: "../../envs/FrEIA_env.yaml"
    log:
        (config["OutPath"] + "/logs/" + ProjDirName +
         "/4_FrEIA/5_FrEIA_score/FrEIA_score.log")
    shell:
        """
        python3 ../../scripts/1_FrEIA/4_FrEIA_score.py \
        -i {params.inPath} \
        -o {params.outPath} \
        -m {params.metaPath} \
        --control_name {params.ctrName} \
        --case_name {params.casName} \
        -p {params.panelMed} \
        --tnc_control {params.tncCtr} \
        --tnc_case {params.tncCas} \
        -b {params.batch} \
        --gini \
        -t {params.threads} \
        2> {log}
        """