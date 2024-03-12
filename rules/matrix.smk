rule contact_matrix:
    input:
        mapped = f'{derived}/alignments/pairtools/{{sample}}/{{sample}}.mapped.pairs',
        ref= hg38,
        tool = "tools/juicertools.jar"
    output:
        hic = f'{derived}/matrix/juicer/{{sample}}/{{sample}}.contact_map.hic'
    log: f'{logs}/{{sample}}/10_contact_matrix.log'
    benchmark: f'{benchmarks}/{{sample}}.contact_matrix.benchmark.txt'
    group:   "preprocess_align"
    threads:  config['ncores']
    priority: 40
    shell:
        """
        
        java -Xmx32000m  -Djava.awt.headless=true -jar \
        {input.tool} pre --threads {threads} \
        {input.mapped} {output.hic} {input.ref}
        
        """ 
  
rule pairix:
    input:
        mapped = f'{derived}/alignments/pairtools/{{sample}}/{{sample}}.mapped.pairs',
        pairix = "tools/pairix"
    output:
        mapped = f'{derived}/matrix/pairix/{{sample}}/{{sample}}.mapped.pairs.gz'
    log: f'{logs}/{{sample}}/11_pairix.log'
    benchmark: f'{benchmarks}/{{sample}}.pairix.benchmark.txt'
    group:   "preprocess_align"
    threads:  1
    priority: 40
    shell:
        """    
        bgzip -c {input.mapped} > {output.mapped}
        
        {input.pairix} {output.mapped}     
        """ 
if k5:  
	rule cooler_5k:
		input:
			mapped = f'{derived}/matrix/pairix/{{sample}}/{{sample}}.mapped.pairs.gz',
			ref= f'{genome_data}/GRCh38.genome'
		output:
			matrix_5kb  = f'{derived}/matrix/cooler/{{sample}}/{{sample}}.matrix_5kb.cool',
			txt_5kb  = f'{derived}/matrix/cooler/{{sample}}/{{sample}}.matrix_5kb.txt'
		log: f'{logs}/{{sample}}/12_cooler_5k.log'
		benchmark: f'{benchmarks}/{{sample}}.cooler_5k.benchmark.txt'
		group:   "preprocess_align"
		threads:  config['ncores']
		priority: 40
		shell:
			"""
			
			cooler cload pairix -p {threads} {input.ref}:5000 {input.mapped} {output.matrix_5kb}         
			cooler dump --join {output.matrix_5kb} > {output.txt_5kb}        
			cooler zoomify --balance -p {threads} {output.matrix_5kb}        

			""" 

rule cooler_10k:
    input:
        mapped = f'{derived}/matrix/pairix/{{sample}}/{{sample}}.mapped.pairs.gz',
        ref= f'{genome_data}/GRCh38.genome'
    output:
        matrix_10kb  = f'{derived}/matrix/cooler/{{sample}}/{{sample}}.matrix_10kb.cool',
        txt_10kb  = f'{derived}/matrix/cooler/{{sample}}/{{sample}}.matrix_10kb.txt'
    log: f'{logs}/{{sample}}/13_cooler_10k.log'
    benchmark: f'{benchmarks}/{{sample}}.cooler_10k.benchmark.txt'
    group:   "preprocess_align"
    threads:  config['ncores']
    priority: 40
    shell:
        """
        
        cooler cload pairix -p {threads} {input.ref}:10000 {input.mapped} {output.matrix_10kb}         
        cooler dump --join {output.matrix_10kb} > {output.txt_10kb}        
        cooler zoomify --balance -p {threads} {output.matrix_10kb}        

        """ 
        
rule cooler_50k:
    input:
        mapped = f'{derived}/matrix/pairix/{{sample}}/{{sample}}.mapped.pairs.gz',
        ref= f'{genome_data}/GRCh38.genome'
    output:
        matrix_50kb  = f'{derived}/matrix/cooler/{{sample}}/{{sample}}.matrix_50kb.cool',
        txt_50kb  = f'{derived}/matrix/cooler/{{sample}}/{{sample}}.matrix_50kb.txt'
    log: f'{logs}/{{sample}}/14_cooler_50k.log'
    benchmark: f'{benchmarks}/{{sample}}.cooler_50k.benchmark.txt'
    group:   "preprocess_align"
    threads:  max(1, int(config['ncores'] // 3))
    priority: 40
    shell:
        """
        
        cooler cload pairix -p {threads} {input.ref}:50000 {input.mapped} {output.matrix_50kb}         
        cooler dump --join {output.matrix_50kb} > {output.txt_50kb}        
        cooler zoomify --balance -p {threads} {output.matrix_50kb}        

        """      
           
rule add_gene_info:
    input:
        txt  = f'{derived}/matrix/cooler/{{sample}}/{{sample}}.matrix_{{bin}}kb.txt',
        ref = get_transcripts,
        script = "scripts/add_gene_info.R"    
    output:
        txt_ID = temp(f'{derived}/matrix/annotation/{{sample}}/{{sample}}.matrix_{{bin}}kb.ID.txt'),
        txt_1  = temp(f'{derived}/matrix/annotation/{{sample}}/{{sample}}.matrix_{{bin}}kb.1.txt'),
        txt_2  = temp(f'{derived}/matrix/annotation/{{sample}}/{{sample}}.matrix_{{bin}}kb.2.txt'),
        txt_A1 = temp(f'{derived}/matrix/annotation/{{sample}}/{{sample}}.matrix_{{bin}}kb.A1.txt'),
        txt_A2 = temp(f'{derived}/matrix/annotation/{{sample}}/{{sample}}.matrix_{{bin}}kb.A2.txt'),
        annot  = f'{derived}/matrix/annotation/{{sample}}/{{sample}}.matrix_{{bin}}kb.ANNOT.txt'
    log: f'{logs}/{{sample}}/15_add_gene_info_{{bin}}kb.log'
    benchmark: f'{benchmarks}/{{sample}}.add_gene_info_{{bin}}kb.benchmark.txt'
    group:   "preprocess_align"
    threads:  1
    priority: 40
    shell:
        """
        
        awk -F"\t" 'NR>0{{$0=$0"\t"NR-1}} 1' \
        {input.txt} > {output.txt_ID}
        
        awk 'OFS="\t" {{print $1,$2,$3,$8}}' \
        {output.txt_ID} \
        > {output.txt_1}
        
        awk 'OFS="\t" {{print $4,$5,$6,$8}}' \
        {output.txt_ID} \
        > {output.txt_2}
        
        bedtools intersect -wao -a {output.txt_1} \
        -b {input.ref} | awk 'OFS="\t" {{print $4,$8,$9}}' > {output.txt_A1}
        
        bedtools intersect -wao -a  {output.txt_2} \
        -b {input.ref} | awk 'OFS="\t" {{print $4,$8,$9}}' > {output.txt_A2}
        
        {input.script} {output.txt_ID} {output.txt_A1} \
        {output.txt_A2} {output.annot} {log}
        
        """
        
rule gather_contact_matrices:
	input:
		txt_5kb  = f'{derived}/matrix/annotation/{{sample}}/{{sample}}.matrix_5kb.ANNOT.txt' if k5 else [],
		txt_10kb  = f'{derived}/matrix/annotation/{{sample}}/{{sample}}.matrix_10kb.ANNOT.txt',
		txt_50kb  = f'{derived}/matrix/annotation/{{sample}}/{{sample}}.matrix_50kb.ANNOT.txt',
		hic = f'{derived}/matrix/juicer/{{sample}}/{{sample}}.contact_map.hic'
	output:
		zip=f'{final}/contact_matrix/{{sample}}.matrix.zip'
	run:
		if config["k5"]:
			cmd = "zip -j '{output.zip}' '{input.txt_5kb}' \
			'{input.txt_10kb}' '{input.txt_50kb}' "
		else:
			cmd = "zip -j '{output.zip}' \
			'{input.txt_10kb}' '{input.txt_50kb}' "
		shell(cmd)	
			
