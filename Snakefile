configfile: "config.json"

rule all:
	input:
		expand("tmp/{genome}_{sample}/exonerate", genome=config["ref_genome"], sample=config["vcf_phased"])

rule phase_homozoygous:
	input:
		"vcf_phased/{genome}_{sample}_phased.bcf"
	output:
		first_nine="vcf_phased/{genome}_{sample}_first_nine.vcf",
		headers="vcf_phased/{genome}_{sample}_headers.vcf",
		tenth="vcf_phased/{genome}_{sample}_tenth.vcf",
		vcf_phased="vcf_phased/{genome}_{sample}_hz_phased.bcf"
	threads: 1
	run:
		shell("bcftools view -Ov {input} | grep $'^SL' | cut -f1,2,3,4,5,6,7,8,9 > {output.first_nine}")
		shell("bcftools view -Ov {input} | grep $'^SL' | cut -f10 | sed 's/^1\/1/1\|1/g' > {output.tenth}")
		shell("bcftools view -Ov {input} | grep $'#' > {output.headers}")
		shell("paste {output.first_nine} {output.tenth} >> {output.headers}")
		shell("bcftools view -Ob {output.headers} > {output.vcf_phased}")
		shell("bcftools index {output.vcf_phased}")

rule identify_phased_genes:
	input:
		work="vcf_phased/{genome}_{sample}_hz_phased.bcf",
		rm1="vcf_phased/{genome}_{sample}_first_nine.vcf",
		rm2="vcf_phased/{genome}_{sample}_headers.vcf",
		rm3="vcf_phased/{genome}_{sample}_tenth.vcf"
	output:
		gen_phased="candidates/{genome}_{sample}_fully_phased.txt",
		single_block="candidates/{genome}_{sample}_single_block_phased.txt",
		one_miss="candidates/{genome}_{sample}_one_phased_missing.txt"
	params:
		names="ref_genome/ITAG4.0_proteins_names.txt",
		cds=expand("ref_genome/{cds_reference}.gff", cds_reference=config["ref_genome_annot"])

	threads: 2
	conda:
		"envs/cds_extraction.yaml"
	shell:
		"cat {params.names} | xargs -I@ -P 2 sh -c 'scripts/phased_genes.sh @ {params.cds} {input.work} {output.gen_phased} {output.single_block} {output.one_miss}'"

rule locus_bed:
	input:
		one_miss="candidates/{genome}_{sample}_one_phased_missing.txt"
	output:
		directory("tmp/{genome}_{sample}/locus_bed")
	params:
		mrna="ref_genome/ITAG4.0_mrna_models.gff"
	threads: 2
	run:
		shell("mkdir {output}")
		shell("cat {input.one_miss} | xargs -I@ -P 2 sh -c 'scripts/awk_gene.sh @ {params.mrna} {output}'")

	
rule phased_cds_bed:
	input:
		"candidates/{genome}_{sample}_one_phased_missing.txt"
	output:
		"vcf_phased/{genome}_{sample}_hz_phased_one_missing.bed"
	params:
		cds=expand("ref_genome/{cds_reference}.gff", cds_reference=config["ref_genome_annot"]),
	threads: 1
	conda:
		"envs/cds_extraction.yaml"
	shell:
		"cat {input} | xargs -I@ sh -c 'grep @ {params.cds} | cut -f1,4,5 >> {output}'"

rule phased_cds_intersect:
	input:
		bed="vcf_phased/{genome}_{sample}_hz_phased_one_missing.bed",
		bcf="vcf_phased/{genome}_{sample}_hz_phased.bcf"
	output:
		vcf="vcf_phased/{genome}_{sample}_hz_phased_one_missing.vcf",
		vcf_gz="vcf_phased/{genome}_{sample}_hz_phased_one_missing.vcf.gz"
	threads: 1
	run:
		shell("bcftools view -Ov {input.bcf} | grep $'#' - > {output.vcf}")
		shell("bcftools view -Oz {input.bcf} | bedtools intersect -a stdin -b {input.bed} >> {output.vcf}")
		shell("bcftools view -Oz {output.vcf} > {output.vcf_gz}")
		shell("bcftools index {output.vcf_gz}")

rule phased_genome_h1:
	input:
		vcf_gz="vcf_phased/{genome}_{sample}_hz_phased_one_missing.vcf.gz",
	output:
		"tmp/{genome}_{sample}/h1_one_missing.fasta"
	params:
		htype=1,
		ref_genome=expand("ref_genome/{genome}.fa", genome=config["ref_genome"])
	threads: 1
	conda:
		"envs/cds_extraction.yaml"
	shell:
		"cat {params.ref_genome} | bcftools consensus -H {params.htype} {input.vcf_gz} > {output}"

rule index_genome:
	input:
		"tmp/{genome}_{sample}/h1_one_missing.fasta"
	output:
		"tmp/{genome}_{sample}/h1_one_missing.fasta.fai"
	threads: 1
	conda:
		"envs/cds_extraction.yaml"
	shell:
		"samtools faidx {input};"
		"touch {output}"


rule mrna_phased_1:
	input:
		one_miss="candidates/{genome}_{sample}_one_phased_missing.txt",
		locus="tmp/{genome}_{sample}/locus_bed",
		genome="tmp/{genome}_{sample}/h1_one_missing.fasta",
		index="tmp/{genome}_{sample}/h1_one_missing.fasta.fai"

	output:
		directory("tmp/{genome}_{sample}/locus_fasta_h1")
	params:	htype=1
	threads: 2
        run:
		shell("mkdir {output}")
		shell("cat {input.one_miss} | xargs -I@ -P 2 sh -c 'scripts/mrna_phased.sh @ {input.genome} {input.locus} {output} {params.htype}'")


rule exonerate:
	input:
		one_miss="candidates/{genome}_{sample}_one_phased_missing.txt",
		h1_fasta="tmp/{genome}_{sample}/locus_fasta_h1"

	output:
		tmp=directory("tmp/{genome}_{sample}/exonerate"),
		result=directory("output_cds/{genome}_{sample}")
	params:
		proteome="proteome"
	threads: 2
	run:
		shell("mkdir {output.tmp}; mkdir {output.result}; cat {input.one_miss} | xargs -I@ -n 1 -P 2 sh -c 'scripts/exonerate.sh @ {input.h1_fasta} {params.proteome} {output.tmp} {output.result}'")

## NEXT: I'm gettin an error in the exonerate search (initial HSP -1???) doesn't happen with many, but worth looking into
## Then: generate second haplotype as well and try pipeline for more than one species

