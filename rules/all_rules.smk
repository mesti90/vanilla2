rule validate_fastq:
	input:
		r1=lambda wc: R1[wc.sample],
		r2=lambda wc: R2[wc.sample]

	output:
		touch("work/validation/{sample}.ok")

	container:
		config["containers"]["seqkit"]

	shell:
		r"""
		set -euo pipefail

		gzip -t {input.r1}
		gzip -t {input.r2}

		n1=$(seqkit stats -T {input.r1} | tail -n1 | cut -f4)
		n2=$(seqkit stats -T {input.r2} | tail -n1 | cut -f4)

		[ "$n1" -eq "$n2" ]

		paste \
			<(seqkit seq -n {input.r1} | cut -d' ' -f1) \
			<(seqkit seq -n {input.r2} | cut -d' ' -f1) \
		| awk '
			$1 != $2 {{
				print "ERROR: read pair mismatch"
				print "R1:", $1
				print "R2:", $2
				exit 1
			}}
		'

		touch {output}
		"""

# Preprocess reads: trimming
rule fastp:
	input:
		r1=lambda wc: R1[wc.sample],
		r2=lambda wc: R2[wc.sample],
		validate=rules.validate_fastq.output
	output:
		r1="work/clean_reads/{sample}_R1.fastq.gz",
		r2="work/clean_reads/{sample}_R2.fastq.gz",
		html="work/clean_reads/{sample}.html",
		json="work/clean_reads/{sample}.json"
	threads:
		config["threads"]["fastp"]
	container:
		config["containers"]["fastp"]
	shell:
		r"""
		fastp \
			-i {input.r1} \
			-I {input.r2} \
			-o {output.r1} \
			-O {output.r2} \
			-h {output.html} \
			-j {output.json} \
			-w {threads}
		"""


rule snippy:
	input:
		r1=rules.fastp.output.r1,
		r2=rules.fastp.output.r2,
		ref=lambda wc: f"References/{SAMPLE_REF_NAME[wc.sample]}.gbk"
	output:
		vcf="work/snippy/{sample}/snps.vcf",
		tab="work/snippy/{sample}/snps.tab",
		snpeff_reference=directory("work/snippy/{sample}/reference"),
		ref_fa="work/snippy/{sample}/reference/ref.fa"
	params:
		outdir="work/snippy/{sample}",
		tmpdir=lambda wc: os.path.abspath(f"work/snippy/{wc.sample}")
	threads:
		config["threads"]["snippy"]
	container:
		config["containers"]["snippy"]
	shell:
		"""
		set -euo pipefail
		mkdir -p {params.tmpdir}
		mkdir -p {params.outdir}
		export TMPDIR={params.tmpdir}	
		snippy --outdir {params.outdir} --tmpdir {params.tmpdir} --ref {input.ref} --R1 {input.r1} --R2 {input.r2} --cpus {threads} --force --quiet 2> /dev/null
		"""

rule breseq:
	input:
		r1=rules.fastp.output.r1,
		r2=rules.fastp.output.r2,
		ref=rules.snippy.input.ref
	output:
		vcf="work/breseq/{sample}/output/output.vcf"
	params:
		outdir="work/breseq/{sample}",
		tmpdir="work/breseq/{sample}/tmp"
	threads:
		config["threads"]["breseq"]
	container:
		config["containers"]["breseq"]

	shell:
		"""
		set -euo pipefail
		mkdir -p {params.tmpdir}
		mkdir -p {params.outdir}
		export TMPDIR={params.tmpdir}
		breseq -j {threads} -n {wildcards.sample} -o {params.outdir} -r {input.ref} {input.r1} {input.r2} > {params.outdir}/breseq.logm
		"""

rule restore_breseq_contigs:
	input:
		vcf=rules.breseq.output.vcf,
		ref=rules.snippy.input.ref
	output:
		vcf="work/breseq/{sample}/output/output.original_contigs.vcf"
	script:
		"../scripts/restore_breseq_contigs.py"

rule snpeff_prepare_from_snippy:
	input:
		refdir=rules.snippy.output.snpeff_reference
	output:
		refdir=directory("work/snpeff_db/{sample}")
	container:
		config["containers"]["snippy"]
	shell:
		r"""
		set -euo pipefail

		rm -rf {output.refdir}
		mkdir -p {output.refdir}

		cp -a {input.refdir}/. {output.refdir}/

		CONFIG={output.refdir}/snpeff.config

		# Remove any existing data.dir
		sed -i '/^data\.dir[[:space:]]*=/d' "$CONFIG"

		# Prepend an absolute data.dir
		sed -i "1idata.dir = $(realpath {output.refdir})" "$CONFIG"
		"""

rule snpeff_for_breseq:
	input:
		vcf=rules.restore_breseq_contigs.output.vcf,
		refdir=rules.snpeff_prepare_from_snippy.output.refdir
	output:
		vcf="work/breseq/{sample}/breseq_out.restored_contigs.annotated.vcf"
	container:
		config["containers"]["snippy"]
	shell:
		r"""
		snpEff \
			-c {input.refdir}/snpeff.config \
			-noLog \
			-noStats \
			-no-downstream \
			-no-upstream \
			-no-utr \
			ref \
			{input.vcf} \
			> {output.vcf}
		"""



rule merge_breseq_and_snippy:
	input:
		breseq_vcf=rules.snpeff_for_breseq.output.vcf,
		snippy_vcf=rules.snippy.output.vcf,
		gbk=lambda wc: f"References/{SAMPLE_REF_NAME[wc.sample]}.gbk"
	output:
		tsv="work/variants/{sample}/snippy_breseq_variants.tsv"
	script:
		"../scripts/breseq_snippy_merge.py"

rule collect_breseq_snippy_tables:
	input:
		expand("work/variants/{sample}/snippy_breseq_variants.tsv", sample=SAMPLE_NAMES ),
	output:
		tsv="results/variant_table.tsv"
	params:
		samples=SAMPLE_NAMES,
		sample_ref_name=SAMPLE_REF_NAME
	script:
		"../scripts/collect_tables.py"


rule variant_stats:
	input:
		tsv=rules.collect_breseq_snippy_tables.output.tsv
	output:
		tsv="results/variant_statistics.tsv"
	params:
		sample_ref_name=SAMPLE_REF_NAME
	script:
		"../scripts/variants_stats.py"

rule summarize_mutants:
	input:
		expand("work/snippy/{sample}/snps.vcf", sample=SAMPLE_NAMES)
	output:
		summary="results/mutant_variant_summary.tsv"
	params:
		samples = SAMPLE_NAMES
	script:
		"../scripts/summarize_mutants.py"
		
		
rule collect_mutants_to_table:
	input:
		snps=expand("work/snippy/{sample}/snps.tab", sample=SAMPLE_NAMES),
	output:
		"results/all_variants.tsv"
	script:
		"../scripts/collect_mutants_to_table.py"