import os

configfile: "config/config.yaml"

include: "rules/common.smk"
include: "rules/all_rules.smk"

rule all:
	input:
		rules.collect_breseq_snippy_tables.output.tsv,
		rules.variant_stats.output.tsv
