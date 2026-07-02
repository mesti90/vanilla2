#!/usr/bin/env python3
# pyright: reportUndefinedVariable=false

import pandas as pd
import pdb
import urllib.parse
from Bio import SeqIO

ANN_FIELDS = [ "allele", "effect", "impact", "gene", "gene_id", "feature_type", "feature_id", "biotype", "rank", "hgvs_c", "hgvs_p", "cdna", "cds", "aa", "distance", "warnings", ]

def parse_info(info):
	d = {}
	for item in info.split(";"):
		if "=" in item:
			k, v = item.split("=", 1)
			d[k] = v
	return d

def decode_snippy_field(x):
	if pd.isna(x):
		return x
	return urllib.parse.unquote(str(x))  # fixes %2C

def parse_best_ann(ann_string):
	if not ann_string:
		return {}
	anns = ann_string.split(",")
	parsed = []
	for ann in anns:
		parts = ann.split("|")
		if len(parts) < len(ANN_FIELDS):
			parts += [""] * (len(ANN_FIELDS) - len(parts))
		record = dict(zip(ANN_FIELDS, parts))
		parsed.append(record)
	# prefer non-MODIFIER
	for r in parsed:
		if r.get("impact") != "MODIFIER":
			return r
	return parsed[0] if parsed else {}


def variant_type(ref, alt):
	if len(ref) == 1 and len(alt) == 1:
		return "snp"
	elif len(ref) < len(alt):
		return "ins"
	elif len(ref) > len(alt):
		return "del"
	else:
		return "complex"

def load_gbk_annotations(gbk_path):
	ann = {}
	for record in SeqIO.parse(gbk_path, "genbank"):
		for feat in record.features:
			if "locus_tag" in feat.qualifiers:
				locus = feat.qualifiers["locus_tag"][0]

				product = feat.qualifiers.get("product", [""])[0]
				function = feat.qualifiers.get("function", [""])[0]

				ann[locus] = {
					"product": product,
					"function": function
				}
	return ann


def read_vcf(path, motor):
	rows = []
	template = {
		"CHROM": None,
		"POS": None,
		"REF": None,
		"ALT": None,
		"TYPE": None,
		f"{motor}_effect": None,
		f"{motor}_impact": None,
		f"{motor}_gene": None,
		f"{motor}_NUC": None,
		f"{motor}_PROT": None,
		f"in_{motor}": None,
	}


	with open(path) as f:
		for line in f:
			if line.startswith("#"):
				continue
			chrom, pos, _id, ref, alt, qual, filt, info = line.rstrip().split("\t")[:8]
			info_dict = parse_info(info)
			ann = parse_best_ann(info_dict.get("ANN", ""))
			
			gene = ann.get("gene")
			if motor == "snippy":
				gene = decode_snippy_field(gene).split(",")[0]

			rows.append({ "CHROM": chrom, "POS": int(pos), "REF": ref, "ALT": alt, "TYPE": variant_type(ref, alt), f"{motor}_effect": ann.get("effect"), f"{motor}_impact": ann.get("impact"), f"{motor}_gene": gene, f"{motor}_NUC": ann.get("hgvs_c"), f"{motor}_PROT": ann.get("hgvs_p"), f"in_{motor}": 1, })
	df = pd.DataFrame(rows)
	if df.empty:
		df = pd.DataFrame(columns=template.keys())
	return df

def resolve(a, b):
	if pd.isna(a) and pd.isna(b):
		return None
	if pd.isna(a):
		return b
	if pd.isna(b):
		return a
	if a == b:
		return a
	if a == b.split("-")[0]:
		return b
	return f"{a}|{b}"

def get_gbk_field(gbk_ann, locus, key):
	if pd.isna(locus):
		return None
	entry = gbk_ann.get(locus)
	if not entry:
		return None
	return entry.get(key)

def main():
	gbk_ann = load_gbk_annotations(snakemake.input.gbk)
	snippy = read_vcf(snakemake.input.snippy_vcf, "snippy")
	breseq = read_vcf(snakemake.input.breseq_vcf, "breseq")
	merged = pd.merge( snippy, breseq, on=["CHROM", "POS", "REF", "ALT", "TYPE"], how="outer", )
	
	merged["in_snippy"] = merged["in_snippy"].fillna(0).astype(int)
	merged["in_breseq"] = merged["in_breseq"].fillna(0).astype(int)
	merged["in_both"] = ((merged["in_snippy"] == 1) & (merged["in_breseq"] == 1)).astype(int)
	


	for field in ["effect", "impact", "gene",  "NUC", "PROT"]:
		s_col = f"snippy_{field}"
		b_col = f"breseq_{field}"
		u_col = field

		merged[u_col] = merged.apply(
			lambda r: resolve(r.get(s_col), r.get(b_col)),
			axis=1
		)
	merged = merged.drop(columns=[c for c in merged.columns if c.startswith("snippy_") or c.startswith("breseq_")])

	merged["product"] = merged["gene"].apply(lambda x: get_gbk_field(gbk_ann, x, "product"))
	merged["function"] = merged["gene"].apply(lambda x: get_gbk_field(gbk_ann, x, "function"))

	column_order = ["in_snippy", "in_breseq", "in_both", "CHROM", "POS", "REF", "ALT", "TYPE", "NUC", "PROT", "effect", "impact", "gene", "product", "function"]
	merged = merged[column_order]

	merged.to_csv(snakemake.output.tsv, sep="\t", index=False)


if __name__ == "__main__":
	main()
