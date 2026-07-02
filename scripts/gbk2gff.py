from Bio import SeqIO
from urllib.parse import quote

def format_attributes(qualifiers):
	attrs = []
	if "locus_tag" in qualifiers:
		attrs.append(f"ID={quote(qualifiers['locus_tag'][0], safe='')}")
	elif "gene" in qualifiers:
		attrs.append(f"ID={quote(qualifiers['gene'][0], safe='')}")
	if "gene" in qualifiers:
		attrs.append(f"Name={quote(qualifiers['gene'][0], safe='')}")
	if "product" in qualifiers:
		attrs.append(f"product={quote(qualifiers['product'][0], safe='')}")
	return ";".join(attrs)


with open(snakemake.output.gff, "w") as out:
	out.write("##gff-version 3\n")
	for record in SeqIO.parse(snakemake.input.gbk, "genbank"):
		for feature in record.features:

			if feature.type != "CDS":
				continue

			start = int(feature.location.start) + 1
			end = int(feature.location.end)

			if feature.location.strand == 1:
				strand = "+"
			elif feature.location.strand == -1:
				strand = "-"
			else:
				strand = "."
			out.write( "\t".join([ record.id, "Biopython", "CDS", str(start), str(end), ".", strand, "0", format_attributes(feature.qualifiers), ]) + "\n")