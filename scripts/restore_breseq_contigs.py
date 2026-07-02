from Bio import SeqIO
import pdb
# snakemake.input.vcf
# snakemake.input.ref
# snakemake.output.vcf

sanitized_to_original = {}

for record in SeqIO.parse(snakemake.input.ref, "genbank"):
	original = record.id
	sanitized = original.replace("|", "_")
	
	if sanitized in sanitized_to_original:
		raise ValueError( f"Contig name collision after sanitization: {sanitized}" )
	sanitized_to_original[sanitized] = original


with open(snakemake.input.vcf) as fin, open(snakemake.output.vcf, "w") as fout:
	for line in fin:

		if line.startswith("##contig=<ID="):
			prefix = "##contig=<ID="
			rest = line[len(prefix):]
			contig = rest.split(",", 1)[0].rstrip(">")

			if contig in sanitized_to_original:
				line = line.replace(
					f"ID={contig}",
					f"ID={sanitized_to_original[contig]}",
					1,
				)

		elif not line.startswith("#"):
			fields = line.rstrip("\n").split("\t")

			if fields[0] in sanitized_to_original:
				fields[0] = sanitized_to_original[fields[0]]

			line = "\t".join(fields) + "\n"

		fout.write(line)