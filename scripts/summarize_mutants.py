from pathlib import Path
import pandas as pd

rows = []

for path in snakemake.input:
	sample = Path(path).parts[-2]
	ref = snakemake.config["sample_table"].loc[sample, "Reference_name"]

	n = 0
	with open(path) as f:
		for line in f:
			if not line.startswith("#"):
				n += 1

	rows.append((sample, ref, n))

df = pd.DataFrame(rows, columns=["Sample", "Reference", "Variants"])
df.to_csv(snakemake.output.summary, sep="\t", index=False)