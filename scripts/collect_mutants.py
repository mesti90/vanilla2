from pathlib import Path
import pandas as pd

dfs = []

for snpfile in snakemake.input.snps:

	sample = Path(snpfile).parts[-2]

	df = pd.read_csv( snpfile, sep="\t", comment="#", dtype=str )

	if df.empty:
		dfs.append(pd.DataFrame([{ "Sample": sample }]))
	else:
		df.insert(0, "Sample", sample)
		dfs.append(df)

final = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
final = final.fillna("-")

final.to_csv(snakemake.output[0], sep="\t", index=False)