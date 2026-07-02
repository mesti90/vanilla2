import pandas as pd

def concat_tables(samples, inputs, output, sample_ref_name):
	dfs = []

	for sample, table in zip(samples, inputs):
		df = pd.read_csv(table, sep="\t")

		df["Sample"] = sample
		df["Reference"] = sample_ref_name[sample]

		dfs.append(df)

	final = pd.concat(dfs, ignore_index=True)

	fixed_cols = ["Sample", "Reference"]
	columns = fixed_cols + [c for c in final.columns if c not in fixed_cols]
	final = final[columns]

	final.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
	concat_tables(
		snakemake.params.samples,
		snakemake.input,
		snakemake.output.tsv,
		snakemake.params.sample_ref_name
	)