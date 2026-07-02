import pandas as pd

df = pd.read_csv(snakemake.input.tsv, sep="\t")

df["Reference"] = df["Sample"].map(snakemake.params.sample_ref_name)

df["snippy"] = df["in_snippy"].eq(1)
df["breseq"] = df["in_breseq"].eq(1)

df["both"] = df["snippy"] & df["breseq"]

# expand all combinations in one shot
df["snippy_snp"]   = df.snippy & df.TYPE.eq("snp")
df["breseq_snp"]   = df.breseq & df.TYPE.eq("snp")
df["both_snp"]     = df.both   & df.TYPE.eq("snp")

df["snippy_complex"]     = df.snippy & df.TYPE.eq("complex")
df["breseq_complex"]     = df.breseq & df.TYPE.eq("complex")
df["both_complex"]       = df.both   & df.TYPE.eq("complex")

df["snippy_ins"]     = df.snippy & df.TYPE.eq("ins")
df["breseq_ins"]     = df.breseq & df.TYPE.eq("ins")
df["both_ins"]       = df.both   & df.TYPE.eq("ins")

df["snippy_del"]     = df.snippy & df.TYPE.eq("del")
df["breseq_del"]     = df.breseq & df.TYPE.eq("del")
df["both_del"]       = df.both   & df.TYPE.eq("del")

agg = df.groupby("Sample")[[
    "snippy_snp","breseq_snp","both_snp",
    "snippy_complex","breseq_complex","both_complex",
    "snippy_ins","breseq_ins","both_ins",
    "snippy_del","breseq_del","both_del",
    "snippy","breseq","both"
]].sum().reset_index()

agg.to_csv(snakemake.output.tsv, sep="\t", index=False)