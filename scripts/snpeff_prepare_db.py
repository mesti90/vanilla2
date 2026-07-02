import os
from pathlib import Path

input_gbk = snakemake.input.gbk
config_path = snakemake.output.config
gbk_out = snakemake.output.gbk
ref = snakemake.params.refname

base = Path(config_path).parent
data_dir = base / "data" / ref
data_dir.mkdir(parents=True, exist_ok=True)

# copy GBK
Path(gbk_out).write_bytes(Path(input_gbk).read_bytes())

# write config (ABSOLUTE PATH = important)
with open(config_path, "w") as f:
	f.write(f"data.dir = {base.resolve()}\n")
	f.write(f"ref.genome : ref\n")