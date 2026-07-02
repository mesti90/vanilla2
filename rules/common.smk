import pandas as pd
import os
from pathlib import Path

def abspath(x):
	return os.path.abspath(x)


#sample layer

SAMPLES = pd.read_csv(config["samples"], sep="\t")
SAMPLES = SAMPLES.set_index("Sample")
SAMPLE_NAMES = list(SAMPLES.index)

R1 = SAMPLES["R1"].to_dict()
R2 = SAMPLES["R2"].to_dict()


#reference layer

REFERENCES =  SAMPLES[ [ "Reference_name", "Reference_fasta", "Reference_gbk", ] ].drop_duplicates().set_index('Reference_name')

REF_NAMES = list(REFERENCES.index)
REF_GBK = REFERENCES["Reference_gbk"].to_dict()
REF_FASTA = REFERENCES["Reference_fasta"].to_dict()

#mixed layer

PAIRS = SAMPLES.reset_index()[["Sample", "Reference_name"]]
SAMPLE_PAIRS = list(zip(PAIRS["Sample"], PAIRS["Reference_name"]))


#helper
SAMPLE_REF_NAME = SAMPLES["Reference_name"].to_dict()

