#Source:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_protein.faa.gz
gunzip GCF_000013425.1_ASM1342v1_protein.faa.gz
diamond makedb --in GCF_000013425.1_ASM1342v1_protein.faa --db GCF_000013425.1.protein.dmnd

