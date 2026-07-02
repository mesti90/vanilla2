#Source:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_protein.faa.gz
gunzip GCF_000240185.1_ASM24018v2_protein.faa.gz
diamond makedb --in GCF_000240185.1_ASM24018v2_protein.faa --db GCF_000240185.1.protein.dmnd
