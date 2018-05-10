# sm_assemblies
Analysis of protein coding exons in single molecule assemblies

# dependencies
* Snakemake
* BioPython
* conda
* splign (https://www.ncbi.nlm.nih.gov/sutils/splign/splign.cgi?textpage=downloads)
* Python 3.5
* BLAT
* wget
* samtools

# clone this repo
```
git clone https://github.com/WatsonLab/sm_assemblies.git
cd sm_assemblies
```

# download genomes
/bin/bash scripts/download.sh
gunzip *.fasta.gz
mkdir genomes
mv *.fasta genomes

# run
```
snakemake --use-conda
```

