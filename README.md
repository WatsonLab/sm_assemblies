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
```
/bin/bash scripts/download.sh
gunzip *.fasta.gz
mkdir genomes
mv *.fasta genomes
```

# run
```
snakemake --use-conda
```

## parse splign output

Perl script alnparse.pl can be used to summarise the splign output.  

The way the pipeline stores splign results is in a "one query, one subject" file i.e. "ENST00000052569.10.OCVW01001666.1.aln" - this would be all of the splign hits from ENST00000052569.10 (query) against OCVW01001666.1 (subject).

To summarise the best hit from a single alignment file, run the script like this:

```sh
perl scripts/alnparse.pl ENST00000052569.10.OCVW01001666.1.aln
```

However, often we want to consider the hits against multiple subjects in order to find the best hit.  In this case, we run it like this:

```sh
perl scripts/alnparse.pl <(cat ENST00000052569.10.*.aln)
```

This will find the best hit from alignments of ENST00000052569.10 against all subjects it has been aligned against.

The output is tab-delimited:
* query name
* hit name
* query start
* length of alignment
* number of mismatch events
* number of bases in mismatch events
* number of insertion events
* number of bases in insertion events
* number of deletion events
* number of bases in deletion events
* protein sequence of aligned bases


