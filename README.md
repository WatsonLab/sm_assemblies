# sm_assemblies
Analysis of protein coding exons in single molecule assemblies

# dependencies
* Snakemake
* Python 3.5
* BioPython
* BLAT
* wget

# what it does
* downloads human primary assembly from Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
* downloads PacBio-only human assembly CA_001013985.1_ASM101398v1
* downloads Nanopore+Illumina human assembly GCA_900232925.1_Nanopore-only_assembly_with_Illumina_polishing
* downloads protein coding exons from BioMart and filters >300bp
* BLAT protein coding exons against each assembly, taking the top hit into a tsv report
