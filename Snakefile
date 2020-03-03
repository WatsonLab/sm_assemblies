shell.executable("/bin/bash")
shell.prefix("source $HOME/.bashrc; ")

from Bio import SeqIO
import sys
import os

IDS, = glob_wildcards("genomes/{id}.fasta")

rule all:
	input: expand("{sample}.splign.finished.txt", sample=IDS), expand("genomes/{sample}.fasta.fai", sample=IDS), expand("{sample}.transcript_ids.txt", sample=IDS), expand("blat_{sample}.tsv", sample=IDS), expand("blat_{sample}.blast", sample=IDS), "unique_transcript_ids.txt"

rule coding_transcripts_download:
	output: "protein_coding_transcripts.fasta"
	shell:
		'''
		wget -O protein_coding_transcripts.fasta.gz ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
		gunzip protein_coding_transcripts.fasta.gz
		samtools faidx protein_coding_transcripts.fasta
		'''

rule exons_download:
        output: "protein_coding_exons.fasta"
        shell: "wget -O protein_coding_exons.fasta 'http://www.ensembl.org/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"FASTA\" header = \"0\" uniqueRows = \"1\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"hsapiens_gene_ensembl\" interface = \"default\" ><Filter name = \"biotype\" value = \"protein_coding\"/><Filter name = \"transcript_biotype\" value = \"protein_coding\"/><Attribute name = \"ensembl_gene_id\" /><Attribute name = \"ensembl_transcript_id\" /><Attribute name = \"gene_exon\" /><Attribute name = \"ensembl_exon_id\" /></Dataset></Query>'"

rule length_filter:
        input: "protein_coding_exons.fasta"
        output: "protein_coding_exons.filtered.fasta"
        params:
                len=300
        run:
                input_seq_iterator = SeqIO.parse(open(input[0], "rU"), "fasta")
                short_seq_iterator = (record for record in input_seq_iterator \
                        if len(record.seq) > int(300))

                output_handle = open(output[0], "w")
                SeqIO.write(short_seq_iterator, output_handle, "fasta")
                output_handle.close()


rule make_faidx:
	input: "genomes/{id}.fasta"
	output: "genomes/{id}.fasta.fai"
	shell: "samtools faidx {input}"

rule splign:
	input:
		idl="trouble_list.txt",
		blo="blat_{id}.blast",
		ass="genomes/{id}.fasta",
		fai="genomes/{id}.fasta.fai"
	output: "{id}.splign.finished.txt"
	params: dir="{id}.splign"
	conda: "envs/bioperl.yml"
	shell:
		"""
		perl scripts/run_splign.pl {input.idl} {input.blo} {input.ass} {params.dir} && touch {output}
		"""

rule make_ooc:
	input: "genomes/{id}.fasta"
	output: "{id}.ooc"
	shell: "blat {input} /dev/null /dev/null -makeOoc={output} -repMatch=1024"

rule blat:
	input:
		gen="genomes/{id}.fasta",
		ooc="{id}.ooc",
		exn="protein_coding_exons.filtered.fasta"
	output: "blat_{id}.psl"
	shell: "blat {input.gen} {input.exn} -out=blast9 -ooc={input.ooc} {output}"

rule blat_cdna:
	input:
		gen="genomes/{id}.fasta",
		ooc="{id}.ooc",
		cdn="protein_coding_transcripts.fasta"
	output: "blat_{id}.blast"
	shell: "blat {input.gen} {input.cdn} -out=blast -ooc={input.ooc} {output}"

rule unique_tid:
	input: "protein_coding_exons.filtered.fasta"
	output: "unique_transcript_ids.txt"
	shell: "cat {input} | grep '>' | awk -F'|' '{{print $2}}' | perl -e 'while(<>) {{print join(\"\\n\", split(\";\"))}}' | sort | uniq > {output}"

rule get_list:
	input: 
		uti="unique_transcript_ids.txt",
		rep="blat_{id}.tsv"
	output: "{id}.transcript_ids.txt"
	params:
		tmp="{id}"
	shell:
		"""
		cat {input.rep} | awk '$6>0' | awk -F'|' '{{print $2}}' | perl -e 'while(<>) {{print join(\"\\n\", split(\";\"))}}' | sort | uniq > {params.tmp}.faulty
		cat {input.rep} | awk -F'|' '{{print $2}}' | perl -e 'while(<>) {{print join(\"\\n\", split(\";\"))}}' | sort | uniq > {params.tmp}.allhits
		comm {params.tmp}.allhits {input.uti} | awk -F\"\\t\" '$2~/ENST/ {{print $2}}' > {params.tmp}.missing
		cat {params.tmp}.faulty {params.tmp}.missing > {output}
		"""

rule combine_lists:
	input: expand("{allids}.transcript_ids.txt", allids=IDS)
	output: "trouble_list.txt"
	shell: "cat {input} | sort | uniq > {output}"

rule report:
	input: "blat_{id}.psl", "protein_coding_exons.filtered.fasta"
	output: "blat_{id}.tsv"
	run:
		seq_length = dict()

 		with open(input[1], "rU") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				seq_length[record.id] = len(record.seq)

		# open the input file
		psl_file = open(input[0], mode="r")

		top_hits = dict()

		# open the output file
		f = open(output[0], 'w')

		# iterate over file
		for row in psl_file:

			if row.startswith("#"):
				continue
  
			# split on whitespace
			arr = row.strip().split()

			if arr[0] in top_hits:
				continue

			top_hits[arr[0]] = 1
			
			if arr[0] in seq_length.keys():
				print("%s\t%s\t%s\t%s\t%s\t%s" % (arr[0], seq_length[arr[0]], arr[2], arr[3], arr[4], arr[5]), file=f)
			else:
				print("This key isn't in seq_length: ", arr[0], end='\n\n')

		f.close()

