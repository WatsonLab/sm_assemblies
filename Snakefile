shell.executable("/bin/bash")
shell.prefix("source $HOME/.bashrc; ")

from Bio import SeqIO
import sys
import os

rule all:
	input: "blat_grch38.tsv", "blat_pacbio.tsv", "blat_nanopore.tsv", "blat_chm1.round2.psl", "blat_chm1.round2.tsv"

rule grch38_download:
	output: "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
	shell: "wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

rule pacbio_download:
        output: "pacbio.fasta.gz"
        shell: "wget -O {output} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/013/985/GCA_001013985.1_ASM101398v1/GCA_001013985.1_ASM101398v1_genomic.fna.gz"

rule nanopore_download:
        output: "nanopore.fasta.gz"
        shell: "wget -O {output} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/232/925/GCA_900232925.1_Nanopore-only_assembly_with_Illumina_polishing/GCA_900232925.1_Nanopore-only_assembly_with_Illumina_polishing_genomic.fna.gz"

rule canu_chm_quiver_download:
	output: "chm1.round2.fasta"
	shell: "wget http://gembox.cbcb.umd.edu/shared/canu/quiver/canu/chm1.round2.fasta"

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


rule make_ooc_nanopore:
	input:"nanopore.fasta.gz"
	output: "nanopore.ooc"
	shell: "blat {input} /dev/null /dev/null -makeOoc={output} -repMatch=1024"

rule make_ooc_pacbio:
	input:"pacbio.fasta.gz"
	output: "pacbio.ooc"
	shell: "blat {input} /dev/null /dev/null -makeOoc={output} -repMatch=1024"

rule make_ooc_grch38:
	input:"Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
	output: "grch38.ooc"
	shell: "blat {input} /dev/null /dev/null -makeOoc={output} -repMatch=1024"

rule make_ooc_canu_chm_quiver:
	input: "chm1.round2.fasta"
	output: "chm1.round2.ooc"
	shell: "blat {input} /dev/null /dev/null -makeOoc={output} -repMatch=1024"

rule blat_grch38:
	input: 
		ooc="grch38.ooc",
		exn="protein_coding_exons.filtered.fasta"
	output: "blat_grch38.psl"
	shell: "blat Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz {input.exn} -out=blast9 -ooc={input.ooc} {output}"

rule blat_nanopore:
	input: 
		ooc="nanopore.ooc",
		exn="protein_coding_exons.filtered.fasta"
	output: "blat_nanopore.psl"
	shell: "blat nanopore.fasta.gz {input.exn} -out=blast9 -ooc={input.ooc} {output}"

rule blat_pacbio:
	input: 
		ooc="pacbio.ooc",
		exn="protein_coding_exons.filtered.fasta"
	output: "blat_pacbio.psl"
	shell: "blat pacbio.fasta.gz {input.exn} -out=blast9 -ooc={input.ooc} {output}"

rule blat_canu_chm_quiver:
	input:
		ooc="chm1.round2.ooc",
		exn="protein_coding_exons.filtered.fasta"
	output: "blat_chm1.round2.psl"
	shell: "blat chm1.round2.fasta {input.exn} -out=blast9 -ooc={input.ooc} {output}"

rule report_grch38:
	input: "blat_grch38.psl", "protein_coding_exons.filtered.fasta"
	output: "blat_grch38.tsv"
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

        		print("%s\t%s\t%s\t%s\t%s\t%s" % (arr[0], seq_length[arr[0]], arr[2], arr[3], arr[4], arr[5]), file=f)
		
		f.close()

rule report_nanopore:
        input: "blat_nanopore.psl", "protein_coding_exons.filtered.fasta"
        output: "blat_nanopore.tsv"
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

                        print("%s\t%s\t%s\t%s\t%s\t%s" % (arr[0], seq_length[arr[0]], arr[2], arr[3], arr[4], arr[5]), file=f)

                f.close()

rule report_pacbio:
        input: "blat_pacbio.psl", "protein_coding_exons.filtered.fasta"
        output: "blat_pacbio.tsv"
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

                        print("%s\t%s\t%s\t%s\t%s\t%s" % (arr[0], seq_length[arr[0]], arr[2], arr[3], arr[4], arr[5]), file=f)

                f.close()

rule report_chm1_canu_quiver:
	input: "blat_chm1.round2.psl", "protein_coding_exons.filtered.fasta"
	output: "blat_chm1.round2.tsv"
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

			print("%s\t%s\t%s\t%s\t%s\t%s" % (arr[0], seq_length[arr[0]], arr[2], arr[3], arr[4], arr[5]), file=f)

		f.close()

