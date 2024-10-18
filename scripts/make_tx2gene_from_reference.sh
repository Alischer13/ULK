#!/bin/bash

# Create tx2gene df from reference used when aligning
# Use csv file created from reference [https://nbisweden.github.io/workshop-RNAseq/2011/lab_kallisto.html]
# Ref (v112) from https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/

ref=../reference/Homo_sapiens.GRCh38.cdna.all.fa.gz
tx2gene=../reference/tx2gene.GRCh38.cdna.csv

gzcat ${ref} | grep '>' | awk '{FS= " "}BEGIN{ print "TXNAME,GENEID"};{print substr($1,2) "," substr($4,6)};' > ${tx2gene}


