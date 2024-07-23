#!/bin/bash

GENOME=../data/hg19_files/hg19.fa

# Download genome FASTA
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv hg19.fa.gz $GENOME.gz
gzip -d $GENOME.gz

# Create genome indexes
samtools faidx $GENOME
bowtie2-build $GENOME $GENOME
bowtie-build -C $GENOME $GENOME.colorspace
bwa index $GENOME

# Download RefSeq annotations