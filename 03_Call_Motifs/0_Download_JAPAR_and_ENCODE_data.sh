#!/bin/bash

# Download JASPAR PWMs (.meme) and ENCODE narrowPeak (.bed) files for building
# "lowly bound" TF Motif RefPTs. Some of these downloads will not be used in
# larger analysis so they have been omitted from being saved in the global
# `../data/` directory.

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/03_Call_Motifs
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/03_Call_Motifs
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/03_Call_Motifs
METADATA=TF_JASPAR_ENCODE_config.txt
###

# Dependencies
# - wget

JDIR=../data/JASPAR

set -exo
[ -d $JDIR ] || mkdir $JDIR
[ -d narrowPeaks ] || mkdir narrowPeaks

while read line; do
	# Format input variables
	TF=`echo $line | awk '{print $1}'`
	JASPAR=`echo $line | awk '{print $2}'`
	ENCFF=`echo $line | awk '{print $3}'`

	# Download JASPAR motif
	wget -c -O $JDIR/$TF\_$JASPAR.meme https://jaspar.elixir.no/api/v1/matrix/$JASPAR.meme

	# Download ENCODE peaks
	wget -c -O narrowPeaks/$TF\_$ENCFF.bed.gz https://www.encodeproject.org/files/$ENCFF/@@download/$ENCFF.bed.gz

done < $METADATA
