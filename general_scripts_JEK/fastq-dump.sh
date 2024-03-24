# purpose - download fastq files

# usage
# qq
#
# example
# purpose - qq

# usage
# qq
#
# example
#
# 'qq'

#set intermediate folder
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/fastq_files

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

mkdir -p $OUTPUT

JOBSTATS="#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=24gb
#PBS -l walltime=10:00:00
#PBS -A open
cd $OUTPUT

source ~/.bashrc #configures shell to use conda activate
conda activate bioinfo"

        sampleID=fastq_dump.pbs
        rm -f $sampleID
	echo "$JOBSTATS" >> $sampleID
	echo "#download each fastq file" >> $sampleID
	echo "fastq-dump --split-files SRR3211679" >> $sampleID
	echo "fastq-dump --split-files SRR3211680" >> $sampleID
	echo "fastq-dump --split-files SRR3211681" >> $sampleID
	echo "fastq-dump --split-files SRR3211682" >> $sampleID
	echo "fastq-dump --split-files SRR16815400" >> $sampleID
