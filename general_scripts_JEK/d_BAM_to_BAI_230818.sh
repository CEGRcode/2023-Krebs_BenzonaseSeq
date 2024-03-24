# purpose - make index files (BAI) for all final BAM files
# usage
# qq
#
# example
#
# 'qq'

#input directory
INPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/final_files

#picard.jar
PICARD=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/picard.jar

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

JOBSTATS="#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=24gb
#PBS -l walltime=1:00:00
#PBS -A open
cd $INPUT

source ~/.bashrc #configures shell to use conda activate
conda activate bioinfo"

for file in $INPUT/*.bam; do
	#Get basename for PBS
	fileID=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
	echo $fileID
	#set output file names
	BAI=$(echo $file | rev | cut -d"/" -f1 | rev | awk '{print $1".bai"}')

        sampleID=$fileID\_index.pbs
        rm -f $sampleID
        echo "$JOBSTATS" >> $sampleID
        echo "#index input bam file" >> $sampleID
	echo "java -jar $PICARD BuildBamIndex --INPUT $file --OUTPUT $INPUT/$BAI" >> $sampleID
	echo "#finish script" >> $sampleID

done
