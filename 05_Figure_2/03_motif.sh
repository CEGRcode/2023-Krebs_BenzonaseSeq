# purpose - make 4-color plot of each bedfile for each motif to confirm position of motif center for Fig. 5.

# usage
# qq
#
# example
#
# 'qq'

#set input and output directory
INPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/final_bedfiles
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/03_output
#set blacklist, genome and .genome file
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed
GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/referenceDATA_Will/GENOMES/hg19.fa
HG19_GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/files/human.hg19.genome

#set scriptmanager and job
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

JOBSTATS="#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=24GB
#SBATCH --time=0:30:00
#SBATCH --partition=open

source ~/.bashrc #configures shell to use conda activate
module load anaconda
conda activate bioinfo"


for file in $INPUT/*_final_1000bp.bed; do
	#Get basename for PBS
	fileID=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
	echo $fileID
	#set output file names
	BEDFILE_3bp=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_3bp.bed"}')
	BEDFILE_20bp=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_20bp.bed"}')
	BEDFILE_3bp_shuffled=$(echo $BEDFILE_3bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled.bed"}')
	BEDFILE_20bp_shuffled=$(echo $BEDFILE_20bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled.bed"}')
	BEDFILE_3bp_50sites=$(echo $BEDFILE_3bp_shuffled | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_50sites.bed"}')
	BEDFILE_20bp_50sites=$(echo $BEDFILE_20bp_shuffled | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_50sites.bed"}')
	FASTA_3bp=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_3bp_50sites.fasta"}')
	FASTA_20bp=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_20bp_50sites.fasta"}')
	PNG_3bp=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_3bp_4color.png"}')
	PNG_20bp=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_20bp_4color.png"}')

	sampleID=$fileID\_motif_script03_v4\.slurm
        rm -f $sampleID
	echo "$JOBSTATS" >> $sampleID
	echo "#set output directory" >> $sampleID
	echo "cd $OUTPUT" >> $sampleID
	echo "#expand bedfile by 3 bp." >> $sampleID
	echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=3 $file -o=$BEDFILE_3bp" >> $sampleID
	echo "#expand bedfile by 20 bp." >> $sampleID
	echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=20 $file -o=$BEDFILE_20bp" >> $sampleID
	echo "#shuffle bedfiles" >> $sampleID
	echo "shuf $BEDFILE_3bp > $BEDFILE_3bp_shuffled" >> $sampleID
	echo "shuf $BEDFILE_20bp > $BEDFILE_20bp_shuffled" >> $sampleID
	echo "#take top 50 sites" >> $sampleID
	echo "cat $BEDFILE_3bp_shuffled | head -50 > $BEDFILE_3bp_50sites" >> $sampleID
	echo "cat $BEDFILE_20bp_shuffled | head -50 > $BEDFILE_20bp_50sites" >> $sampleID
	echo "#extract fasta. Default is bedfile header." >> $sampleID
	echo "java -jar $SCRIPTMANAGER sequence-analysis fasta-extract -o=$FASTA_3bp $GENOME $BEDFILE_3bp_50sites" >> $sampleID
	echo "java -jar $SCRIPTMANAGER sequence-analysis fasta-extract -o=$FASTA_20bp $GENOME $BEDFILE_20bp_50sites" >> $sampleID
	echo "#make four color plot" >> $sampleID
	echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation four-color --output=$PNG_3bp $FASTA_3bp" >> $sampleID
	echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation four-color --output=$PNG_20bp $FASTA_20bp" >> $sampleID
	echo "#remove intermediate files" >> $sampleID
	echo "rm $BEDFILE_3bp" >> $sampleID
	echo "rm $BEDFILE_20bp" >> $sampleID
	echo "rm $BEDFILE_3bp_shuffled" >> $sampleID
	echo "rm $BEDFILE_20bp_shuffled" >> $sampleID
	echo "rm $BEDFILE_3bp_50sites" >> $sampleID
	echo "rm $BEDFILE_20bp_50sites" >> $sampleID
	echo "#rm $FASTA_3bp" >> $sampleID
	echo "#rm $FASTA_20bp" >> $sampleID
	echo "# finish script" >> $sampleID
done
