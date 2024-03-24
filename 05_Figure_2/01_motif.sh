# purpose - take PWMs to make gff files and then to bedfiles. Remove rows from blacklist filter.

# usage
# qq
#
# example
#
# 'qq'

#set input and output directory
MOTIF_DATABASE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/JASPAR2022_CORE_verebrates_non-redundant_pfms_meme
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/final_bedfiles

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

for file in $MOTIF_DATABASE/*.meme; do
	#Get basename for PBS
	fileID_1=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
	fileID_2=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $2}')	
	fileID=$(echo "$fileID_1""_""$fileID_2" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
	echo $fileID
	#set output file names
	PWMa=$(echo $fileID | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
	EPS1=$(echo "$PWMa" | awk -F. '{print $1".eps"}')
	PWM_FIMO=$(echo "$OUTPUT""/""$PWMa""_FIMO")
	PWM_FIMO_gff=$(echo "$OUTPUT""/""$PWMa""_FIMO""/""fimo.gff")
	PWM_FIMO_gff2=$(echo "$OUTPUT""/""$PWMa""_FIMO""/""fimo_noHeader.gff")
	PWM_bedfile=$(echo "$PWMa" | awk -F. '{print $1".bed"}')
	PWM_bedfile2=$(echo "$PWMa" | awk -F. '{print $1"_1000bp.bed"}')
	FINAL_BEDFILE=$(echo "$PWMa" | awk -F. '{print $1"_final_1000bp.bed"}')

	sampleID=$fileID\_motif_script01\.slurm
        rm -f $sampleID
	echo "$JOBSTATS" >> $sampleID
	echo "#set output directory" >> $sampleID
	echo "cd $OUTPUT" >> $sampleID
	echo "#make logos of original meme file" >> $sampleID
	echo "ceqlogo -i $file -m 1 > $OUTPUT/$EPS1" >> $sampleID
	echo "#get sites in hg19 with original motifswith FIMO" >> $sampleID
	echo "fimo --o $PWM_FIMO --verbosity 1 --bfile --motif-- --thresh 1.0E-4 $file $GENOME" >> $sampleID
	echo "#remove first line of gff file" >> $sampleID
	echo "sed '1d' $PWM_FIMO_gff > $PWM_FIMO_gff2" >> $sampleID
	echo "#convert gff to bedfile" >> $sampleID
	echo "java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed -o=$OUTPUT/$PWM_bedfile $PWM_FIMO_gff2" >> $sampleID
	echo "#expand bedfile by 1000 bp. This matches the max legnth of peaks in script 02" >> $sampleID
	echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $OUTPUT/$PWM_bedfile -o=$OUTPUT/$PWM_bedfile2" >> $sampleID
	echo "#remove blacklisted regions" >> $sampleID
	echo "bedtools intersect -v -a $OUTPUT/$PWM_bedfile2 -b $BLACKLIST -bed > $OUTPUT/$FINAL_BEDFILE" >> $sampleID
	echo "rmdir $OUTPUT/""$PWMa""_FIMO""" >> $sampleID
	echo "rm $PWM_bedfile" >> $sampleID
	echo "rm $PWM_bedfile2" >> $sampleID
	echo "# finish script" >> $sampleID
done
