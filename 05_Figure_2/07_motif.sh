# purpose - run MEME on lowly-bound ZKSCAN1 sites to determine if there is an extnded motif. If there is, get motif.
# usage
# qq
#
# example
#
# 'qq'

#set bedfiles
LOWLY_BOUND=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/04_DNAshape_output/MA1585_1_final_1000bp_intersected_lowlyBound_1000bp.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/extended_ZKSCAN1_motif_231211/01_output

#set scriptmanager
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar

#set genome and human.hg19.genome file
GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/referenceDATA_Will/GENOMES/hg19.fa
HG19_GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/files/human.hg19.genome

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

mkdir -p $OUTPUT

JOBSTATS="#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=24GB
#SBATCH --time=2:00:00
#SBATCH --partition=open

source ~/.bashrc #configures shell to use conda activate
conda activate bioinfo"

#Get basename for PBS
fileID=$(echo $LOWLY_BOUND | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
LOWLY_BOUND_FASTA=$(echo $LOWLY_BOUND | rev | cut -d"/" -f1 | rev | awk -F. '{print $1".fa"}')
MEME=$(echo "meme.txt")
MEME_motif1=$(echo $LOWLY_BOUND | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_motif1.meme"}')


	sampleID=$fileID\_script01_extMOTIF_ZKSCAN1_v1_final_231212\.slurm
        rm -f $sampleID
	echo "$JOBSTATS" >> $sampleID
	echo "#set output file names" >> $sampleID
	echo "cd $OUTPUT" >> $sampleID
	echo "#get fasta seq of each DNA fragment." >> $sampleID 
	echo "java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME $LOWLY_BOUND -o=$LOWLY_BOUND_FASTA" >> $sampleID
	echo "#run MEME" >> $sampleID
	echo "meme $LOWLY_BOUND_FASTA -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 6 -maxw 100 -objfun classic -revcomp -markov_order 0" >> $sampleID
	echo "#extract MEME file from MEME.txt file" >> $sampleID
	echo "meme-get-motif -id RVACARTGCCTGGCACAYAGTAGGTGCTCARTAAATRT $OUTPUT/$MEME > $MEME_motif1" >> $sampleID
	echo "# finish script" >> $sampleID
