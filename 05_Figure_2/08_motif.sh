# purpose - run MEME on lowly bound fasta file
# usage
# qq
#
# example
#
# 'qq'

#set bedfiles
LOWLY_BOUND_FASTA=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/extended_ZKSCAN1_motif_231211/01_output/MA1585_1_final_1000bp_intersected_lowlyBound_1000bp.fa
MEME_motif1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/extended_ZKSCAN1_motif_231211/01_output/MA1585_1_final_1000bp_intersected_lowlyBound_1000bp_motif1.meme

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/extended_ZKSCAN1_motif_231211/02_output

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
#SBATCH --time=0:15:00
#SBATCH --partition=open

source ~/.bashrc #configures shell to use conda activate
conda activate bioinfo"

#Get basename for PBS
fileID=$(echo $LOWLY_BOUND_FASTA | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')

	sampleID=$fileID\_script02_extMOTIF_ZKSCAN1_v1_same_sites_231212\.slurm
        rm -f $sampleID
	echo "$JOBSTATS" >> $sampleID
	echo "#set output file names" >> $sampleID
	echo "cd $OUTPUT" >> $sampleID
	echo "#run FIMO with extended motif in MEME file on LOWLy bound fasta" >> $sampleID
	echo "fimo --oc . --verbosity 1 --bgfile --nrdb-- --thresh 1.0E-4 --norc $MEME_motif1 $LOWLY_BOUND_FASTA" >> $sampleID
	echo "# finish script" >> $sampleID
