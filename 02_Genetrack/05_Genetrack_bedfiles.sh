# purpose - make redundant bedfile
# usage
# qq
#
# example
#
# 'qq'

#input nucleosomal peak bedfiles
NUCLEOSOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES/K562_benzonase-seq_master_MIDPOINT_128-164_NUC_s40e80F5_nucleosomes_midpoint_164bp_final.bed
HEXASOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES/K562_benzonase-seq_master_MIDPOINT_92-127_NUC_s30e60F6_hexasomes_midpoint_127bp_final.bed
TETRASOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES/K562_benzonase-seq_master_MIDPOINT_55-91_NUC_s20e40F5_tetrasomes_midpoint_91bp_final.bed
SUBTETRASOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES/K562_benzonase-seq_master_MIDPOINT_0-54_NUC_s10e20_subtetrasomes_midpoint_54bp_final.bed
SUPRAOCTASOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES/K562_benzonase-seq_master_MIDPOINT_165-Inf_NUC_s50e100F3_supraoctasomes_midpoint_200bp_final.bed

#output directory
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/02_intermediate_FILES

#set scriptmanager
SCRIPTMANAGER=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar

#set blacklist
BLACKLIST=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/hg19_Blacklist.bed

#set human.hg19.genome file
HG19_GENOME=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230705_masterBedfile/files/human.hg19.genome

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

#set output file names
SUBTETRASOMES_BEDFILE=$(echo $SUBTETRASOMES | rev | cut -d"/" -f1 | rev | awk '{gsub(/.gff/,"_subtetrasomes.bed"); print}')
BEDFILE1=$(echo $SUBTETRASOMES_BEDFILE_10K | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$BEDFILE1" | awk -F. '{print $1"_midpoint.out"}')

#set output directory
cd $OUTPUT

#bedtools intersect to get unique hexasomes (those that do NOT overlap nucleosomes)
bedtools intersect -v -a $HEXASOMES -b $NUCLEOSOMES > K562_uHex.bed

#concatenate nucleosomes and unique hexasomes (not nucleosomes) files together
cat $NUCLEOSOMES K562_uHex.bed > K562_nuc_uHex.bed

#bedtools intersect to get unique tetrasomes (those that do NOT overlap nucleosomes or hexasomes)
bedtools intersect -v -a $TETRASOMES -b K562_nuc_uHex.bed > K562_uTetra.bed

#concatenate nucleosomes, unique hexasomes (not nucleosomes), and unique tetrasomes (not nucleosomes, hexasomes) files together
cat K562_nuc_uHex.bed K562_uTetra.bed > K562_nuc_uHex_uTetra.bed

#bedtools intersect to get unique subtetrasomes (those that do NOT overlap nucleosomes, hexasomes, or tetrasomes)
bedtools intersect -v -a $SUBTETRASOMES -b K562_nuc_uHex_uTetra.bed > K562_uSubtetra.bed

#concatenate nucleosomes, unique hexasomes (not nucleosomes), unique tetrasomes (not nucleosomes, hexasomes), and unique subtetrasomes (not nucleosomes, hexasomes, tetrasomes) files together
cat K562_nuc_uHex_uTetra.bed K562_uSubtetra.bed > K562_nuc_uHex_uTetra_uSubtetra.bed

#bedtools intersect to get unique supraoctasomes (those that do NOT overlap nucleosomes, hexasomes, tetrasomes, or subtetrasomes)
bedtools intersect -v -a $SUPRAOCTASOMES -b K562_nuc_uHex_uTetra_uSubtetra.bed > K562_uSupraoct.bed

#concatenate nucleosomes, unique hexasomes (not nucleosomes), unique tetrasomes (not nucleosomes, hexasomes), unique subtetrasomes (not nucleosomes, hexasomes, tetrasomes), and unique supraoctasomes (not nucleosomes, hexasomes, tetrasomes, or subtetrasomes) files together
cat K562_nuc_uHex_uTetra_uSubtetra.bed K562_uSupraoct.bed > K562_master_nonredundant_particles.bed
cat K562_master_nonredundant_particles.bed | wc -l



# finish script
echo "DONE"
