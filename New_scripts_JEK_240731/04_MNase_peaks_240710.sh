# purpose - prepare bedfile for all nucleosome peak gf files and merge into master bedfile; NOW (as of 1b) bedfiles are expanded by maximum length of peak
# usage
# qq
#
# example
#
# 'qq'

#input gff files of replicate peak calls
MNase=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/02_scIDX-NUC_output/genetrack_s40e80F32/final_shifted_scIDX_MNase_NUC_s40e80F32.gff

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/03_Genetrack_peaks_v2

#SHIFT_CHECK=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES/shift_check

#set scriptmanager
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar

#set blacklist and .genome file
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed
HG19_GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/files/human.hg19.genome

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

#set output file names
MNase_BEDFILE=$(echo $MNase | rev | cut -d"/" -f1 | rev | awk '{gsub(/.gff/,"_MNase.bed"); print}')
MNase_BEDFILE_1bp=$(echo $MNase_BEDFILE | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_1bp.bed"); print}')
MNase_BEDFILE_MIDPOINT=$(echo $MNase_BEDFILE_1bp | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_midpoint.bed"); print}')
MNase_BEDFILE_EXPANDED=$(echo $MNase_BEDFILE_MIDPOINT | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_164bp.bed"); print}')
MNase_BEDFILE_EXPANDED_FIX=$(echo $MNase_BEDFILE_EXPANDED | rev | cut -d"/" -f1 | rev | awk '{gsub(/_NUC_s40e80F32_MNase_1bp_midpoint_164bp.bed/,"_164bp_final.bed"); print}')

#set output directory
cd $OUTPUT

#convert gff to bedfiles
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed -s $MNase > $MNase_BEDFILE

#set output directory
cd $OUTPUT

#expand all bedfile to 1 bp so that midpoint of each fragment is designated
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1 $MNase_BEDFILE

#put midpoint of fragment in column 4 with st. dev.
cat $MNase_BEDFILE_1bp | awk '{if ($6=="+") print $1"\t"$2"\t"$3"\t"$4",midpoint="$2",MNase\t"$5"\t"$6; else print ""}' > $MNase_BEDFILE_MIDPOINT

#ensure that final start / end match previous bedfiles

#expand bedfile by maximum DNA length for each indicated nucleosomal peak; supraoctasomes here will be expanded to 200 bp, there is no maximum for them. 
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=164 $MNase_BEDFILE_MIDPOINT

#Then remove any lines in expanded bedfile with a region expanded past a chromosome end (as seen by a "-" in column 2 or 3.
cat $MNase_BEDFILE_EXPANDED | awk '{if ($2>="1") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' | awk '{if ($3>="1") print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $MNase_BEDFILE_EXPANDED_FIX

#concatenate above files to master bedfile
cat $MNase_BEDFILE_EXPANDED_FIX | wc -l

# finish script
echo "DONE"
