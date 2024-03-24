# purpose - prepare bedfile for all nucleosome peak gf files and merge into master bedfile; NOW (as of 1b) bedfiles are expanded by maximum length of peak
# usage
# qq
#
# example
#
# 'qq'

#input gff files of replicate peak calls
SUBTETRASOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230719_BI_Genetrack/scIDX-Nuc/genetrack_s10e20/K562_benzonase-seq_master_MIDPOINT_0-54_NUC_s10e20.gff
TETRASOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230719_BI_Genetrack/scIDX-Nuc/genetrack_s20e40F5/K562_benzonase-seq_master_MIDPOINT_55-91_NUC_s20e40F5.gff
HEXASOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230719_BI_Genetrack/scIDX-Nuc/genetrack_s30e60F6/K562_benzonase-seq_master_MIDPOINT_92-127_NUC_s30e60F6.gff
NUCLEOSOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230719_BI_Genetrack/scIDX-Nuc/genetrack_s40e80F5/K562_benzonase-seq_master_MIDPOINT_128-164_NUC_s40e80F5.gff
SUPRAOCTASOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230719_BI_Genetrack/scIDX-Nuc/genetrack_s50e100F3/K562_benzonase-seq_master_MIDPOINT_165-Inf_NUC_s50e100F3.gff
#output directory
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES
SHIFT_CHECK=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES/shift_check

#set scriptmanager
SCRIPTMANAGER=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar

#set bam library file to BI_rep1
BAM1=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/MEP_Project/01_BAM/Chromatin/Input_Benzonase/25769-25976_Input_-_K562_-_-_50Unuclease10min-0cycSonic_BI_hg19.bam

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
TETRASOMES_BEDFILE=$(echo $TETRASOMES | rev | cut -d"/" -f1 | rev | awk '{gsub(/.gff/,"_tetrasomes.bed"); print}')
HEXASOMES_BEDFILE=$(echo $HEXASOMES | rev | cut -d"/" -f1 | rev | awk '{gsub(/.gff/,"_hexasomes.bed"); print}')
NUCLEOSOMES_BEDFILE=$(echo $NUCLEOSOMES | rev | cut -d"/" -f1 | rev | awk '{gsub(/.gff/,"_nucleosomes.bed"); print}')
SUPRAOCTASOMES_BEDFILE=$(echo $SUPRAOCTASOMES| rev | cut -d"/" -f1 | rev | awk '{gsub(/.gff/,"_supraoctasomes.bed"); print}')
SUBTETRASOMES_BEDFILE_10K=$(echo $SUBTETRASOMES_BEDFILE | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_10K.bed"); print}')
TETRASOMES_BEDFILE_10K=$(echo $TETRASOMES_BEDFILE | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_10K.bed"); print}')
HEXASOMES_BEDFILE_10K=$(echo $HEXASOMES_BEDFILE | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_10K.bed"); print}')
NUCLEOSOMES_BEDFILE_10K=$(echo $NUCLEOSOMES_BEDFILE | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_10K.bed"); print}')
SUPRAOCTASOMES_BEDFILE_10K=$(echo $SUPRAOCTASOMES_BEDFILE | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_10K.bed"); print}')
SUBTETRASOMES_BEDFILE_10K_1000bp=$(echo $SUBTETRASOMES_BEDFILE_10K | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_1000bp.bed"); print}')
TETRASOMES_BEDFILE_10K_1000bp=$(echo $TETRASOMES_BEDFILE_10K | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_1000bp.bed"); print}')
HEXASOMES_BEDFILE_10K_1000bp=$(echo $HEXASOMES_BEDFILE_10K | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_1000bp.bed"); print}')
NUCLEOSOMES_BEDFILE_10K_1000bp=$(echo $NUCLEOSOMES_BEDFILE_10K | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_1000bp.bed"); print}')
SUPRAOCTASOMES_BEDFILE_10K_1000bp=$(echo $SUPRAOCTASOMES_BEDFILE_10K | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_1000bp.bed"); print}')
SUBTETRASOMES_BEDFILE_1bp=$(echo $SUBTETRASOMES_BEDFILE | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_1bp.bed"); print}')
TETRASOMES_BEDFILE_1bp=$(echo $TETRASOMES_BEDFILE | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_1bp.bed"); print}')
HEXASOMES_BEDFILE_1bp=$(echo $HEXASOMES_BEDFILE | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_1bp.bed"); print}')
NUCLEOSOMES_BEDFILE_1bp=$(echo $NUCLEOSOMES_BEDFILE | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_1bp.bed"); print}')
SUPRAOCTASOMES_BEDFILE_1bp=$(echo $SUPRAOCTASOMES_BEDFILE | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_1bp.bed"); print}')
SUBTETRASOMES_BEDFILE_MIDPOINT=$(echo $SUBTETRASOMES_BEDFILE_1bp | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_midpoint.bed"); print}')
TETRASOMES_BEDFILE_MIDPOINT=$(echo $TETRASOMES_BEDFILE_1bp | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_midpoint.bed"); print}')
HEXASOMES_BEDFILE_MIDPOINT=$(echo $HEXASOMES_BEDFILE_1bp | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_midpoint.bed"); print}')
NUCLEOSOMES_BEDFILE_MIDPOINT=$(echo $NUCLEOSOMES_BEDFILE_1bp | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_midpoint.bed"); print}')
SUPRAOCTASOMES_BEDFILE_MIDPOINT=$(echo $SUPRAOCTASOMES_BEDFILE_1bp | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_midpoint.bed"); print}')
SUBTETRASOMES_BEDFILE_EXPANDED=$(echo $SUBTETRASOMES_BEDFILE_MIDPOINT | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_54bp.bed"); print}')
TETRASOMES_BEDFILE_EXPANDED=$(echo $TETRASOMES_BEDFILE_MIDPOINT | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_91bp.bed"); print}')
HEXASOMES_BEDFILE_EXPANDED=$(echo $HEXASOMES_BEDFILE_MIDPOINT | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_127bp.bed"); print}')
NUCLEOSOMES_BEDFILE_EXPANDED=$(echo $NUCLEOSOMES_BEDFILE_MIDPOINT | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_164bp.bed"); print}')
SUPRAOCTASOMES_BEDFILE_EXPANDED=$(echo $SUPRAOCTASOMES_BEDFILE_MIDPOINT | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_200bp.bed"); print}')
SUBTETRASOMES_BEDFILE_EXPANDED_FIX=$(echo $SUBTETRASOMES_BEDFILE_EXPANDED | rev | cut -d"/" -f1 | rev | awk '{gsub(/_1bp_midpoint_54bp.bed/,"_midpoint_54bp_final.bed"); print}')
TETRASOMES_BEDFILE_EXPANDED_FIX=$(echo $TETRASOMES_BEDFILE_EXPANDED | rev | cut -d"/" -f1 | rev | awk '{gsub(/_1bp_midpoint_91bp.bed/,"_midpoint_91bp_final.bed"); print}')
HEXASOMES_BEDFILE_EXPANDED_FIX=$(echo $HEXASOMES_BEDFILE_EXPANDED | rev | cut -d"/" -f1 | rev | awk '{gsub(/_1bp_midpoint_127bp.bed/,"_midpoint_127bp_final.bed"); print}')
NUCLEOSOMES_BEDFILE_EXPANDED_FIX=$(echo $NUCLEOSOMES_BEDFILE_EXPANDED | rev | cut -d"/" -f1 | rev | awk '{gsub(/_1bp_midpoint_164bp.bed/,"_midpoint_164bp_final.bed"); print}')
SUPRAOCTASOMES_BEDFILE_EXPANDED_FIX=$(echo $SUPRAOCTASOMES_BEDFILE_EXPANDED | rev | cut -d"/" -f1 | rev | awk '{gsub(/_1bp_midpoint_200bp.bed/,"_midpoint_200bp_final.bed"); print}')
SUBTETRASOMES_BEDFILE_EXPANDED_FIX_LABEL=$(echo $SUBTETRASOMES_BEDFILE_EXPANDED_FIX | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_label.bed"); print}')
TETRASOMES_BEDFILE_EXPANDED_FIX_LABEL=$(echo $TETRASOMES_BEDFILE_EXPANDED_FIX | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_label.bed"); print}')
HEXASOMES_BEDFILE_EXPANDED_FIX_LABEL=$(echo $HEXASOMES_BEDFILE_EXPANDED_FIX | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_label.bed"); print}')
NUCLEOSOMES_BEDFILE_EXPANDED_FIX_LABEL=$(echo $NUCLEOSOMES_BEDFILE_EXPANDED_FIX | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_label.bed"); print}')
SUPRAOCTASOMES_BEDFILE_EXPANDED_FIX_LABEL=$(echo $SUPRAOCTASOMES_BEDFILE_EXPANDED_FIX | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_label.bed"); print}')
BEDFILE1=$(echo $SUBTETRASOMES_BEDFILE_10K | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BEDFILE2=$(echo $TETRASOMES_BEDFILE_10K | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BEDFILE3=$(echo $HEXASOMES_BEDFILE_10K | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BEDFILE4=$(echo $NUCLEOSOMES_BEDFILE_10K | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BEDFILE5=$(echo $SUPRAOCTASOMES_BEDFILE_10K | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$BEDFILE1" | awk -F. '{print $1"_midpoint.out"}')
CDT1=$(echo "$BAM1a""_""$BEDFILE1" | awk -F. '{print $1}')
OUT2=$(echo "$BAM1a""_""$BEDFILE2" | awk -F. '{print $1"_midpoint.out"}')
CDT2=$(echo "$BAM1a""_""$BEDFILE2" | awk -F. '{print $1}')
OUT3=$(echo "$BAM1a""_""$BEDFILE3" | awk -F. '{print $1"_midpoint.out"}')
CDT3=$(echo "$BAM1a""_""$BEDFILE3" | awk -F. '{print $1}')
OUT4=$(echo "$BAM1a""_""$BEDFILE4" | awk -F. '{print $1"_midpoint.out"}')
CDT4=$(echo "$BAM1a""_""$BEDFILE4" | awk -F. '{print $1}')
OUT5=$(echo "$BAM1a""_""$BEDFILE5" | awk -F. '{print $1"_midpoint.out"}')
CDT5=$(echo "$BAM1a""_""$BEDFILE5" | awk -F. '{print $1}')

#set output directory
cd $OUTPUT

#convert gff to bedfiles
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed -s $SUBTETRASOMES > $SUBTETRASOMES_BEDFILE
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed -s $TETRASOMES > $TETRASOMES_BEDFILE
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed -s $HEXASOMES > $HEXASOMES_BEDFILE
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed -s $NUCLEOSOMES > $NUCLEOSOMES_BEDFILE
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed -s $SUPRAOCTASOMES > $SUPRAOCTASOMES_BEDFILE

#confirm files are shifted correctly
#set to directory for checking shift of bedfiles
cd $SHIFT_CHECK

#shuffle and take 10K sites from each bedfile
cat $OUTPUT/$SUBTETRASOMES_BEDFILE | shuf | head -10000 > $SUBTETRASOMES_BEDFILE_10K
cat $OUTPUT/$TETRASOMES_BEDFILE | shuf | head -10000 > $TETRASOMES_BEDFILE_10K
cat $OUTPUT/$HEXASOMES_BEDFILE | shuf | head -10000 > $HEXASOMES_BEDFILE_10K
cat $OUTPUT/$NUCLEOSOMES_BEDFILE | shuf | head -10000 > $NUCLEOSOMES_BEDFILE_10K
cat $OUTPUT/$SUPRAOCTASOMES_BEDFILE | shuf | head -10000 > $SUPRAOCTASOMES_BEDFILE_10K

#expand bedfiles by 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $SUBTETRASOMES_BEDFILE_10K
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $TETRASOMES_BEDFILE_10K
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $HEXASOMES_BEDFILE_10K
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $NUCLEOSOMES_BEDFILE_10K
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $SUPRAOCTASOMES_BEDFILE_10K

##do tag pile-up and confirm center of peak is at 0 bp. Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT1 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 $SHIFT_CHECK/$SUBTETRASOMES_BEDFILE_10K_1000bp $BAM1
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT2 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 $SHIFT_CHECK/$TETRASOMES_BEDFILE_10K_1000bp $BAM1
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3 $SHIFT_CHECK/$HEXASOMES_BEDFILE_10K_1000bp $BAM1
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT4 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4 $SHIFT_CHECK/$NUCLEOSOMES_BEDFILE_10K_1000bp $BAM1
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT5 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5 $SHIFT_CHECK/$SUPRAOCTASOMES_BEDFILE_10K_1000bp $BAM1

#make figures for composite plots
java -jar $SCRIPTMANAGER figure-generation composite-plot --title="center of subtetrasomes" $OUT1
java -jar $SCRIPTMANAGER figure-generation composite-plot --title="center of tetrasomes" $OUT2
java -jar $SCRIPTMANAGER figure-generation composite-plot --title="center of hexasomes" $OUT3
java -jar $SCRIPTMANAGER figure-generation composite-plot --title="center of nucleosomes" $OUT4
java -jar $SCRIPTMANAGER figure-generation composite-plot --title="center of supraoctasomes" $OUT5

#set output directory
cd $OUTPUT

#expand all bedfile to 1 bp so that midpoint of each fragment is designated
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1 $SUBTETRASOMES_BEDFILE
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1 $TETRASOMES_BEDFILE
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1 $HEXASOMES_BEDFILE
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1 $NUCLEOSOMES_BEDFILE
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1 $SUPRAOCTASOMES_BEDFILE

#put midpoint of fragment in column 4 with st. dev.
cat $SUBTETRASOMES_BEDFILE_1bp | awk '{if ($6=="+") print $1"\t"$2"\t"$3"\t"$4",midpoint="$2"\t"$5"\t"$6; else print ""}' > $SUBTETRASOMES_BEDFILE_MIDPOINT
cat $TETRASOMES_BEDFILE_1bp | awk '{if ($6=="+") print $1"\t"$2"\t"$3"\t"$4",midpoint="$2"\t"$5"\t"$6; else print ""}' > $TETRASOMES_BEDFILE_MIDPOINT
cat $HEXASOMES_BEDFILE_1bp | awk '{if ($6=="+") print $1"\t"$2"\t"$3"\t"$4",midpoint="$2"\t"$5"\t"$6; else print ""}' > $HEXASOMES_BEDFILE_MIDPOINT
cat $NUCLEOSOMES_BEDFILE_1bp | awk '{if ($6=="+") print $1"\t"$2"\t"$3"\t"$4",midpoint="$2"\t"$5"\t"$6; else print ""}' > $NUCLEOSOMES_BEDFILE_MIDPOINT
cat $SUPRAOCTASOMES_BEDFILE_1bp | awk '{if ($6=="+") print $1"\t"$2"\t"$3"\t"$4",midpoint="$2"\t"$5"\t"$6; else print ""}' > $SUPRAOCTASOMES_BEDFILE_MIDPOINT

#ensure that final start / end match previous bedfiles

#expand bedfile by maximum DNA length for each indicated nucleosomal peak; supraoctasomes here will be expanded to 200 bp, there is no maximum for them. 
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=54 $SUBTETRASOMES_BEDFILE_MIDPOINT
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=91 $TETRASOMES_BEDFILE_MIDPOINT
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=127 $HEXASOMES_BEDFILE_MIDPOINT
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=164 $NUCLEOSOMES_BEDFILE_MIDPOINT
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=200 $SUPRAOCTASOMES_BEDFILE_MIDPOINT

#Then remove any lines in expanded bedfile with a region expanded past a chromosome end (as seen by a "-" in column 2 or 3.
cat $SUBTETRASOMES_BEDFILE_EXPANDED | awk 'BEGIN {OFS=FS="\t"} $2 !~ /^(-)/' | awk 'BEGIN {OFS=FS="\t"} $3 !~ /^(-)/' > $SUBTETRASOMES_BEDFILE_EXPANDED_FIX
cat $TETRASOMES_BEDFILE_EXPANDED | awk 'BEGIN {OFS=FS="\t"} $2 !~ /^(-)/' | awk 'BEGIN {OFS=FS="\t"} $3 !~ /^(-)/' > $TETRASOMES_BEDFILE_EXPANDED_FIX
cat $HEXASOMES_BEDFILE_EXPANDED | awk 'BEGIN {OFS=FS="\t"} $2 !~ /^(-)/' | awk 'BEGIN {OFS=FS="\t"} $3 !~ /^(-)/' > $HEXASOMES_BEDFILE_EXPANDED_FIX
cat $NUCLEOSOMES_BEDFILE_EXPANDED | awk 'BEGIN {OFS=FS="\t"} $2 !~ /^(-)/' | awk 'BEGIN {OFS=FS="\t"} $3 !~ /^(-)/' > $NUCLEOSOMES_BEDFILE_EXPANDED_FIX
cat $SUPRAOCTASOMES_BEDFILE_EXPANDED | awk 'BEGIN {OFS=FS="\t"} $2 !~ /^(-)/' | awk 'BEGIN {OFS=FS="\t"} $3 !~ /^(-)/' > $SUPRAOCTASOMES_BEDFILE_EXPANDED_FIX

#add 7th column with particle type indicated; add numbr to category for later sorting
cat $SUBTETRASOMES_BEDFILE_EXPANDED_FIX | awk '{if ($6=="+") $7="<:"; else $7=""; print}' > $SUBTETRASOMES_BEDFILE_EXPANDED_FIX_LABEL
cat $TETRASOMES_BEDFILE_EXPANDED_FIX | awk '{if ($6=="+") $7="T:"; else $7=""; print}' > $TETRASOMES_BEDFILE_EXPANDED_FIX_LABEL
cat $HEXASOMES_BEDFILE_EXPANDED_FIX | awk '{if ($6=="+") $7="H:"; else $7=""; print}' > $HEXASOMES_BEDFILE_EXPANDED_FIX_LABEL
cat $NUCLEOSOMES_BEDFILE_EXPANDED_FIX | awk '{if ($6=="+") $7="N:"; else $7=""; print}' > $NUCLEOSOMES_BEDFILE_EXPANDED_FIX_LABEL
cat $SUPRAOCTASOMES_BEDFILE_EXPANDED_FIX | awk '{if ($6=="+") $7="<:"; else $7=""; print}' > $SUPRAOCTASOMES_BEDFILE_EXPANDED_FIX_LABEL

#concatenate above files to master bedfile
cat $SUBTETRASOMES_BEDFILE_EXPANDED_FIX_LABEL $TETRASOMES_BEDFILE_EXPANDED_FIX_LABEL $HEXASOMES_BEDFILE_EXPANDED_FIX_LABEL $NUCLEOSOMES_BEDFILE_EXPANDED_FIX_LABEL $SUPRAOCTASOMES_BEDFILE_EXPANDED_FIX_LABEL > K562_redundant_master_nucleosomal_particles.bed
cat K562_redundant_master_nucleosomal_particles.bed | wc -l

#remove intermediate files
rm $SHIFT_CHECK/$SUBTETRASOMES_BEDFILE_10K
rm $SHIFT_CHECK/$TETRASOMES_BEDFILE_10K
rm $SHIFT_CHECK/$HEXASOMES_BEDFILE_10K
rm $SHIFT_CHECK/$NUCLEOSOMES_BEDFILE_10K
rm $SHIFT_CHECK/$SUPRAOCTASOMES_BEDFILE_10K
rm $SHIFT_CHECK/$SUBTETRASOMES_BEDFILE_10K_1000bp
rm $SHIFT_CHECK/$TETRASOMES_BEDFILE_10K_1000bp
rm $SHIFT_CHECK/$HEXASOMES_BEDFILE_10K_1000bp
rm $SHIFT_CHECK/$NUCLEOSOMES_BEDFILE_10K_1000bp
rm $SHIFT_CHECK/$SUPRAOCTASOMES_BEDFILE_10K_1000bp
rm $SUBTETRASOMES_BEDFILE
rm $TETRASOMES_BEDFILE
rm $HEXASOMES_BEDFILE
rm $NUCLEOSOMES_BEDFILE
rm $SUPRAOCTASOMES_BEDFILE
rm $SUBTETRASOMES_BEDFILE_1bp
rm $TETRASOMES_BEDFILE_1bp
rm $HEXASOMES_BEDFILE_1bp
rm $NUCLEOSOMES_BEDFILE_1bp
rm $SUPRAOCTASOMES_BEDFILE_1bp
rm $SUBTETRASOMES_BEDFILE_MIDPOINT
rm $TETRASOMES_BEDFILE_MIDPOINT
rm $HEXASOMES_BEDFILE_MIDPOINT
rm $NUCLEOSOMES_BEDFILE_MIDPOINT
rm $SUPRAOCTASOMES_BEDFILE_MIDPOINT
rm $SUBTETRASOMES_BEDFILE_EXPANDED
rm $TETRASOMES_BEDFILE_EXPANDED
rm $HEXASOMES_BEDFILE_EXPANDED
rm $NUCLEOSOMES_BEDFILE_EXPANDED
rm $SUPRAOCTASOMES_BEDFILE_EXPANDED
#rm $SUBTETRASOMES_BEDFILE_EXPANDED_LABEL
#rm $TETRASOMES_BEDFILE_EXPANDED_LABEL
#rm $HEXASOMES_BEDFILE_EXPANDED_LABEL
#rm $NUCLEOSOMES_BEDFILE_EXPANDED_LABEL
#rm $SUPRAOCTASOMES_BEDFILE_EXPANDED_LABEL

# finish script
echo "DONE"
