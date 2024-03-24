# purpose - this version script intersects all redundant (smaller) peaks with reundant peaks. When multiple smaller peaks intersect with a single larger particle, all are kept and tagged. 
# usage
# qq
#
# example
#
# 'qq'

#input nucleosomal peak bedfiles. All "NUCLEOSOMES" are unqiue.
NUCLEOSOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES/K562_benzonase-seq_master_MIDPOINT_128-164_NUC_s40e80F5_nucleosomes_midpoint_164bp_final.bed
HEXASOMES_REDUNDANT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES/K562_benzonase-seq_master_MIDPOINT_92-127_NUC_s30e60F6_hexasomes_midpoint_127bp_final.bed
TETRASOMES_REDUNDANT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES/K562_benzonase-seq_master_MIDPOINT_55-91_NUC_s20e40F5_tetrasomes_midpoint_91bp_final.bed
SUBTETRASOMES_REDUNDANT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES/K562_benzonase-seq_master_MIDPOINT_0-54_NUC_s10e20_subtetrasomes_midpoint_54bp_final.bed
SUPRAOCTASOMES_REDUNDANT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/01_C_intermediate_FILES/K562_benzonase-seq_master_MIDPOINT_165-Inf_NUC_s50e100F3_supraoctasomes_midpoint_200bp_final.bed
UNIQUE_HEXASOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/02_intermediate_FILES/K562_uHex.bed
UNIQUE_TETRASOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/02_intermediate_FILES/K562_uTetra.bed
UNIQUE_SUBTETRASOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/02_intermediate_FILES/K562_uSubtetra.bed
UNIQUE_SUPRAOCTASOMES=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/02_intermediate_FILES/K562_uSupraoct.bed

#output directory
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/03_C_intermediate_FILES

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
NUC_UniqueIDs=$(echo $NUCLEOSOMES | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_uniqueIDs.bed"); print}')
RedHEX_UniqueIDs=$(echo $HEXASOMES_REDUNDANT | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_uniqueIDs.bed"); print}')
RedTETRA_UniqueIDs=$(echo $TETRASOMES_REDUNDANT | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_uniqueIDs.bed"); print}')
RedSUBTETRA_UniqueIDs=$(echo $SUBTETRASOMES_REDUNDANT | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_uniqueIDs.bed"); print}')
RedSUPRAOCT_UniqueIDs=$(echo $SUPRAOCTASOMES_REDUNDANT | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_uniqueIDs.bed"); print}')
uHEX_UniqueIDs=$(echo $UNIQUE_HEXASOMES | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_uniqueIDs.bed"); print}')
uTetra_UniqueIDs=$(echo $UNIQUE_TETRASOMES | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_uniqueIDs.bed"); print}')
uSubtetra_UniqueIDs=$(echo $UNIQUE_SUBTETRASOMES | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_uniqueIDs.bed"); print}')
uSupraoct_UniqueIDs=$(echo $UNIQUE_SUPRAOCTASOMES | rev | cut -d"/" -f1 | rev | awk '{gsub(/.bed/,"_uniqueIDs.bed"); print}')

#set output directory
cd $OUTPUT

#add unique IDs to every bedfile (nucleosomes and redundant hexasomes for now). add unique ID by ("Nuc"), update column 4 with chr#,start,end,stdev,Nuc#; Nuc# = row# -> make final file
cat $NUCLEOSOMES | awk '{print $0"\t"NR}' | awk '{print $0"\t"$1","$2","$3","$4",""Nuc"$7}' | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$5"\t"$6}' > $NUC_UniqueIDs
cat $HEXASOMES_REDUNDANT | awk '{print $0"\t"NR}' | awk '{print $0"\t"$1","$2","$3","$4",""Hex"$7}' | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$5"\t"$6}' > $RedHEX_UniqueIDs
cat $TETRASOMES_REDUNDANT | awk '{print $0"\t"NR}' | awk '{print $0"\t"$1","$2","$3","$4",""Tetra"$7}' | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$5"\t"$6}' > $RedTETRA_UniqueIDs
cat $SUBTETRASOMES_REDUNDANT | awk '{print $0"\t"NR}' | awk '{print $0"\t"$1","$2","$3","$4",""Subtetra"$7}' | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$5"\t"$6}' > $RedSUBTETRA_UniqueIDs
cat $SUPRAOCTASOMES_REDUNDANT | awk '{print $0"\t"NR}' | awk '{print $0"\t"$1","$2","$3","$4",""Supraoct"$7}' | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$5"\t"$6}' > $RedSUPRAOCT_UniqueIDs
cat $UNIQUE_HEXASOMES | awk '{print $0"\t"NR}' | awk '{print $0"\t"$1","$2","$3","$4",""uHex"$7}' | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$5"\t"$6}' > $uHEX_UniqueIDs
cat $UNIQUE_TETRASOMES | awk '{print $0"\t"NR}' | awk '{print $0"\t"$1","$2","$3","$4",""uTetra"$7}' | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$5"\t"$6}' > $uTetra_UniqueIDs
cat $UNIQUE_SUBTETRASOMES | awk '{print $0"\t"NR}' | awk '{print $0"\t"$1","$2","$3","$4",""uSubtetra"$7}' | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$5"\t"$6}' > $uSubtetra_UniqueIDs
cat $UNIQUE_SUPRAOCTASOMES | awk '{print $0"\t"NR}' | awk '{print $0"\t"$1","$2","$3","$4",""uSupraoct"$7}' | awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$5"\t"$6}' > $uSupraoct_UniqueIDs

#intersect all particles with full-length nucleosomes
#intersect all redundant hexasomes with nucleosomes. Use min (-F in reference to fraction of B) of 1.0 to ensure that an any intersected hexamer fully overlaps with a nucleosome. All others are not considered. Maximum overlap is 126 bp as seen in data. THEN substitute "-1" with 0 in columns 8, 9, and 11. Make final bedfile with chr# / start / end /HEX label;intersecting Nuc label; bp overlap/tag count / strand of hexasome  ****if this file does not work for tag pile-up take chr# and strand from nucleosome section -> unique lines without "."
bedtools intersect -wao -F 1.0 -a $NUC_UniqueIDs -b $RedHEX_UniqueIDs | awk '{if ($8=="-1") $8="0"; print}' | awk '{if ($9=="-1") $9="0"; print}' | awk '{if ($11=="-1") $11="0"; print}' | awk '{print $7"\t"$8"\t"$9"\t"$10";"$4";"$13"\t"$11"\t"$12}' > Nucleosomes_intersect_redundantHex.bed

#intersect all redundant tetrasomes with nucleosomes. Use min (-F in reference to fraction of B) of 1.0 to ensure that an any intersected tetramer fully overlaps with a nucleosome. All others are not considered. Maximum overlap is 90 bp as seen in data. THEN substitute "-1" with 0 in columns 8, 9, and 11. Make final bedfile with chr# / start / end /TETRA label;intersecting Nuc label; bp overlap/tag count / strand of tetrasome  ****if this file does not work for tag pile-up take chr# and strand from nucleosome section -> unique lines without "."
bedtools intersect -wao -F 1.0 -a $NUC_UniqueIDs -b $RedTETRA_UniqueIDs | awk '{if ($8=="-1") $8="0"; print}' | awk '{if ($9=="-1") $9="0"; print}' | awk '{if ($11=="-1") $11="0"; print}' | awk '{print $7"\t"$8"\t"$9"\t"$10";"$4";"$13"\t"$11"\t"$12}' > Nucleosomes_intersect_redundantTetra.bed

#intersect all redundant subtetrasomes with nucleosomes. Use min (-F in reference to fraction of B) of 1.0 to ensure that an any intersected subtetramer fully overlaps with a nucleosome. All others are not considered. THEN substitute "-1" with 0 in columns 8, 9, and 11. Make final bedfile with chr# / start / end /SUBTETRA label;intersecting Nuc label; bp overlap/tag count / strand of subtetrasome      ****if this file does not work for tag pile-up take chr# and strand from nucleosome section -> unique lines without "."
bedtools intersect -wao -F 1.0 -a $NUC_UniqueIDs -b $RedSUBTETRA_UniqueIDs | awk '{if ($8=="-1") $8="0"; print}' | awk '{if ($9=="-1") $9="0"; print}' | awk '{if ($11=="-1") $11="0"; print}' | awk '{print $7"\t"$8"\t"$9"\t"$10";"$4";"$13"\t"$11"\t"$12}' > Nucleosomes_intersect_redundantSubtetra.bed

#intersect all redundant supraoctasomes with nucleosomes. Minimum (-f) of 1.0 used so that only supraoctasomes that overlap nucleosomes by 164 bp are considered. THEN substitute "-1" with 0 in columns 8, 9, and 11. Make final bedfile with chr# / start / end /SUPRAOCT label;intersecting Nuc label; bp overlap/tag count / strand of supraoctasome      ****if this file does not work for tag pile-up take chr# and strand from nucleosome section -> unique lines without "."
bedtools intersect -wao -f 1.0 -a $NUC_UniqueIDs -b $RedSUPRAOCT_UniqueIDs | awk '{if ($8=="-1") $8="0"; print}' | awk '{if ($9=="-1") $9="0"; print}' | awk '{if ($11=="-1") $11="0"; print}' | awk '{print $7"\t"$8"\t"$9"\t"$10";"$4";"$13"\t"$11"\t"$12}' > Nucleosomes_intersect_redundantSupraoct.bed

#intersect all particles with non-redundant hexasomes
#intersect all redundant tetrasomes with non-redundant hexasomes. Use min (-F in reference to fraction of B) of 1.0 to ensure that an any intersected tetramer fully overlaps with a hexasome. All others are not considered. THEN substitute "-1" with 0 in columns 8, 9, and 11. Make final bedfile with chr# / start / end /TETRA label;intersecting Nuc label; bp overlap/tag count / strand of tetrasome  
bedtools intersect -wao -F 1.0 -a $uHEX_UniqueIDs -b $RedTETRA_UniqueIDs | awk '{if ($8=="-1") $8="0"; print}' | awk '{if ($9=="-1") $9="0"; print}' | awk '{if ($11=="-1") $11="0"; print}' | awk '{print $7"\t"$8"\t"$9"\t"$10";"$4";"$13"\t"$11"\t"$12}' > uHex_intersect_redundantTetra.bed

#intersect all redundant subtetrasomes with non-redundant hexasomes. Use min (-F in reference to fraction of B) of 1.0 to ensure that an any intersected subtetramer fully overlaps with a hexasome. All others are not considered. THEN substitute "-1" with 0 in columns 8, 9, and 11. Make final bedfile with chr# / start / end /SUBTETRA label;intersecting Nuc label; bp overlap/tag count / strand of subtetrasome
bedtools intersect -wao -F 1.0 -a $uHEX_UniqueIDs -b $RedSUBTETRA_UniqueIDs | awk '{if ($8=="-1") $8="0"; print}' | awk '{if ($9=="-1") $9="0"; print}' | awk '{if ($11=="-1") $11="0"; print}' | awk '{print $7"\t"$8"\t"$9"\t"$10";"$4";"$13"\t"$11"\t"$12}' > uHex_intersect_redundantSubtetra.bed

#intersect all redundant supraoctasomes with non-redundant hexasomes. Minimum (-f) of 1.0 used so that only supraoctasomes that overlap nucleosomes by 164 bp are considered. THEN substitute "-1" with 0 in columns 8, 9, and 11. Make final bedfile with chr# / start / end /SUPRAOCT label;intersecting Nuc label; bp overlap/tag count / strand of supraoctasome
bedtools intersect -wao -f 1.0 -a $uHEX_UniqueIDs -b $RedSUPRAOCT_UniqueIDs | awk '{if ($8=="-1") $8="0"; print}' | awk '{if ($9=="-1") $9="0"; print}' | awk '{if ($11=="-1") $11="0"; print}' | awk '{print $7"\t"$8"\t"$9"\t"$10";"$4";"$13"\t"$11"\t"$12}' > uHex_intersect_redundantSupraoct.bed

#intersect all particles with non-redundant tetrasomes
#intersect all redundant subtetrasomes with non-redundant tetrasomes. Use min (-F in reference to fraction of B) of 1.0 to ensure that an any intersected subtetramer fully overlaps with a tetrasome. All others are not considered. THEN substitute "-1" with 0 in columns 8, 9, and 11. Make final bedfile with chr# / start / end /SUBTETRA label;intersecting Nuc label; bp overlap/tag count / strand of subtetrasome
bedtools intersect -wao -F 1.0 -a $uTetra_UniqueIDs -b $RedSUBTETRA_UniqueIDs | awk '{if ($8=="-1") $8="0"; print}' | awk '{if ($9=="-1") $9="0"; print}' | awk '{if ($11=="-1") $11="0"; print}' | awk '{print $7"\t"$8"\t"$9"\t"$10";"$4";"$13"\t"$11"\t"$12}' > uTetra_intersect_redundantSubtetra.bed

#intersect all redundant supraoctasomes with non-redundant tetrasomes. Minimum (-f) of 1.0 used so that only supraoctasomes that overlap nucleosomes by 164 bp are considered. THEN substitute "-1" with 0 in columns 8, 9, and 11. Make final bedfile with chr# / start / end /SUPRAOCT label;intersecting Nuc label; bp overlap/tag count / strand of supraoctasome
bedtools intersect -wao -f 1.0 -a $uTetra_UniqueIDs -b $RedSUPRAOCT_UniqueIDs | awk '{if ($8=="-1") $8="0"; print}' | awk '{if ($9=="-1") $9="0"; print}' | awk '{if ($11=="-1") $11="0"; print}' | awk '{print $7"\t"$8"\t"$9"\t"$10";"$4";"$13"\t"$11"\t"$12}' > uTetra_intersect_redundantSupraoct.bed

#intersect all particles with non-redundant subtetrasomes
#intersect all redundant supraoctasomes with non-redundant subtetrasomes. Minimum (-f) of 1.0 used so that only supraoctasomes that overlap nucleosomes by 164 bp are considered. THEN substitute "-1" with 0 in columns 8, 9, and 11. Make final bedfile with chr# / start / end /SUPRAOCT label;intersecting Nuc label; bp overlap/tag count / strand of supraoctasome
bedtools intersect -wao -f 1.0 -a $uSubtetra_UniqueIDs -b $RedSUPRAOCT_UniqueIDs | awk '{if ($8=="-1") $8="0"; print}' | awk '{if ($9=="-1") $9="0"; print}' | awk '{if ($11=="-1") $11="0"; print}' | awk '{print $7"\t"$8"\t"$9"\t"$10";"$4";"$13"\t"$11"\t"$12}' > uSubtetra_intersect_redundantSupraoct.bed

# finish script
echo "DONE"
