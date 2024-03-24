# purpose - determine cloest inferrred nucleosome to CoPRO-determined TSS

# usage
# qq
#
# example
#
# 'qq'

#output directory
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_plus1_minus1/output_v2_NonRed_OCT_Hex_Tet_Unexpressed_230905

#set scriptmanager
SCRIPTMANAGER=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar

#set bedfiles and *.genome file
BEDFILE=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/02_intermediate_FILES/K562_nuc_uHex_uTetra.bed
TSS=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/figures/fig1_atTSS_CpGsort/bedfiles/K562_CoPRO-unexpressed_Gene-TSS_2000bp.bed
HG19_GENOME=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_master_bedfile/files/human.hg19.genome

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

#set output directory
cd $OUTPUT

#get -1 nucleosomes in reference to CoPRO TSS
#expand TSS bedfile by 1 bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1 $TSS -o=$OUTPUT/K562_CoPRO-unexpressed_Gene-TSS_1bp.bed

#SORT TSS bedfile HERE so that all files are SORTed later
bedtools sort -i K562_CoPRO-unexpressed_Gene-TSS_1bp.bed > K562_CoPRO-unexpressed_Gene-TSS_1bp_SORTed.bed

#expand TSS upstream in length of 5 kb in stranded fashion
bedtools slop -i K562_CoPRO-unexpressed_Gene-TSS_1bp_SORTed.bed -g $HG19_GENOME -s -r 0 -l 5000 > K562_CoPRO-1bpTSS-upstream5kb.bed

#intersect 5kb region upstream of TSS with all inferred nucleosomes
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b $BEDFILE -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect.bed

#create 14th column that determines distance from TSS to nucleosome THEN sort gene-specific ID (column 4) by increasing bp distance; sort bedfile by increasing distance from TSS to +1 nucleosomal particle; make bedfile for +1 nucleosomal particles [can't use bedtools closest here; features are too close] ****check Math here
#fixed so that midpoint of dyad is taken or (($8+$9)/2) instead of $8
cat K562_CoPRO_TSS_upstreamOctamers_Intersect.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); else $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k14,14n -r | awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$6}' > K562_Minus1_SORTdistanceToTSS_Unexpressed_nonRedOct_Hex_Tet.bed

#same as above but with TSS ID in 7th column
cat K562_CoPRO_TSS_upstreamOctamers_Intersect.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); else $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k14,14n -r | awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$6"\t"$4}' > K562_Minus1_SORTdistanceToTSS_Unexpressed_nonRedOct_Hex_Tet.tsv

#now sort output nuclesomes by #_TSS to revert back to decreasing RNA_expression
cat K562_CoPRO_TSS_upstreamOctamers_Intersect.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); else $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k4,4n -r | awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$6}' > K562_Minus1_NOsort_Unexpressed_nonRedOct_Hex_Tet.bed

#***should we shuffle before doing SORT so that chr. regions are not nearby
#keep the same sort as in upstream file but get +1 nucleosome

#not get closest  +1 nucleosome to CoPRO TSS via downstream region
#expand TSS downstream in length of 5 kb in stranded fashion ********add plus 50 bp flanking side????
bedtools slop -i K562_CoPRO-unexpressed_Gene-TSS_1bp_SORTed.bed -g $HG19_GENOME -s -r 5000 -l 0 > K562_CoPRO-1bpTSS-downstream5kb.bed

#intersect 5kb region downstream of TSS with all inferred nucleosomes
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b $BEDFILE -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect.bed

#create 14th column that determines distance from TSS to nucleosome THEN sort gene-specific ID (column 4) by increasing bp distance; sort bedfile by increasing distance from TSS to +1 nucleosomal particle; make bedfile for +1 nucleosomal particles [can't use bedtools closest here; features are too close]
#fixed so that midpoint of dyad is taken or (($8+$9)/2) instead of $8
cat K562_CoPRO_TSS_downstreamOctamers_Intersect.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); else $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k14,14n -r | awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$6}' > K562_Plus1_SORTdistanceToTSS_Unexpressed_nonRedOct_Hex_Tet.bed

#same as above but with TSS ID in 7th column
cat K562_CoPRO_TSS_downstreamOctamers_Intersect.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); else $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k14,14n -r | awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$6"\t"$4}' > K562_Plus1_SORTdistanceToTSS_Unexpressed_nonRedOct_Hex_Tet.tsv

#now sort output nuclesomes by #_TSS to revert back to decreasing RNA_expression
cat K562_CoPRO_TSS_downstreamOctamers_Intersect.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); else $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k4,4n -r | awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$6}' > K562_Plus1_NOsort_Unexpressed_nonRedOct_Hex_Tet.bed

# finish script
echo "DONE"
