# purpose - determine closest inferrred nucleosome to CoPRO-determined TSS.**this code now keeps the RNA value in the bedfile for later correaltion. 240708 v1 - combines MNase and BNase peaks -> takes closest peak to TSS agonist to method

# usage
# qq
#
# example
#
# 'qq'

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/04_combined_plus1_240730

#set scriptmanager
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar

#set bedfiles and *.genome file
BEDFILE1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/01_C_intermediate_FILES/K562_benzonase-seq_master_MIDPOINT_128-164_NUC_s40e80F5_nucleosomes_midpoint_164bp_final.bed
BEDFILE2=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/02_intermediate_FILES/K562_uHex.bed
BEDFILE3=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/02_intermediate_FILES/K562_uTetra.bed
BEDFILE4=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/03_Genetrack_peaks_v2/final_shifted_scIDX_MNase_164bp_final.bed
TSS=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_plus1_minus1/K562_CoPRO-expressed_Gene-refSeqTSS_2000bp.bed
TSS_RNAsortOrder=$OUTPUT/K562_CoPRO-expressed_Gene-refSeqTSS_2000bp_RNAsortOrder.bed
HG19_GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/files/human.hg19.genome

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

#set output directory
cd $OUTPUT

#assume above TSS bedfile was sorted by RNA expression. Update column 4 to include RNA sort. This can ne used to revert back to this order later.
cat $TSS | awk '{print $0"\t"NR}' | awk '{print $1"\t"$2"\t"$3"\t"$7"_"$4"\t"$5"\t"$6}' > $TSS_RNAsortOrder

#get -1 nucleosomes in reference to CoPRO TSS
#expand TSS bedfile by 1 bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1 $TSS_RNAsortOrder -o=$OUTPUT/K562_CoPRO-expressed_Gene-refSeqTSS_1bp.bed

#SORT TSS bedfile HERE so that all files are SORTed later
bedtools sort -i K562_CoPRO-expressed_Gene-refSeqTSS_1bp.bed > K562_CoPRO-expressed_Gene-refSeqTSS_1bp_SORTed.bed

#expand TSS upstream in length of 5 kb in stranded fashion
bedtools slop -i K562_CoPRO-expressed_Gene-refSeqTSS_1bp_SORTed.bed -g $HG19_GENOME -s -r 0 -l 5000 > K562_CoPRO-1bpTSS-upstream5kb.bed

#split full-length nucleosome bedfile into up to 0.5 million rows each
cat $BEDFILE1 | head -1000000 > Nuclesomes_fullLength_first1M.bed
cat $BEDFILE1 | head -2000000 | tail -1000000 > Nuclesomes_fullLength_second1M.bed
cat $BEDFILE1 | head -3000000 | tail -1000000 > Nuclesomes_fullLength_third1M.bed
cat $BEDFILE1 | head -4000000 | tail -1000000 > Nuclesomes_fullLength_fourth1M.bed
cat $BEDFILE1 | head -5000000 | tail -1000000 > Nuclesomes_fullLength_fifth1M.bed
cat $BEDFILE1 | tail -76544 > Nuclesomes_fullLength_lastset.bed

#split full-length nucleosome bedfile into up to 0.5 million rows each
cat $BEDFILE4 | head -1000000 > MNase_fullLength_first1M.bed
cat $BEDFILE4 | head -2000000 | tail -1000000 > MNase_fullLength_second1M.bed
cat $BEDFILE4 | head -3000000 | tail -1000000 > MNase_fullLength_third1M.bed
cat $BEDFILE4 | head -4000000 | tail -1000000 > MNase_fullLength_fourth1M.bed
cat $BEDFILE4 | head -5000000 | tail -1000000 > MNase_fullLength_fifth1M.bed
cat $BEDFILE4 | head -6000000 | tail -1000000 > MNase_fullLength_sixth1M.bed
cat $BEDFILE4 | head -7000000 | tail -1000000 > MNase_fullLength_seventh1M.bed
cat $BEDFILE4 | head -8000000 | tail -1000000 > MNase_fullLength_eighth1M.bed
cat $BEDFILE4 | head -9000000 | tail -1000000 > MNase_fullLength_ninth1M.bed
cat $BEDFILE4 | head -10000000 | tail -1000000 > MNase_fullLength_tenth1M.bed
cat $BEDFILE4 | tail -808810 > MNase_fullLength_lastset.bed

#intersect 5kb region upstream of TSS with all inferred nucleosomes. This code was changed to intersect each category separately (vs all 3 together).
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b Nuclesomes_fullLength_first1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_first1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b Nuclesomes_fullLength_second1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_second1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b Nuclesomes_fullLength_third1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_third1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b Nuclesomes_fullLength_fourth1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_fourth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b Nuclesomes_fullLength_fifth1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_fifth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b Nuclesomes_fullLength_lastset.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_lastset.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b MNase_fullLength_first1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_first1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b MNase_fullLength_second1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_second1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b MNase_fullLength_third1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_third1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b MNase_fullLength_fourth1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_fourth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b MNase_fullLength_fifth1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_fifth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b MNase_fullLength_sixth1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_sixth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b MNase_fullLength_seventh1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_seventh1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b MNase_fullLength_eighth1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_eighth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b MNase_fullLength_ninth1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_ninth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b MNase_fullLength_tenth1M.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_tenth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-upstream5kb.bed -b MNase_fullLength_lastset.bed -bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_lastset.bed

#cat all intersected regions together
#cat K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_first1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_second1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_third1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_fourth1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_fifth1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_lastset.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_uHex.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_uTetra.bed > K562_CoPRO_TSS_upstreamOctamers_Intersect.bed
cat K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_first1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_second1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_third1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_fourth1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_fifth1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_lastset.bed > K562_CoPRO_TSS_upstreamOctamers_Intersect.bed

cat K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_first1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_second1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_third1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_fourth1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_fifth1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_sixth1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_seventh1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_eighth1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_ninth1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_tenth1M.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_lastset.bed > K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase.bed

#combine files (benz without Hex or Tetra AND MNase
cat K562_CoPRO_TSS_upstreamOctamers_Intersect.bed K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase.bed | shuf > K562_CoPRO_TSS_upstreamOctamers_Intersect_BNase_MNase.bed

#create 14th column that determines distance from TSS to nucleosome THEN sort gene-specific ID (column 4) by increasing bp distance; sort bedfile by increasing distance from TSS to +1 nucleosomal particle; make bedfile for +1 nucleosomal particles [can't use bedtools closest here; features are too close] ****check Math here
#fixed so that midpoint of dyad is taken or (($8+$9)/2) instead of $8
##240708- added TSS ID to output (column4)
cat K562_CoPRO_TSS_upstreamOctamers_Intersect_BNase_MNase.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); else $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k14,14n -r | awk '{print $7"\t"$8"\t"$9"\t"$4","$10"\t"$11"\t"$6}' > K562_Minus1_SORTdistanceToTSS_BNase_MNase.bed

#same as above but with TSS ID in 7th column
cat K562_CoPRO_TSS_upstreamOctamers_Intersect_BNase_MNase.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); else $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k14,14n -r | awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$6"\t"$4}' > K562_Minus1_SORTdistanceToTSS_BNase_MNase.tsv

#now sort output nuclesomes by #_TSS to revert back to decreasing RNA_expression
cat K562_CoPRO_TSS_upstreamOctamers_Intersect_BNase_MNase.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); else $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k4,4n -r | awk '{print $7"\t"$8"\t"$9"\t"$4","$10","$5"\t"$11"\t"$6}' > K562_Minus1_SORTbyRNAexp_BNase_MNase.bed

#***should we shuffle before doing SORT so that chr. regions are not nearby
#keep the same sort as in upstream file but get +1 nucleosome

#not get closest  +1 nucleosome to CoPRO TSS via downstream region
#expand TSS downstream in length of 5 kb in stranded fashion ********add plus 50 bp flanking side????
bedtools slop -i K562_CoPRO-expressed_Gene-refSeqTSS_1bp_SORTed.bed -g $HG19_GENOME -s -r 5000 -l 0 > K562_CoPRO-1bpTSS-downstream5kb.bed

#intersect 5kb region downstream of TSS with all inferred nucleosomes
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b Nuclesomes_fullLength_first1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_first1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b Nuclesomes_fullLength_second1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_second1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b Nuclesomes_fullLength_third1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_third1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b Nuclesomes_fullLength_fourth1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_fourth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b Nuclesomes_fullLength_fifth1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_fifth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b Nuclesomes_fullLength_lastset.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_lastset.bed
#bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b $BEDFILE2 -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_uHex.bed
#bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b $BEDFILE3 -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_uTetra.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b MNase_fullLength_first1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_first1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b MNase_fullLength_second1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_second1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b MNase_fullLength_third1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_third1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b MNase_fullLength_fourth1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_fourth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b MNase_fullLength_fifth1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_fifth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b MNase_fullLength_sixth1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_sixth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b MNase_fullLength_seventh1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_seventh1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b MNase_fullLength_eighth1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_eighth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b MNase_fullLength_ninth1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_ninth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b MNase_fullLength_tenth1M.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_tenth1M.bed
bedtools intersect -wo -a K562_CoPRO-1bpTSS-downstream5kb.bed -b MNase_fullLength_lastset.bed -bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_lastset.bed

#cat all intersected regions together
#cat K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_first1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_second1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_third1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_fourth1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_fifth1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_lastset.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_uHex.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_uTetra.bed > K562_CoPRO_TSS_downstreamOctamers_Intersect.bed

cat K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_first1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_second1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_third1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_fourth1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_fifth1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_lastset.bed > K562_CoPRO_TSS_downstreamOctamers_Intersect.bed

cat K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_first1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_second1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_third1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_fourth1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_fifth1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_sixth1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_seventh1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_eighth1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_ninth1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_tenth1M.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_lastset.bed > K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase.bed

#combine files (benz without Hex or Tetra AND MNase
cat K562_CoPRO_TSS_downstreamOctamers_Intersect.bed K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase.bed | shuf > K562_CoPRO_TSS_downstreamOctamers_Intersect_BNase_MNase.bed

#create 14th column that determines distance from TSS to nucleosome THEN sort gene-specific ID (column 4) by increasing bp distance; sort bedfile by increasing distance from TSS to +1 nucleosomal particle; make bedfile for +1 nucleosomal particles [can't use bedtools closest here; features are too close]
#fixed so that midpoint of dyad is taken or (($8+$9)/2) instead of $8
cat K562_CoPRO_TSS_downstreamOctamers_Intersect_BNase_MNase.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); else $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k14,14n -r | awk '{print $7"\t"$8"\t"$9"\t"$4","$10"\t"$11"\t"$6}' > K562_Plus1_SORTdistanceToTSS_BNase_MNase.bed

#print lines for original code (singletons)
cat K562_CoPRO_TSS_downstreamOctamers_Intersect_BNase_MNase.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); else $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k4,4n | awk -F"\t" 'BEGIN{OFS=FS}{print $0}' > K562_Plus1_singletons.tab

#code for pairs
cat K562_CoPRO_TSS_downstreamOctamers_Intersect_BNase_MNase.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); else $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); print}' | sort -k4,4d -k14,14n -r | awk '!seen[$0]++' | awk 'NR == 1 {prev_line = $0; prev_col4 = $4; prev_col14 = $14; next} $4 == prev_col4 && ($14 - prev_col14 <= 50) {print prev_line; print $0} {prev_line = $0; prev_col4 = $4; prev_col14 = $14}' | awk '$10 ~ /MNase/ {print $0}' | sort -k4,4d -u  | sort -k4,4n  | awk -F"\t" 'BEGIN{OFS=FS}{print $0}' > K562_Plus1_pairs.tab

#combined and sort > take pair if within 50 bp of singleton; 
cat K562_Plus1_singletons.tab K562_Plus1_pairs.tab | sort -k4,4n -k14,14n | awk '!seen[$0]++' | awk 'NR == 1 {prev_line = $0; prev_col4 = $4; prev_col14 = $14; next} $4 == prev_col4 && ($14 - prev_col14 <= 50) {print $0; next} {print prev_line; prev_line = $0; prev_col4 = $4; prev_col14 = $14}' | sort -k4,4d -u  |  sort -k14,14n | awk '{print $7"\t"$8"\t"$9"\t"$4","$10"\t"$11"\t"$6}' > K562_Plus1_SORTdistanceToTSS_final_singletons_pairs.bed

##initial code
#do math, then sort by TSS ID and increased TSS distance; delete repeated lines; check that column 4 matches previous line AND than column 14 is no greater than 50 than previous column 14 (to see if they are in a pair) -> pairs of benzonase and MNase peaks now -> chose MNase peak within pair

#same as above but with TSS ID in 7th column
#cat K562_CoPRO_TSS_downstreamOctamers_Intersect_BNase_MNase.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); else $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k14,14n -r | awk '{print $7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$6"\t"$4}' > K562_Plus1_SORTdistanceToTSS_BNase_MNase.tsv

#now sort output nuclesomes by #_TSS to revert back to decreasing RNA_expression
cat K562_CoPRO_TSS_downstreamOctamers_Intersect_BNase_MNase.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); else $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k4,4n -r | awk '{print $7"\t"$8"\t"$9"\t"$4","$10","$5"\t"$11"\t"$6}' > K562_Plus1_SORTbyRNAexp_BNase_MNase.bed

#remove intermediate files
rm Nuclesomes_fullLength_first1M.bed
rm Nuclesomes_fullLength_second1M.bed
rm Nuclesomes_fullLength_third1M.bed
rm Nuclesomes_fullLength_fourth1M.bed
rm Nuclesomes_fullLength_fifth1M.bed
rm Nuclesomes_fullLength_lastset.bed
rm MNase_fullLength_first1M.bed
rm MNase_fullLength_second1M.bed
rm MNase_fullLength_third1M.bed
rm MNase_fullLength_fourth1M.bed
rm MNase_fullLength_fifth1M.bed
rm MNase_fullLength_sixth1M.bed
rm MNase_fullLength_seventh1M.bed
rm MNase_fullLength_eighth1M.bed
rm MNase_fullLength_ninth1M.bed
rm MNase_fullLength_tenth1M.bed
rm MNase_fullLength_lastset.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_first1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_second1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_third1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_fourth1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_fifth1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_Nuc_lastset.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_first1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_second1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_third1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_fourth1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_fifth1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_sixth1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_seventh1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_eighth1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_ninth1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_tenth1M.bed
rm K562_CoPRO_TSS_upstreamOctamers_Intersect_MNase_lastset.bed
#rm K562_CoPRO_TSS_upstreamOctamers_Intersect_uHex.bed
#rm K562_CoPRO_TSS_upstreamOctamers_Intersect_uTetra.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_first1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_second1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_third1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_fourth1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_fifth1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_Nuc_lastset.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_first1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_second1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_third1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_fourth1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_fifth1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_sixth1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_seventh1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_eighth1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_ninth1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_tenth1M.bed
rm K562_CoPRO_TSS_downstreamOctamers_Intersect_MNase_lastset.bed
#rm K562_CoPRO_TSS_downstreamOctamers_Intersect_uHex.bed
#rm K562_CoPRO_TSS_downstreamOctamers_Intersect_uTetra.bed

# finish script
echo "DONE"
