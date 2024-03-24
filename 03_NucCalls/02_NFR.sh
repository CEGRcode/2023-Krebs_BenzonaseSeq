# purpose - determine NFRs based on +1 and -1 Nuc calls in input folder.

# usage
# qq
#
# example
#
# 'qq'

#output directory
INPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_plus1_minus1/output_v2_NonRed_Oct_Hex_Tet_230825
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_plus1_minus1/02_NFR_output_231128

#set scriptmanager
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

#set output directory
cd $OUTPUT

#same sort by TSS distance but output has more columns: column 7 = sortOrder_TSSid; column 8 = distance from TSS to Minus1
cat $INPUT/K562_CoPRO_TSS_upstreamOctamers_Intersect.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); else $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k14,14n -r | awk '{print $4"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$6"\t"$14}' > K562_Minus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures.tsv

#***should we shuffle before doing SORT so that chr. regions are not nearby
#same sort by TSS distance but output has more columns: column 7 = sortOrder_TSSid; column 8 = distance from TSS to Minus1
cat $INPUT/K562_CoPRO_TSS_downstreamOctamers_Intersect.bed | awk '{if ($6=="+") $14=sqrt(((($8+$9)/2)-$2)*((($8+$9)/2)-$2)); else $14=sqrt(((($8+$9)/2)-$3)*((($8+$9)/2)-$3)); print}' | sort -k4,4d -k14,14n -r | sort -k4,4d -u | sort -k14,14n -r | awk '{print $4"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$6"\t"$14}' > K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures.tsv

#pull TSS ID(s) per each file
cat K562_Minus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures.tsv | sort -k1,1n | awk '{print $1}' > K562_Minus1_TSSid_sortOrder.tsv
cat K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures.tsv | sort -k1,1n | awk '{print $1}' > K562_Plus1_TSSid_sortOrder.tsv

#make file with all shared TSS(s) between the two files.
awk 'FNR==NR{a[$1];next} $1 in a' K562_Minus1_TSSid_sortOrder.tsv K562_Plus1_TSSid_sortOrder.tsv > TSS_IDs_shared_bothFiles.tsv

#reSORT files with shared bedfile 
##if number of inferred -1 and +1 positions do not match, remove those that are not shared. A NFR can not be calculated for a TSS without both inferred positions anyway.
awk -F"\t" 'FNR==NR{a[$1]=$0; next} $1 in a {print a[$1]}' K562_Minus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures.tsv TSS_IDs_shared_bothFiles.tsv > K562_Minus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv
#11,693 sites
awk -F"\t" 'FNR==NR{a[$1]=$0; next} $1 in a {print a[$1]}' K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures.tsv TSS_IDs_shared_bothFiles.tsv > K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv
#11,693 sites

#calculate NFRs. NOTE: this code restricts NFR calculate where both +1 and -1 were inferred per a single TSS. **May need to remove *same sites.
#double check that TSS IDs are the same for each row as seen in column 17; if inferred -1 and +1 nucleosomes are the same indicate that in column 18; final file are the following columns: 1 (TSSid), 2 (RNA exp. sort order), 3 (-1 to TSS), 4 (TSS to +1), 5 (NFR), 6 (match if TSS IDs match for both files), 7 (same if inferred -1 and +1 match).
paste K562_Minus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv |  awk '{if ($1==$9) $17="match"; else $17="-"; print}' | awk '{if ($5==$13) $18="same"; else $18="-"; print}' | awk '{print $1"\t"$1"\t"$8"\t"$16"\t"$8+$16"\t"$17"\t"$18}' | awk -F_ '{print $2"\t"$0}' | awk '{print $1"\t"$2"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' >  TSSid_RNAsortOrder_DistanceFeatures.tsv

#make bedfile of NFR. If TSS IDs match, then calculate NFR in 17th column. Next sort by increasing NFR length. Final befile has 11,693 sites.
paste K562_Minus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv |  awk '{if ($1==$9) $17=$8+$16; else $17="-"; print}' | sort -k17,17n | awk '{if ($5==$13) $18="same"; else $18="-"; print}' | awk '{if ($18="same") print $2"\t"$3"\t"$4"\t"$1"\t"$6"\t"$7; else {if ($18="-" || $7="+") print $2"\t"$4+1"\t"$11-1"\t"$1"\t"$6"\t"$7; else {if ($18="-" || $7="-") print $2"\t"$11-1"\t"$4+1"\t"$1"\t"$6"\t"$7; else ""}}}' > K562_NFR_bedfile_increasingLength.bed

#make bedfile of NFR when -1 and +1 nucleosomes are NOT the same AND that NFR has at least 1 bp. Output is bedfile for true NFRs. Final output is: 1 (chr #), 2 (NFR start), 3 (NFR end), 4( St dev, midpoint, TSS ID), 5 (number; ? number of reads per this nucleosomal peak), 6 (+ or - strand)
paste K562_Minus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv |  awk '{if ($1==$9) $17=$8+$16; else $17="-"; print}' | sort -k17,17n | awk '{if ($5==$13) $18="same"; else $18="-"; print}' | awk '$18 == "-"' | awk '{if ($18=="-") $19=$11-$4; else $19="-"; print}' | awk '($19 != 0)' | awk '{if ($7 == "+") print $2"\t"($4+1)"\t"($11-1)"\t"$5","$1"\t"$6"\t"$7; else print $2"\t"($11-1)"\t"($4+1)"\t"$5","$1"\t"$6"\t"$7}' > K562_trueNFR.bed

#**the following are unused and have not been checked since Sept.**
#make bedfile of NFR when -1 and +1 nucleosomes are NOT the same. Output is bedfile -1 flanking true NFR.
paste K562_Minus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv |  awk '{if ($1==$9) $17=$8+$16; else $17="-"; print}' | sort -k17,17n | awk '{if ($5==$13) $18="same"; else $18="-"; print}' | awk '$18 == "-"' | awk '{if ($18=="-") $19=$11-$4; else $19="-"; print}' | awk '($19 != 0)' | awk '{print $2"\t"$3"\t"$4"\t"$1","$5"\t"$6"\t"$7}' > K562_trueNFR_Minus1region.bed

#make bedfile of NFR when -1 and +1 nucleosomes are NOT the same. Output is bedfile +1 flanking true NFR.
paste K562_Minus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv |  awk '{if ($1==$9) $17=$8+$16; else $17="-"; print}' | sort -k17,17n | awk '{if ($5==$13) $18="same"; else $18="-"; print}' | awk '$18 == "-"' | awk '{if ($18=="-") $19=$11-$4; else $19="-"; print}' | awk '($19 != 0)' | awk '{print $10"\t"$11"\t"$12"\t"$9","$13"\t"$14"\t"$15}' > K562_trueNFR_Plus1region.bed

#make bedfile with -1 and +1 dyad with NFR sort **not sure if these are correct
#reference point as -1 dyad with same sort as above (increasing length of NFR)
paste K562_Minus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv |  awk '{if ($1==$9) $17=$8+$16; else $17="-"; print}' | sort -k17,17n | awk '{if ($5==$13) $18="same"; else $18="-"; print}' | awk '{if ($3>0) $19=sqrt(((($3+$4)/2)-$2)*((($3+$4)/2)-$2)); else $19=""; print}' | awk '{print $2"\t"$19"\t"$19"\t"$1"\t"$6"\t"$7}' > K562_Minus1_IncreasingNFR.bed

#reference point as -1 dyad with same sort as above (increasing length of NFR)
paste K562_Minus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet_distanceFeatures-reSORTed.tsv |  awk '{if ($1==$9) $17=$8+$16; else $17="-"; print}' | sort -k17,17n | awk '{if ($5==$13) $18="same"; else $18="-"; print}' | awk '{if ($3>0) $19=sqrt(((($11+$12)/2)-$2)*((($11+$12)/2)-$2)); else $19=""; print}' | awk '{print $2"\t"$19"\t"$19"\t"$1"\t"$6"\t"$7}' > K562_Minus1_IncreasingNFR.bed

# finish script
echo "DONE"
