#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=24GB
#SBATCH --time=0:30:00
#SBATCH --partition=open

# (Not tested, awaiting Jordan's response to questions)

# Determines the total number of full-length nucleosomes at the +1 position, 
# makes a bedfile and then uses this bedfile to calculate the number of tags 
# (representing SNs) in the proximal and distal SNs regions (based on 
# midpoints in Fig. 4). The final files gives the number of rows with 
# proximal SNs, distal SNs, or both at +1 positions that are based on full-
# length nucleosomes.

set -eux
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

PLUSONE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_plus1_minus1/output_v2_NonRed_Oct_Hex_Tet_230825/K562_Plus1_SORTbyRNAexp_nonRedOct_Hex_Tet.bed
OTABLE=STable2-2.tsv
GENOME=../data/hg19_files/hg19.fa
BLACKLIST=../data/hg19_files/hg19_blacklist.bed
BAMFILE=Merged_BX_H3K4me3.bam
BAM=`basename $BAMFILE ".bam"`

SUMROW=../bin/sum_Row_CDT.pl
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar

TEMP=temp-STable2-FullLength
[ -d $TEMP ] || mkdir $TEMP

# Create ID & Particle size text file
awk 'BEGIN{OFS="\t";FS="\t"}{if ($2>0) print $4,$3-$2; else print $4,""}' $PLUSONE > $TEMP/Plus1_Sizes.txt

# Filter +1 Nuc RefPT to 164bp particle sized annotations
perl $FILTERL $PLUSONE <(awk '{if ($2==164) print $1}' $TEMP/Plus1_Sizes.txt) 3 keep $TEMP/Plus1_FILTER-164bp.bed

##determine the number of full-length nucleosome that have tags for SNs
# Expand 2bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2 $TEMP/Plus1_FILTER-164bp.bed -o $TEMP/Plus1_FILTER-164bp_2bp.bed

# Build RefPT: proximal (upstream) or distal (downstream) half of nucleosome relative to +1 midpoint (dyad)
# math is based on 80 bp SNs for shifts: proximal SN is ~2 bp upstream of dyad so 80 + 2 = 82 for start, -2 + -2 = -4 for end; distal SN is ~4 bp downstream of dyad so -2 + -4 = -6 for start, 80 + 4 = 84 for end
bedtools slop -i $TEMP/Plus1_FILTER-164bp_2bp.bed -g $GENOME -l 82 -r -4 -s > $TEMP/Plus1_FILTER-164bp_proximal.bed
bedtools slop -i $TEMP/Plus1_FILTER-164bp_2bp.bed -g $GENOME -l -6 -r 84 -s > $TEMP/Plus1_FILTER-164bp_distal.bed
# TagPileup
# Settings: midpoint(m), No smoothing (N), max insert size 80bp, load blacklist
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m --combined -x 80 -N -p \
    --cpu 4 -f $BLACKLIST \
    -M $TEMP/Plus1_FILTER-164bp_proximal_$BAM.cdt -o $TEMP/Plus1_FILTER-164bp_proximal_$BAM.out \
    $TEMP/Plus1_FILTER-164bp_proximal.bed $BAMFILE
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m --combined -x 80 -N -p \
    --cpu 4 -f $BLACKLIST \
    -M $TEMP/Plus1_FILTER-164bp_distal_$BAM.cdt -o $TEMP/Plus1_FILTER-164bp_distal_$BAM.out \
    $TEMP/Plus1_FILTER-164bp_distal.bed $BAMFILE

# Row-wise sum tags in CDT
perl $SUMROW $TEMP/Plus1_FILTER-164bp_proximal_$BAM.cdt $TEMP/Plus1_FILTER-164bp_proximal_$BAM.tab
perl $SUMROW $TEMP/Plus1_FILTER-164bp_distal_$BAM.cdt $TEMP/Plus1_FILTER-164bp_distal_$BAM.tab


#make a file showing how many of above full-length nucleosomes at +1 position have at least 1 tag in either proximal SN region. First only rows with matching IDs are kept.
paste $TEMP/Plus1_FILTER-164bp.bed $TEMP/Plus1_FILTER-164bp_proximal_$BAM.cdt $TEMP/Plus1_FILTER-164bp_distal_$BAM.cdt > $TEMP/Plus1_BED-CDT-proximal-distal.tab

# Quality-check that ids match (should output empty file)
awk '{if($4!=$7 || $4!=$9) print}' $TEMP/Plus1_BED-CDT-proximal-distal.tab > $TEMP/QC_proximal_IDmatch.tab

# ======Build STable======

# Count +1 Nuc rows/BED records
wc -l $PLUSONE | awk 'BEGIN{OFS="\t"}{print "total",$1}' > $OTABLE
# Count annotations for each particle size
cut -f 2 $TEMP/Plus1_Sizes.txt | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2"bp",$1}' >> $OTABLE

awk 'BEGIN{FS="\t"}{if($8!=0) print}' $TEMP/Plus1_BED-CDT-proximal-distal.tab | wc -l | awk 'BEGIN{OFS="\t"}{print "proximal_SN_region_withSNs",$1}' > $OTABLE
awk 'BEGIN{FS="\t"}{if($8==0) print}' $TEMP/Plus1_BED-CDT-proximal-distal.tab | wc -l | awk 'BEGIN{OFS="\t"}{print "proximal_SN_region_withoutSNs",$1}' >> $OTABLE
awk 'BEGIN{FS="\t"}{if($10!=0) print}' $TEMP/Plus1_BED-CDT-proximal-distal.tab | wc -l | awk 'BEGIN{OFS="\t"}{print "distal_SN_region_withSNs",$1}' >> $OTABLE
awk 'BEGIN{FS="\t"}{if($10==0) print}' $TEMP/Plus1_BED-CDT-proximal-distal.tab | wc -l | awk 'BEGIN{OFS="\t"}{print "distal_SN_region_withoutSNs",$1}' >> $OTABLE
awk 'BEGIN{FS="\t"}{if($8!=0 && $10!=0) print}' $TEMP/Plus1_BED-CDT-proximal-distal.tab | wc -l | awk 'BEGIN{OFS="\t"}{print "both_SNs_present",$1}' >> $OTABLE
wc -l $TEMP/Plus1_BED-CDT-proximal-distal.tab | awk 'BEGIN{OFS="\t"}{print "total_possible_for_combined",$1}' >> $OTABLE

# clean-up
# rm -r $TEMP