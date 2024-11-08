#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 20:00:00
#SBATCH -A open
#SBATCH -o logs/6_Match_Dyads_for_NFR.log.out
#SBATCH -e logs/6_Match_Dyads_for_NFR.log.err

# Determine NFRs based on +1 and -1 Nuc calls from previous script.
# see 03_NucCalls/02_NFR.sh

# data/RefPT-Krebs
#   |--NFR_SORT-NFRLength.bed

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
###

# Dependencies
# - java
# - perl
# - python

set -exo
module load anaconda3
module load bedtools
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
KREBS=$WRK/../data/RefPT-Krebs/
NUCLEOSOME=$KREBS/Nucleosome_uHex_uTetra.bed
EXPRESSED=$KREBS/TSS_GROUP-Expressed_SORT-Expression.bed
UNEXPRESSED=$KREBS/TSS_GROUP-Unexpressed.bed

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
FILTERL=../bin/filter_BED_by_list_ColumnSelect.pl

TEMP=MakePlusMinus
[ -d $TEMP ] || mkdir $TEMP

## =====Match PlusOne and MinusOne Dyads=====

# Get shared Rank IDs (found in both plus an minus dyad files)
cut -f5 $TEMP/MinusOneDyad_SORT-DistToExpressedTSS.tsv $TEMP/PlusOneDyad_SORT-DistToExpressedTSS.tsv \
	| sort | uniq -c | awk '{if($1==2) print $2}' > $TEMP/Shared_RankIDs.ids

# Filter to keep only shared IDs
perl $FILTERL $TEMP/MinusOneDyad_SORT-DistToExpressedTSS.tsv $TEMP/Shared_RankIDs.ids 4 keep $TEMP/Matched-MinusOneDyad_SORT-DistToExpressedTSS.tsv
perl $FILTERL $TEMP/PlusOneDyad_SORT-DistToExpressedTSS.tsv  $TEMP/Shared_RankIDs.ids 4 keep $TEMP/Matched-PlusOneDyad_SORT-DistToExpressedTSS.tsv

# Quality Check for TSS loss
wc -l $TEMP/MinusOneDyad_SORT-DistToExpressedTSS.tsv
wc -l $TEMP/PlusOneDyad_SORT-DistToExpressedTSS.tsv
wc -l $TEMP/Shared_RankIDs.ids
wc -l $TEMP/Matched-MinusOneDyad_SORT-DistToExpressedTSS.tsv # should match Shared_RankIDs
wc -l $TEMP/Matched-PlusOneDyad_SORT-DistToExpressedTSS.tsv  # should match Shared_RankIDs

# Sort TSV by RankExpression
sort -nk5,5 $TEMP/Matched-MinusOneDyad_SORT-DistToExpressedTSS.tsv > $TEMP/Matched-MinusOneDyad_SORT-RankExpression.tsv
sort -nk5,5 $TEMP/Matched-PlusOneDyad_SORT-DistToExpressedTSS.tsv > $TEMP/Matched-PlusOneDyad_SORT-RankExpression.tsv

# Create Matched Dyad file
paste $TEMP/Matched-MinusOneDyad_SORT-RankExpression.tsv $TEMP/Matched-PlusOneDyad_SORT-RankExpression.tsv > $TEMP/MatchedDyads_SORT-RankExpression.tsv

# Quality check for proper matching (expect 1 line: "<NLINES> True True True" without any "False" tokens)
awk '{FS="\t"}{
		CHRMATCH="False";
		IDMATCH="False";
		SCOREMATCH="False";
		STRANDMATCH="False";
		if ($1==$14 && $1==$7 && $7==$20) CHRMATCH="True";
		if ($4==$17) IDMATCH="True";
		if ($5==$18) SCOREMATCH="True";
		if ($6==$19) STRANDMATCH="True";
		print CHRMATCH,IDMATCH,SCOREMATCH,STRANDMATCH;
	}' $TEMP/MatchedDyads_SORT-RankExpression.tsv | sort | uniq -c

# Add distance metrics and reduce redundant info in TSV and sort by NFR length
# col 1-2 --> Global info with col1-chr and col2-Strand
# col 1-6 --> TSS info with col3-Start, col4-End, col5-GeneName, col6-TSSRankExpressionID,
# col 7-11 --> MinusOneDyad info with col7-Start, col8-End, col9-NucID, col10-Nuc Genetrack score, col11-Dist TSS to MinusOneDyad
# col 12-16 --> PlusOneDyad info with col12-Start, col13-End, col14-NucID, col15-Nuc Genetrack score, col16-Dist TSS to PlusOneDyad
# col 17-20 --> NFR info with col17-Start, col18-End, col19-TSSRankExpressionID, col20-length (outer dyad edge to outer dyad edge)
awk '{OFS="\t"}{FS="\t"}{
		NFRSTART=-1;
		NFREND=-1;
		NFRLENGTH=0;
		if ($6=="-") {
			NFRSTART=$21;
			NFREND=$9;
		} else {
			NFRSTART=$8;
			NFREND=$22;
		}
		NFRLENGTH=NFREND-NFRSTART;
		print $1,$6, $2,$3,$4,$5, $8,$9,$10,$11,$13, $21,$22,$23,$24,$26, NFRSTART,NFREND,$5,NFRLENGTH;
	}' $TEMP/MatchedDyads_SORT-RankExpression.tsv \
	> $TEMP/MatchedDyads_SORT-RankExpression_WithNFRInfo.tsv

# Quality check for same nucleosome assigned to +1 and -1 (expect all to be False if cases handled correctly)
awk '{FS="\t"}{
		IDMATCH="False";
		if ($9==$14) IDMATCH="True";
		print IDMATCH;
	}' $TEMP/MatchedDyads_SORT-RankExpression_WithNFRInfo.tsv \
	| sort | uniq -c


## =====Slice Master Info TSV into RefPT=====

# NFR sorted by NFR length: Format +1 and -1 coords for proper BED style NFR
sort -nk20,20 $TEMP/MatchedDyads_SORT-RankExpression_WithNFRInfo.tsv \
	| awk '{OFS="\t"}{FS="\t"}{if ($9!=$14) print $1,$17,$18,$19,$20,$2;}' \
	> $KREBS/NFR_SORT-NFRLength.bed
# Filter by +1 and -1 distance?

# Quality check for how many "NFR" have same nucleosome assigned to +1 and -1 (expect mostly 1s, and low count of 2s)
cut -f1-3,6 $KREBS/NFR_SORT-NFRLength.bed \
	| sort | uniq -c | awk '{print $1}' | sort | uniq -c
# If too many duplicates, we should attempt to reduce the TSS count

# Expand 2000bp
# java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/NFR_SORT-NFRLength.bed -o $KREBS/2000bp/NFR_SORT-NFRLength_2000bp.bed
