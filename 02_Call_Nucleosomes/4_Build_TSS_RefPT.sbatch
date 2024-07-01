#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 2:00:00
#SBATCH -A open
#SBATCH -o logs/4_Build_TSS_RefPT.log.out
#SBATCH -e logs/4_Build_TSS_RefPT.log.err

# Call TSS using UCSC knownGene download shifted to mode CoPRO signal and
# renamed to a human-readable gene name.

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/storage/home/owl5022/scratch/2023-Krebs_Benzonase-seq/02_Call_Nucleosomes
###

# Dependencies
# - bedtools
# - java

set -exo
module load anaconda3
module load bedtools
source activate bx

# Inputs and outputs
UCSCDB=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/
COPRO=../data/BAM/CoPRO_-_merge_hg19.bam
CAPPED=../data/BAM/CoPRO_Capped_merge_hg19.bam
UNCAPPED=../data/BAM/CoPRO_Uncapped_merge_hg19.bam
CPG=../data/RefPT-Other/CpG_Islands.bed
KREBS=../data/RefPT-Krebs

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar
GENETRACK=../bin/genetrack_v2.py
UNION=../bin/union_peak_by_ID.py
FILTERL=../bin/filter_BED_by_list_ColumnSelect.pl
UPDATES=../bin/update_BED_score_with_TAB_score.pl
MODE=../bin/call_TSS-localMax_CDT.py

TEMP=MakeTSS
[ -d $TEMP ] || mkdir $TEMP

## =====Format RefSeq annotations=====

# Download annotations
wget -c $UCSCDB/bigZips/genes/hg19.knownGene.gtf.gz
wget -c $UCSCDB/database/knownCanonical.txt.gz
wget -c $UCSCDB/database/knownGenePep.txt.gz
wget -c $UCSCDB/database/knownToLynx.txt.gz
wget -c $UCSCDB/database/knownToMalacards.txt.gz
# wget -c $UCSCDB/gc5Base/hg19.gc5Base.txt.gz

# Build CanonicalPep ID set
gzip -dc knownCanonical.txt.gz | cut -f5 > $TEMP/knownCanonical.ids
gzip -dc knownGenePep.txt.gz | cut -f1 > $TEMP/knownGenePep.ids
grep -f $TEMP/knownGenePep.ids $TEMP/knownCanonical.ids > $TEMP/knownCanonicalPep.ids

# Build master ID to name set (use Malacard name map, then Lynx, and add self-map for remaining Canonical Pep)
gzip -dc knownToLynx.txt.gz > $TEMP/knownToLynx.txt
gzip -dc knownToMalacards.txt.gz | cut -f1 > $TEMP/knownToMalacards.ids
perl $FILTERL $TEMP/knownToLynx.txt $TEMP/knownToMalacards.ids 0 remove $TEMP/knownToLynx_FILTER-RemoveMalacards.txt
gzip -dc knownToMalacards.txt.gz | cat - $TEMP/knownToLynx_FILTER-RemoveMalacards.txt | sort | uniq > $TEMP/knownToMalacards-wLynx.txt
perl $FILTERL $TEMP/knownCanonicalPep.ids <(cut -f1 $TEMP/knownToMalacards-wLynx.txt) 0 remove $TEMP/knownCanonicalPep-NoNames.txt
awk '{OFS="\t"}{print $1,$1}' $TEMP/knownCanonicalPep-NoNames.txt | cat $TEMP/knownToMalacards-wLynx.txt - > $TEMP/knownTo_NameMap.txt

# Filter annotations to keep "transcript" types w/ canonical chr
gzip -dc hg19.knownGene.gtf.gz \
	| awk -F $'\t' '{OFS="\t"}{if($3=="transcript") print;}' \
	| grep $'^chr[0-9XYM]\+\t' \
	> $TEMP/hg19.knownGene.transcripts.gtf
# rm hg19.knownGene.gtf.gz

# Parse out gene ID info (knownGene)
cut -f9 $TEMP/hg19.knownGene.transcripts.gtf \
	| awk -F ";" '{
		gsub("transcript_id","",$2);
		gsub("[ \"]+","",$2);
		print $2
	}' > $TEMP/hg19.knownGene.transcripts.ids

# Append annotation info with new gene ID info
paste <(cut -f1-8 $TEMP/hg19.knownGene.transcripts.gtf) $TEMP/hg19.knownGene.transcripts.ids > $TEMP/hg19.knownGene.id-transcripts.gtf

# Subset all knownGene IDs to knownCanonical IDs
perl $FILTERL $TEMP/hg19.knownGene.id-transcripts.gtf $TEMP/knownCanonicalPep.ids 8 keep $TEMP/hg19.knownCanonicalPep.id-transcripts.gtf

# Convert ORF --> TSS and GTF --> BED (dedup for isoforms w same start site and ID)
awk '{FS="\t"}{OFS="\t"}{
		if($7=="-") {
			$4=$5-1;
		} else {
			$5=$4+1;
		}
		print $1,$4,$5,$9,".",$7
	}' $TEMP/hg19.knownCanonicalPep.id-transcripts.gtf \
	> $TEMP/hg19_knownCanonicalPep-TSS.bed

# Expand 200bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 200 $TEMP/hg19_knownCanonicalPep-TSS.bed -o $TEMP/hg19_knownCanonicalPep-TSS_200bp.bed

## =====Call CoPRO TSS (mode)=====

# Pileup CoPRO data (Capped)
java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -2 $TEMP/hg19_knownCanonicalPep-TSS_200bp.bed $CAPPED -M $TEMP/Capped_READ2
# java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -2 $TEMP/hg19_knownCanonicalPep-TSS_200bp.bed $COPRO -M $TEMP/CoPRO_READ2
# java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -2 $TEMP/hg19_knownCanonicalPep-TSS_200bp.bed $UNCAPPED -M $TEMP/Uncapped_READ2

# Shift TSS to position of local max (mode)
python $MODE -i $TEMP/Capped_READ2_sense.cdt -r $TEMP/hg19_knownCanonicalPep-TSS_200bp.bed -o $TEMP/maxPeak.bed

# Map ID to gene name
perl $UPDATES $TEMP/maxPeak.bed $TEMP/knownTo_NameMap.txt $TEMP/TSS.bed

## =====Add TSS expression scores=====
[ -d $KREBS/2000bp ] || mkdir -p $KREBS/2000bp

# Re-expand 200bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 200 $TEMP/TSS.bed -o $TEMP/TSS_200bp.bed

# Re-pileup CoPRO data
java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -2 $TEMP/TSS_200bp.bed ../data/BAM/CoPRO_Capped_merge_hg19.bam -M $TEMP/Capped_READ2_TSS_200bp

# Row-wise sum tags to get expression score
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum --merge -o $TEMP/CappedExpression.out $TEMP/Capped_READ2_TSS_200bp_sense.cdt

# Update TSS score col-5 with expression score
perl $UPDATES $TEMP/TSS.bed $TEMP/CappedExpression.out $TEMP/TSS_SCORE-CappedExpression.bed

# Sort by expression value
sort -rnk5,5 $TEMP/TSS_SCORE-CappedExpression.bed > $KREBS/TSS_GROUP-All_SORT-CappedExpression.bed

# Group by Expressed/Unexpressed
awk '{FS="\t"}{if($5>0) print}' $KREBS/TSS_GROUP-All_SORT-CappedExpression.bed > $KREBS/TSS_GROUP-Expressed_SORT-Expression.bed
awk '{FS="\t"}{if($5==0) print}' $KREBS/TSS_GROUP-All_SORT-CappedExpression.bed > $KREBS/TSS_GROUP-Unexpressed.bed

# Expand 2000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/TSS_GROUP-All_SORT-CappedExpression.bed   -o $KREBS/2000bp/TSS_GROUP-All_SORT-CappedExpression_2000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/TSS_GROUP-Expressed_SORT-Expression.bed   -o $KREBS/2000bp/TSS_GROUP-Expressed_SORT-Expression_2000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/TSS_GROUP-Unexpressed.bed -o $KREBS/2000bp/TSS_GROUP-Unexpressed_2000bp.bed

## =====Sort TSS by CpG island length=====

# Peak-align CpG islands on TSS
java -jar $SCRIPTMANAGER peak-analysis peak-align-ref $CPG $KREBS/2000bp/TSS_GROUP-Expressed_SORT-Expression_2000bp.bed -o $TEMP/TSS_CpG.cdt

# Sort BED
java -jar $SCRIPTMANAGER coordinate-manipulation sort-bed -x 900 2000 $KREBS/TSS_GROUP-Expressed_SORT-Expression.bed $TEMP/TSS_CpG.cdt -o $TEMP/TSS_CpG_SORTED

# Move to global RefPT directory and format name
mv $TEMP/TSS_CpG_SORTED.bed $KREBS/TSS_GROUP-Expressed_SORT-CpG.bed

# Expand 2000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/TSS_GROUP-Expressed_SORT-CpG.bed -o $KREBS/2000bp/TSS_GROUP-Expressed_SORT-CpG_2000bp.bed
