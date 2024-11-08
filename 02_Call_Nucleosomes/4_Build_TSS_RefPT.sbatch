#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 2:00:00
#SBATCH -A open
#SBATCH -o logs/4_Build_TSS_RefPT.log.out
#SBATCH -e logs/4_Build_TSS_RefPT.log.err

# Call TSS using UCSC knownGene download shifted to mode CoPRO signal and
# renamed to a human-readable gene name.

# data/RefPT-Krebs
#   |--TSS_GROUP-All_SORT-CappedExpression.bed
#   |--TSS_GROUP-Expressed_SORT-Expression.bed
#   |--TSS_GROUP-Unexpressed.bed
#   |--TSS_GROUP-Expressed_SORT-CpG.bed
#   |--2000bp
#     |--TSS_GROUP-All_SORT-CappedExpression_2000bp.bed
#     |--TSS_GROUP-Expressed_SORT-Expression_2000bp.bed
#     |--TSS_GROUP-Unexpressed_2000bp.bed
#     |--TSS_GROUP-Expressed_SORT-CpG_2000bp.bed

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
###

# Dependencies
# - bedtools
# - java
# - perl
# - wget

set -exo
module load anaconda3
module load bedtools
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
UCSCDB=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/
COPRO=../data/BAM/CoPRO_-_merge_hg38.bam
CAPPED=../data/BAM/CoPRO_Capped_merge_hg38.bam
UNCAPPED=../data/BAM/CoPRO_Uncapped_merge_hg38.bam
CPG=../data/RefPT-Other/CpGIslands.bed
KREBS=../data/RefPT-Krebs

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
GENETRACK=../bin/genetrack_v2.py
UNION=../bin/union_peak_by_ID.py
FILTERL=../bin/filter_BED_by_list_ColumnSelect.pl
UPDATES=../bin/update_BED_score_with_TAB_score.pl
MODE=../bin/call_TSS-localMax_CDT.py
GTFCOL=../bin/parse_gtf_column9.py

TEMP=temp-4_Build_TSS_RefPT
[ -d $TEMP ] || mkdir $TEMP

## =====Call de novo CoPRO TSS=====

# BAM to ScIdx
java -jar $SCRIPTMANAGER bam-format-converter bam-to-scidx -2 $CAPPED -o $TEMP/CoPRO-capped.tab
# ~8 min

# Genetrack call peaks
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/genetrack
python $GENETRACK -s 20 -e 40 -F 1 $TEMP/CoPRO-capped.tab
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx
# ~3 min

# Convert to BED
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed $TEMP/genetrack_s20e40F1/CoPRO-capped_s20e40F1.gff -o $TEMP/CoPRO-capped.bed

## =====Format GENCODE annotations (~11 sec)=====

# Download ENCODE GENCODE annotations
wget -c https://www.encodeproject.org/files/gencode.v29.primary_assembly.annotation_UCSC_names/@@download/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz

# Pull out protein coding full transcript annotations
gzip -dc gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz \
	| grep 'gene_type \"protein_coding\";' \
	| awk '{if ($3=="transcript")  print}' \
	> $TEMP/gencode.pcoding_transcripts.gff

# Reduce coordinate windows to unique transcript starts (TSS)
#   - update score col with transcript len
#   - reduce coordinate to strand-aware start end
#   - keep longest transcript for each unique TSS
awk 'BEGIN{OFS="\t";FS="\t"}{
		$6=$5-$4;
		if ($7=="-") {
			$4=$5-1;
			$5=$5-1;
		} else {
			$5=$4-1;
			$4=$4-1;
		}
		print $0;
	}' $TEMP/gencode.pcoding_transcripts.gff \
	| sort -nk6,6 \
	| sort -u -k1,1 -k4,4 \
	> $TEMP/gencode.pcoding_tss.gff

# Update TSS 9th column with identifier and bump 9th column gene info to 10th
awk 'BEGIN{OFS="\t";FS="\t"}{
		$10=$9;
		$9="TSS-"NR"-"$1"-"$4;
		print $0;
	}' $TEMP/gencode.pcoding_tss.gff \
	> $TEMP/gencode.pcoding_tss_uid.tsv

# Parse out gene ID info (set aside this reference for later use...)
# 	TSS-NID-chrZ-coord	gene_id "ENSG00000187634.11"; transcript_id "ENST00000618779.4"; gene_type "protein_coding"; gene_name "SAMD11"; transcript_type "protein_coding"; transcript_name "SAMD11-215"; level 3; protein_id "ENSP00000484256.1"; transcript_support_level "5"; tag "basic"; tag "appris_alternative_2"; havana_gene "OTTHUMG00000040719.10";
# 	gene_id	transcript_id	gene_type	gene_name	transcript_type	transcript_name	level	protein_id	transcript_support_level	tag	tag	havana_gene
python $GTFCOL -i <(cut -f9,10 $TEMP/gencode.pcoding_tss_uid.tsv) -o $TEMP/gencode.pcoding_tss_uid.geneinfo

# Reformat TSV back to GFF
cut -f1-9 $TEMP/gencode.pcoding_tss_uid.tsv > $TEMP/gencode.pcoding_tss_uid.gff

# Convert to BED
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed $TEMP/gencode.pcoding_tss_uid.gff -o $TEMP/gencode.pcoding_tss_uid.bed


## =====Associate CoPRO TSS with annotations (~7 sec)=====
# TSS annotations have multiple per gene (ENSG), CoPRO includes TSS agnostic of gene annotations. Link the two

# Sort annotations and de novo TSS
bedtools sort -i $TEMP/CoPRO-capped.bed > $TEMP/CoPRO-capped_SORT.bed
bedtools sort -i $TEMP/gencode.pcoding_tss_uid.bed > $TEMP/gencode_SORT.bed

# Get closest for each CoPRO peak
bedtools closest -s -D a -a $TEMP/CoPRO-capped_SORT.bed -b $TEMP/gencode_SORT.bed > $TEMP/CoPRO-gencode_CLOSEST.bed

# Filter out distances > 1000bp away or if no association at all (edges of chrom)
awk 'BEGIN{OFS="\t";FS="\t"}{if ($13 < 1000 && $13 > -1000 && $7!=".") print}' $TEMP/CoPRO-gencode_CLOSEST.bed > $TEMP/CoPRO-gencode_CLOSEST_FILTER-Distance.bed

# Restructure columns to BED6 with extra info in right columns
#	ID    - uid
#	SCORE - "SCORE" (to update next w/ Name)
#	7     - stdev
#   8     - tag count
#   9     - gene len
#	10    - dist CoPRO to Annotation
awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3,$10,"SCORE",$6, $4,$5,$11,$13}' $TEMP/CoPRO-gencode_CLOSEST_FILTER-Distance.bed > $TEMP/CoPRO_uid-SCORE_info.bed

# Build map for uid <--> gname from geneinfo
cut -f 1,5 $TEMP/gencode.pcoding_tss_uid.geneinfo > $TEMP/uid_to_gname.txt

# Update uid --> gname
perl $UPDATES $TEMP/CoPRO_uid-SCORE_info.bed $TEMP/uid_to_gname.txt $TEMP/CoPRO_uid-Name_info.bed

# TODO: Filter minimum tag count? - No, take care of that in Genetrack

# Dedup to keep largest tag count
sort -rnk8,8 $TEMP/CoPRO_uid-Name_info.bed | sort -u -k5,5 > $TEMP/CoPRO_uid-Name_info_DEDUP.bed

# Reorganize columns to a BED6
awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3,$5,$7,$6}' $TEMP/CoPRO_uid-Name_info_DEDUP.bed > $TEMP/CoPRO_NonZeroTags.bed


## =====Organize remaining TSS annotations (~5 min)=====

# Filter gene info to exclude genes already called in NonZeroTags for candidate "ZeroTag" TSS annotations
perl $FILTERL $TEMP/gencode.pcoding_tss_uid.geneinfo <(cut -f5 $TEMP/CoPRO_uid-Name_info_DEDUP.bed) 4 remove $TEMP/ZeroTagCandidates.geneinfo

# Filter candidate "ZeroTag" TSS annotations further to exclude genes with more than one ENSG using the duplicate filter (ENSG_Filter.txt)
perl $FILTERL $TEMP/ZeroTagCandidates.geneinfo <(cut -f1 ENSG_Filter.txt) 3 remove $TEMP/ZeroTagCandidates_wFilter.geneinfo

# Filter original annotations to keep candidate "ZeroTag" TSSs uids
perl $FILTERL $TEMP/gencode.pcoding_tss_uid.bed <(cut -f1 $TEMP/ZeroTagCandidates_wFilter.geneinfo) 3 keep $TEMP/ZeroTagCandidates.bed

# Paste transcript length to 7th column
paste $TEMP/ZeroTagCandidates.bed <(cut -f5 $TEMP/ZeroTagCandidates.bed) > $TEMP/ZeroTagCandidates_tlen.bed

# Update uid --> gene name
perl $UPDATES $TEMP/ZeroTagCandidates_tlen.bed $TEMP/uid_to_gname.txt $TEMP/ZeroTagCandidates_uid-Name_tlen.bed

# Keep longest transcript annotations per gene name
sort -rnk7,7 $TEMP/ZeroTagCandidates_uid-Name_tlen.bed \
	| sort -u -k5,5 \
	> $TEMP/ZeroTag_uid-Name_tlen.bed

# Reorganize columns to a BED6
awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3,$5,"0",$6}' $TEMP/ZeroTag_uid-Name_tlen.bed > $TEMP/ZeroTags.bed


## =====Add TSS expression scores (~45 sec)=====

[ -d $KREBS/2000bp ] || mkdir -p $KREBS/2000bp

# Merge zero and non-zero tags with a random shuffle
cat $TEMP/CoPRO_NonZeroTags.bed $TEMP/ZeroTags.bed | shuf > $TEMP/TSS.bed

# Re-expand 200bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 200 $TEMP/TSS.bed -o $TEMP/TSS_200bp.bed

# Re-pileup CoPRO data
java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -2 $TEMP/TSS_200bp.bed $CAPPED -M $TEMP/Capped_READ2_TSS_200bp

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


## =====Sort TSS by CpG island length (~13 sec)=====

# Peak-align CpG islands on TSS
java -jar $SCRIPTMANAGER peak-analysis peak-align-ref $CPG $KREBS/2000bp/TSS_GROUP-Expressed_SORT-Expression_2000bp.bed -o $TEMP/TSS_CpG

# Sort BED
java -jar $SCRIPTMANAGER coordinate-manipulation sort-bed -x 900 2000 $KREBS/TSS_GROUP-Expressed_SORT-Expression.bed $TEMP/TSS_CpG_combined.cdt -o $TEMP/TSS_CpG_SORTED

# Move to global RefPT directory and format name
mv $TEMP/TSS_CpG_SORTED.bed $KREBS/TSS_GROUP-Expressed_SORT-CpG.bed

# Expand 2000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/TSS_GROUP-Expressed_SORT-CpG.bed -o $KREBS/2000bp/TSS_GROUP-Expressed_SORT-CpG_2000bp.bed