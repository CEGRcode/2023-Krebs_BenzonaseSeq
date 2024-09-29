#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 20:00:00
#SBATCH -A open
#SBATCH -o logs/3_Aggregate_Nucleosome_Peaks.log.out
#SBATCH -e logs/3_Aggregate_Nucleosome_Peaks.log.err

# Intersect all redundant (smaller) peaks with redundant peaks. When multiple
# smaller peaks intersect with a single larger particle, all are kept and
# tagged.
# see 02_Genetrack/06_Genetrack_bedfiles.sh

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
###

# Dependencies
# - bedtools

set -exo
module load anaconda3
module load bedtools
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
ALL=AllParticles
UNIQUE=UniqueParticles
INTERSECT=Intersect

SUB=$ALL/Subtetra.bed
TET=$ALL/Tetra.bed
HEX=$ALL/Hex.bed
NUC=$ALL/Nucleosome.bed
SUP=$ALL/Supraoct.bed

USUB=$UNIQUE/uSubtetra.bed
UTET=$UNIQUE/uTetra.bed
UHEX=$UNIQUE/uHex.bed
USUP=$UNIQUE/uSupraoct.bed

# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d $INTERSECT ] || mkdir $INTERSECT

## Intersect all particles with full-length nucleosomes
# - Use -F 1.0 (min fraction of B) to ensure complete overlap w/ nucleosome.
# - All others are not considered.
# - Max overlap is 126/90/ bp as seen in data.
# - THEN substitute "-1" with 0 in columns 8, 9, and 11.
# Make final bedfile with chr# / start / end /HEX or TET or SUB or NUC label;intersecting Nuc label; bp overlap/tag count / strand of hexasome/tet/supraoctasome
# ****if this file does not work for tag pile-up take chr# and strand from nucleosome section -> unique lines without "."

# Redundant hexasomes w/ nucleosomes.
bedtools intersect -wao -F 1.0 -a $NUC -b $HEX \
	| awk '{OFS="\t"}{FS="\t"}{
		if ($8=="-1") $8="0";
		if ($9=="-1") $9="0";
		if ($11=="-1") $11="0";
		print $7,$8,$9,$10";"$4";"$13,$11,$12
	}' > $INTERSECT/AllHex_NestedIn-Nucleosomes.bed

# Redundant tetrasomes w/ nucleosomes.
bedtools intersect -wao -F 1.0 -a $NUC -b $TET \
	| awk '{OFS="\t"}{FS="\t"}{
		if ($8=="-1") $8="0";
		if ($9=="-1") $9="0";
		if ($11=="-1") $11="0";
		print $7,$8,$9,$10";"$4";"$13,$11,$12
	}' > $INTERSECT/AllTet_NestedIn-Nucleosomes.bed

# Redundant subtetrasomes w/ nucleosomes
bedtools intersect -wao -F 1.0 -a $NUC -b $SUB \
	| awk '{OFS="\t"}{FS="\t"}{
		if ($8=="-1") $8="0";
		if ($9=="-1") $9="0";
		if ($11=="-1") $11="0";
		print $7,$8,$9,$10";"$4";"$13,$11,$12
	}' > $INTERSECT/AllSubtetra_NestedIn-Nucleosomes.bed

# Redundant supraoctasomes with nucleosomes.
# Minimum (-f) of 1.0 used so that only supraoctasomes that overlap nucleosomes by 164 bp are considered. (Ask Jordan?)
bedtools intersect -wao -f 1.0 -a $NUC -b $SUP \
	| awk '{OFS="\t"}{FS="\t"}{
		if ($8=="-1") $8="0";
		if ($9=="-1") $9="0";
		if ($11=="-1") $11="0";
		print $7,$8,$9,$10";"$4";"$13,$11,$12
	}' > $INTERSECT/AllSubtetra_WithNested-Nucleosome.bed

## Intersect all particles with non-redundant hexasomes

# Redundant tetrasomes w/ non-redundant hexasomes.
bedtools intersect -wao -F 1.0 -a $UHEX -b $TET \
	| awk '{OFS="\t"}{FS="\t"}{
		if ($8=="-1") $8="0";
		if ($9=="-1") $9="0";
		if ($11=="-1") $11="0";
		print $7,$8,$9,$10";"$4";"$13,$11,$12
	}' > $INTERSECT/AllTetra_NestedIn-uHex.bed

# Redundant subtetrasomes w/ non-redundant hexasomes.
bedtools intersect -wao -F 1.0 -a $UHEX -b $SUB \
	| awk '{OFS="\t"}{FS="\t"}{
		if ($8=="-1") $8="0";
		if ($9=="-1") $9="0";
		if ($11=="-1") $11="0";
		print $7,$8,$9,$10";"$4";"$13,$11,$12
	}' > $INTERSECT/AllSubtetra_NestedIn-uHex.bed

# Redundant supraoct w/ non-redundant hexasomes.
bedtools intersect -wao -f 1.0 -a $UHEX -b $SUP \
	| awk '{OFS="\t"}{FS="\t"}{
		if ($8=="-1") $8="0";
		if ($9=="-1") $9="0";
		if ($11=="-1") $11="0";
		print $7,$8,$9,$10";"$4";"$13,$11,$12
	}' > $INTERSECT/AllSupraoct_WithNested-uHex.bed

## Intersect all particles with non-redundant tetrasomes.

# Redundant subtetra w/ non-redundant tetrasomes.
bedtools intersect -wao -F 1.0 -a $UTET -b $SUB \
	| awk '{OFS="\t"}{FS="\t"}{
		if ($8=="-1") $8="0";
		if ($9=="-1") $9="0";
		if ($11=="-1") $11="0";
		print $7,$8,$9,$10";"$4";"$13,$11,$12
	}' > $INTERSECT/AllSubtetra_NestedIn-uTet.bed

# Redundant supraoctasomes w/ non-redundant tetrasomes.
bedtools intersect -wao -f 1.0 -a $UTET -b $SUP \
	| awk '{OFS="\t"}{FS="\t"}{
		if ($8=="-1") $8="0";
		if ($9=="-1") $9="0";
		if ($11=="-1") $11="0";
		print $7,$8,$9,$10";"$4";"$13,$11,$12
	}' > $INTERSECT/AllSupraoct_WithNested-uTet.bed

## Intersect all particles with non-redundant subtetrasomes.

# Redundant supraoctasomes w/ non-redundant subtetrasomes.
bedtools intersect -wao -f 1.0 -a $USUB -b $SUP \
	| awk '{OFS="\t"}{FS="\t"}{
		if ($8=="-1") $8="0";
		if ($9=="-1") $9="0";
		if ($11=="-1") $11="0";
		print $7,$8,$9,$10";"$4";"$13,$11,$12
	}' > $INTERSECT/AllSupraoct_WithNested-uSub.bed
