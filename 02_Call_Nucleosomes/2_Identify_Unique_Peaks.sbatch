#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 0:10:00
#SBATCH -A open
#SBATCH -o logs/2_Identify_Unique_Peaks.log.out
#SBATCH -e logs/2_Identify_Unique_Peaks.log.err

# Match sub-nucleosome sized peaks to nucleosome sized ones
# see 02_Genetrack/05_Genetrack_bedfiles.sh
# see 02_Genetrack/06_Genetrack_bedfiles.sh (partial - assign unique ids)

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
###

# Dependencies
# - bedtools
# - python

set -exo
module load anaconda3
module load bedtools
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
BAMFILE=$WRK/../data/BAM/BNase-seq_50U-10min_merge_hg38.bam
BLACKLIST=$WRK/../data/hg38_files/ENCFF356LFX_hg38_exclude.bed
KREBS=../data/RefPT-Krebs/
UNIQUE=UniqueParticles

SUBTETRASOMES=AllParticles/Subtetra.bed
TETRASOMES=AllParticles/Tetra.bed
HEXASOMES=AllParticles/Hex.bed
NUCLEOSOMES=AllParticles/Nucleosome.bed
SUPRAOCTASOMES=AllParticles/Supraoct.bed

# Script shortcuts
GENETRACK=../bin/genetrack_v2.py

# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d $KREBS ] || mkdir $KREBS
[ -d $UNIQUE ] || mkdir $UNIQUE

# Get unique Hexasomes (those that do NOT overlap nucleosomes)
bedtools intersect -v -a $HEXASOMES -b $NUCLEOSOMES \
	| awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,"uHex-" NR "-" $1"_"$2"_"$3"_"$4,$5,$6}' \
	> $UNIQUE/uHex.bed
# Nucleosomes + unique Hexasomes
cat $NUCLEOSOMES $UNIQUE/uHex.bed > $UNIQUE/Nucleosome_uHex.bed

# Get unique Tetrasomes (those that do NOT overlap nucleosomes or hexasomes)
bedtools intersect -v -a $TETRASOMES -b $UNIQUE/Nucleosome_uHex.bed \
	| awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,"uTetra-" NR "-" $1"_"$2"_"$3"_"$4,$5,$6}' \
	> $UNIQUE/uTetra.bed
# Nucleosomes + unique Hexasomes + unique Tetrasomes (not nucleosomes, hexasomes)
cat $UNIQUE/Nucleosome_uHex.bed $UNIQUE/uTetra.bed > $UNIQUE/Nucleosome_uHex_uTetra.bed

# Get unique subtetrasomes (those that do NOT overlap nucleosomes, hexasomes, or tetrasomes)
bedtools intersect -v -a $SUBTETRASOMES -b $UNIQUE/Nucleosome_uHex_uTetra.bed \
	| awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,"uSubtetra-" NR "-" $1"_"$2"_"$3"_"$4,$5,$6}' \
	> $UNIQUE/uSubtetra.bed
# Nucleosomes + unique Hexasomes + unique Tetrasomes + unique subtetrasomes (not nucleosomes, hexasomes, tetrasomes)
cat $UNIQUE/Nucleosome_uHex_uTetra.bed $UNIQUE/uSubtetra.bed > $UNIQUE/Nucleosome_uHex_uTetra_uSubtetra.bed

# Get unique supraoctasomes (those that do NOT overlap nucleosomes, hexasomes, tetrasomes, or subtetrasomes)
bedtools intersect -v -a $SUPRAOCTASOMES -b $UNIQUE/Nucleosome_uHex_uTetra_uSubtetra.bed \
	| awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,"uSupraoct-" NR "-" $1"_"$2"_"$3"_"$4,$5,$6}' \
	> $UNIQUE/uSupraoct.bed
# Nucleosomes + unique Hexasomes + unique Tetrasomes + unique subtetrasomes + unique supraoctasomes (not nucleosomes, hexasomes, tetrasomes, or subtetrasomes)
cat $UNIQUE/Nucleosome_uHex_uTetra_uSubtetra.bed $UNIQUE/uSupraoct.bed > $UNIQUE/Merged_Nonredundant_particles.bed

# Count merged nonredundant particles
NSITES=`wc -l $UNIQUE/Merged_Nonredundant_particles.bed | awk '{print $1}'`
echo "Complete overlap merge of all particles: $NSITES sites"
