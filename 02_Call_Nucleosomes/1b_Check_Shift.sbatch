#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 2:00:00
#SBATCH -A open
#SBATCH -o logs/1b_Check_Shift.log.out
#SBATCH -e logs/1b_Check_Shift.log.err

# Confirm files are shifted correctly
# see 02_Genetrack/04_Genetrack_bedfiles.sh (partial)

### CHANGE ME
THREADS=4
WRK=/path/to/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
###

# Dependencies
# - java

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
GENOME=$WRK/../data/hg38_files/hg38.fa
BAMFILE=$WRK/../data/BAM/BNase-seq_50U-10min_merge_hg38.bam
BLACKLIST=$WRK/../data/hg38_files/ENCFF356LFX_hg38_exclude.bed
ODIR=ShiftCheck

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

[ -d ShiftCheck ] || mkdir ShiftCheck

for NAME in "Subtetra" "Tetra" "Hex" "Nucleosome" "Supraoct";
do
	# Shuffle and take 10K sites from each bedfile
	shuf AllParticles/$NAME.bed | head -10000 > $ODIR/$NAME\_10k.bed
	# Expand 1000bp
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $ODIR/$NAME\_10k.bed -o $ODIR/$NAME\_10k_1000bp.bed
	# Tag Pileup and confirm center of peak is at 0 bp
	# Settings: midpoint(m), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist
	java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined -N -p --cpu $THREADS \
		-M $ODIR/$NAME\_10k_1000bp_midpoint -o $ODIR/$NAME\_10k_1000bp_midpoint.out \
		$ODIR/$NAME\_10k_1000bp.bed $BAMFILE
done

# Generate composite plots
java -jar $SCRIPTMANAGER figure-generation composite-plot --title="center of subtetrasomes" $ODIR/Subtetra_10k_1000bp_midpoint.out
java -jar $SCRIPTMANAGER figure-generation composite-plot --title="center of tetrasomes" $ODIR/Tetra_10k_1000bp_midpoint.out
java -jar $SCRIPTMANAGER figure-generation composite-plot --title="center of hexasomes" $ODIR/Hex_10k_1000bp_midpoint.out
java -jar $SCRIPTMANAGER figure-generation composite-plot --title="center of nucleosomes" $ODIR/Nucleosome_10k_1000bp_midpoint.out
java -jar $SCRIPTMANAGER figure-generation composite-plot --title="center of supraoctasomes" $ODIR/Supraoct_10k_1000bp_midpoint.out
