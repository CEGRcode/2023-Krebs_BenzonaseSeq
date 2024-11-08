#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 3:00:00
#SBATCH -A open
#SBATCH -o logs/1_Benzonase_Peak_Calling.log.out
#SBATCH -e logs/1_Benzonase_Peak_Calling.log.err

# Perform peak calling of various nucleosome particle sized fragment populations.
# see 02_Genetrack/01_BAM_to_scIDX.sh
# see 02_Genetrack/02_convert_scIDX_for_GeneTrack.sh
# see 02_Genetrack/03_scIDX_to_Genetrack.sh
# see 02_Genetrack/04_Genetrack_bedfiles.sh (partial)
# see 02_Genetrack/06_Genetrack_bedfiles.sh (partial - assign unique ids)

# data/RefPT-Krebs
#   |--BNase-Nucleosomes.bed

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
###

# Dependencies
# - java
# - python2
# - numpy (py2)

set -exo
module load anaconda3
module load bedtools
# conda create -p /storage/group/bfp2/default/owl5022-OliviaLang/conda/genetrack -c conda-forge -c bioconda python=2.7 numpy
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/genetrack

# Inputs and outputs
BAMFILE=$WRK/../data/BAM/BNase-seq_50U-10min_merge_hg38.bam
BLACKLIST=$WRK/../data/hg38_files/ENCFF356LFX_hg38_exclude.bed

# Script shortcuts
GENETRACK=../bin/genetrack_v2.py
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d SCIDX ] || mkdir SCIDX
[ -d sub ] || mkdir sub
[ -d tet ] || mkdir tet
[ -d hex ] || mkdir hex
[ -d nuc ] || mkdir nuc
[ -d sup ] || mkdir sup
[ -d AllParticles ] || mkdir AllParticles

# BAM to scIDX
echo "BAM to scIDX..."
time1=$(date +%s)
java -jar $SCRIPTMANAGER bam-format-converter bam-to-scidx -m -p -x 54         $BAMFILE -o SCIDX/sub.tab &
java -jar $SCRIPTMANAGER bam-format-converter bam-to-scidx -m -p -n 55  -x 91  $BAMFILE -o SCIDX/tet.tab &
java -jar $SCRIPTMANAGER bam-format-converter bam-to-scidx -m -p -n 92  -x 127 $BAMFILE -o SCIDX/hex.tab &
java -jar $SCRIPTMANAGER bam-format-converter bam-to-scidx -m -p -n 128 -x 164 $BAMFILE -o SCIDX/nuc.tab &
java -jar $SCRIPTMANAGER bam-format-converter bam-to-scidx -m -p -n 165        $BAMFILE -o SCIDX/sup.tab &
wait
time2=$(date +%s)
echo "Elapsed Time (sub): $(($time2-$time1)) seconds"

# Reformat combined value to "forward" and zero out "reverse"
echo "Reformat..."
time1=$(date +%s)
for NAME in "sub" "tet" "hex" "nuc" "sup";
do
	[ -d $NAME ] || mkdir $NAME
	awk 'BEGIN {getline; print $0; OFS="\t"; FS="\t";}{
		print $1, $2, $3+$4, "0", $5
	}' SCIDX/$NAME.tab > $NAME/formatted.tab &
done
wait
time2=$(date +%s)
echo "Elapsed Time (sub): $(($time2-$time1)) seconds"
# 73 sec

# GeneTrack
echo "GeneTrack..."
time1=$(date +%s)
python2 $GENETRACK -s 10 -e 20 -F 0  sub/formatted.tab &
python2 $GENETRACK -s 20 -e 40 -F 5  tet/formatted.tab &
python2 $GENETRACK -s 30 -e 60 -F 6  hex/formatted.tab &
python2 $GENETRACK -s 40 -e 80 -F 5  nuc/formatted.tab &
python2 $GENETRACK -s 50 -e 100 -F 3 sup/formatted.tab &
wait
time2=$(date +%s)
echo "Elapsed Time (sub): $(($time2-$time1)) seconds"
# <1 hour (3562 sec)

# Convert to BED
echo "Convert GFF to BED..."
time1=$(date +%s)
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed sub/genetrack_s10e20/formatted_s10e20.gff       -o sub/genetrack_output.bed &
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed tet/genetrack_s20e40F5/formatted_s20e40F5.gff   -o tet/genetrack_output.bed &
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed hex/genetrack_s30e60F6/formatted_s30e60F6.gff   -o hex/genetrack_output.bed &
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed nuc/genetrack_s40e80F5/formatted_s40e80F5.gff   -o nuc/genetrack_output.bed &
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed sup/genetrack_s50e100F3/formatted_s50e100F3.gff -o sup/genetrack_output.bed &
wait
time2=$(date +%s)
echo "Elapsed Time (sub): $(($time2-$time1)) seconds"
# 6 min (377 sec)

# Expand BED (max DNA frag size for indicated nuc peak except supraoctasomes expanded to 200bp)
echo "Expand BED..."
time1=$(date +%s)
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 54  sub/genetrack_output.bed -o sub/expanded.bed &
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 91  tet/genetrack_output.bed -o tet/expanded.bed &
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 127 hex/genetrack_output.bed -o hex/expanded.bed &
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 164 nuc/genetrack_output.bed -o nuc/expanded.bed &
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 200 sup/genetrack_output.bed -o sup/expanded.bed &
wait
time2=$(date +%s)
echo "Elapsed Time (sub): $(($time2-$time1)) seconds"
# 49 sec

# Filter to exclude lines expanded into negative value coordinate indexes & add 7th column with particle type indicated
echo "Filter out negative coordinates and build ID..."
time1=$(date +%s)
awk '{OFS="\t"}{FS="\t"}{if ($2>0) {$4="Subtetra-" NR "-" $1"_"$2"_"$3"_"$4; print;}}' sub/expanded.bed > AllParticles/Subtetra.bed &
awk '{OFS="\t"}{FS="\t"}{if ($2>0) {$4="Tetra-"    NR "-" $1"_"$2"_"$3"_"$4; print;}}' tet/expanded.bed > AllParticles/Tetra.bed &
awk '{OFS="\t"}{FS="\t"}{if ($2>0) {$4="Hex-"      NR "-" $1"_"$2"_"$3"_"$4; print;}}' hex/expanded.bed > AllParticles/Hex.bed &
awk '{OFS="\t"}{FS="\t"}{if ($2>0) {$4="Nuc-"      NR "-" $1"_"$2"_"$3"_"$4; print;}}' nuc/expanded.bed > AllParticles/Nucleosome.bed &
awk '{OFS="\t"}{FS="\t"}{if ($2>0) {$4="Supraoct-" NR "-" $1"_"$2"_"$3"_"$4; print;}}' sup/expanded.bed > AllParticles/Supraoct.bed &
wait
time2=$(date +%s)
echo "Elapsed Time (sub): $(($time2-$time1)) seconds"
# 41 sec

# Merge all labeled nucleosome particles
cat AllParticles/Supraoct.bed AllParticles/Nucleosome.bed \
	AllParticles/Hex.bed AllParticles/Subtetra.bed \
	AllParticles/Tetra.bed | bedtools sort -i \
	> AllParticles/Merged_Redundant_Nucleosome-Particles.bed #K562_redundant_master_nucleosomal_particles.bed
wc -l AllParticles/Merged_Redundant_Nucleosome-Particles.bed

# Create a peak file for just the full-length nucleosomes
KREBS=../data/RefPT-Krebs
[ -d $KREBS ] || mkdir $KREBS
bedtools sort -i AllParticles/Nucleosome.bed > $KREBS/BNase-Nucleosomes.bed
wc -l $KREBS/BNase-Nucleosomes.bed