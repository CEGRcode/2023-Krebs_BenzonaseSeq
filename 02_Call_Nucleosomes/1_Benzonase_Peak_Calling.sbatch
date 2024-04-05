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

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/storage/home/owl5022/scratch/2023-Krebs_Benzonase-seq/02_Call_Nucleosomes
###

# Dependencies
# - java
# - python2
# - numpy (py2)

set -exo
module load anaconda3
module load bedtools
# conda create -n genetrack -c conda-forge -c bioconda python=2.7 numpy
source activate genetrack

# Inputs and outputs
BAMFILE=../data/BAM/BNase-seq_50U_merge_hg19.bam
BLACKLIST=../data/hg19_files/hg19_exclude.bed

# Script shortcuts
GENETRACK=../bin/genetrack_v2.py
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

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
#

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
