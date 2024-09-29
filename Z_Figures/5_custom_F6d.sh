#!/bin/bash

# (F6d) Calculate density violin plots for active histone mod(s) w/ respepct to their base histone.
# see 04_Figures/Fig_6d.sh

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/Z_Figures
###

# Dependencies
# - java
# - python
# - pandas
# - seaborn

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
BED=PlusOneDyad_SORT-Expression_2000bp
CDIR=../X_Bulk_Processing/Library/$BED/CDT/
FDIR=../data/BAM/NormalizationFactors
TEMP=temp-F6d

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar
VIOLIN=$WRK/../bin/make_violin_plot.py

[ -d F6/d ] || mkdir -p F6/d
[ -d $TEMP ] || mkdir $TEMP


for TARGET in "H3" "H3K4me3" "H3K9ac" "H3K27ac";
do
	BAM=BNase-ChIP_${TARGET}_merge_hg19
	FACTOR=`grep '^Scaling factor' $FDIR/$BAM\_TotalTag_ScalingFactors.out | awk '{print $3}'`

	# Calculate per-site Proximal
	gzip -dc $CDIR/$BAM\_$BED\_5read1-MAX80_sense.cdt.gz \
		| cut -f1,2,920-999 > $TEMP/$TARGET\_Proximal_NOSCALE.cdt
	java -jar $SCRIPTMANAGER read-analysis scale-matrix -s $FACTOR $TEMP/$TARGET\_Proximal_NOSCALE.cdt -o $TEMP/$TARGET\_Proximal.cdt

	# Calculate per-site Distal
	gzip -dc $CDIR/$BAM\_$BED\_5read1-MAX80_anti.cdt.gz \
		| cut -f1,2,1003-1082 > $TEMP/$TARGET\_Distal_NOSCALE.cdt
	java -jar $SCRIPTMANAGER read-analysis scale-matrix -s $FACTOR $TEMP/$TARGET\_Distal_NOSCALE.cdt -o $TEMP/$TARGET\_Distal.cdt
done

# Merge each (Proximal/Distal) occupancy pair into one tab-delimited file
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K4me3-H3_Proximal.tab $TEMP/H3K4me3_Proximal.cdt $TEMP/H3_Proximal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K4me3-H3_Distal.tab   $TEMP/H3K4me3_Distal.cdt   $TEMP/H3_Distal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K9ac-H3_Proximal.tab  $TEMP/H3K9ac_Proximal.cdt  $TEMP/H3_Proximal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K9ac-H3_Distal.tab    $TEMP/H3K9ac_Distal.cdt    $TEMP/H3_Distal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K27ac-H3_Proximal.tab $TEMP/H3K27ac_Proximal.cdt $TEMP/H3_Proximal.cdt
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -o $TEMP/H3K27ac-H3_Distal.tab   $TEMP/H3K27ac_Distal.cdt   $TEMP/H3_Distal.cdt

# Calculate log2 density at each coordinate given the pair of values
sed '1d' $TEMP/H3K4me3-H3_Proximal.tab | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K4me3-H3_Proximal"}' > $TEMP/H3K4me3-H3_Proximal.density
sed '1d' $TEMP/H3K4me3-H3_Distal.tab   | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K4me3-H3_Distal"}'   > $TEMP/H3K4me3-H3_Distal.density
sed '1d' $TEMP/H3K9ac-H3_Proximal.tab  | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K9ac-H3_Proximal"}'  > $TEMP/H3K9ac-H3_Proximal.density
sed '1d' $TEMP/H3K9ac-H3_Distal.tab    | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K9ac-H3_Distal"}'    > $TEMP/H3K9ac-H3_Distal.density
sed '1d' $TEMP/H3K27ac-H3_Proximal.tab | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K27ac-H3_Proximal"}' > $TEMP/H3K27ac-H3_Proximal.density
sed '1d' $TEMP/H3K27ac-H3_Distal.tab   | awk 'BEGIN {OFS="\t"}{z = (log(($2+1)/($3+1))/log(2)); print $1,z,"H3K27ac-H3_Distal"}'   > $TEMP/H3K27ac-H3_Distal.density

# Compile density info
cat $TEMP/*.density | gzip > F6/d/DensityInfo.tab.gz

# Generate violin plot
python $VIOLIN -i <(gzip -dc F6/d/DensityInfo.tab | cut -f2,3) -o F6/d/DensityInfo.svg \
	--width 8 --height 4 --preset2 \
	--title "Density at +1 nucleosome" \
	--xlabel "modification" --ylabel "Density (log2)"

# Clean-up
rm -r $TEMP
