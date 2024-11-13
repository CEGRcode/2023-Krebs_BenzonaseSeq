#!/bin/bash

# Make weblogos for every JASPAR MEME file in `data/JASPAR/*.meme`

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
WRK=/storage/group/bfp2/default/owl5022-OliviaLang/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
###

# Dependencies
# - ceqlogo

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Fill in placeholder constants with your directories
PWM=$WRK/../data/JASPAR
OUTDIR=$WRK/Library/WebLogos

# Set up output directories
[ -d $OUTDIR ] || mkdir -p $OUTDIR

# Loop through the PWM files
for PWMFILE in $PWM/*.meme;
do
    BASE=`basename $PWMFILE ".meme"`
    # Generate logo
    ceqlogo -i $PWMFILE -m 1 -o $OUTDIR/$BASE\_logo.eps -f EPS
done
