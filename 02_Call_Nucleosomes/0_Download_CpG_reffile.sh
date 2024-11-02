#!/bin/bash

REFPT=../data/RefPT-Other
CPGTSV=$REFPT/CpGIslands.tsv.gz
CPGBED=$REFPT/CpG_Islands.bed

# =====UCSC Downloads=====

# Download ../data/RefPT-Other/CpGIslands.tsv
#  1. Navigate to https://genome.ucsc.edu/cgi-bin/hgTables
#  2. Select dataset
#    - Clade: "Mammal"
#    - Genome: "Human"
#    - Assembly: "Dec. 2013 (GRCh38/hg38)"
#    - Group: "Regulation"
#    - Track: "CpG Islands"
#    - Table: "cpgIslandExt"
#  3. Define region of interest
#    - Region: "Genome"
#  4. Click "Get output"
#  5. Save file to ../data/RefPT-Other/CpGIslands.tsv

# Download BigWig for more dynamic signal
# wget -c -O ../data/RefPT-Other/gc5Base.bw https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/gc5BaseBw/gc5Base.bw

# =====Reformat CpG Islands into BED format=====

# Reorganize columns into BED6
gzip -dc $CPGTSV | sed '1d' | awk 'BEGIN{OFS="\t";FS="\t"}{print $2,$3,$4,$5,$6,"."}' > $REFPT/CpGIslands.bed


