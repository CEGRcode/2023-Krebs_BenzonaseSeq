# purpose - take MNase or Benzonase genetrack called peak nearest to each TSS, calculate dyad to dyad distance, output to violin plot; overlap with CpG vs. non-CpG?

# usage
# qq
#
# example
#
# 'qq'

#output
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/07_dyad_distance_v2

#set files with nearest peak
TSS_IDs=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/04_BNase_plus1_05a_240711/K562_CoPRO-expressed_Gene-refSeqTSS_2000bp_RNAsortOrder.bed

Plus1_BN_peaks=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/04_BNase_plus1_05a_240711/K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet.bed
Minus1_BN_peaks=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/04_BNase_plus1_05a_240711/K562_Minus1_SORTdistanceToTSS_nonRedOct_Hex_Tet.bed
Plus1_MN_peaks=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/04_MNase_plus1_05b_240710/K562_Plus1_SORTdistanceToTSS_MNase.bed
Minus1_MN_peaks=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/04_MNase_plus1_05b_240710/K562_Minus1_SORTdistanceToTSS_MNase.bed

#set genome and human.hg19.genome file
GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/referenceDATA_Will/GENOMES/hg19.fa
HG19_GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/files/human.hg19.genome

#set blacklist
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed

#set scriptmanager and job
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
JOB=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/jobs/sum_Col_CDT.pl
PLOT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/jobs/violin_plots_dyad_distance_240711.py

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

mkdir -p $OUTPUT

JOBSTATS="#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=24GB
#SBATCH --time=4:00:00
#SBATCH --partition=open

module load anaconda #configures shell to use conda activate
conda activate plot"

#set output file names

#set output directory
cd $OUTPUT

#prep initial bedfiles so that columns can be matched up (both plus1)
cat /storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/04_BNase_plus1_05a_240711/K562_Plus1_SORTdistanceToTSS_nonRedOct_Hex_Tet.bed | awk -F "," '{split($0,a, ""); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | awk '{print $4"\t"$6}' > Benz_peaks.tab
cat /storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/04_MNase_plus1_05b_240710/K562_Plus1_SORTdistanceToTSS_MNase.bed | awk -F "," '{split($0,a, ""); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | awk '{print $4"\t"$6}' > MNase_peaks.tab

#take column 4 from TSS_IDs file
cat '/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/04_BNase_plus1_05a_240711/K562_CoPRO-expressed_Gene-refSeqTSS_2000bp_RNAsortOrder.bed' | cut -f4 > TSS_IDs.tab

#match column1 in "Benz_peaks.tab" to column 1 of TSS IDs
awk 'NR==FNR{A[$1]=$0; next} ($1 in A){print A[$1], $2}; !($1 in A){print $0, "NA"};' FS="\t" OFS="\t" Benz_peaks.tab TSS_IDs.tab > TSS_Benz.tab

#match column1 in "Benz_peaks.tab" to column 1 of TSS IDs
awk 'NR==FNR{A[$1]=$0; next} ($1 in A){print A[$1], $2}; !($1 in A){print $0, "NA"};' FS="\t" OFS="\t" MNase_peaks.tab TSS_IDs.tab > TSS_MNase.tab

#bring both files, only include rows with no "NA" in column 2 or 4, then do dyad to dyad calculation (benz dyad - MNase dyad)
paste TSS_Benz.tab TSS_MNase.tab | awk '{if ($2 !="NA" && $4 !="NA") print $1"\t"$2"\t"$3"\t"$4}' | awk -F"=" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' | awk '{print $0"\t"(sqrt(($3-$6)*($3-$6)))}' | awk '{print $7"\t""dyad_distance"}' | awk '{if ($1 <=82) print $1"\t"$2}' > distance.tab

#make violin plot wity modified violin_plots_mod.py
python $PLOT -i distance.tab -o distance.svg
# finish script
