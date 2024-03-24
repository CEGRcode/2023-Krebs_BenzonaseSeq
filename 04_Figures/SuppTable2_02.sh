# purpose - this code determines the total number of full-length nucleosomes at the +1 position, makesa bedfile and then uses this bedfile to calculate the number of tags (representing SNs) in the proximal and distal SNs regions (based on midpoints in Fig. 4). The final files gives the number of rows with proximal SNs, distal SNs, or both at +1 positions that are based on full-length nucleosomes.

# usage
# qq
#
# example
# purpose - qq

# usage
# qq
#
# example
#
# 'qq'

#set bedfiles with Plus1 (downstream) nucleosome positions based on non-redundant bedfile-based calls, sorted by decreasing RNA-expression
Plus1_BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_plus1_minus1/output_v2_NonRed_Oct_Hex_Tet_230825/K562_Plus1_SORTbyRNAexp_nonRedOct_Hex_Tet.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/Fig4_transcription/H3K4me3_SNs_Plus1_fullLength_JEK_240121

#set bam library file to BI_rep1 **testing with subsampled master BAM file
BAM1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/25861_25869_25965_25972_28805_28809_Benz_0sonicCycles_BX_H3K4me3_master.bam

#set blacklist
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed
HG19_GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/files/human.hg19.genome

#set scriptmanager and job
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
JOB=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/jobs/sum_Col_CDT.pl
JOB_ROW=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig6_Subnucleosomes/job/sum_Row_CDT.pl

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
#SBATCH --time=0:30:00
#SBATCH --partition=open

source ~/.bashrc #configures shell to use conda activate
module load anaconda
conda activate bioinfo"

#set output file names
Plus1_labeled=$(echo $Plus1_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_particles_labeled.bed"}')
NUMBERS=$(echo "$Plus1_BEDFILE" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_rowsNumber.tab"}')
BEDFILE=$(echo $Plus1_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_FullLengthNuc_Only.bed"}')
BEDFILE_2bp=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_2bp.bed"}')
Plus1_proximal_SN=$(echo $BEDFILE_2bp| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_proximal_SN.bed"}')
Plus1_distal_SN=$(echo $BEDFILE_2bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_distal_SN.bed"}')
Plus1_proximal_SN_a=$(echo $Plus1_proximal_SN | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
Plus1_distal_SN_a=$(echo $Plus1_distal_SN | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$Plus1_proximal_SN_a" | awk -F. '{print $1"_proximal_midpoint.out"}')
CDT1=$(echo "$BAM1a""_""$Plus1_proximal_SN_a" | awk -F. '{print $1"_proximal"}')
CDT1_gz=$(echo "$BAM1a""_""$Plus1_proximal_SN_a" | awk -F. '{print $1"_proximal_combined.cdt.gz"}')
CDT1_unzipped=$(echo "$BAM1a""_""$Plus1_proximal_SN_a" | awk -F. '{print $1"_proximal_combined.cdt"}')
OUT2=$(echo "$BAM1a""_""$Plus1_distal_SN_a" | awk -F. '{print $1"_distal_midpoint.out"}')
CDT2=$(echo "$BAM1a""_""$Plus1_distal_SN_a" | awk -F. '{print $1"_distal"}')
CDT2_gz=$(echo "$BAM1a""_""$Plus1_distal_SN_a" | awk -F. '{print $1"_distal_combined.cdt.gz"}')
CDT2_unzipped=$(echo "$BAM1a""_""$Plus1_distal_SN_a" | awk -F. '{print $1"_distal_combined.cdt"}')
Target_proximal=$(echo "$BAM1a""_""$Plus1_proximal_SN_a" | awk -F. '{print $1"_RowCount.tab"}')
Target_distal=$(echo "$BAM1a""_""$Plus1_distal_SN_a" | awk -F. '{print $1"_RowCount.tab"}')
NUMBERS_SNs=$(echo "$Plus1_BEDFILE" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_rowsNumber_proximal_distal_both.tab"}')

sampleID=Fig4_H3K4me3_SNs_+1_full_length_JEK_240120.slurm
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#make new file a new variable" >> $sampleID
echo "cat $Plus1_BEDFILE | awk '{if (\$2>0) \$7=\$3-\$2; else \$7=\"\"; print}' | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6\"\t\"\$7}' > $Plus1_labeled" >> $sampleID
echo "#get full number of rows of bedfile" >> $sampleID
echo "cat $Plus1_BEDFILE | wc -l | awk '{print \"total\"\"\t\"\$1}' > $NUMBERS" >> $sampleID
echo "#get number of Hexasomes (126 bp), full-length nuclesosomes (164 bp), and tetrasome (90 bp)." >> $sampleID
echo "cat $Plus1_labeled | cut -f 7 | sort | uniq -c | awk '{print \$2\"bp\"\"\t\"\$1}' >> $NUMBERS" >> $sampleID
echo "#make bedfile of NUC only at +1 bp" >> $sampleID
echo "cat $Plus1_labeled | awk '{if (\$7==\"164\") print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6\"\t\"\$7; else \"\"}' > $BEDFILE" >> $sampleID
echo "##determine the number of full-length nucleosome that have tags for SNs" >> $sampleID
echo "#expand bedfile" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=2 $BEDFILE -o=$BEDFILE_2bp" >> $sampleID
echo "#get proximal (upstream) or distal (downstream) half of nucleosome relative to +1 midpoint (dyad)" >> $sampleID
echo "#math is based on 80 bp SNs for shifts: proximal SN is ~2 bp upstream of dyad so 80 + 2 = 82 for start, -2 + -2 = -4 for end; distal SN is ~4 bp downstream of dyad so -2 + -4 = -6 for start, 80 + 4 = 84 for end" >> $sampleID
echo "bedtools slop -i $BEDFILE_2bp -g $HG19_GENOME -l 82 -r -4 -s > $Plus1_proximal_SN" >> $sampleID
echo "bedtools slop -i $BEDFILE_2bp -g $HG19_GENOME -l -6 -r 84 -s > $Plus1_distal_SN" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT1 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 --max-insert=80 $Plus1_proximal_SN $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT2 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 --max-insert=80 $Plus1_distal_SN $BAM1" >> $sampleID
echo "#unzip files" >> $sampleID
echo "gunzip -c $CDT1_gz > $CDT1_unzipped" >> $sampleID
echo "gunzip -c $CDT2_gz > $CDT2_unzipped" >> $sampleID
echo "#sum tags in each above CDT file" >> $sampleID
echo "perl $JOB_ROW $CDT1_unzipped $Target_proximal" >> $sampleID
echo "perl $JOB_ROW $CDT2_unzipped $Target_distal" >> $sampleID
echo "#make a file showing how many of above full-length nucleosomes at +1 position have at least 1 tag in either proximal SN region. First only rows with matching IDs are kept." >> $sampleID
echo "paste $Plus1_proximal_SN $Target_proximal  |  awk '(\$4==\$8){print \$0}' | awk '(\$9!=0){print \$0}' | wc -l | awk '{print \"proximal_SN_region_withSNs\"\"\t\"\$1}' > $NUMBERS_SNs" >> $sampleID
echo "paste $Plus1_proximal_SN $Target_proximal  |  awk '(\$4==\$8){print \$0}' | awk '(\$9==0){print \$0}' | wc -l | awk '{print \"proximal_SN_region_withoutSNs\"\"\t\"\$1}' >> $NUMBERS_SNs" >> $sampleID
echo "paste $Plus1_distal_SN $Target_distal  |  awk '(\$4==\$8){print \$0}' | awk '(\$9!=0){print \$0}' | wc -l | awk '{print \"distal_SN_region_withSNs\"\"\t\"\$1}' >> $NUMBERS_SNs" >> $sampleID
echo "paste $Plus1_distal_SN $Target_distal  |  awk '(\$4==\$8){print \$0}' | awk '(\$9==0){print \$0}' | wc -l | awk '{print \"distal_SN_region_withoutSNs\"\"\t\"\$1}' >> $NUMBERS_SNs" >> $sampleID
echo "paste $BEDFILE_2bp $Target_proximal $Target_distal  |  awk '(\$4==\$8 && \$4==\$10){print \$0}' | awk '(\$8!=0 && \$10!=0){print \$0}' | wc -l | awk '{print \"both_SNs_present\"\"\t\"\$1}' >> $NUMBERS_SNs" >> $sampleID
echo "paste $BEDFILE_2bp $Target_proximal $Target_distal  |  awk '(\$4==\$8 && \$4==\$10){print \$0}' | wc -l | awk '{print \"total_possible_for_combined\"\"\t\"\$1}' >> $NUMBERS_SNs" >> $sampleID
echo "#remove intermediate files" >> $sampleID
echo "rm $OUT1" >> $sampleID
echo "rm $CDT1_gz" >> $sampleID
echo "rm $CDT1_unzipped" >> $sampleID
echo "rm $OUT2" >> $sampleID
echo "rm $CDT2_gz" >> $sampleID
echo "rm $CDT2_unzipped" >> $sampleID
echo "# finish script" >> $sampleID
