# purpose - code for Fig. 3c showing exo stop sites at core histones for all DNA fragments (128-164 bp insert size selection). Now this version (v8) does NOT make heatmaps and uses all sites.

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

#set bedfiles with  Plus1 (downstream) nucleosome positions based on non-redundant bedfile-based calls
BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_plus1_minus1/output_v2_NonRed_Oct_Hex_Tet_230825/K562_Plus1_SORTbyRNAexp_nonRedOct_Hex_Tet.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/nre_fig3_240219/fig3c_occupancy_128to164bp_output_240304

#set bam library file to BI_rep1
BAM1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/28452_28460_Benz_0sonicCycles_BX_H2A_master.bam
BAM2=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/28453_28794_28461_28799_Benz_0sonicCycles_BX_H2A_Z_master.bam
BAM3=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/28454_28795_28462_28800_Benz_0sonicCycles_BX_H2B_master.bam
BAM4=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/25860_25868_25964_25971_28804_28808_Benz_0sonicCycles_BX_H3_master.bam
BAM5=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/28457_28797_28465_28802_Benz_0sonicCycles_BX_H4_master.bam
BAM6=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230718_MERGE/K562_benzonase-seq_master.bam


#set blacklist
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed

#set scriptmanager and job
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
JOB=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/jobs/sum_Col_CDT.pl

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

source ~/.bashrc #configures shell to use conda activate
module load anaconda
conda activate bioinfo"

#set output file names
BEDFILE_1200=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_1200bp.bed"}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1a_500=$(echo "$BAM1a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1.out"}')
CDT1d_500=$(echo "$BAM1a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1"}')
CDT1_sense_gz_500=$(echo "$BAM1a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_sense.cdt.gz"}')
CDT1_anti_gz_500=$(echo "$BAM1a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_anti.cdt.gz"}')
CDT1_sense_500=$(echo "$BAM1a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_sense.cdt"}')
CDT1_anti_500=$(echo "$BAM1a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_anti.cdt"}')
SCALE1_OUT=$(echo "$BAM1a" | awk -F. '{print $1"_ForComposite_read1"}')
SCALE1a_OUT=$(echo "$BAM1a" | awk -F. '{print $1"_ForComposite_read1_ScalingFactors.out"}')
CDT1_SCALED_COMP_sense_500=$(echo "$BAM1a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_sense.cdt"}')
CDT1_SCALED_COMP_anti_500=$(echo "$BAM1a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_anti.cdt"}')
SCALED_OUT1_sense_500=$(echo "$BAM1a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_sense.tab"}')
SCALED_OUT1_anti_500=$(echo "$BAM1a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_anti.tab"}')
SCALED_OUT1_final_500=$(echo "$BAM1a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_final.tab"}')
BAM2a=$(echo $BAM2 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT2a_500=$(echo "$BAM2a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1.out"}')
CDT2d_500=$(echo "$BAM2a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1"}')
CDT2_sense_gz_500=$(echo "$BAM2a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_sense.cdt.gz"}')
CDT2_anti_gz_500=$(echo "$BAM2a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_anti.cdt.gz"}')
CDT2_sense_500=$(echo "$BAM2a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_sense.cdt"}')
CDT2_anti_500=$(echo "$BAM2a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_anti.cdt"}')
SCALE2_OUT=$(echo "$BAM2a" | awk -F. '{print $1"_ForComposite_read1"}')
SCALE2a_OUT=$(echo "$BAM2a" | awk -F. '{print $1"_ForComposite_read1_ScalingFactors.out"}')
CDT2_SCALED_COMP_sense_500=$(echo "$BAM2a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_sense.cdt"}')
CDT2_SCALED_COMP_anti_500=$(echo "$BAM2a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_anti.cdt"}')
SCALED_OUT2_sense_500=$(echo "$BAM2a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_sense.tab"}')
SCALED_OUT2_anti_500=$(echo "$BAM2a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_anti.tab"}')
SCALED_OUT2_final_500=$(echo "$BAM2a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_final.tab"}')
BAM3a=$(echo $BAM3 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT3a_500=$(echo "$BAM3a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1.out"}')
CDT3d_500=$(echo "$BAM3a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1"}')
CDT3_sense_gz_500=$(echo "$BAM3a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_sense.cdt.gz"}')
CDT3_anti_gz_500=$(echo "$BAM3a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_anti.cdt.gz"}')
CDT3_sense_500=$(echo "$BAM3a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_sense.cdt"}')
CDT3_anti_500=$(echo "$BAM3a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_anti.cdt"}')
SCALE3_OUT=$(echo "$BAM3a" | awk -F. '{print $1"_ForComposite_read1"}')
SCALE3a_OUT=$(echo "$BAM3a" | awk -F. '{print $1"_ForComposite_read1_ScalingFactors.out"}')
CDT3_SCALED_COMP_sense_500=$(echo "$BAM3a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_sense.cdt"}')
CDT3_SCALED_COMP_anti_500=$(echo "$BAM3a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_anti.cdt"}')
SCALED_OUT3_sense_500=$(echo "$BAM3a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_sense.tab"}')
SCALED_OUT3_anti_500=$(echo "$BAM3a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_anti.tab"}')
SCALED_OUT3_final_500=$(echo "$BAM3a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_final.tab"}')
BAM4a=$(echo $BAM4 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT4a_500=$(echo "$BAM4a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1.out"}')
CDT4d_500=$(echo "$BAM4a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1"}')
CDT4_sense_gz_500=$(echo "$BAM4a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_sense.cdt.gz"}')
CDT4_anti_gz_500=$(echo "$BAM4a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_anti.cdt.gz"}')
CDT4_sense_500=$(echo "$BAM4a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_sense.cdt"}')
CDT4_anti_500=$(echo "$BAM4a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_anti.cdt"}')
SCALE4_OUT=$(echo "$BAM4a" | awk -F. '{print $1"_ForComposite_read1"}')
SCALE4a_OUT=$(echo "$BAM4a" | awk -F. '{print $1"_ForComposite_read1_ScalingFactors.out"}')
CDT4_SCALED_COMP_sense_500=$(echo "$BAM4a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_sense.cdt"}')
CDT4_SCALED_COMP_anti_500=$(echo "$BAM4a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_anti.cdt"}')
SCALED_OUT4_sense_500=$(echo "$BAM4a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_sense.tab"}')
SCALED_OUT4_anti_500=$(echo "$BAM4a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_anti.tab"}')
SCALED_OUT4_final_500=$(echo "$BAM4a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_final.tab"}')
BAM5a=$(echo $BAM5 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT5a_500=$(echo "$BAM5a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1.out"}')
CDT5d_500=$(echo "$BAM5a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1"}')
CDT5_sense_gz_500=$(echo "$BAM5a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_sense.cdt.gz"}')
CDT5_anti_gz_500=$(echo "$BAM5a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_anti.cdt.gz"}')
CDT5_sense_500=$(echo "$BAM5a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_sense.cdt"}')
CDT5_anti_500=$(echo "$BAM5a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_read1_anti.cdt"}')
SCALE5_OUT=$(echo "$BAM5a" | awk -F. '{print $1"_ForComposite_read1"}')
SCALE5a_OUT=$(echo "$BAM5a" | awk -F. '{print $1"_ForComposite_read1_ScalingFactors.out"}')
CDT5_SCALED_COMP_sense_500=$(echo "$BAM5a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_sense.cdt"}')
CDT5_SCALED_COMP_anti_500=$(echo "$BAM5a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_anti.cdt"}')
SCALED_OUT5_sense_500=$(echo "$BAM5a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_sense.tab"}')
SCALED_OUT5_anti_500=$(echo "$BAM5a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_anti.tab"}')
SCALED_OUT5_final_500=$(echo "$BAM5a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled_read1_final.tab"}')
BAM6a=$(echo $BAM6 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT6a_500=$(echo "$BAM6a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT6d_500=$(echo "$BAM6a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite"}')
CDT6e_500=$(echo "$BAM6a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT6f_500=$(echo "$BAM6a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
SCALE6_OUT=$(echo "$BAM6a" | awk -F. '{print $1"_ForComposite"}')
SCALE6a_OUT=$(echo "$BAM6a" | awk -F. '{print $1"_ForComposite_ScalingFactors.out"}')
CDT6_SCALED_COMP_500=$(echo "$BAM6a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT6_500=$(echo "$BAM6a""_""$BEDFILE_1200" | awk -F. '{print $1"_ForComposite_scaled.tab"}')

sampleID=fig3c_occupancy_128to164bp_v8_240304.slurm
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#expand bedfiles by 1200 bp" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1200 $BEDFILE -o=$BEDFILE_1200" >> $sampleID
echo "#scale output files: options - total tag scaling -t; *this is same regardless of bedfile" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE1_OUT $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE2_OUT $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE3_OUT $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE4_OUT $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE5_OUT $BAM5" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE6_OUT $BAM6" >> $sampleID
echo "#prep files for plotter" >> $sampleID
echo "#do another tag-pileUp (output is input directory). Settings: midpoint(m), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT1d_500 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1a_500 --min-insert=128 --max-insert=164 $BEDFILE_1200 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT2d_500 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2a_500 --min-insert=128 --max-insert=164 $BEDFILE_1200 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT3d_500 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3a_500 --min-insert=128 --max-insert=164 $BEDFILE_1200 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT4d_500 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4a_500 --min-insert=128 --max-insert=164 $BEDFILE_1200 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT5d_500 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5a_500 --min-insert=128 --max-insert=164 $BEDFILE_1200 $BAM5" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT6d_500 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT6a_500 $BEDFILE_1200 $BAM6" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT1_sense_gz_500 > $CDT1_sense_500" >> $sampleID
echo "gunzip -c $CDT1_anti_gz_500 > $CDT1_anti_500" >> $sampleID
echo "gunzip -c $CDT2_sense_gz_500 > $CDT2_sense_500" >> $sampleID
echo "gunzip -c $CDT2_anti_gz_500 > $CDT2_anti_500" >> $sampleID
echo "gunzip -c $CDT3_sense_gz_500 > $CDT3_sense_500" >> $sampleID
echo "gunzip -c $CDT3_anti_gz_500 > $CDT3_anti_500" >> $sampleID
echo "gunzip -c $CDT4_sense_gz_500 > $CDT4_sense_500" >> $sampleID
echo "gunzip -c $CDT4_anti_gz_500 > $CDT4_anti_500" >> $sampleID
echo "gunzip -c $CDT5_sense_gz_500 > $CDT5_sense_500" >> $sampleID
echo "gunzip -c $CDT5_anti_gz_500 > $CDT5_anti_500" >> $sampleID
echo "gunzip -c $CDT6e_500 > $OUTPUT/$CDT6f_500" >> $sampleID
echo "#scale CDT files data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_COMP_sense_500 --scaling-factor=\$(cat $SCALE1a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT1_sense_500" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_COMP_anti_500 --scaling-factor=\$(cat $SCALE1a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT1_anti_500" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_COMP_sense_500 --scaling-factor=\$(cat $SCALE2a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT2_sense_500" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_COMP_anti_500 --scaling-factor=\$(cat $SCALE2a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT2_anti_500" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED_COMP_sense_500 --scaling-factor=\$(cat $SCALE3a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT3_sense_500" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED_COMP_anti_500 --scaling-factor=\$(cat $SCALE3a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT3_anti_500" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED_COMP_sense_500 --scaling-factor=\$(cat $SCALE4a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT4_sense_500" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED_COMP_anti_500 --scaling-factor=\$(cat $SCALE4a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT4_anti_500" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT5_SCALED_COMP_sense_500 --scaling-factor=\$(cat $SCALE5a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT5_sense_500" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT5_SCALED_COMP_anti_500 --scaling-factor=\$(cat $SCALE5a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT5_anti_500" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT6_SCALED_COMP_500 --scaling-factor=\$(cat $SCALE6a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT6f_500" >> $sampleID
echo "#make scaled OUT file for each strand" >> $sampleID
echo "perl $JOB $CDT1_SCALED_COMP_sense_500 $SCALED_OUT1_sense_500" >> $sampleID
echo "perl $JOB $CDT1_SCALED_COMP_anti_500 $SCALED_OUT1_anti_500" >> $sampleID
echo "perl $JOB $CDT2_SCALED_COMP_sense_500 $SCALED_OUT2_sense_500" >> $sampleID
echo "perl $JOB $CDT2_SCALED_COMP_anti_500 $SCALED_OUT2_anti_500" >> $sampleID
echo "perl $JOB $CDT3_SCALED_COMP_sense_500 $SCALED_OUT3_sense_500" >> $sampleID
echo "perl $JOB $CDT3_SCALED_COMP_anti_500 $SCALED_OUT3_anti_500" >> $sampleID
echo "perl $JOB $CDT4_SCALED_COMP_sense_500 $SCALED_OUT4_sense_500" >> $sampleID
echo "perl $JOB $CDT4_SCALED_COMP_anti_500 $SCALED_OUT4_anti_500" >> $sampleID
echo "perl $JOB $CDT5_SCALED_COMP_sense_500 $SCALED_OUT5_sense_500" >> $sampleID
echo "perl $JOB $CDT5_SCALED_COMP_anti_500 $SCALED_OUT5_anti_500" >> $sampleID
echo "perl $JOB $CDT6_SCALED_COMP_500 $SCALED_OUT6_500" >> $sampleID
echo "#concatenate OUT fles and take lines 1,2,4 to final composite files for each library" >> $sampleID
echo "cat $SCALED_OUT1_sense_500 $SCALED_OUT1_anti_500 | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT1_final_500" >> $sampleID
echo "cat $SCALED_OUT2_sense_500 $SCALED_OUT2_anti_500 | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT2_final_500" >> $sampleID
echo "cat $SCALED_OUT3_sense_500 $SCALED_OUT3_anti_500 | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT3_final_500" >> $sampleID
echo "cat $SCALED_OUT4_sense_500 $SCALED_OUT4_anti_500 | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT4_final_500" >> $sampleID
echo "cat $SCALED_OUT5_sense_500 $SCALED_OUT5_anti_500 | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT5_final_500" >> $sampleID
echo "#remove intermediate files" >> $sampleID
echo "rm $OUT1a_500" >> $sampleID
echo "rm $CDT1_sense_gz_500" >> $sampleID
echo "rm $CDT1_anti_gz_500" >> $sampleID
echo "rm $CDT1_sense_500" >> $sampleID
echo "rm $CDT1_anti_500" >> $sampleID
echo "rm $SCALE1a_OUT" >> $sampleID
echo "rm $CDT1_SCALED_COMP_sense_500" >> $sampleID
echo "rm $CDT1_SCALED_COMP_anti_500" >> $sampleID
echo "rm $SCALED_OUT1_sense_500" >> $sampleID
echo "rm $SCALED_OUT1_anti_500" >> $sampleID
echo "rm $OUT2a_500" >> $sampleID
echo "rm $CDT2_sense_gz_500" >> $sampleID
echo "rm $CDT2_anti_gz_500" >> $sampleID
echo "rm $CDT2_sense_500" >> $sampleID
echo "rm $CDT2_anti_500" >> $sampleID
echo "rm $SCALE2a_OUT" >> $sampleID
echo "rm $CDT2_SCALED_COMP_sense_500" >> $sampleID
echo "rm $CDT2_SCALED_COMP_anti_500" >> $sampleID
echo "rm $SCALED_OUT2_sense_500" >> $sampleID
echo "rm $SCALED_OUT2_anti_500" >> $sampleID
echo "rm $OUT3a_500" >> $sampleID
echo "rm $CDT3_sense_gz_500" >> $sampleID
echo "rm $CDT3_anti_gz_500" >> $sampleID
echo "rm $CDT3_sense_500" >> $sampleID
echo "rm $CDT3_anti_500" >> $sampleID
echo "rm $SCALE3a_OUT" >> $sampleID
echo "rm $CDT3_SCALED_COMP_sense_500" >> $sampleID
echo "rm $CDT3_SCALED_COMP_anti_500" >> $sampleID
echo "rm $SCALED_OUT3_sense_500" >> $sampleID
echo "rm $SCALED_OUT3_anti_500" >> $sampleID
echo "rm $OUT4a_500" >> $sampleID
echo "rm $CDT4_sense_gz_500" >> $sampleID
echo "rm $CDT4_anti_gz_500" >> $sampleID
echo "rm $CDT4_sense_500" >> $sampleID
echo "rm $CDT4_anti_500" >> $sampleID
echo "rm $SCALE4a_OUT" >> $sampleID
echo "rm $CDT4_SCALED_COMP_sense_500" >> $sampleID
echo "rm $CDT4_SCALED_COMP_anti_500" >> $sampleID
echo "rm $SCALED_OUT4_sense_500" >> $sampleID
echo "rm $SCALED_OUT4_anti_500" >> $sampleID
echo "rm $OUT5a_500" >> $sampleID
echo "rm $CDT5_sense_gz_500" >> $sampleID
echo "rm $CDT5_anti_gz_500" >> $sampleID
echo "rm $CDT5_sense_500" >> $sampleID
echo "rm $CDT5_anti_500" >> $sampleID
echo "rm $SCALE5a_OUT" >> $sampleID
echo "rm $CDT5_SCALED_COMP_sense_500" >> $sampleID
echo "rm $CDT5_SCALED_COMP_anti_500" >> $sampleID
echo "rm $SCALED_OUT5_sense_500" >> $sampleID
echo "rm $SCALED_OUT5_anti_500" >> $sampleID
echo "rm $OUT6a_500" >> $sampleID
echo "rm $CDT6e_500" >> $sampleID
echo "rm $CDT6f_500" >> $sampleID
echo "rm $SCALE6a_OUT" >> $sampleID
echo "rm $CDT6_SCALED_COMP_500" >> $sampleID
echo "# finish script" >> $sampleID
