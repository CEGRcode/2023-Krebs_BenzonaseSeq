SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
INPUT_BAM=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/01_scIDX_output/K562_MNase_read1.tab
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_MNase_NUC_calls/test_run/01_scIDX_output/MNase

#picard.jar
PICARD=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/picard.jar

mkdir -p $OUTPUT

JOBSTATS="#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --partition=open

module load anaconda #configures shell to use conda activate
conda activate bioinfo"

	#Get basename for Slurm file
	fileID=$(echo $INPUT_BAM | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
	echo $fileID
	# output file1
	OUTNAME1=$(echo $INPUT_BAM | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_read1.tab"}')

        sampleID=$fileID\_shift.slurm
        rm -f $sampleID
        echo "$JOBSTATS" >> $sampleID
	echo "#set output directory" >> $sampleID
	echo "cd $OUTPUT" >> $sampleID
	echo "#remove first 2 rows, if column 3 (forward) has a value greater than or equal to 1 AND column 4 is equal to 0 THEN take those columns, add 82 bp to index" >> $sampleID
	echo "cat $INPUT_BAM | awk 'NR>2' | awk '{if (\$3>=\"1\" && \$4==\"0\") print \$1\"\t\"(\$2+82)\"\t\"\$3\"\t\"\$4\"\t\"\$5}' > file1.tab" >> $sampleID
	echo "#remove first 2 rows, if column 4 (reverse) has a value greater than or equal to 1 AND column 3 equal to 0 THEN take those columns, subtract 82 bp from index" >> $sampleID
	echo "cat $INPUT_BAM | awk 'NR>2' | awk '{if (\$3==\"0\" && \$4>=\"1\") print \$1\"\t\"(\$2-82)\"\t\"\$3\"\t\"\$4\"\t\"\$5}' > file2.tab" >> $sampleID
	echo "#collect rows that have values in forward and reverse; duplicate this by during this twice -> when adding 82 to forward, make reverse column =0; total values = forward column" >> $sampleID
	echo "cat $INPUT_BAM | awk 'NR>2' | awk '{if (\$3>=\"1\" && \$4>=\"1\") print \$1\"\t\"(\$2+82)\"\t\"\$3\"\t\"\"0\"\"\t\"\$3}' > file3.tab" >> $sampleID
	echo "#collect rows that have values in forward and reverse; duplicate this by during this twice -> when subtracting 82 to reverse, make forward column =0; total values = forward column" >> $sampleID
	echo "cat $INPUT_BAM | awk 'NR>2' | awk '{if (\$3>=\"1\" && \$4>=\"1\") print \$1\"\t\"(\$2-82)\"\t\"\"0\"\"\t\"\$4\"\t\"\$4}' > file4.tab" >> $sampleID
	echo "#make header - add header (2 lines) back" >> $sampleID
	echo "cat $INPUT_BAM | head -2 | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5}' > file_header.tab" >> $sampleID
	echo "#combine, remove rows that are now off the edge of the chr. (\$2>=\"1\" on index), and sort file" >> $sampleID
	echo "cat file1.tab file2.tab file3.tab file4.tab | awk '{if (\$2>=\"1\") print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5}' | sort -k1,1 -k2,2n > shifted_scIDX_MNase.tab" >> $sampleID
	echo "cat file_header.tab shifted_scIDX_MNase.tab > final_shifted_scIDX_MNase.tab" >> $sampleID

