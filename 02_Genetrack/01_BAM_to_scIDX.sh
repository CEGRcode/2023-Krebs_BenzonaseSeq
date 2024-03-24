SCRIPTMANAGER=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
INPUT_BAM=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE/K562_benzonase-seq_master.bam
INPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230719_BI_Genetrack/scIDX

#picard.jar
PICARD=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/picard.jar

mkdir -p $OUTPUT

JOBSTATS="#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=24gb
#PBS -l walltime=6:00:00
#PBS -A open
cd $OUTPUT

source ~/.bashrc #configures shell to use conda activate
conda activate bioinfo"

	#Get basename for PBS
	fileID=$(echo $INPUT_BAM | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
	echo $fileID
	BAI=$(echo $INPUT_BAM | rev | cut -d"/" -f1 | rev | awk '{print $1".bai"}')
	# output file1 for <55bp insert size
	OUTNAME1=$(echo $INPUT_BAM | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_MIDPOINT_0-54.tab"}')
	# output file2 for 55-91 insert size
	OUTNAME2=$(echo $INPUT_BAM | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_MIDPOINT_55-91.tab"}')
	# output file3 for 92-127 insert size
	OUTNAME3=$(echo $INPUT_BAM | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_MIDPOINT_92-127.tab"}')
	# output file4 for 128-164 insert size
	OUTNAME4=$(echo $INPUT_BAM | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_MIDPOINT_128-164.tab"}')
	# output file5 for >164 insert size
	OUTNAME5=$(echo $INPUT_BAM | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_MIDPOINT_165-Inf.tab"}')

        sampleID=$fileID\_split.pbs
        rm -f $sampleID
        echo "$JOBSTATS" >> $sampleID
        echo "#index input bam file" >> $sampleID
	echo "java -jar $PICARD BuildBamIndex --INPUT $INPUT_BAM --OUTPUT $INPUT/$BAI" >> $sampleID
	echo "#convert BAM to scIDX via command line version of scriptmanager with insert size <55bp" >> $sampleID
        echo "java -jar $SCRIPTMANAGER bam-format-converter bam-to-scidx -m -p -x=54 -o $OUTPUT/$OUTNAME1 $INPUT_BAM" >> $sampleID
        echo "#convert BAM to scIDX via command line version of scriptmanager with insert sizes 55-91" >> $sampleID
        echo "java -jar $SCRIPTMANAGER bam-format-converter bam-to-scidx -m -p -n=55 -x=91 -o $OUTPUT/$OUTNAME2 $INPUT_BAM" >> $sampleID
        echo "#convert BAM to scIDX via command line version of scriptmanager with insert sizes 92-127" >> $sampleID
        echo "java -jar $SCRIPTMANAGER bam-format-converter bam-to-scidx -m -p -n=92 -x=127 -o $OUTPUT/$OUTNAME3 $INPUT_BAM" >> $sampleID
        echo "#convert BAM to scIDX via command line version of scriptmanager with insert sizes 128-164" >> $sampleID
        echo "java -jar $SCRIPTMANAGER bam-format-converter bam-to-scidx -m -p -n=128 -x=164 -o $OUTPUT/$OUTNAME4 $INPUT_BAM" >> $sampleID
        echo "#convert BAM to scIDX via command line version of scriptmanager with insert sizes >164" >> $sampleID
        echo "java -jar $SCRIPTMANAGER bam-format-converter bam-to-scidx -m -p -n=165 -o $OUTPUT/$OUTNAME5 $INPUT_BAM" >> $sampleID

