SCRIPTMANAGER=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
GENETRACK=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/220803_Nuc-Genetrack/job/genetrack_v2.py

INPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230719_BI_Genetrack/scIDX-Nuc
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230719_BI_Genetrack/Genetrack_peaks

#mkdir -p $OUTPUT

JOBSTATS="#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=24gb
#PBS -l walltime=4:00:00
#PBS -A open
"

for file in $INPUT/*0-54*_NUC.tab; do
        #Get basename for PBS
        fileID=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
        echo $fileID

        sampleID=$fileID\.pbs
        rm -f $sampleID
        echo "$JOBSTATS" >> $sampleID
        echo "module load anaconda3/2020.07" >> $sampleID
        echo "python2 $GENETRACK -s 10 -e 20 -F 0 $file" >> $sampleID
done

for file in $INPUT/*55-91*_NUC.tab; do
        #Get basename for PBS
        fileID=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
        echo $fileID

        sampleID=$fileID\.pbs
        rm -f $sampleID
        echo "$JOBSTATS" >> $sampleID
        echo "module load anaconda3/2020.07" >> $sampleID
        echo "python2 $GENETRACK -s 20 -e 40 -F 5 $file" >> $sampleID
done

for file in $INPUT/*92-127*_NUC.tab; do
        #Get basename for PBS
        fileID=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
        echo $fileID

        sampleID=$fileID\.pbs
        rm -f $sampleID
        echo "$JOBSTATS" >> $sampleID
        echo "module load anaconda3/2020.07" >> $sampleID
        echo "python2 $GENETRACK -s 30 -e 60 -F 6 $file" >> $sampleID
done

for file in $INPUT/*128-164*_NUC.tab; do
	#Get basename for PBS
	fileID=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
	echo $fileID

        sampleID=$fileID\.pbs
        rm -f $sampleID
        echo "$JOBSTATS" >> $sampleID
        echo "module load anaconda3/2020.07" >> $sampleID
        echo "python2 $GENETRACK -s 40 -e 80 -F 5 $file" >> $sampleID
done

for file in $INPUT/*165-Inf*_NUC.tab; do
        #Get basename for PBS
        fileID=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
        echo $fileID

        sampleID=$fileID\.pbs
        rm -f $sampleID
        echo "$JOBSTATS" >> $sampleID
        echo "module load anaconda3/2020.07" >> $sampleID
        echo "python2 $GENETRACK -s 50 -e 100 -F 3 $file" >> $sampleID
done


for file in *.pbs; do qsub $file; done
