#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 10:00:00
#SBATCH -A open
#SBATCH -o logs/0_Download_SRA.log.out-%a
#SBATCH -e logs/0_Download_SRA.log.err-%a

# Download external dataset (raw FASTQ files)

# OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/fastq_files
OUTPUT=../data/FASTQ

[ -d $OUTPUT ] || mkdir -p $OUTPUT
cd $OUTPUT

set -eux
module load anaconda
conda activate bioinfo

# MNase titrations
[ -d MNase ] || mkdir MNase
fasterq-dump -O MNase --split-files SRR3211679
fasterq-dump -O MNase --split-files SRR3211680
fasterq-dump -O MNase --split-files SRR3211681
fasterq-dump -O MNase --split-files SRR3211682

# DNase-seq
[ -d DNase ] || mkdir DNase
fasterq-dump -O DNase --split-files SRR16815400

# H3K4me3 MNase ChIP-seq
[ -d MNaseChIP ] || mkdir MNaseChIP
fasterq-dump -O MNaseChIP --split-files SRR6010175
fasterq-dump -O MNaseChIP --split-files SRR6010177
fasterq-dump -O MNaseChIP --split-files SRR6010180
fasterq-dump -O MNaseChIP --split-files SRR7441419
fasterq-dump -O MNaseChIP --split-files SRR7441420
