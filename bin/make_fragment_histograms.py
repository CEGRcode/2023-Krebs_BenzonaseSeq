#!/bin/python
import os, sys, argparse
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator)
import seaborn as sns
sns.set()

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Build insert-size histograms from ScriptManager bam-statistics pe-stats output.')

	parser.add_argument('-i','--input', metavar='input_fn', required=True, help='the "*InsertHistogram.out" file listing insert size counts')
	parser.add_argument('-o','--output', metavar='png_fn', required=True, help='the output figure image (use .svg or .png)')

	parser.add_argument('--reflines', type=int, nargs="+", default=[], help="A list of integers for vertical line marks (optional).")
	parser.add_argument('--ymax', default=None, type=int, help='set a maximum to the y-axis range')

	args = parser.parse_args()
	return(args)

'''
# 2024-10-30 16:05:34.442
# BNase-seq_50U-3min_1_hg38.bam
# Chromosome_ID	Chromosome_Size	Aligned_Reads	Unaligned_Reads
# chr1	248956422	1113330.0	0.0
# chr10	133797422	447050.0	0.0
# chr11	135086622	449232.0	0.0
# ...
# chrX	156040895	288140.0	0.0
# chrY	57227415	8466.0	0.0
# chrY_KI270740v1_random	37240	0.0	0.0
# Total Genome Size: 3.209286105E9	Total Aligned Tags: 1.0103658E7
# bwa	# 0.7.17-r1188
# bwa mem -t 8 -v 1 /storage/group/bfp2/default/00_pughlab/tool_data/hg38/bwa_mem_index/hg38/hg38.fa /storage/group/bfp2/default/00_pughlab/pulsar/files/staging/801231/inputs/459-28385_S59_R1_001.fastq.gz /storage/group/bfp2/default/00_pughlab/pulsar/files/staging/801231/inputs/459-28385_S59_R2_001.fastq.gz
# MarkDuplicates	# 2.18.2-SNAPSHOT
# MarkDuplicates TAGGING_POLICY=All INPUT=[Map_with_BWA-MEM_on_data_4_and_data_3__mapped_reads_in_BAM_format_] OUTPUT=/storage/group/bfp2/default/00_pughlab/pulsar/files/staging/801235/outputs/dataset_69c47ed8-8812-49cc-aecd-af17264bcf22.dat METRICS_FILE=/storage/group/bfp2/default/00_pughlab/pulsar/files/staging/801235/outputs/dataset_3d2f72e7-2fa0-4837-9cab-ed60bab3bb94.dat REMOVE_DUPLICATES=false ASSUME_SORTED=true DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*. OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 VERBOSITY=ERROR QUIET=true VALIDATION_STRINGENCY=LENIENT    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 TAG_DUPLICATE_SET_MEMBERS=false REMOVE_SEQUENCING_DUPLICATES=false CLEAR_DT=true ADD_PG_TAG_TO_READS=true PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
# Average Insert Size: 131.72716328945938
# Median Insert Size: 130.0
# Mode Insert Size: 91.0
# Std deviation of Insert Size: 41.40447742645101
# Number of ReadPairs: 4875330.0
# Histogram
# Size (bp)	Frequency
0	0.0
1	0.0
2	17.0
3	25.0
4	25.0
...
'''

# Main program which takes in input parameters
if __name__ == '__main__':
	'''Collect metadata and EpitopeID results to get detection stats on the YEP data'''

	args = getParams()

	# Populate dataframe with tab file data
	filedata = pd.read_table(args.input, sep='\t', skiprows=[0,1], names=['InsertSize','ReadCount'], comment='#')

	# Initialize plot
	sns.set_style("ticks")
	fig, ax = plt.subplots()

	# Plot the filedata
	ax = sns.barplot(data=filedata, x='InsertSize', y='ReadCount', color='black', linewidth = 0)

	# Format axes and tickmarks
	ax.xaxis.grid(False)
	ax.set_xlim(0, 340)
	plt.xticks([0, 73, 146, 219, 292])
	ax.yaxis.grid(True)

	# Set y-axis limit if user provided
	if (args.ymax != None):
		ax.set_ylim(0, args.ymax)

	# Add vertical reference lines
	for x_pos in args.reflines:
		plt.axvline(x=x_pos, color='gray', linestyle='--', linewidth=0.8)

	# Title
	ax.set_title(os.path.basename(args.input))

	# Fill line
	plt.fill_between(filedata.InsertSize.values, filedata.ReadCount.values, color='black')
	plt.xticks()

	# Save figure as image
	sfig = ax.get_figure()
	sfig.savefig(args.output)
