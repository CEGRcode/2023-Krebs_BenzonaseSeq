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

	args = parser.parse_args()
	return(args)

'''
2023-03-23 18:16:24.561
InsertSize (bp) K562_CTCF_BX_hg19_rep1.bam
0       0.0
1       0.0
2       179.0
3       162.0
4       186.0
5       160.0
6       175.0
7       135.0
8       153.0
9       180.0
10      168.0
11      180.0
12      183.0
13      141.0
'''

# Main program which takes in input parameters
if __name__ == '__main__':
	'''Collect metadata and EpitopeID results to get detection stats on the YEP data'''

	args = getParams()

	# Populate dataframe with tab file data
	filedata = pd.read_table(args.input, sep='\t', skiprows=[0,1], names=['InsertSize','ReadCount'])

	# Initialize plot
	sns.set_style("ticks")
	fig, ax = plt.subplots()
	
	# Plot the filedata
	ax = sns.barplot(data=filedata, x='InsertSize', y='ReadCount', color='black', linewidth = 0)
	
	# Format axes and tickmarks
	ax.xaxis.grid(False)
	ax.xaxis.set_major_locator(MultipleLocator(20))
	ax.xaxis.set_minor_locator(MultipleLocator(10))
	ax.set_xlim(0, 340)
	ax.yaxis.grid(True)
	
	# Title
	ax.set_title(os.path.basename(args.input))

	# Fill line
	plt.fill_between(filedata.InsertSize.values, filedata.ReadCount.values, color='black')
	plt.xticks()

	# Save figure as image
	sfig = ax.get_figure()
	sfig.savefig(args.output)
