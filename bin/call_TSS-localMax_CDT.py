import sys
import argparse
import numpy as np
from scipy.stats import poisson

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Calls max signal within (sense) CDT pileup and returns BED-formatted coordinates for the max signal positions')
	parser.add_argument('-i','--input', metavar='cdt_fn', required=True, help='A CDT-formatted file of the data\'s tag pileup')
	parser.add_argument('-r','--reference', metavar='bed_fn', required=True, help='A BED-formatted file that the pileup was generated off of (assumes unique ID col-4)')
	parser.add_argument('-o','--output', metavar='outfile', required=True, help='a BED file of max signal position within reference BED intervals')

	parser.add_argument('--tag-count', action='store_true', help='Replace score column with max tag count if flagged')
	parser.add_argument('-s','--smoothing', metavar='num_bp', default=1, type=int, help='The smoothing window size for determining max position, must be odd. (default=1)')
	parser.add_argument('-t','--threshold', metavar='tag_count', default=5, type=float, help='The threshold data signal value minimum for calling a max. (default=5)')
	#parser.add_argument('-t','--title', metavar='figure_title', dest='title', required=True, help='')

	args = parser.parse_args()

	# validate smoothing input
	if(args.smoothing % 2 != 1):
		raise Exception("Invalid input: smoothing value must be odd ("+args.smoothing+")\n")

	return(args)


def validateBED(bed_fn):
	list_BED = []
	reader = open(bed_fn,'r')
	for line in reader:
		tokens = line.strip().split('\t')
		id = tokens[3]
		list_BED.append((id,tokens))
	reader.close()
	return(list_BED)

# tagPileup_CDT_fn	string		filename of sorted scidx file read counts
# ref_features_bed	string		filename of bed file containing coordinates of most 3' side of window (exclusive)
# 									->written to take BED file of ORFs, where most 5' coordinate is the A in ATG
# window				int			specifies window size as the number of basepairs to check on each side of the given annotations
# threshold			numeric		the threshold value for defining a TSS (inclusive)
# header				boolean		set to true if you wish for the output file to have column names
if __name__ == "__main__":
	'''Calls max signal within CDT pileup and returns BED-formatted coordinates for the max signal positions.'''

	# Load reference BED coordinates
	args = getParams()
	coords = validateBED(args.reference)

	c_index = 0

	writer = open(args.output,'w')
	# Parse pileup
	cdt_reader = open(args.input,'r')
	for line in cdt_reader:
		# Skip first line
		if(line.find("YORF")==0):
			continue
		tokens = line.strip().split('\t')
		# Check CDT id matches BED id
		if(tokens[0]!=coords[c_index][0]):
			c_index += 1
			raise Exception("Unmatched ID error: "+tokens[0]+"!="+str(coords[c_index][0])+" @line"+str(c_index)+"\n")
		# Create smoothed value array
		raw_values = [float(i) for i in tokens[2:]]
		smoothed_values = [sum(raw_values[i:i+args.smoothing]) for i in range(len(raw_values)-int(args.smoothing/2))]

		if(len(smoothed_values)==0):
			print(tokens[:2])
			c_index += 1
			continue
		# Get index of maximum
		max_index = smoothed_values.index(max(smoothed_values))
		max_val = smoothed_values[max_index]

		# Filter by threshold
		if (args.threshold > max_val):
			c_index += 1
			continue

		# Create new coordinate adjusted to max signal
		new_coord = coords[c_index][1]
		# Only modify coordinate if non-zero tally
		if (max(smoothed_values)>0):
			if(new_coord[5]=='-'):
				new_coord[1] = str(int(new_coord[2])-1 - max_index)
			else:
				new_coord[1] = str(int(new_coord[1]) + max_index)
			new_coord[2] = new_coord[1]

		# Update score column with tag count if flagged
		if (args.tag_count):
			new_coord[4] = str(max_val)

		# Write new coordinate to BED file output
		writer.write("\t".join(new_coord) + "\n")

		# Increment to next BED entry
		c_index += 1
	cdt_reader.close()
	writer.close()
