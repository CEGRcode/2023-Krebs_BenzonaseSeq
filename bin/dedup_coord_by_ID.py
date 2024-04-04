import sys, argparse

NT = ["A","T","C","G"]

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description = """
============
Deduplicate coordinates in a BED file by the ID column (4th). Keep the first instance.
============
""", formatter_class = argparse.RawTextHelpFormatter)

	parser.add_argument('-i','--input', metavar='bed_fn', required=True, help='the BED file to deduplicate')
	parser.add_argument('-o','--output', metavar='tsv_fn', required=True, help='the output BED file of deduplicated coordinates')

	# parser.add_argument('--first', required=False, help='keep the first instance (default)')
	# parser.add_argument('--greatest-score', required=False, help='keep the instance with the greater score (first if equivalent)')
	# parser.add_argument('--least-score', required=False, help='keep the instance with the least score (first if equivalent)')

	args = parser.parse_args()
	return(args)

if __name__ == '__main__':
	'''Main program which takes in input parameters'''
	args = getParams()

	print("BED file: ", args.input)
	print("Output file: ", args.output)

	# Initialize file reader and writer
	reader = open(args.input, 'r')
	writer = open(args.output, 'w')

	id_list = []
	# Parse BED
	for line in reader:
		tokens = line.strip().split('\t')
		if (line.find("#")==0):
			writer.write(line)
			continue
		if (tokens[3] not in id_list):
			writer.write(line)
			id_list.append(tokens[3])
	# Close files
	reader.close()
	writer.close()
