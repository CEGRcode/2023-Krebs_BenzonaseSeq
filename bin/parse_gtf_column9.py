import sys, argparse


def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description = """
============
Parse a two column file (col1=UID, col2=GTF column 9 string) to extract features into a tsv file.
============
""", formatter_class = argparse.RawTextHelpFormatter)

	parser.add_argument('-i','--input', metavar='bed_fn', required=True, help='the two-column input file')
	parser.add_argument('-o','--output', metavar='tsv_fn', required=True, help='the output TSV file of field values')

	args = parser.parse_args()
	return(args)

if __name__ == '__main__':
	'''Main program which takes in input parameters'''
	args = getParams()

	print("Input file: ", args.input)
	print("Output file: ", args.output)

	# Initialize file reader
	reader = open(args.input, 'r')

	id2field = {}
	field_names = set()
	# Parse two-col
	for line in reader:
		tokens = line.strip().split('\t')
		uid = tokens[0]
		if (uid in id2field.keys()):
			raise Exception("First column contains non-unique value (%s). First tab-delimited column needs to contain unique ids" % uid)
		id2field.update({uid : {}})
		for field_str in tokens[1].strip().split(";"):
			# Clean field string
			field_str = field_str.strip()
			# Split field string by space
			split_idx = field_str.find(' ')
			fn = field_str[:split_idx]
			fv = field_str[split_idx+1:]
			# Update parsed field name
			field_names.add(fn)
			# Update parsed field value (if fn is unique or not)
			if (fn in id2field[uid].keys()):
				id2field[uid][fn] += ";" + fv
			else:
				id2field[uid].update({fn:fv})

	# Close file
	reader.close()

	# Cast field names as list
	field_names = sorted(list(field_names))

	# Initialize file writer
	writer = open(args.output, 'w')

	# Write header
	writer.write('UID\t%s\n' % ('\t'.join(field_names)))
	for uid in id2field.keys():
		fields = [""] * len(field_names)
		for i,fn in enumerate(field_names):
			fields[i] = id2field[uid].get(fn,"")
		writer.write('%s\t%s\n' % (uid, '\t'.join(fields)))

	# Close file
	writer.close()
