import sys, argparse

NT = ["A","T","C","G"]

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description = """
============
Get 0/1 matrix (CDT format) of dinucleotides
============
""", formatter_class = argparse.RawTextHelpFormatter)

	parser.add_argument('-i','--input', metavar='fasta_fn', required=True, help='the FAStA file to analyze')
	parser.add_argument('-s','--seq', metavar='dinucleotides_str', required=True, help='the \"-\" delimited set of dinucleotides to check for')
	parser.add_argument('-o','--output', metavar='tsv_fn', required=True, help='the output CDT formatted 0/1 matrix of dinucleotide matches')

	args = parser.parse_args()
	return(args)

def dint_str_to_list(dint_str):
	'''Check dinucleotide string format and return list'''
	# Parse dinucleotide string and force uppercase chars
	dints = [dnt.upper() for dnt in dint_str.split("-")]
	# Check dinucleotides for length and char usage (ATCG)
	for d in dints:
		if (len(d)!=2):
			raise Exception("Dinucleotides must be *two* nucleotides long: %s" % d)
		if (d[0] not in NT and d[1] not in NT):
			raise Exception("Dinucleotides must be made up of \"A\", \"T\", \"C\", or \"G\":" % d)
	return(dints)

if __name__ == '__main__':
	'''Main program which takes in input parameters'''
	args = getParams()

	print("FASTA file: ", args.input)
	print("Dinucleotides to check: ", args.seq)
	print("Output file: ", args.output)

	# Parse dinucleotides
	dints = dint_str_to_list(args.seq)

	# Initialize file reader
	reader = open(args.input, 'r')

	LINES = []
	NAME = None
	SEQ = ""
	# Parse FASTA
	for line in reader:
		# FASTA header
		if (line.find(">")==0):
			# If not the first entry...
			if (NAME != None):
				SEQ = SEQ.upper()
				# Substring dinucleotide at current position and set 0/1 as appropriate
				list = [ "1" if SEQ[i:i+2] in dints else "0" for i in range(len(SEQ)-1)]
				# Write to CDT (pad 0 at right)
				LINES.append("%s\t%s\t%s\t0\n" % (NAME,NAME,"\t".join(list)))
			# Reset vars
			NAME = line.strip().split(" ")[0][1:]
			SEQ = ""
			continue
		# FASTA sequence lines
		SEQ += line.strip()
	# Close files
	reader.close()

	# Initialize writer
	writer = open(args.output, 'w')
	# Write header
	writer.write("NAME\tYORF\t%s\n" % ("\t".join( [str(i) for i in range(len(SEQ))])))
	# Write lines
	[ writer.write(line) for line in LINES ]
	# Close writer
	writer.close()
