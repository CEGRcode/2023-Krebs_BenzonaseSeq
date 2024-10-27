import sys

def process_file(input_file, output_file, divisor):
    # Open the input file for reading
    with open(input_file, 'r') as infile:
        # Read lines from the file
        lines = infile.readlines()
        
        # First line (positions), split into columns
        header = lines[0].strip().split('\t')
        
        # Second line (values), split into columns
        values = lines[1].strip().split('\t')
        
        # Name column remains the same (first column)
        name_column = values[0]
        
        # Convert the values (skip the first column) and divide by the divisor
        divided_values = [str(float(value) / divisor) for value in values[1:]]
        
        # Create the output (concatenate name column with divided values)
        divided_row = [name_column] + divided_values
        
    # Open the output file for writing
    with open(output_file, 'w') as outfile:
        # Write the header (positions row)
        outfile.write('\t'.join(header) + '\n')
        # Write the divided values row
        outfile.write('\t'.join(divided_row) + '\n')

if __name__ == "__main__":
    # Command line arguments: input file, output file, divisor
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    divisor = float(sys.argv[3])
    
    # Call the process function
    process_file(input_file, output_file, divisor)

