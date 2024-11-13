import sys

# Substitution mapping based on user input
substitutions = {
    '9': '1',
    '8': '2',
    '7': '3',
    '6': '4',
    '4': '6',
    '3': '7',
    '2': '8',
    '1': '9'
}

# Read the input and output file paths from command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the input file
with open(input_file, 'r') as infile:
    content = infile.read().strip()

# Split by tab and apply the substitutions
columns = content.split('\t')
for i in range(len(columns)):
    columns[i] = substitutions.get(columns[i], columns[i])  # Replace if in the dictionary

# Join the columns back together with a tab and write to the output file
with open(output_file, 'w') as outfile:
    outfile.write('\t'.join(columns))

print(f"Substitutions completed. Output saved to {output_file}.")
