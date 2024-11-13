import pandas as pd
import sys

# Get the input file path from command line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Load the tab-delimited file
data = pd.read_csv(input_file, delimiter='\t', header=None)

# Extract the relevant columns (201-326), using 1-based indexing
start_col = 201 - 1  # Convert to 0-based indexing
end_col = 326        # End column is inclusive in 1-based, exclusive in 0-based

# Extract data from row 2 (index 1) for the specified columns
row_2_data = data.iloc[1, start_col:end_col].astype(float)

# Find the maximum value and its column position within row 2
max_value = row_2_data.max()
max_column = row_2_data.idxmax()  # Column index of the maximum value in row 2

# Extract the corresponding value from row 1 (header row)
max_position = data.iloc[0, max_column]

# Adjust the column index by -500 (convert to 1-based indexing)
adjusted_max_position = (max_column) - 500

# Prepare the results for output
results = {
    'Filename': input_file,
    'Max Value': max_value,
    'Max Position (Header Row)': max_position,
    'Adjusted Column Index': adjusted_max_position
}

# Create a DataFrame for the results
results_df = pd.DataFrame([results])

# Save the results to the specified output file
results_df.to_csv(output_file, sep='\t', index=False)

print(f'Results saved to {output_file}')
