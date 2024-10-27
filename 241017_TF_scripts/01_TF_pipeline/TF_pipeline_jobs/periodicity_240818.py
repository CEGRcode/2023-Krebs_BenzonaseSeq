import pandas as pd
import sys

def find_local_maxima(input_file, output_file):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(input_file, sep='\t', header=None)
    
    # Extract the columns of interest
    values_col3 = df[2]  # Column 3 (index 2)
    values_col2 = df[1]  # Column 2 (index 1)
    
    # Initialize a list to store output rows
    output_rows = []
    
    # Iterate through the values to find local maxima
    for i in range(1, len(values_col3) - 1):
        if values_col3[i] > values_col3[i - 1] and values_col3[i] > values_col3[i + 1]:
            # Append a row with "periodicity" in column 1 and the value from column 2 in column 2
            output_rows.append(["periodicity", values_col2[i]])
    
    # Convert the list of output rows to a DataFrame
    output_df = pd.DataFrame(output_rows)
    
    # Write the DataFrame to an output file in TSV format
    output_df.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: py JOB.py input.tab output.tab")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    find_local_maxima(input_file, output_file)