import pandas as pd
import argparse

def count_rows_based_on_column(input_file, output_file):
    # Read the tab-delimited file into a DataFrame
    df = pd.read_csv(input_file, sep='\t')
    
    # Column 2 (1-based) in pandas is index 1 (0-based)
    column_index = 1
    
    # Count the number of rows where the value in column 2 is less than 501
    count = (df.iloc[:, column_index] < 501).sum()
    
    # Write the count to the output file with the updated header message
    with open(output_file, 'w') as f:
        f.write(f"Unique, significant peaks 5' to motif: {count}\n")

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Count rows where values in column 2 are less than 501.")
    parser.add_argument('input_file', type=str, help='Path to the input tab-delimited file')
    parser.add_argument('output_file', type=str, help='Path to the output file')

    args = parser.parse_args()

    # Call the function with command-line arguments
    count_rows_based_on_column(args.input_file, args.output_file)