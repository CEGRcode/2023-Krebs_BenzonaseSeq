import pandas as pd
import argparse

def process_file(input_file, output_file):
    # Load the data into a DataFrame
    df = pd.read_csv(input_file, delimiter='\t')
    
    # Find the maximum value in column 2 (excluding the header)
    max_value = df.iloc[1:, 1].max()
    
    # Add a new column (column 5) with the maximum value of column 2 divided by the value of column 2 in that row
    df['Column 5'] = (max_value / df.iloc[:, 1]).round(2)
    
    # Save the updated DataFrame to a new file (or overwrite the existing one)
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    # Setup command-line argument parsing
    parser = argparse.ArgumentParser(description='Process a tab-delimited file.')
    parser.add_argument('input_file', type=str, help='Input tab-delimited file')
    parser.add_argument('output_file', type=str, help='Output tab-delimited file')
    
    args = parser.parse_args()
    
    # Call the processing function with command-line arguments
    process_file(args.input_file, args.output_file)
