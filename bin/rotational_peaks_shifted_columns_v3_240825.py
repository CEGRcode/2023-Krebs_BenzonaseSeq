import pandas as pd
import sys

def calculate_bin_statistics(input1_file, input2_file, output_file):
    # Read the input files
    input1_df = pd.read_csv(input1_file, sep='\t', header=None)
    input2_df = pd.read_csv(input2_file, sep='\t', header=None)

    # Convert positions and values to numeric types
    positions = pd.to_numeric(input1_df.iloc[0].values, errors='coerce')
    values = pd.to_numeric(input1_df.iloc[1].values, errors='coerce')

    results = []

    # Process each bin in input2_df
    for _, row in input2_df.iterrows():
        # Read start and end values from input2_df, which are 1-based
        start = pd.to_numeric(row[0], errors='coerce')
        end = pd.to_numeric(row[1], errors='coerce')

        # Adjust start and end for 0-based indexing
        start_0_based = start - 1
        end_0_based = end - 1

        # Get indices for the current bin (adjusted for 0-based indexing)
        bin_indices = (positions >= start_0_based) & (positions <= end_0_based)
        
        # Extract bin data
        bin_positions = positions[bin_indices]
        bin_values = values[bin_indices]
        
        if len(bin_values) == 0:
            continue
        
        # Calculate statistics
        max_value = bin_values.max()
        min_value = bin_values.min()
        max_position = bin_positions[bin_values.argmax()] # Convert to 1-based and add 2
        min_position = bin_positions[bin_values.argmin()] # Convert to 1-based and add 2

        # Adjust positions (subtracting 499)
        adjusted_max_position = max_position - 499
        adjusted_min_position = min_position - 499
        
        # Compute range (max - min)
        range_value = max_value - min_value
        
        # Append result, including both adjusted and indexed positions
        results.append([
            start, end, max_value, min_value, range_value, 
            adjusted_max_position, adjusted_min_position,
            max_position, min_position
        ])

    # Create a DataFrame for results and save to output file
    results_df = pd.DataFrame(results, columns=[
        'Start', 'End', 'Max_Value', 'Min_Value', 'Range', 
        'Adjusted_Max_Position', 'Adjusted_Min_Position',
        'Max_Position', 'Min_Position'
    ])
    results_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    # Ensure correct number of arguments are provided
    if len(sys.argv) != 4:
        print("Usage: python JOB.py input1.tab input2.tab output.tab")
        sys.exit(1)
    
    # Read arguments from command line
    input1_file = sys.argv[1]
    input2_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Run the function
    calculate_bin_statistics(input1_file, input2_file, output_file)