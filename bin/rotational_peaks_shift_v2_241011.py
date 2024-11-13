import pandas as pd
import sys

def shift_bins(input1_path, input2_path, output_path):
    try:
        # Load input files
        input1 = pd.read_csv(input1_path, sep='\t', header=None, names=['Start_Column', 'End_Column'])
        input2 = pd.read_csv(input2_path, sep='\t', header=None)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Extract the mode from input2 (assuming mode is in the 2nd column of the first row)
    mode_value = input2.iloc[0, 1]
    
    # Calculate the shift amount as 'mode-5'
    shift_amount = mode_value - 5
    
    # Create a copy of the input1 DataFrame for shifting
    shifted_nucleosome_df = input1.copy()
    
    # Apply the shift
    shifted_nucleosome_df[['Start_Column', 'End_Column']] += shift_amount
    
    # Save the shifted bins to output file
    shifted_nucleosome_df.to_csv(output_path, sep='\t', index=False, header=False)
    print(f"Successfully saved shifted bins to {output_path}")

# Example usage
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script_name.py input1.tab input2.tab output.tab")
        sys.exit(1)

    script_name, input1_path, input2_path, output_path = sys.argv
    shift_bins(input1_path, input2_path, output_path)
