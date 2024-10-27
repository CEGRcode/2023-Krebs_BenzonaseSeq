import pandas as pd
import sys

def process_files(input1_path, input2_path, output_path):
    # Read the input files
    input1 = pd.read_csv(input1_path, sep='\t', header=None)
    input2 = pd.read_csv(input2_path, sep='\t', header=None)

    # Check for the negative value after "Average Value:"
    average_value_row = input1[input1[0].str.contains("Average Value:")]
    if not average_value_row.empty:
        # Extract the value and convert it to float
        average_value_str = average_value_row.iloc[0, 0].split(":")[1].strip()
        average_value = float(average_value_str)

        if average_value < 0:
            # Change signs for the second row of input2 (except the first column)
            input2.iloc[1, 1:] = -input2.iloc[1, 1:]

    # Add 10 to all values in the second row of input2 (except the first column)
    input2.iloc[1, 1:] += 10

    # Write the modified or original data to output
    input2.to_csv(output_path, sep='\t', index=False, header=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python JOB.py input1.tab input2.tab output.tab")
        sys.exit(1)

    input1_file = sys.argv[1]
    input2_file = sys.argv[2]
    output_file = sys.argv[3]

    process_files(input1_file, input2_file, output_file)
