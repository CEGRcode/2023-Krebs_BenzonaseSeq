import sys

# List of labels corresponding to each input file
labels = ["AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG", "CA", "CT", "CC", "CG", "GA", "GT", "GC", "GG"]

def get_masked_region(masked_file):
    with open(masked_file, 'r') as file:
        lines = file.readlines()
        # col1 and col2 are in the second row, 1-based indices
        col1, col2 = map(int, lines[1].strip().split("\t"))
        return col1, col2

def calculate_average(values):
    avg_value = sum(values) / len(values)
    return avg_value

def get_avg_and_max(input_file, col1, col2):
    with open(input_file, 'r') as file:
        lines = file.readlines()
        # The second row contains the values, skip the first column (name column)
        values = list(map(float, lines[1].strip().split("\t")[1:]))
        # Adjust for 0-based indexing in Python by subtracting 1 from col1 and col2
        masked_values = values[:col1-1] + values[col2:]  # Exclude the range col1 to col2
        # Calculate average and maximum
        avg_value = calculate_average(masked_values)
        max_value = max(masked_values)
        return avg_value, max_value

def write_output(output_file, results):
    with open(output_file, 'w') as file:
        # Write each result (average, maximum, and relative difference) with corresponding label, separated by tabs, rounded to 3 decimal places
        for i, (avg_value, max_value) in enumerate(results):
            if avg_value != 0:
                relative_difference = (max_value - avg_value) / avg_value  # Calculate (max - avg) / avg
            else:
                relative_difference = "N/A"  # Handle division by zero case
            file.write(f"{labels[i]}\t{round(avg_value, 3)}\t{round(max_value, 3)}\t{round(relative_difference, 3) if relative_difference != 'N/A' else relative_difference}\n")

if __name__ == "__main__":
    if len(sys.argv) != 19:  # 16 input files + MASKED_region.tab + 1 output file + script name
        print("Usage: Python JOB.py input1.tab input2.tab ... input16.tab MASKED_region.tab output.tab")
        sys.exit(1)

    input_files = sys.argv[1:17]  # First 16 arguments are input files
    masked_file = sys.argv[17]    # The 17th argument is the MASKED_region.tab file
    output_file = sys.argv[18]    # The 18th argument is the output file

    # Get the masked region (col1, col2) from the MASKED_region.tab file
    col1, col2 = get_masked_region(masked_file)

    # Get the average and max values from each file, excluding the masked region
    results = [get_avg_and_max(input_file, col1, col2) for input_file in input_files]

    # Write the results to the output file
    write_output(output_file, results)

    print(f"Average, maximum values, and (max-avg)/avg ratio (excluding columns {col1} to {col2}) extracted and saved to {output_file}")
