import sys

# List of labels corresponding to each input file
labels = ["AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG", "CA", "CT", "CC", "CG", "GA", "GT", "GC", "GG"]

def get_max_value(input_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()
        # The second row contains the values, skip the first column (name column)
        values = list(map(float, lines[1].strip().split("\t")[1:]))
        max_value = max(values)
        return max_value

def write_output(output_file, max_values):
    with open(output_file, 'w') as file:
        # Write each maximum value with its corresponding label, rounded to 3 decimal places
        for i, max_value in enumerate(max_values):
            file.write(f"Maximum value for {labels[i]}: {round(max_value, 3)}\n")

if __name__ == "__main__":
    if len(sys.argv) != 18:  # 16 input files + 1 output file + script name
        print("Usage: Python JOB.py input1.tab input2.tab ... input16.tab output.tab")
        sys.exit(1)

    input_files = sys.argv[1:17]  # First 16 arguments are input files
    output_file = sys.argv[17]    # The 17th argument is the output file

    # Get the max values from each file
    max_values = [get_max_value(input_file) for input_file in input_files]

    # Write the results to the output file
    write_output(output_file, max_values)

    print(f"Maximum values extracted and saved to {output_file}")
