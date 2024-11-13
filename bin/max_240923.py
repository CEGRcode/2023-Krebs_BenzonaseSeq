import sys

def get_max_value(input_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()
        # The second row contains the values, skip the first column (name column)
        values = list(map(float, lines[1].strip().split("\t")[1:]))
        max_value = max(values)
        return max_value

def write_output(output_file, max_values):
    with open(output_file, 'w') as file:
        # Writing labeled maximum values, rounded to 1 decimal place
        file.write(f"Maximum value for WW: {round(max_values[0], 3)}\n")
        file.write(f"Maximum value for SS: {round(max_values[1], 3)}\n")
        file.write(f"Maximum value for YY: {round(max_values[2], 3)}\n")
        file.write(f"Maximum value for RR: {round(max_values[3], 3)}\n")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: Python JOB.py input1.tab input2.tab input3.tab input4.tab output.tab")
        sys.exit(1)

    input_files = sys.argv[1:5]
    output_file = sys.argv[5]

    # Get the max values from each file
    max_values = [get_max_value(input_file) for input_file in input_files]

    # Write the results to the output file
    write_output(output_file, max_values)

    print(f"Maximum values extracted and saved to {output_file}")
