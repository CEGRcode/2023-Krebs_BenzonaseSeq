import sys
import statistics

def get_col1_value(input2_file):
    with open(input2_file, 'r') as infile:
        lines = infile.readlines()
        # Get the second row's first column value (1-based index, so index 0)
        col1_value = int(lines[1].strip().split('\t')[0])
    return col1_value

def process_column6(input_file, output_file, col1, mode):
    column6_values = []

    with open(input_file, 'r') as infile:
        for line in infile:
            columns = line.strip().split('\t')

            if (mode=="sense"):
                # Check if the value in column 8 is less than col1
                if int(columns[7]) < col1:  # Column 8 is index 7 (0-based)
                    last_char = columns[5][-1]  # Get the last character of column 6 (index 5, 0-based)
                    if last_char == '0':
                        last_char = '10'
                    column6_values.append(last_char)
            elif (mode=="anti"):
                # Check if the value in column 8 is greater than col1
                if int(columns[7]) > col1:  # Column 8 is index 7 (0-based)
                    last_char = columns[5][-1]  # Get the last character of column 6 (index 5, 0-based)
                    if last_char == '0':
                        last_char = '10'
                    column6_values.append(last_char)

    if column6_values:
        # Calculate the mode
        mode_value = statistics.mode(column6_values)
    else:
        mode_value = "No valid rows to calculate mode"

    # Write the mode to the output file
    with open(output_file, 'w') as outfile:
        outfile.write(f"Mode\t{mode_value}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py mode input1.tab input2.tab output.tab")
    else:
        mode = sys.argv[1]
        input1_file = sys.argv[2]
        input2_file = sys.argv[3]
        output_file = sys.argv[4]

        try:
            col1 = get_col1_value(input2_file)
            process_column6(input1_file, output_file, col1, mode)
        except ValueError:
            print("Error: col1 value must be an integer.")
        except IndexError:
            print("Error: input2.tab does not have a second row or column 1 is missing.")