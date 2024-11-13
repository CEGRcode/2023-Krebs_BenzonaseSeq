import sys

def calculate_average(file1, file2, output_file):
    def get_value_from_file(file):
        with open(file, 'r') as f:
            for i, line in enumerate(f):
                if i == 5:  # Row 6 (0-based index)
                    columns = line.strip().split('\t')
                    if len(columns) > 0:  # Ensure there is at least 1 column
                        return float(columns[0])  # Column 1
        return None
    
    # Get the values from both files
    value1 = get_value_from_file(file1)
    value2 = get_value_from_file(file2)
    
    if value1 is None or value2 is None:
        print("Error: Could not find row 6 or column 1 in one or both files.")
        return
    
    # Calculate the average
    average = (value1 + value2) / 2
    
    # Round the average to 3 decimal places
    average_rounded = round(average, 3)
    
    # Write the rounded average to the output file
    with open(output_file, 'w') as f:
        f.write(f"{average_rounded}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py file1 file2 output_file")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_file = sys.argv[3]
    
    calculate_average(file1, file2, output_file)