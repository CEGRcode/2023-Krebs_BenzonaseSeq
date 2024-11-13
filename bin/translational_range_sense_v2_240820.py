import sys

def find_max_min(input_file, output_file):
    with open(input_file, 'r') as f:
        rows = [line.strip().split('\t') for line in f.readlines()]
    
    # Extract relevant portion
    start_index = 151 - 1  # convert 1-based to 0-based
    end_index = 351  # already inclusive in Python slice

    positions = rows[0][start_index:end_index]
    values = list(map(float, rows[1][start_index:end_index]))

    # Find max value and its position
    max_value = max(values)
    max_index = values.index(max_value)
    max_position = positions[max_index]

    # Find min value before max position
    min_value = float('inf')
    min_position = None
    
    for i in range(max_index):
        if values[i] < min_value:
            min_value = values[i]
            min_position = positions[i]

    # Calculate the range (max - min)
    value_range = max_value - min_value if min_position is not None else 'N/A'

    # Output the results
    with open(output_file, 'w') as f:
        f.write(f"Max_Value\tMax_Position\n")
        f.write(f"{max_value}\t{max_position}\n")
        if min_position is not None:
            f.write(f"Min_Value\tMin_Position\n")
            f.write(f"{min_value}\t{min_position}\n")
        f.write(f"Range\n")
        f.write(f"{value_range}\n")

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    find_max_min(input_file, output_file)