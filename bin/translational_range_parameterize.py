import sys

if __name__ == "__main__":
    mode = sys.argv[1]
    input_file = sys.argv[2]
    output_file = sys.argv[3]

    # Extract relevant portion
    if (mode=="sense"):
        start_index = 151 - 1  # convert 1-based to 0-based
        end_index = 351  # already inclusive in Python slice
    elif (mode=="anti"):
        start_index = 651 - 1  # convert 1-based to 0-based
        end_index = 851  # already inclusive in Python slice
    else:
        raise Exception("Unexpected mode value. Please use one of ['sense','anti']")

    with open(input_file, 'r') as f:
        rows = [line.strip().split('\t') for line in f.readlines()]


    positions = rows[0][start_index:end_index]
    values = list(map(float, rows[1][start_index:end_index]))

    # Find max value and its position
    max_value = max(values)
    max_index = values.index(max_value)
    max_position = positions[max_index]

    # Find min value before (sense)/after (anti) max position
    min_value = float('inf')
    min_position = None

    if (mode=="sense"):
        for i in range(max_index):
            if values[i] < min_value:
                min_value = values[i]
                min_position = positions[i]
    elif (mode=="anti"):
        for i in range(max_index + 1, len(values)):
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