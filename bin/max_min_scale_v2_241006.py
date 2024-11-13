import sys

def find_scale_based_on_diff(input_file, masked_region_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Extract values from the second row (ignoring the first column)
    values = list(map(float, lines[1].strip().split('\t')[1:]))  # 0-based indexing

    # Slice to get the values from columns 251 to 751 (1-based -> 0-based = 250 to 750)
    values_subset = values[250:751]

    # Calculate the average of the specified range (for all values, including the masked region)
    average_value = sum(values_subset) / len(values_subset)

    # Reading the masked region from the masked_region_file
    with open(masked_region_file, 'r') as masked_f:
        masked_lines = masked_f.readlines()
        col1, col2 = map(int, masked_lines[1].strip().split('\t'))  # Reading col1 and col2 from 2nd row

    # Adjust col1 and col2 to 0-based indexing
    col1 -= 1
    col2 -= 1

    # Exclude the masked region from the subset
    values_excluded = values_subset[:col1-250] + values_subset[col2-250+1:]

    # Find the max and min values in the remaining columns (excluding the masked region)
    max_value = max(values_excluded)
    min_value = min(values_excluded)

    # Calculate the absolute difference between the average and max/min
    diff_max = abs(average_value - max_value)
    diff_min = abs(average_value - min_value)

    # Calculate the scales: 2 divided by each of the differences
    if diff_max != 0:
        scale_max = 2 / diff_max
    else:
        scale_max = float('inf')  # Handle division by zero

    if diff_min != 0:
        scale_min = 2 / diff_min
    else:
        scale_min = float('inf')  # Handle division by zero

    # The final scale is the smaller of the two scales
    final_scale = min(scale_max, scale_min)

    # Round the scales to 2 decimal places
    scale_max = round(scale_max, 2)
    scale_min = round(scale_min, 2)
    final_scale = round(final_scale, 2)

    # Write the result to the output file
    with open(output_file, 'w') as out_f:
        out_f.write(f"Average Value: {average_value}\n")
        out_f.write(f"Max Value: {max_value}\n")
        out_f.write(f"Min Value: {min_value}\n")
        out_f.write(f"Scale based on max: {scale_max}\n")
        out_f.write(f"Scale based on min: {scale_min}\n")
        out_f.write(f"Final Scale: {final_scale}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python job.py input_file.tab MASKED_region.tab output_file.tab")
    else:
        input_file = sys.argv[1]
        masked_region_file = sys.argv[2]
        output_file = sys.argv[3]
        find_scale_based_on_diff(input_file, masked_region_file, output_file)
