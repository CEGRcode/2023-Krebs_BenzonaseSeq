import sys

def calculate_boundaries(numeric_row_count, center_column=256):
    if numeric_row_count % 2 == 0:
        # Even case
        distance_5_prime = (numeric_row_count // 2) - 1
        distance_3_prime = numeric_row_count // 2
    else:
        # Odd case
        distance_from_center = numeric_row_count // 2

        distance_5_prime = distance_from_center
        distance_3_prime = distance_from_center

    # Calculate boundaries
    boundary_5_prime = center_column - distance_5_prime - 5
    boundary_3_prime = center_column + distance_3_prime + 5

    return boundary_5_prime, boundary_3_prime

if __name__ == "__main__":
    # Expecting command line arguments
    if len(sys.argv) != 3:
        print("Usage: py JOB.py input.tab output.tab")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Read numeric_row_count from the input file
    try:
        with open(input_file, 'r') as file:
            # Assuming the first line contains the numeric_row_count
            numeric_row_count = int(file.readline().strip())
    except ValueError:
        print("Error: The first line of the input file must be an integer representing numeric_row_count.")
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: File {input_file} not found.")
        sys.exit(1)

    # Calculate the boundaries
    boundary_5_prime, boundary_3_prime = calculate_boundaries(numeric_row_count)

    # Write the result to the output file
    with open(output_file, 'w') as file:
        file.write("5'_boundary\t3'_boundary\n")
        file.write(f"{boundary_5_prime}\t{boundary_3_prime}\n")

    print(f"Masked motif region boundaries saved to {output_file}")
