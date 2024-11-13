import sys

def count_numeric_rows_in_meme(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    start_keyword = "letter-probability matrix"
    end_keyword = "URL"
    count_numeric_rows = False
    numeric_row_count = 0
    for line in lines:
        if start_keyword in line:
            count_numeric_rows = True
            continue
        elif end_keyword in line:
            count_numeric_rows = False
            break
        if count_numeric_rows:
            if any(char.isdigit() for char in line):
                numeric_row_count += 1
    return numeric_row_count

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: py JOB.py <input.meme> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    numeric_row_count = count_numeric_rows_in_meme(input_file)

    # Write the numeric_row_count to the output file
    with open(output_file, 'w') as f:
        f.write(f"{numeric_row_count}\n")

    # Print the result for confirmation
    print(f"Numeric row count saved to {output_file}")