import sys
import csv

def filter_rows(input1_file, input2_file, output_file):
    print(f"Reading values from {input2_file}")
    # Read col1 and col2 from the 2nd row of input2.tab
    with open(input2_file, 'r') as file2:
        reader = csv.reader(file2, delimiter='\t')
        rows = list(reader)
        print(f"Rows in input2.tab: {rows}")
        if len(rows) < 2 or len(rows[1]) < 2:
            raise ValueError("input2.tab must contain at least two rows with at least two columns.")
        try:
            col1, col2 = int(rows[1][0]), int(rows[1][1])
            print(f"Range values from input2.tab: col1={col1}, col2={col2}")
        except ValueError:
            raise ValueError("The values in input2.tab must be integers.")

    print(f"Filtering rows from {input1_file}")
    # Filter rows from input1.tab
    with open(input1_file, 'r') as file1, open(output_file, 'w', newline='') as output:
        reader = csv.reader(file1, delimiter='\t')
        writer = csv.writer(output, delimiter='\t')

        # Write the header to the output file
        header = next(reader)
        writer.writerow(header)
        print(f"Header written to output.tab: {header}")

        row_count = 0
        included_row_count = 0
        filtered_row_count = 0

        for row in reader:
            row_count += 1
            if len(row) < 8:
                print(f"Row {row_count}: Skipping row due to insufficient columns.")
                continue
            
            # Debugging output for the raw row data
            print(f"Debug: Raw data for Row {row_count}: {row}")
            
            try:
                # Attempt to convert to integer; if it's a float, convert to integer
                col8_value_str = row[7].strip()
                if '.' in col8_value_str:
                    col8_value = int(float(col8_value_str))
                else:
                    col8_value = int(col8_value_str)
                print(f"Row {row_count}: Column 8 value (after strip and conversion) = {col8_value}")
            except ValueError:
                print(f"Row {row_count}: Skipping row due to invalid data in column 8. Raw value: '{row[7]}'")
                continue  # Skip rows with invalid data in column 8

            # Debugging statements to show values being used for filtering
            print(f"Debug: Processing Row {row_count} with Column 8 value = {col8_value}")
            
            if col8_value == col1:
                print(f"Row {row_count}: Filtering out row. Column 8 value {col8_value} is equal to col1 {col1}.")
            elif col8_value == col2:
                print(f"Row {row_count}: Filtering out row. Column 8 value {col8_value} is equal to col2 {col2}.")
            elif col1 <= col8_value <= col2:
                print(f"Row {row_count}: Filtering out row. Column 8 value {col8_value} is within range {col1}-{col2}.")
            else:
                writer.writerow(row)  # Include rows that are not filtered
                included_row_count += 1
                print(f"Row {row_count}: Including row.")
                continue
            
            filtered_row_count += 1

        print(f"Total rows processed: {row_count}")
        print(f"Total rows included: {included_row_count}")
        print(f"Total rows filtered out: {filtered_row_count}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python JOB.py input1.tab input2.tab output.tab")
        sys.exit(1)

    input1_file = sys.argv[1]
    input2_file = sys.argv[2]
    output_file = sys.argv[3]

    try:
        filter_rows(input1_file, input2_file, output_file)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)