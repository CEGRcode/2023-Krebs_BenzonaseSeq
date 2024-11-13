import sys

def process_file(file_path, filter_condition):
    values = []
    with open(file_path, 'r') as file:
        headers = file.readline().strip().split('\t')
        for line in file:
            row = list(map(float, line.strip().split('\t')))
            if filter_condition(row[1]):  # Column 2 is index 1 (1-based index)
                values.append(row[4])  # Column 5 is index 4 (1-based index)
    return values

def calculate_average(values):
    if not values:
        return None
    return sum(values) / len(values)

def main(input1, input2, output):
    # Process input1.tab for rows with column 2 < 501
    values1 = process_file(input1, lambda x: x < 501)
    
    # Process input2.tab for rows with column 2 > 501
    values2 = process_file(input2, lambda x: x > 501)
    
    # Combine all values into a single list
    all_values = values1 + values2
    
    # Calculate the overall average for column 5
    overall_avg = calculate_average(all_values)
    
    # Output the result to output.tab
    with open(output, 'w') as out_file:
        out_file.write('average_range\n')
        if overall_avg is not None:
            out_file.write(f'{overall_avg:.3f}\n')
        else:
            out_file.write('No data met the criteria.\n')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python JOB.py input1.tab input2.tab output.tab")
    else:
        input1 = sys.argv[1]
        input2 = sys.argv[2]
        output = sys.argv[3]
        main(input1, input2, output)