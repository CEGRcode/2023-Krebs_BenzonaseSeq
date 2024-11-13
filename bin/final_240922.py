import os
import sys

def process_file(input_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
        if lines:
            first_line = lines[0].strip()
            last_line = lines[-1].strip()
            return first_line, last_line
    return None, None

def main(input_files, output_file):
    headers = ["HelT", "MGW", "PropT", "Roll"]
    output_lines = []

    for i, input_file in enumerate(input_files):
        if os.path.exists(input_file):
            first_line, last_line = process_file(input_file)
            if first_line and last_line:
                output_lines.append(f"{headers[i]}\n{first_line}")
                output_lines.append(last_line)
                output_lines.append("***********")

    # Remove the last "***********"
    if output_lines and output_lines[-1] == "***********":
        output_lines.pop()

    # Save results to output file
    with open(output_file, 'w') as out:
        out.write('\n'.join(output_lines))

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python job.py input1.tab input2.tab input3.tab input4.tab output.tab")
    else:
        input_files = sys.argv[1:5]
        output_file = sys.argv[5]
        main(input_files, output_file)
