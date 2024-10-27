import sys

def process_and_write_files(output_file, substrings_to_remove, input_files):
    with open(output_file, 'w') as outfile:
        mode_count = 0  # Track occurrences of "Mode" across all files
        
        for i, input_file in enumerate(input_files):
            print(f"Processing file: {input_file}")  # Debug print to see which file is being processed
            
            with open(input_file, 'r') as infile:
                for line in infile:
                    # Modify specific headers only in the output file
                    if "Max Value" in line:
                        line = line.replace("Max Value", "Max_Value")
                    if "Column 5" in line:
                        line = line.replace("Column 5", "Scale")
                    if "Filename" in line:
                        line = line.replace("Filename", "Motif#_Quartile")
                    if "Max Position (Header Row)" in line:
                        line = line.replace("Max Position (Header Row)", "Position_Max_Value(Index)")
                    if "Adjusted Column Index" in line:
                        line = line.replace("Adjusted Column Index", "Position_Max_Value(bp)")
                    if "Unique, significant peaks 5' to motif" in line:
                        line = line.replace("Unique, significant peaks 5' to motif", "Unique, significant peaks 5' to motif (motif strand)")
                    if "Unique, significant peaks 3' to motif" in line:
                        line = line.replace("Unique, significant peaks 3' to motif", "Unique, significant peaks 3' to motif (opposite strand)")
                    if "periodicity" in line:
                        line = line.replace("periodicity", "Periodicity:")
                    if "category" in line:
                        line = line.replace("category", "quartile")
                    if "Phase (0,9):" in line and "\t10" in line:
                        line = line.replace("\t10", "\t0")

                    # Handle "Mode" replacements based on count
                    if "Mode" in line:
                        mode_count += 1
                        if mode_count == 1:
                            line = line.replace("Mode", "Motif strand phase (0,9):")
                        elif mode_count == 2:
                            line = line.replace("Mode", "Opposite strand phase (0,9):")

                    # Remove specified substrings
                    for substring in substrings_to_remove:
                        line = line.replace(substring, "")

                    # Print final line for debugging
                    print(f"Final line for {input_file}: {line.strip()}")
                    
                    # Write the processed line to the output file
                    outfile.write(line)
            
            # Add a blank line only between input4.tab and input5.tab
            if i == 3:
                outfile.write("\n")  # Blank line between input4 and input5
            
            # Add separator between other files (but not after input4.tab)
            elif i < len(input_files) - 1:
                outfile.write("****************\n")

if __name__ == "__main__":
    if len(sys.argv) != 10:
        print("Usage: python JOB.py input1.tab input2.tab input3.tab input4.tab input5.tab input6.tab input7.tab input8.tab output.tab")
    else:
        # Substrings to remove from the output (only for input1.tab)
        substrings_to_remove = [
            "BNase-seq_50U-10min_merge_hg38_", "_final_1000bp_intersected_164bp",
            "_1000bp_ForComposite_allReads_sense_smoothed_20bp.tab"
        ]

        # Process all input files (from input1 to input8)
        input_files = sys.argv[1:9]
        process_and_write_files(sys.argv[-1], substrings_to_remove, input_files)
