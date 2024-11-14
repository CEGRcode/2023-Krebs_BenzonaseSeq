import sys
import pandas as pd
import numpy as np
from scipy import stats

def calculate_bin_statistics(df, start_col, end_col, adjust_position=False, recenter_bins=False):
    bin_stats = []
    bin_ratios = []

    # Loop through the 10 bp bins
    for start in range(start_col, end_col + 1, 10):
        stop = min(start + 9, end_col)

        bin_values = df.iloc[:, start-1:stop]
        bin_values_row2 = bin_values.iloc[1]

        bin_max = bin_values_row2.max()
        bin_min = bin_values_row2.min()
        bin_ratio = bin_max / bin_min if bin_min != 0 else np.nan  # Updated to max/min (Ratio)

        # Calculate the 1-indexed positions of max and min values
        max_pos = bin_values_row2.idxmax() + 1
        min_pos = bin_values_row2.idxmin() + 1

        # Calculate the index of the max and min values
        index_max = bin_values_row2.idxmax()
        index_min = bin_values_row2.idxmin()

        # Adjust the positions by subtracting 501
        if adjust_position:
            max_pos_adjusted = max_pos - 501
            min_pos_adjusted = min_pos - 501
        else:
            max_pos_adjusted = max_pos
            min_pos_adjusted = min_pos

        # Store the ratio for later calculation
        bin_ratios.append(bin_ratio)

        # Only retain the required columns
        bin_stats.append([start, stop, bin_max, bin_min, bin_ratio, max_pos_adjusted, min_pos_adjusted, index_max + 1, index_min + 1])

    if recenter_bins:
        # Calculate mode based on bins within columns 301-500
        mode_start_col = 301
        mode_end_col = 500
        mode_bin_stats = []

        for start in range(mode_start_col, mode_end_col + 1, 10):
            stop = min(start + 9, mode_end_col)

            bin_values = df.iloc[:, start-1:stop]
            bin_values_row2 = bin_values.iloc[1]

            max_pos = bin_values_row2.idxmax() + 1
            mode_bin_stats.append(int(str(max_pos - 501)[-1]) if int(str(max_pos - 501)[-1]) != 0 else 10)

        mode_last_char = stats.mode(mode_bin_stats)[0][0]

        shifted_bin_stats = []
        for stats_row in bin_stats:
            shift = mode_last_char - 5
            shifted_start = stats_row[0] + shift
            shifted_stop = stats_row[1] + shift

            # Add shifted bins to the result
            shifted_bin_stats.append([shifted_start, shifted_stop] + stats_row[2:])

        return pd.DataFrame(shifted_bin_stats, columns=['Start_Column', 'End_Column', 'Max', 'Min', 'Ratio', 'Max_Position_Adjusted', 'Min_Position_Adjusted', 'Index_Max', 'Index_Min']), bin_ratios
    else:
        return pd.DataFrame(bin_stats, columns=['Start_Column', 'End_Column', 'Max', 'Min', 'Ratio', 'Max_Position_Adjusted', 'Min_Position_Adjusted', 'Index_Max', 'Index_Min']), bin_ratios

def calculate_avg_plus_2sd(bin_ratios):
    avg_ratio = np.mean(bin_ratios)
    std_ratio = np.std(bin_ratios, ddof=1)  # Sample standard deviation
    avg_plus_2sd = avg_ratio + 2 * std_ratio
    return avg_plus_2sd

if __name__ == "__main__":
    mode = sys.argv[1]
    input_file = sys.argv[2]
    output_file_nucleosome = sys.argv[3]

    df = pd.read_csv(input_file, sep="\t", header=None)

    # Process nucleosome region
    nucleosome_df, _ = calculate_bin_statistics(df, 301, 701, adjust_position=True)

    # Process flank region to calculate 'avg + 2SD'
    if (mode=="sense"):
        _, bin_ratios = calculate_bin_statistics(df, 11, 110, adjust_position=False)
    elif (mode=="anti"):
        _, bin_ratios = calculate_bin_statistics(df, 892, 991, adjust_position=False)
    else:
        raise Exception("Unexpected mode value. Please use one of ['sense','anti']")
    avg_plus_2sd = calculate_avg_plus_2sd(bin_ratios)

    # Filter nucleosome bins by 'avg + 2SD'
    filtered_nucleosome_df = nucleosome_df[nucleosome_df['Ratio'] > avg_plus_2sd]

    # Save results (including 'Ratio' instead of 'Range')
    filtered_nucleosome_df.to_csv(output_file_nucleosome, sep="\t", header=False, index=False)