import pandas as pd
import numpy as np
import sys

def smooth_data(data, window):
    # Apply {window} bp sized smoothing using convolution
    window = np.ones(window) / window
    return np.convolve(data, window, mode='same')

# Get the input and output file paths from command line arguments
window_size = int(sys.argv[1])
input_file = sys.argv[2]
output_file = sys.argv[3]

# Load the tab-delimited file
data = pd.read_csv(input_file, delimiter='\t', header=None)

# Extract the header and the first column
header = data.iloc[0]
first_column = data.iloc[1:, 0].values  # Get the first column values
numeric_data = data.iloc[1, 1:].astype(float).values

# Apply smoothing
smoothed_data = smooth_data(numeric_data, window_size)

# Convert smoothed data back to a DataFrame
smoothed_df = pd.DataFrame(smoothed_data.reshape(1, -1), columns=data.columns[1:])
# Add the first column as the first column in the output DataFrame
smoothed_df.insert(0, header[0], first_column)

# Add the header to the DataFrame
smoothed_df.columns = [header[0]] + list(header[1:])

# Save the smoothed data to the specified output file with headers
smoothed_df.to_csv(output_file, sep='\t', index=False, header=True)

print(f'Smoothed data saved to {output_file}')
