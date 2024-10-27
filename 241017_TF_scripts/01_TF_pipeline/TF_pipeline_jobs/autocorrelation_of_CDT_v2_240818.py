import sys, argparse
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

def getParams():
    '''Parse parameters from the command line'''
    parser = argparse.ArgumentParser(description = """
============

============
""", formatter_class = argparse.RawTextHelpFormatter)

    parser.add_argument('-i','--input', metavar='composite_fn', required=True, help='the tab-delimited matrix file to analyze (skip first row as header)')
    parser.add_argument('-o','--output', metavar='base_fn', required=True, help='the output tab-delimited correlation across lag values')
    
    parser.add_argument('--cdt', help='first two columns are index')

    parser.add_argument('--start', type=int, default=0, help='the start column (0-indexed) of the range to autocorrelate (default=1 for ScriptManager composite-style input)')
    parser.add_argument('--end', type=int, default=-1, help='the end column (0-indexed, exclusive) of the range to autocorrelate (default=max)')

    parser.add_argument('--lag-min', type=int, default=1, help='the shift value start range (default=1)')
    parser.add_argument('--lag-max', type=int, default=20, help='the shift value end range (inclusive, default=20)')

    args = parser.parse_args()
    return(args)

if __name__ == '__main__':
    '''Main program which takes in input parameters'''
    args = getParams()

    cidx = [0]
    if (args.cdt):
        cidx = [0,1]

    # Load data
    df = pd.read_csv(args.input, sep="\t", header=0, index_col=cidx)
    print(df)

    # Set default end column range to last column if default
    if (args.end==-1):
        args.end = df.shape[1]

    # Validate params
    if (args.start < 0 or args.end < 1):
        raise Exception("Cannot retrieve data range values from negative indexes (remember end range value is exclusive and also cannot be 0)")
    if (args.start > df.shape[1] or args.end > df.shape[1]):
        raise Exception("Cannot retrieve data range values larger than max column index in this data input file (max_col=" + df.shape[1] + ")")
    if (args.start >= args.end):
        raise Exception("Invalid column index range")
    if (args.lag_min < 1 or args.lag_max < 1 or args.lag_min > args.lag_max):
        raise Exception("Invalid lag range")

    # Initialize results
    row_results = []
    
    print(args.start)
    print(args.end)

    for row in range(df.shape[0]):

        # Parse composite id        
        id = df.index[row]

        # Slice range
        sliced = df.iloc[row,args.start:args.end]

        print('====(%i)====' % row)    
        print(sliced)

        for lag in range(args.lag_min, args.lag_max+1):
            corr = sliced.autocorr(lag)
            row_results.append({'index':id, 'lag':lag, 'autocorrelation':corr})

    results = pd.DataFrame(row_results)

    # Write output
    results.to_csv(f'{args.output}.tsv', sep="\t", index=False)

    # Initialize plot
    sns.set_style("ticks")
    fig, ax = plt.subplots()

    # Plot results
    sns.lineplot(data=results, x='lag', y='autocorrelation', hue='index')
    
    # Format axes and tickmarks
    # ax.xaxis.grid(False)
    # ax.set_xlim(0, 340)
    # plt.xticks([0, 73, 146, 219, 292])
    # ax.yaxis.grid(True)

    # Title
    # ax.set_title(os.path.basename(args.input))

    # Fill line
    # plt.fill_between(filedata.InsertSize.values, filedata.ReadCount.values, color='black')
    # plt.xticks()

    # Save figure as image
    sfig = ax.get_figure()
    sfig.savefig(f'{args.output}.svg')