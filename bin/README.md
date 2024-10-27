
This directory stores generally used scripts and executables.

### ScriptManager-v0.15.jar

Download this is the Java binary executable for ScriptManager that includes a collection of tools including TagPileup which is used to count tags and calculate coverage of samples around reference points.

```
wget https://github.com/CEGRcode/scriptmanager/releases/download/v0.15/ScriptManager-v0.15.jar
```

### bedGraphToBigWig
Download the appropriate binary for your OS from UCSC.

### bigBedToBed
Download the appropriate binary for your OS from UCSC.
```
# Linux
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
```

### dinucleotide_CDT_from_FASTA.py
Scan a FASTA file for positional dinucleotide content
```
python dinucleotide_CDT_from_FASTA.py  -h
usage: dinucleotide_CDT_from_FASTA.py [-h] -i fasta_fn -s dinucleotides_str -o tsv_fn

============
Get 0/1 matrix (CDT format) of dinucleotides
============

optional arguments:
 -h, --help            show this help message and exit
 -i fasta_fn, --input fasta_fn
                       the FASTA file to analyze
 -s dinucleotides_str, --seq dinucleotides_str
                       the "-" delimited set of dinucleotides to check for
 -o tsv_fn, --output tsv_fn
                       the output CDT formatted 0/1 matrix of dinucleotide matches
```

### sum_Col_CDT.pl
This script sums the columns of a CDT matrix file by column values (CDT to composite).
```
usage:		perl sum_Col_CDT.pl	Input_CDT_File	Output_TAB_File
Example:	perl sum_Col_CDT.pl input.cdt composite.out
```

### generate_BAM_file_from_PEGR.py
Download this script from `EGC_utility_scripts` repository.

### generate_FQ_file_from_PEGR.py
Download this script from `EGC_utility_scripts` repository.

### kmer_tally_to_pwm.py
Can feed this script the tally output from `upstream_seq_tally.py` to summarize the kmer tallies' positional nucleotide content.
```
usage: kmer_tally_to_pwm.py [-h] -i bam_fn -o tsv_fn [-c COLUMN]
kmer_tally_to_pwm.py: error: the following arguments are required: -i/--input, -o/--output
```

### make_excel_composite_v2.py
Concatenate all composite files from an input directory into a formatted Excel file to use with Excel plotting tool.
```
usage: make_excel_composite_v2.py [-h] -i composite-dir -o outfile

This script takes a directory of composite (*.out) files and combines them into an excel spreadsheet.

optional arguments:
  -h, --help                                    show this help message and exit
  -i <composite-dir>, --input <composite-dir>   directory with all the composite data files
  -o <outfile>, --output <outfile>              output name to save workbook to
```

### make_violin_plot.py
Make violin plots with presets formatted for this manuscript.
```
usage: make_violin_plot.py [-h] [-i two_col_file] [--width width] [--height height] [--title title] [--xlabel xlabel] [--ylabel ylabel] [--preset1] [--preset2] [-o output_svg]

optional arguments:
  -h, --help            show this help message and exit
  -i two_col_file, --input two_col_file
                        tab-delimited file made of two columns: first column y values to plot (must all be numeric values), second column is the grouping (which violin group along x-axis to contribute to)
  --width width         width of figure
  --height height       height of figure
  --title title         title of figure
  --xlabel xlabel       x-axis label
  --ylabel ylabel       y-axis label
  --preset1             use proximal/distal presets for nucleosome intervals (Fig 4d)
  --preset2             use proximal/distal presets for half-nucleosome intervals (Fig 4d)
  -o output_svg, --output output_svg
                        name of SVG filepath to save figure to (if none provided, figure pops up in new window)
```

### upstream_seq_tally.py
Tally up kmers upstream of each read (can use proper-pair flag for paired-end data). Reverse-stand mapped reads are reverse complemented to orient kmers with the read on the right.
```
usage: upstream_seq_tally.py [-h] -i bam_fn -g fasta_fn -o tsv_fn [-p]
                             [-k KMER]
upstream_seq_tally.py: error: the following arguments are required: -i/--input, -g/--genome, -o/--output
```
