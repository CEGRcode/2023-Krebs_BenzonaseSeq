
This directory stores generally used scripts and executables.

### ScriptManager-v0.14.jar
Download this is the Java binary executable for ScriptManager that includes a collection of tools including TagPileup which is used to count tags and calculate coverage of samples around reference points.
```
wget https://github.com/CEGRcode/scriptmanager/releases/download/v0.14/ScriptManager-v0.14.jar
```

### bedGraphToBigWig
Download the appropriate binary for your OS from UCSC.

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


### upstream_seq_tally.py
Tally up kmers upstream of each read (can use proper-pair flag for paired-end data). Reverse-stand mapped reads are reverse complemented to orient kmers with the read on the right.
```
usage: upstream_seq_tally.py [-h] -i bam_fn -g fasta_fn -o tsv_fn [-p]
                             [-k KMER]
upstream_seq_tally.py: error: the following arguments are required: -i/--input, -g/--genome, -o/--output
```

### kmer_tally_to_pwm.py
Can feed this script the tally output from `upstream_seq_tally.py` to summarize the kmer tallies' positional nucleotide content.
```
usage: kmer_tally_to_pwm.py [-h] -i bam_fn -o tsv_fn [-c COLUMN]
kmer_tally_to_pwm.py: error: the following arguments are required: -i/--input, -o/--output
```