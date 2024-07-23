# Translational and rotational setting of nucleosomes across a human genome

### Jordan E. Krebs <sup>1,2</sup>, Haining Chen <sup>2</sup>, Olivia W. Lang <sup>2</sup>, William K.M. Lai <sup>2</sup>, B. Franklin Pugh <sup>2</sup>\*

<sup>1</sup>MD/PhD Medical Scientist Training Program, Penn State College of Medicine, Hershey, PA, USA
<sup>2</sup>Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York, 14853, USA

### Correspondence: fp265@cornell.edu

### PMID : [XXXXXXXX](https://pubmed.ncbi.nlm.nih.gov/XXXXXXXX/)
### GEO ID : [GSE266547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266547)

## Abstract
Eukaryotic DNA is wrapped around a complex of histones, such that one side of the helix is accessible, and the other is buried. Sliding and rotating DNA by 5 bp in either direction reverses this polarity. This has substantial ramifications for regulating transcription factor (TF) binding. Yet, there does not exist a means to measure the rotational setting of DNA on nucleosomes in vivo and on a genomic scale. We developed Benzonase-seq, where Benzonase cleaves and marks the accessible rotationally exposed DNA surface in addition to marking linker regions between nucleosomes. Further, Benzonase readily maps nucleosomes in CpG-rich mammalian promoters, which tend to be more resistant to A/T24 biased micrococcal nuclease (MNase). When Benzonase cleavages are analyzed in the context of TF binding sites, we determine whether the TF motif has a preferred rotational setting on the nucleosome surface. When coupled to chromatin immunoprecipitation (ChIP-exo) of histones, histone variants, and modifications, we find evidence for transcription-linked subnucleosomal structures. Together, this study reveals the translational and rotational setting of nucleosomes in K562 cells, along with altered or subnucleosomal structures.

## Directions
To recreate the figures for this manuscript, please execute the scripts in each directory in numerical order. Each directory's README includes more specific details on execution. To be more explicit, run the scripts in each directory in the following order: `00_Download_and_Preprocessing`, `01_Run_GenoPipe`, `02_Call_Nucleosomes`, `03_Call_Motifs`, `X_Bulk_Processing`, and then finally `Z_Figures`.

## Dependencies
Use the following [anaconda](https://anaconda.org/) environment initialization for setting up dependencies

```
conda create -n bx -c bioconda -c conda-forge bedtools bowtie2 bwa cutadapt meme opencv pandas samtools scipy sra-tools wget
```

For genetrack-executing script, a python2 environment needed to be created. The create command for that env is as follows:

```
conda create -n genetrack -c conda-forge -c bioconda python=2.7 numpy
```

## Table of Contents

### 00_Download_and_Preprocessing
Perform the preprocessing steps including alignment of raw sequencing data from both novel and previously published data

### 01_Run_GenoPipe
Perform quality control for genetic background on these data by running GenoPipe on the aligned BAMs.

### 02_Call_Nucleosomes
Call nucleosome positions and identify TSS and +1 nucleosome reference points with different sorts.

### 03_Call_Motifs
Build the sequence-specific transcription factor (ssTF) motif reference points.

### X_Bulk_Processing
With the BAM and BED files built from the scripts in the above directories, perform bulk read pileups for heatmaps and composites.

### Z_Figures
Copy/organize results from bulk processing into figure-specific directories corresponding to subfigures in the manuscript. Also includes custom/one-off scripts for analysis that didn't need bulk-style execution.

### data
Store large files to be globally accessed by the scripts in each directory

### bin
Generalized scripts and executables for global access by each of the numbered directories.
