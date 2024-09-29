
Call "lowly bound" motif reference points for motif-centered figures.

<details>
<summary> Full execution summary
</summary>

```
data
  |--JASPAR
    |--ATF7_MA0834-1.meme
    |--BACH1_MA1633-1.meme
    |--CTCF_MA1930-1.meme
    |--ELF1_MA0473-3.meme
    |--MAX_MA0058-3.meme
    |--MEIS2_MA1640-1.meme
    |--NFIC_MA1527-1.meme
    |--REST_MA0138-2.meme
    |--SP1_MA0079-5.meme
    |--SPI1_MA0080-6.meme
    |--ZKSCAN1_MA1585-1.meme
  |--RefPT-Motifs
    |--ATF7_Bound-LowerHalf.bed
    |--ATF7_Bound-NonPromoter.bed
    |--ATF7_Bound-Promoter.bed
    |--BACH1_Bound-LowerHalf.bed
    |--CTCF_Bound-LowerHalf.bed
    |--ELF1_Bound-LowerHalf.bed
    |--MAX_Bound-LowerHalf.bed
    |--MEIS2_Bound-LowerHalf.bed
    |--NFIC_Bound-LowerHalf.bed
    |--REST_Bound-LowerHalf.bed
    |--SP1_Bound-LowerHalf.bed
    |--SPI1_Bound-LowerHalf.bed
    |--ZKSCAN1_Bound-LowerHalf.bed
    |--1000bp/
      |--...<see above with _1000bp.bed suffix>
03_Call_Motifs
  |--narrowPeak
    |--ATF7_ENCFF868QLL.bed.gz
    |--BACH1_ENCFF423EMU.bed.gz
    |--CTCF_ENCFF738TKN.bed.gz
    |--ELF1_ENCFF392MUM.bed.gz
    |--MAX_ENCFF422NGZ.bed.gz
    |--MEIS2_ENCFF613RNG.bed.gz
    |--NFIC_ENCFF370ENX.bed.gz
    |--REST_ENCFF895QLA.bed.gz
    |--SP1_ENCFF300XUA.bed.gz
    |--SPI1_ENCFF664XPS.bed.gz
    |--ZKSCAN1_ENCFF163VUK.bed.gz
```

</details>

### 0_Download_JAPAR_and_ENCODE_data.sh
Download PWMs from JASPAR (`.meme`) and ChIP binding peaks from ENCODE (`.bed.gz`)
```
sh 0_Download_JAPAR_and_ENCODE_data.sh
```
For each TF, there should be one of each the following files:
```
data/JASPAR/TF_MAXXXX-X.meme
03_Call_Motifs/narrowPeaks/TF_ENCFFXXXXXX.meme
```

### 1_FIMO_Motifs_from_Genome.sbatch
Run FIMO for each JASPAR motif to get all instances in the genome and filter to keep motifs that are at least 500bp from all blacklisted regions.
```
sbatch 1_FIMO_Motifs_from_Genome.sbatch
```
Results in four intermediate files for each TF:
```
FIMO/TF/fimo.gff
FIMO/TF/fimo.bed
FIMO/TF/fimo_1000bp.bed
FIMO/TF/filtered.bed
```

### 2_Intersect_Motifs_wENCODE_ChIP-seq_peaks.sbatch
Get "lowly bound" motifs by intersecting the motif occurrences with ENCODE's ChIP-seq peaks and filtering for the bottom half of motifs (sorted by ENCODE binding score),
```
sbatch 2_Intersect_Motifs_wENCODE_ChIP-seq_peaks.sbatch
```
Along with intermediate files, two final motif RefPT files are generated for each TF:
```
../data/RefPT-Motif/TF.bed
../data/RefPT-Motif/1000bp/TF_1000bp.bed
```

### 3_Filter-ATF7_PromoterProximal.sh
Intersect `BoundMotifs.bed` instances with NFRs to get `Promoter` and `NonPromoter` motif sets.
```
sh 3_Filter-ATF7_PromoterProximal.sh
```
```
../data/RefPT-Motif/ATF7_PromoterProximal.bed
../data/RefPT-Motif/ATF7_NotPromoterProximal.bed
../data/RefPT-Motif/1000bp/ATF7_PromoterProximal_1000bp.bed
../data/RefPT-Motif/1000bp/ATF7_NotPromoterProximal_1000bp.bed
```

### 4_Re-MEME_ZKSCAN1_LowlyBoundMotif.sbatch
Re-MEME the lowly bound motif instances to get an extended motif and tally the bound motif instances (full set vs lowly bound half) for how many contain the extended motif.
```
sbatch 4_Re-MEME_ZKSCAN1_LowlyBoundMotif.sbatch
```

```
ZKSCAN1_ExtendedMotif/ZKSCAN1_LowerBound_Re-MEME-M1.meme.txt
ZKSCAN1_ExtendedMotif/FIMO-Full/fimo.gff
ZKSCAN1_ExtendedMotif/FIMO-Lower/fimo.gff
```
