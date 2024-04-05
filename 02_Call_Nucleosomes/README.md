Perform peak-calling on Benzonase-seq data to infer nucleosomes and build RefPT for various nucleosome particles and nucleosome related RefPT (also TSS RefPT).

<details>
<summary> Full execution summary
</summary>

```
data
  |--RefPT-Krebs
02_Call_Nucleosomes
  |--Merged_Redundant_Nucleosome-Particles.bed
  |--AllParticles
    |--Subtetra.bed
    |--Tetra.bed
    |--Hex.bed
    |--Nucleosome.bed
    |--Supraoct.bed
  |--SCIDX
    |--sub.tab
    |--tet.tab
    |--hex.tab
    |--nuc.tab
    |--sup.tab
  |--ShiftCheck
    |--Subtetra_10k.bed
    |--Subtetra_10k_1000bp.bed
    |--Subtetra_10k_1000bp_midpoint.out
    |--Subtetra_10k_1000bp_midpoint_combined.cdt
    |--Tetra_10k.bed
    |--Tetra_10k_1000bp.bed
    |--Tetra_10k_1000bp_midpoint.out
    |--Tetra_1000bp_midpoint_combined.cdt
    |--Hex_10k.bed
    |--Hex_10k_1000bp.bed
    |--Hex_10k_1000bp_midpoint.out
    |--Hex_1000bp_midpoint_combined.cdt
    |--Nucleosome_10k.bed
    |--Nucleosome_10k_1000bp.bed
    |--Nucleosome_10k_1000bp_midpoint.out
    |--Nucleosome_1000bp_midpoint_combined.cdt
    |--Supraoct_10k.bed
    |--Supraoct_10k_1000bp.bed
    |--Supraoct_10k_1000bp_midpoint.out
    |--Supraoct_1000bp_midpoint_combined.cdt
  |--sub
    |--expanded.bed
    |--formatted.tab
    |--genetrack_output.bed
    |--genetrack_s10e20/formatted_s10e20.gff
  |--tet
    |--expanded.bed
    |--formatted.tab
    |--genetrack_output.bed
    |--genetrack_s20e40F5/formatted_s20e40F5.gff
  |--hex
    |--expanded.bed
    |--formatted.tab
    |--genetrack_output.bed
    |--genetrack_s30e60F6/formatted_s30e60F6.gff
  |--nuc
    |--expanded.bed
    |--formatted.tab
    |--genetrack_output.bed
    |--genetrack_s40e80F5/formatted_s40e80F5.gff
  |--sup
    |--expanded.bed
    |--formatted.tab
    |--genetrack_output.bed
    |--genetrack_s50e100F3/formatted_s50e100F3.gff
```

</details>

### 1_Benzonase_Peak_Calling.sbatch
Write scIDX formatted pileups of merged Benzonase-seq data filtered by various sub-nucleosome sized fragments. Format as a BED file expanded to the corresponding particle size and with a uniqueID.
```
sbatch 1_Benzonase_Peak_Calling.sbatch
```

### 1b_Check_Shift.sbatch
Perform a positioning check with a subsampled Tag Pileup composite.
```
sbatch 1b_Check_Shift.sbatch
```
