Perform peak-calling on Benzonase-seq data to infer nucleosomes and build RefPT for various nucleosome particles and nucleosome related RefPT (also TSS RefPT).

<details>
<summary> Full execution summary
</summary>

```
data
  |--RefPT-Krebs
    |--TSS_GROUP-All_SORT-CappedExpression.bed
    |--TSS_GROUP-Expressed_SORT-Expression.bed
    |--TSS_GROUP-Unexpressed.bed
02_Call_Nucleosomes
  |--Merged_Redundant_Nucleosome-Particles.bed
  |--AllParticles
    |--Subtetra.bed
    |--Tetra.bed
    |--Hex.bed
    |--Nucleosome.bed
    |--Supraoct.bed
  |--UniqueParticles
    |--uHex.bed
    |--uTetra.bed
    |--uSubtetra.bed
    |--uSupraoct.bed
    |--Nucleosome_uHex.bed
    |--Nucleosome_uHex_uTetra.bed
    |--Nucleosome_uHex_uTetra_uSubtetra.bed
    |--Merged_Nonredundant_particles.bed
  |--Intersect
    |--Nucleosomes_intersect_redundantHex.bed
    |--Nucleosomes_intersect_redundantTetra.bed
    |--Nucleosomes_intersect_redundantSubtetra.bed
    |--Nucleosomes_intersect_redundantSupraoct.bed
    |--uHex_intersect_redundantTetra.bed
    |--uHex_intersect_redundantSubtetra.bed
    |--uHex_intersect_redundantSupraoct.bed
    |--uTetra_intersect_redundantSubtetra.bed
    |--uTetra_intersect_redundantSupraoct.bed
    |--uSubtetra_intersect_redundantSupraoct.bed
  |--MakeTSS
    |--CappedExpression.out
    |--Capped_READ2_anti.cdt
    |--Capped_READ2_sense.cdt
    |--Capped_READ2_TSS_40bp_anti.cdt
    |--Capped_READ2_TSS_40bp_sense.cdt
    |--hg19.knownCanonicalPep.id-transcripts.gtf
    |--hg19_knownCanonicalPep-TSS_100bp.bed
    |--hg19_knownCanonicalPep-TSS.bed
    |--hg19.knownGene.id-transcripts.gtf
    |--hg19.knownGene.transcripts.gtf
    |--hg19.knownGene.transcripts.ids
    |--knownCanonical.ids
    |--knownCanonicalPep.ids
    |--knownCanonicalPep-NoNames.txt
    |--knownGenePep.ids
    |--knownToLynx_FILTER-RemoveMalacards.txt
    |--knownToLynx.txt
    |--knownToMalacards.ids
    |--knownToMalacards-wLynx.txt
    |--knownTo_NameMap.txt
    |--maxPeak.bed
    |--TSS.bed
    |--TSS_200bp.bed
    |--TSS_SCORE-CappedExpression.bed
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

### 2_Identify_Unique_Peaks.sbatch
Create a non-redundant set of non-overlapping (complete overlap) particle peaks (favoring peaks from larger fragments).
```
sbatch 2_Identify_Unique_Peaks.sbatch
```

### 3_Aggregate_Nucleosome_Peaks.sbatch

```
sbatch 3_Aggregate_Nucleosome_Peaks.sbatch
```

### 4_Build_TSS_RefPT.sbatch
Call TSS reference points, trued up by CoPRO mode signal.
```
sbatch 4_Build_TSS_RefPT.sbatch
```
