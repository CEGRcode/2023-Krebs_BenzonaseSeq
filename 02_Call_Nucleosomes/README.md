Perform peak-calling on Benzonase-seq data to infer nucleosomes and build RefPT for various nucleosome particles and nucleosome related RefPT (also TSS RefPT).

<details>
<summary> Full execution summary
</summary>

```
data
  |--RefPT-Krebs
    |--MinusOneDyad_SORT-DistToExpressedTSS.bed
    |--MinusOneDyad_SORT-DistToUnexpressedTSS.bed
    |--MinusOneDyad_SORT-Expression.bed
    |--NFR_SORT-NFRLength.bed
    |--PlusOneDyad_SORT-DistToExpressedTSS.bed
    |--PlusOneDyad_SORT-DistToUnexpressedTSS.bed
    |--PlusOneDyad_SORT-Expression.bed
    |--PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad.bed
    |--PlusOneDyad_SORT-Expression_WithUnexpressed.bed
    |--TSS_GROUP-All_SORT-CappedExpression.bed
    |--TSS_GROUP-Expressed_SORT-CpG.bed
    |--TSS_GROUP-Expressed_SORT-Expression.bed
    |--TSS_GROUP-Unexpressed.bed
  |--RefPT-Other
    |--CpG_Islands.bed
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
    |--Capped_READ2_TSS_200bp_anti.cdt
    |--Capped_READ2_TSS_200bp_sense.cdt
    |--hg19.knownCanonicalPep.id-transcripts.gtf
    |--hg19_knownCanonicalPep-TSS_200bp.bed
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
    |--TSS_200bp.bed
    |--TSS.bed
    |--TSS_CpG.cdt
    |--TSS_CpG_SORTED.cdt
    |--TSS_SCORE-CappedExpression.bed
  |--MakePlusMinus
    |--MatchedDyads_SORT-RankExpression.tsv
    |--MatchedDyads_SORT-RankExpression_WithNFRInfo.tsv
    |--Matched-MinusOneDyad_SORT-DistToExpressedTSS.tsv
    |--Matched-MinusOneDyad_SORT-RankExpression.tsv
    |--Matched-PlusOneDyad_SORT-DistToExpressedTSS.tsv
    |--Matched-PlusOneDyad_SORT-RankExpression.tsv
    |--MinusOneDyad_SORT-DistToExpressedTSS.tsv
    |--MinusOneDyad_SORT-DistToUnexpressedTSS.tsv
    |--Nuc-Dyad.ids
    |--Nucleosomes.bed
    |--PlusOneDyad_SORT-DistToExpressedTSS.tsv
    |--PlusOneDyad_SORT-DistToUnexpressedTSS.tsv
    |--Shared_RankIDs.ids
    |--TSS_downstream_Octomers.bed
    |--TSS_SORT-Genomic.bed
    |--TSS_SORT-RankExpression_1bp.bed
    |--TSS_SORT-RankExpression.bed
    |--TSS_upstream_Octomers.bed
    |--uTSS_downstream_Octomers.bed
    |--uTSS_SORT-Genomic.bed
    |--uTSS_SORT-RankSort_1bp.bed
    |--uTSS_SORT-RankSort.bed
    |--uTSS_upstream_Octomers.bed
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

### 0_Download_CpG_reffile.sh
Download CpG reference file and format into BED.
```
sh 0_Download_CpG_reffile.sh
```

You can download the reference file from [here](https://genome.ucsc.edu/cgi-bin/hgTables) with the following options selected:

![hgTables-img](/img/hgTables-CpG-download.png)

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

```
../data/RefPT-Krebs/TSS_Sort-CpG.bed
../data/RefPT-Krebs/2000bp/TSS_Sort-CpG_2000bp.bed
```


### 5_Determine_PlusOne-MinusOne-Dyads.sh
```
sh 5_Determine_PlusOne-MinusOne-Dyads.sh
```

### 6_Match_Dyads_for_NFR.sh
```
sh 6_Match_Dyads_for_NFR.sh
```

### 7_SortDyad_pHex-dHex.sh
```
sh 7_SortDyad_pHex-dHex.sh
```

### 8_SortDyad_pHN-dHN.sh
```
sh 8_SortDyad_pHN-dHN.sh
```
