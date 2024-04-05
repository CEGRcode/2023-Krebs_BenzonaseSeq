Perform peak-calling on Benzonase-seq data to infer nucleosomes and build RefPT for various nucleosome particles and nucleosome related RefPT (also TSS RefPT).

<details>
<summary> Full execution summary
</summary>

```
02_Call_Nucleosomes
  |--SCIDX
    |--sub.tab
    |--tet.tab
    |--hex.tab
    |--nuc.tab
    |--sup.tab
  |--sub
    |--formatted.tab
  |--tet
    |--formatted.tab
  |--hex
    |--formatted.tab
  |--nuc
    |--formatted.tab
  |--sup
    |--formatted.tab
```

</details>

### 1_Benzonase_Peak_Calling.sbatch
Write scIDX formatted pileups of merged Benzonase-seq data filtered by various sub-nucleosome sized fragments.
