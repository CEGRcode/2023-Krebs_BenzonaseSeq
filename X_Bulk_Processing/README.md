This directory contains scripts to pileup select BAM files against select RefPT BED file and saves the analysis to `2023-Krebs_BenzonaseSeq/Library`.

Please make sure your BAM files are downloaded and indexed in the appropriate directory (`data/BAM/`) with normalization factors calculated (`data/NormalizationFactors/`) and reference points built (`data/RefPT-XXX`).

See `00_Download_and_preprocessing` to read more about how BAM files were built.
See `02_Call_Nucleosomes` and `03_Call_Motifs` to read more about how RefPT files were built.

Scripts were designed to run on linux systems with a Slurm scheduler. Simply update `WRK` to match the directory path where you cloned this repo and set the `#SBATCH --array 1,2-10` to the range you wish to execute.

- [Read more on Slurm Job arrays](https://slurm.schedmd.com/job_array.html)
- No Slurm scheduler? `SLURM_ARRAY_TASK_ID` can be hardcoded and run as a regular shell script

## 1_make_weblogos.sh

```
sh 1_make_weblogos.sh
```

## 2_MotifAnalyses.sbatch

```
sbatch 2_MotifAnalyses.sbatch
```

## 3_Midpoint_Pileups.sbatch

```
sbatch 3_Midpoint_Pileups.sbatch
```

## 4_Five_Read1_Pileups.sbatch

```
sbatch 4_Five_Read1_Pileups.sbatch
```

## 5_CoPRO_TSS_ActiveSite_Pileups.sbatch

```
sbatch 5_CoPRO_TSS_ActiveSite_Pileups.sbatch
```
