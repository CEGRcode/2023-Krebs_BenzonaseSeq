
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
