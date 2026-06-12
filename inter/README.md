# `inter` dataset
This is a bacterial interaction dataset, as described in [*Luzzatto-Knaan et al.*](https://pubs.acs.org/doi/10.1021/acschembio.8b01120) Study involved *Paenibacillus dendritiformis* (Pd), *Bacillus subtilis* NCIB 3610 (Bs, WT) and the NRPS-mutated dsrf, dpps, and dd that are Surfactin or/and Plipastatin Synthetase deficient. Metabolites were extracted by MeOH directly from the agar cut into small pieces and sample were analyzed together with Media blanks on RP column in DDA positive mode on Bruker TIMS-TOF Pro 2. Raw mzXML files were then processed in [`MS-Dial`](https://systemsomicslab.github.io/compms/msdial/main.html). <br/>

Dataset was used to estimate the applicability of the MetaboCensoR App on Statistical Hypothesis Testing.<br/>
Code script and all relevant data are available in the folder. 
- Main script for analyzing data: [`Process inter.R`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/inter/Process%20inter.R)
- List of known target compounds: [`annot table.csv`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/inter/annot%20table.csv)
- Raw Data: [`MSV000100949`](https://doi.org/doi:10.25345/C5930P825)
<table>
  <tr>
    <td>
      <img src="https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/inter/volcano.png" height="400">
    </td>
    <td>
      <img src="https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/inter/lollipop.png" height="400">
    </td>
  </tr>
</table>
