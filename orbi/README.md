# `orbi` dataset
The LC-MS plant profiling dataset, which is described in [*Houriet et al.*](https://pubs.acs.org/doi/10.1021/acs.analchem.4c05577) <br/>
Methanol extracts from the plant ashwagandha [*Withania somnifera (L.) Dunal*] together with blanks were analyzed on RP column in DDA positive mode on Thermo Q-Exactive Plus Orbitrap. Raw mzML files were then processed in [`mzMine`](mzmine.github.io) in default pre-settings for UPLC-DDA. <br/>
Obtained spectral data were subjected to feature-based molecular networking using the online workflow at [`GNPS 2`](https://gnps2.org/homepage) platform.

Dataset was used to estimate the applicability of the MetaboCensoR App on Molecular Networking.<br/>
Code scripts and all relevant data are available in the folder. 
- Main script for analyzing data: [`Process orbi.R`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/orbi/Process%20orbi.R)
- List of known target compounds: [`annotation.csv`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/orbi/annotation.csv)
- FBMN on Raw Data: [`GNPS Project`](https://gnps2.org/status?task=43da3f82c33f43198f40cdc82837b87d); [`Cytoscape Project`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/orbi/RAW%20FBMN.cys)
- FBMN on Data After MetaboCensoR: [`GNPS Project`](https://gnps2.org/status?task=f686537881184bf08d45e499acf18c3d); [`Cytoscape Project`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/orbi/APP%20FBMN.cys)
- IIMN on Raw Data: [`GNPS Project`](https://gnps2.org/status?task=fde43183ef3940aa94985576117b7854); [`Cytoscape Project`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/orbi/RAW%20IIMN.cys)
- IIMN on Data After MetaboCensoR: [`GNPS Project`](https://gnps2.org/status?task=0c58cb476c80455f8b4c2f82472eb683); [`Cytoscape Project`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/orbi/APP%20IIMN.cys)
<table>
  <tr>
    <td>
      <img src="https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/orbi/mn.png" height="600">
    </td>
    <td>
      <img src="https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/orbi/mn2.png" height="600">
    </td>
  </tr>
</table>

