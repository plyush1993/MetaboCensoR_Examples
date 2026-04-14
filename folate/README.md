# `folate` dataset
The human cell lines profiling dataset, which is provided by *Dr. Mariam Fokra, Dr. Nikita Sarvin, Dr. Tomer Shlomi (Technion, Israel)*.<br/>
The Reh (CRL-8286, ATCC, USA) human cell line that was isolated from tissue from an acute lymphocytic leukemia (ALL) patient was used as a wild type (WT) control cell line. Two types of mutations were introduced to the WT in FPGS gene by CRISPR/Cas9 system. One type is the FPGS gene point mutation (PM) at the position E115K identified in clinical samples diagnosed with relapsed leukemia, which is characterized by reduced FPGS activity. The other type is knockout (KO) of FPGS encoded gene. The cell line was cultured in RPMI-1640 with a physiological level of folic acid. For harvesting, after washing two times with ice-cold PBS, cell pellets were resuspended in an appropriate volume of methanol-acetonitrile-water (5:3:2 vol/vol/vol) extraction mixture. Cell lines were analyzed together with pooled QC samples on HILIC column in SCAN negative mode on Thermo Q-Exactive Orbitrap. Raw mzXML files were then processed in [`xcms`](https://bioconductor.org/packages/devel/bioc/html/xcms.html).<br/>
One specific group comparison (KO-6/WT) was subjected to functional analysis using the [`MetaboAnalyst`](https://www.metaboanalyst.ca/). 

Dataset was used to estimate the applicability of the MetaboCensoR App on Functional Analysis.<br/>
Code script and all relevant data are available in the folder. 
- Main script for analyzing data: [`Process folate.R`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/folate/Process%20folate.R)
- List of known target compounds: [`All_affected_cycles_NEG upd.csv`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/folate/All_affected_cycles_NEG%20upd.csv)
- Raw Data: [`MSV000100951`](...)
<table>
  <tr>
    <td>
      <img src="https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/folate/heatmap.png" height="400">
    </td>
    <td>
      <img src="https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/folate/fa.png" height="400">
    </td>
  </tr>
</table>
