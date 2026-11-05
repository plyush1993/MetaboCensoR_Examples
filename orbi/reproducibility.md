## Reproducing the *MetaboCensoR* filtering for *orbi* dataset

The analysis can be reproduced using either interface:
- **Shiny app:** Upload the dataset and required auxiliary files, import the supplied [`JSON`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/orbi/MetaboCensoR_settings_orbi_dataset.json) settings, and apply the required filters.
- **R:** Run the supplied R script or R Markdown workflow, which specifies the inputs, parameters, and processing order:
```r
setwd(...)
library(MetaboCensoR)

# Input
# from local path
imported <- mc_prepare("orbi_iimn_gnps_quant.csv",type="mzmine",sample_keywords=c(".mzML","Peak area"))
spectra <- mc_read_mgf("orbi_iimn_gnps.mgf")

imported$labels <- mc_labels_from_metadata(
  path = "metadata.csv",
  sample_names = imported$sample_cols,
  sample_col = "Sample name",
  label_col = "Label"
)

# from GitHub
imported <- mc_prepare(read.csv("https://raw.githubusercontent.com/plyush1993/MetaboCensoR_Examples/main/orbi/orbi_iimn_gnps_quant.csv"),type="mzmine",sample_keywords=c(".mzML","Peak area"))
mgf_file <- tempfile(fileext = ".mgf")
download.file(
  "https://raw.githubusercontent.com/plyush1993/MetaboCensoR_Examples/main/orbi/orbi_iimn_gnps.mgf",
  destfile = mgf_file, mode = "wb")
spectra <- mc_read_mgf(mgf_file)

imported$labels <- mc_labels_from_csv(
  path = "https://raw.githubusercontent.com/plyush1993/MetaboCensoR_Examples/main/orbi/Label.csv",
  sample_names = imported$sample_cols
)

# Blank filter
after_blank <- mc_filter_blank(imported,groups="Blank",cutoff=.10, mode = "mean_cutoff")
print(after_blank)
mc_plot_blank(after_blank,interactive=FALSE) 

# MS filter
# ms params
ms_options <- mc_ms_parameters(
  # Isotopes
  iso_rt2 = 0.01,          
  iso_cor2 = 0.80,        
  iso_tol_type2 = "da",
  iso_mz_tol2 = 0.005,
  # Adducts
  add_rt2 = 0.005,
  add_cor2 = 0.80,
  add_tol_type2 = "da",
  add_mz_tol2 = 0.005,
  # In-source fragments
  isf_rt2 = 0.01,
  isf_cor2 = 0.95,
  isf_ms2_check2 = TRUE,  
  isf_mgf_id_field = "acquisitionNum",
  isf_ms2_mz_tol2 = 0.01,        
  isf_ms2_rel2 = 1
)

after_ms <- mc_filter_ms(after_blank,filters=c("isotopes","adducts","isf"), parameters = ms_options, spectra = spectra, adducts = "adducts (nontarget + 2xH2O).csv")
after_ms$ms$counts 
print(after_ms) 

# Output
output <- mc_export(
  after_ms,
  csv = "orbi filtered package.csv",
  mgf = spectra,
  mgf_file = "orbi filtered package.mgf"
)
```
