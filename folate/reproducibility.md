## Reproducing the *MetaboCensoR* filtering for *orbi* dataset

The analysis can be reproduced using either interface:
- **Shiny app:** Upload the dataset and required auxiliary files, import the supplied [`JSON`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/folate/MetaboCensoR_settings_folate_dataset.json) settings, and apply the required filters.
- **R:** Run the supplied R script or R Markdown workflow, which specifies the inputs, parameters, and processing order:
```r
setwd(...)
library(MetaboCensoR)

# Input
# from local path
imported <- mc_prepare("xcms_table.csv",type="xcms",sample_keywords=c(".mzXML"))
# from GitHub
imported <- mc_prepare(read.csv("https://raw.githubusercontent.com/plyush1993/MetaboCensoR_Examples/main/folate/xcms_table.csv"),type="xcms",sample_keywords=c(".mzXML"))

imported$labels <- mc_labels_from_names(
  sample_names = imported$sample_cols,
  separator = "_",
  token = 3
)

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
  add_pol2 = "negative",
  # Neutral losses
  nl_rt2 = 0.002,
  nl_cor2 = 0.95,
  nl_tol_type2 = "da",
  nl_mz_tol2 = 0.005,
  nl_pol2 = "negative"
)

after_ms <- mc_filter_ms(imported,parameters = ms_options,filters=c("isotopes","adducts","neutral_losses"))
after_ms$ms$counts
print(after_ms) 

# QC filter
after_qc <- mc_filter_qc(after_ms,filters=c("rsd"),
                      groups=c("QC"),rsd_cutoff=30, mode = "group_any")
length(unique(with(after_qc$qc_tbl,Feature[metric == "rsd" & !keep])))
print(after_qc)
mc_plot_qc(after_qc,interactive=FALSE) 

# Output
output <- mc_export(
  after_qc,
  csv = "folate filtered package.csv"
)
```
