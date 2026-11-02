## Reproducing the *MetaboCensoR* filtering for *orbi* dataset

The analysis can be reproduced using either interface:
- **Shiny app:** Upload the dataset and required auxiliary files, import the supplied [`JSON`](https://github.com/plyush1993/MetaboCensoR_Examples/blob/main/inter/MetaboCensoR_settings_inter_dataset.json) settings, and apply the required filters.
- **R:** Run the supplied R script or R Markdown workflow, which specifies the inputs, parameters, and processing order:
```r
setwd(...)
library(MetaboCensoR)

# Input
imported <- mc_prepare("Area_5_2026_03_10_18_08_43.csv",type="msdial",sample_keywords=c(".mzXML"),sample_cols = 36:63)

imported$labels <- mc_labels_from_names(
  sample_names = imported$sample_cols,
  separator = "_",
  token = 2
)

# Blank filter
after_blank <- mc_filter_blank(imported,groups="Peptone",cutoff=.10, mode = "mean_cutoff")
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
  # Neutral losses
  nl_rt2 = 0.002,
  nl_cor2 = 0.95,
  nl_tol_type2 = "da",
  nl_mz_tol2 = 0.005,
  # In-source fragments
  isf_rt2 = 0.002,
  isf_cor2 = 0.8,
  isf_ms2_check2 = FALSE 
)

after_ms <- mc_filter_ms(after_blank,parameters = ms_options,filters=c("isotopes","adducts","neutral_losses", "isf"))
after_ms$ms$counts
print(after_ms) 

# QC filter
after_qc <- mc_filter_qc(after_ms,filters=c("rsd"),
                      groups=c("001srf","Bs","dd","dPPS","dSRF","Pd", "Peptone"),rsd_cutoff=15, mode = "group_any")
length(unique(with(after_qc$qc_tbl,Feature[metric == "rsd" & !keep])))
print(after_qc)

# Peak filter
after_peak <- mc_filter_peak(after_qc,filters=c("rt"),
                       parameters=mc_peak_parameters(peak_rt_min=1))
length(with(after_peak$peak_tbl, Feature[!pass_rt]))
print(after_peak)

# Output
output <- mc_export(
  after_peak,
  csv = "inter filtered package.csv"
)
```
