#................................................................
#### MetaboCensoR vs mzMine & MS1FA (orbi dataset) ----
#................................................................

library(tidyverse)
setwd("C:/...")

# Peak Table (mzMine)
data_mzmine <- read_csv("orbi_iimn_gnps_quant_FILT_BLANK.csv")
data_mzmine <- subset(data_mzmine, !is.na(data_mzmine$`best ion`))
nrow(data_mzmine)

# App data (mzMine)
add_app <- read_csv("orbi_iimn_gnps_quant adducts (rt for add 005).csv")
add_app_filt <- add_app[-which(add_app$adducts=="none"),]
unique(add_app_filt$peak_id) %>% length()

# MS1FA
ms1fa <- read_csv("Feature_table_output2026-03-10 MS1FA.csv")
which(!is.na(ms1fa$neutral_loss_annotation)) %>% length()

app <- read_csv("orbi_iimn_gnps_quant neutral loses.csv")
which(!is.na(app$loss_name)) %>% length()

#................................................................
#### MetaboCensoR vs nontarget (inter dataset) ----
#................................................................

library(tidyverse)
setwd("C:/...")

# App data
add_app <- read_csv("Area_5_2026_03_10_18_08_43 adducts.csv")
add_app_filt <- add_app[-which(add_app$adducts=="none"),]
nrow(add_app_filt)

# nontarget
library(nontarget)

data_msdial <- read.csv("Area_5_2026_03_10_18_08_43_FILT BLANK ISO.csv") %>% 
        slice(-c(1:3))
colnames(data_msdial) <- make.unique(as.character(data_msdial[1,])) 
data_msdial <- data_msdial %>% slice(-1)

peaklist <- cbind(mass = as.numeric(data_msdial$`Average Mz`), intensity = as.numeric(data_msdial$`1`), rt = as.numeric(data_msdial$`Average Rt(min)`)) %>% as.data.frame()
data(adducts) 
adducts_pos <- subset(adducts, adducts$Ion_mode == "positive")

adduct<-adduct.search(
  peaklist,
  adducts,
  rttol=0.005,
  mztol=0.005,
  ppm=F,
  use_adducts=adducts_pos$Name,
  ion_mode="positive"
)

which(adduct$adducts$`adduct(s)` != "none") %>% length()

#................................................................
#### MetaboCensoR vs CAMERA (folate dataset) ----
#................................................................

library(tidyverse)
setwd("C:/...")

# CAMERA
library(CAMERA)

ppm_diff <- 10 
mz_diff <- 0.005  
rt_diff <- 3

# load xcms object after peak filling
load("low_folate/xcms obj pk_fil NEG.RData") 
pk_fil <- as(pk_fil, "xcmsSet") 
pk_fil@filepaths <- normalizePath(list.files("Raw Data NEG Low Folate/", full.names = TRUE)[-1])

# perform
xsa <- xsAnnotate(pk_fil, polarity = "negative") 
xsaF <- groupFWHM(xsa, sigma = 6 , perfwhm = 0.6, intval = "maxo") 
xsaC <- groupCorr(xsaF, cor_eic_th=0.8, cor_exp_th = 0.8) 

# Annotate isotopes
xsaFI <- findIsotopes(xsaC, ppm = ppm_diff, mzabs=mz_diff, minfrac=0.1) 

# Annotate adducts
xsaFA <- findAdducts(xsaFI, polarity= "negative", ppm=ppm_diff, mzabs=mz_diff) # set polarity and params according to your data

# Get annotation info
annot_camera <- getPeaklist(xsaFA)

cam2 <- annot_camera %>%
  mutate(
    pcgroup = as.character(pcgroup),
    isotopes_chr = str_trim(as.character(isotopes)),
    adduct_chr   = str_trim(as.character(adduct)),
    has_iso = !is.na(isotopes_chr) & isotopes_chr != "",
    has_add = !is.na(adduct_chr)   & adduct_chr   != ""
  )

# Overall counts
overall_counts <- cam2 %>%
  dplyr::summarise(
    n_features = n(),
    n_isotope_features = sum(has_iso),
    n_adduct_features  = sum(has_add),
    n_pcgroups = n_distinct(pcgroup),
    n_pcgroups_with_isotopes = n_distinct(pcgroup[has_iso]),
    n_pcgroups_with_adducts  = n_distinct(pcgroup[has_add])
  )

overall_counts

# App results
iso_dim_app <- read_csv("xcms_table isotopes-dimers.csv")
which(!is.na(iso_dim_app$iso_group)) %>% length()

add_app <- read_csv("xcms_table adducts.csv")
which(add_app$adducts != "none") %>% length()
#.................................................................