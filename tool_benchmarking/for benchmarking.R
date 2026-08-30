#................................................................
#### Libraries ----
#................................................................

library(tidyverse)
library(MetaboAnnotation)
library(MetaboCoreUtils)
library(ggplot2)
library(data.table)
library(mpactr)
library(mscleanr)

#................................................................
#### MetaboCensoR vs MS-CleanR (orbi dataset) ----
#................................................................

setwd("C:/")

target_df <- read.csv("annotation.csv")
target_df <- target_df[,c(1,2,4)]

mass <- calculateMass(target_df$Formula)
target_df$exactmass <- mass
target_df <- target_df %>%
  distinct()
target_df <- target_df %>%
  mutate(Compound = ifelse(
    duplicated(Compound) | duplicated(Compound, fromLast = TRUE), 
    paste0(Compound, "_RT", rt), 
    Compound
  ))

list_to_remove <- c("p-coumaraldehyde (Q27103652)|Phenylacrylic acid", "Wulignan A1", "6-Aminocaproic acid")
cat("p-coumaraldehyde (Q27103652)|Phenylacrylic acid, Wulignan A1, 6-Aminocaproic acid are detected in blank")
target_df <- subset(target_df, !target_df$Compound %in% list_to_remove)

parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+",
                                   "[M+C2H3N+Na]+", "[M+C2H3N+H]+", "[M+H-CH2O2]+"), 
                       tolerance = 0, ppm = 10, toleranceRt = 0.05)

# raw
data_ms2_only <- read.csv("data metabocensor area.csv") %>% 
  slice(-c(1:4))
colnames(data_ms2_only) <- make.unique(as.character(data_ms2_only[1,])) 
data_ms2_only <- data_ms2_only %>% slice(-1)
data_ms2_only <- subset(data_ms2_only, data_ms2_only$`MS/MS assigned` == "True")

data <- read.csv("data metabocensor area_standard_peak_table.csv") 
data <- subset(data, data$Feature %in% data_ms2_only$`Alignment ID`)
peakIn <- as.data.frame(cbind(mz = data$mz, rt = data$rt, id = data$Feature))
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_Compound)
unique(md$target_Compound) %>% length()

# MetaboCensoR
data_ms2_only <- read.csv("data metabocensor area_filtered (iso add +2H2O nl rsd rmd).csv") %>% 
  slice(-c(1:4))
colnames(data_ms2_only) <- make.unique(as.character(data_ms2_only[1,])) 
data_ms2_only <- data_ms2_only %>% slice(-1)
data_ms2_only <- subset(data_ms2_only, data_ms2_only$`MS/MS assigned` == "True") # keep only ms2

data <- read.csv("data metabocensor area_filtered (iso add +2H2O nl rsd rmd)_standard_peak_tab (1).csv") # blank filt iso add nl rsd rmd only ms2
data <- subset(data, data$Feature %in% data_ms2_only$`Alignment ID`)
peakIn <- as.data.frame(cbind(mz = data$mz, rt = data$rt, id = data$Feature))
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

matched_features <- matchValues(peakIn, target_df, parm)
md2 <- matchedData(matched_features)
md2 <- as.data.frame(md2)
md2 <- na.omit(md2)
unique(md2$target_Compound)
unique(md2$target_Compound) %>% length()

intersect(md$target_Compound, md2$target_Compound)
cat("Differences between raw and app:")
setdiff(md2$target_Compound, md$target_Compound)
setdiff(md$target_Compound, md2$target_Compound)

# MS-CleanR
data <- read.csv("MS_peaks-final_selection - mzrt.csv") # blank filt iso add nl rsd rmd only ms2
peakIn <- as.data.frame(cbind(mz = data$mz, rt = data$rt, id = data$Feature))
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

matched_features <- matchValues(peakIn, target_df, parm)
md3 <- matchedData(matched_features)
md3 <- as.data.frame(md3)
md3 <- na.omit(md3)
unique(md3$target_Compound)
unique(md3$target_Compound) %>% length()

intersect(md2$target_Compound, md3$target_Compound)
cat("Differences between 2 apps:")
setdiff(md2$target_Compound, md3$target_Compound)
setdiff(md3$target_Compound, md2$target_Compound)

#................................................................
#### MetaboCensoR vs MPACT/mpactR (orbi dataset) ----
#................................................................

setwd("C:/")

df <- read.csv("orbi_iimn_gnps_quant_standard_peak_table.csv")
colnames(df) <- str_remove(colnames(df), ".mzML")

# peak table
df <- df %>%
  select(
    Compound = id_number,
    "m/z" = mz,
    "Retention time (min)" = rtime,
    everything()[4:ncol(df)]
  ) 

df <- rbind(colnames(df), df) %>% 
  add_row(.before = 1) %>%
  add_row(.before = 1) 

write_csv(df, "peak table for mpactr.csv", col_names = F)

# metadata
df <- read.csv("orbi_iimn_gnps_quant_standard_peak_table.csv")
lab <- read.csv("Label.csv", header = F)

cndf <- colnames(df[,-c(1:3)])
cndf <- str_remove(cndf, ".mzML") %>% as.data.frame()
metadata <- data.frame(Injection = cndf[,1], Sample_Code=lab,	Biological_Group = lab) 
colnames(metadata) <- c("Injection", "Sample_Code", "Biological_Group")

write_csv(metadata, "metadata for mpactr.csv", col_names = T)

# run 
data <- import_data(peak_table ="peak table for mpactr.csv",
                    meta_data ="metadata for mpactr.csv", format = 'Progenesis')

data_filtered <- data |>
  filter_mispicked_ions(ringwin = 0.5,
                        isowin = 0.01,
                        trwin = 0.005,
                        max_iso_shift = 3,
                        merge_peaks = TRUE, 
                        merge_method = "sum") |>
  filter_group(group_to_remove = "blank", group_threshold = 0.1) |>
  filter_cv(cv_threshold = 0.3, cv_param = "mean") |> filter_insource_ions(cluster_threshold = 0.95)

# status
status <- get_raw_data(data_filtered) %>%
  mutate(Compound = as.character(Compound)) %>%
  select(Compound, mz, rt) %>%
  left_join(qc_summary(data_filtered),
            by = join_by("Compound" == "compounds")
  ) 

qc_sum <- qc_summary(data_filtered) 
plot_qc_tree(data_filtered)

# get filtered
ptf <- get_peak_table(data_filtered)
#write_csv(ptf, "mpactr mp gr cv isf.csv", col_names = T)

# annotation 
target_df <- read.csv("annotation.csv")
target_df <- target_df[,c(1,2,4)]

mass <- calculateMass(target_df$Formula)
target_df$exactmass <- mass
target_df <- target_df %>%
  distinct()
target_df <- target_df %>%
  mutate(Compound = ifelse(
    duplicated(Compound) | duplicated(Compound, fromLast = TRUE), 
    paste0(Compound, "_RT", rt), 
    Compound
  ))

df <- read_csv("orbi_iimn_gnps_quant.csv") %>% as.data.frame()
peakIn <- as.data.frame(cbind(mz = df$`row m/z`, rt = df$`row retention time`, id = df$`row ID`)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+", "[M+H-Hexose-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.01)

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_Compound)
unique(md$target_Compound) %>% length()
md <- subset(md, !md$target_Compound %in% c("Wulignan A1", "p-coumaraldehyde (Q27103652)|Phenylacrylic acid"))
cat("Wulignan A1 & p-coumaraldehyde (Q27103652)|Phenylacrylic acid are detected in blank")
unique(md$target_Compound) %>% length()

df <- read_csv("orbi_iimn_gnps_quant_filtered (2x H2O adducts isf rsd)_filtered.csv") %>% as.data.frame() # blank iso add isf rsd
peakIn <- as.data.frame(cbind(mz = df$`row m/z`, rt = df$`row retention time`, id = df$`row ID`)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

#parm <- Mass2MzParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+"), tolerance = 0, ppm = 5) 
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+", "[M+H-Hexose-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.01)

matched_features <- matchValues(peakIn, target_df, parm)
md2 <- matchedData(matched_features)
md2 <- as.data.frame(md2)
md2 <- na.omit(md2)
unique(md2$target_Compound)
unique(md2$target_Compound) %>% length()

# Compare
#unique(md2$target_compound) %>% as.data.frame() %>% View()
intersect(md$target_Compound, md2$target_Compound)
cat("Differences between raw and app:")
setdiff(md2$target_Compound, md$target_Compound)
setdiff(md$target_Compound, md2$target_Compound)

# mpact
mpact <- read_csv("mpactr mp gr cv isf.csv") 
peakIn <- as.data.frame(cbind(mz = mpact$mz, rt = mpact$rt, id = mpact$Compound)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+", "[M+H-Hexose-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.01)

matched_features <- matchValues(peakIn, target_df, parm)
md3 <- matchedData(matched_features)
md3 <- as.data.frame(md3)
md3 <- na.omit(md3)
md3 <- subset(md3, !md3$target_Compound %in% c("Wulignan A1", "p-coumaraldehyde (Q27103652)|Phenylacrylic acid"))
unique(md3$target_Compound)
unique(md3$target_Compound) %>% length()

# Compare
#unique(md2$target_compound) %>% as.data.frame() %>% View()
intersect(md$target_Compound, md3$target_Compound)
cat("Differences between raw and app:")
setdiff(md3$target_Compound, md$target_Compound)
setdiff(md$target_Compound, md3$target_Compound)

# Compare2
#unique(md2$target_compound) %>% as.data.frame() %>% View()
intersect(md2$target_Compound, md3$target_Compound)
cat("Differences between raw and app:")
setdiff(md3$target_Compound, md2$target_Compound)
setdiff(md2$target_Compound, md3$target_Compound)

#................................................................
#### MetaboCensoR vs khipu (orbi dataset) ----
#................................................................

setwd("C:/")

# convert
df <- read.csv("orbi_iimn_gnps_quant_FILT_BLANK_standard_peak_table.csv")
write_tsv(df, "orbi_iimn_gnps_quant_FILT_BLANK_standard_peak_table.tsv")

# annotate of khipu 
target_df <- read.csv("annotation.csv")
target_df <- target_df[,c(1,2,4)]

mass <- calculateMass(target_df$Formula)
target_df$exactmass <- mass
target_df <- target_df %>%
  distinct()
target_df <- target_df %>%
  mutate(Compound = ifelse(
    duplicated(Compound) | duplicated(Compound, fromLast = TRUE), 
    paste0(Compound, "_RT", rt), 
    Compound
  ))

df <- read_csv("orbi_iimn_gnps_quant_FILT_BLANK.csv") %>% as.data.frame()
peakIn <- as.data.frame(cbind(mz = df$`row m/z`, rt = df$`row retention time`, id = df$`row ID`)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+", "[M+H-Hexose-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.01)

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_Compound)
unique(md$target_Compound) %>% length()
md <- subset(md, !md$target_Compound %in% c("Wulignan A1", "p-coumaraldehyde (Q27103652)|Phenylacrylic acid"))
cat("Wulignan A1 & p-coumaraldehyde (Q27103652)|Phenylacrylic acid are detected in blank")
unique(md$target_Compound) %>% length()

df <- read_csv("orbi_iimn_gnps_quant_filtered (2x H2O adducts).csv") %>% as.data.frame()
peakIn <- as.data.frame(cbind(mz = df$`row m/z`, rt = df$`row retention time`, id = df$`row ID`)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

#parm <- Mass2MzParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+"), tolerance = 0, ppm = 5) 
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+", "[M+H-Hexose-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.01)

matched_features <- matchValues(peakIn, target_df, parm)
md2 <- matchedData(matched_features)
md2 <- as.data.frame(md2)
md2 <- na.omit(md2)
unique(md2$target_Compound)
unique(md2$target_Compound) %>% length()

# Compare
#unique(md2$target_compound) %>% as.data.frame() %>% View()
intersect(md$target_Compound, md2$target_Compound)
cat("Differences between raw and app:")
setdiff(md2$target_Compound, md$target_Compound)
setdiff(md$target_Compound, md2$target_Compound)

# khipu
khipu <- read_tsv("annotated_khipu_jIZhbjFfGAxa8XJjDufhScq1TsnlQNb9RVo_LwcPh_0.tsv") # 10 ppm 0.3 s 
peakIn <- as.data.frame(cbind(mz = khipu$mz, rt = khipu$rtime, id = rownames(khipu))) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+", "[M+H-Hexose-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.01)

matched_features <- matchValues(peakIn, target_df, parm)
md3 <- matchedData(matched_features)
md3 <- as.data.frame(md3)
md3 <- na.omit(md3)
md3 <- subset(md3, !md3$target_Compound %in% c("Wulignan A1", "p-coumaraldehyde (Q27103652)|Phenylacrylic acid"))
unique(md3$target_Compound)
unique(md3$target_Compound) %>% length()

# Compare
#unique(md2$target_compound) %>% as.data.frame() %>% View()
intersect(md$target_Compound, md3$target_Compound)
cat("Differences between raw and app:")
setdiff(md3$target_Compound, md$target_Compound)
setdiff(md$target_Compound, md3$target_Compound)

# Compare2
#unique(md2$target_compound) %>% as.data.frame() %>% View()
intersect(md2$target_Compound, md3$target_Compound)
cat("Differences between raw and app:")
setdiff(md3$target_Compound, md2$target_Compound)
setdiff(md2$target_Compound, md3$target_Compound)

#................................................................
#### MetaboCensoR vs Binner (orbi dataset) ----
#................................................................

setwd("C:/")

# convert
df <- read.csv("orbi_iimn_gnps_quant_FILT_BLANK_standard_peak_table.csv")
write.table(df, "orbi_iimn_gnps_quant_FILT_BLANK_standard_peak_table.txt", sep = "\t", row.names = FALSE)

# annotate of Binner
target_df <- read.csv("annotation.csv")
target_df <- target_df[,c(1,2,4)]

mass <- calculateMass(target_df$Formula)
target_df$exactmass <- mass
target_df <- target_df %>%
  distinct()
target_df <- target_df %>%
  mutate(Compound = ifelse(
    duplicated(Compound) | duplicated(Compound, fromLast = TRUE), 
    paste0(Compound, "_RT", rt), 
    Compound
  ))

df <- read_csv("orbi_iimn_gnps_quant_FILT_BLANK.csv") %>% as.data.frame()
peakIn <- as.data.frame(cbind(mz = df$`row m/z`, rt = df$`row retention time`, id = df$`row ID`)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+", "[M+H-Hexose-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.01)

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_Compound)
unique(md$target_Compound) %>% length()
md <- subset(md, !md$target_Compound %in% c("Wulignan A1", "p-coumaraldehyde (Q27103652)|Phenylacrylic acid"))
cat("Wulignan A1 & p-coumaraldehyde (Q27103652)|Phenylacrylic acid are detected in blank")
unique(md$target_Compound) %>% length()

df <- read_csv("orbi_iimn_gnps_quant_filtered (2x H2O adducts) rsd zeros_filtered.csv") %>% as.data.frame() # 30 % zeros, 400 RSD
peakIn <- as.data.frame(cbind(mz = df$`row m/z`, rt = df$`row retention time`, id = df$`row ID`)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

#parm <- Mass2MzParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+"), tolerance = 0, ppm = 5) 
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+", "[M+H-Hexose-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.01)

matched_features <- matchValues(peakIn, target_df, parm)
md2 <- matchedData(matched_features)
md2 <- as.data.frame(md2)
md2 <- na.omit(md2)
unique(md2$target_Compound)
unique(md2$target_Compound) %>% length()

# Compare
#unique(md2$target_compound) %>% as.data.frame() %>% View()
intersect(md$target_Compound, md2$target_Compound)
cat("Differences between raw and app:")
setdiff(md2$target_Compound, md$target_Compound)
setdiff(md$target_Compound, md2$target_Compound)

# binner
binner <- read_csv("Principal Ions + Unannotated.csv") # deiso mz 0.005 rt 0.01 grouping 0.005 annotation mz 0.005 rt 0.005
peakIn <- as.data.frame(cbind(mz = binner$`m/z`, rt = binner$RT, id = binner$Feature)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+", "[M+H-Hexose-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.01)

matched_features <- matchValues(peakIn, target_df, parm)
md3 <- matchedData(matched_features)
md3 <- as.data.frame(md3)
md3 <- na.omit(md3)
md3 <- subset(md3, !md3$target_Compound %in% c("Wulignan A1", "p-coumaraldehyde (Q27103652)|Phenylacrylic acid"))
unique(md3$target_Compound)
unique(md3$target_Compound) %>% length()

# Compare
#unique(md2$target_compound) %>% as.data.frame() %>% View()
intersect(md$target_Compound, md3$target_Compound)
cat("Differences between raw and app:")
setdiff(md3$target_Compound, md$target_Compound)
setdiff(md$target_Compound, md3$target_Compound)

# Compare2
#unique(md2$target_compound) %>% as.data.frame() %>% View()
intersect(md2$target_Compound, md3$target_Compound)
cat("Differences between raw and app:")
setdiff(md3$target_Compound, md2$target_Compound)
setdiff(md2$target_Compound, md3$target_Compound)

#............................................................................................................
