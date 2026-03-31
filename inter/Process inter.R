#................................................................
# Upload -----
#................................................................

library(tidyverse)
library(MetaboAnnotation)
library(MetaboCoreUtils)
library(ggplot2)
library(data.table)
library(ggsci)
library(batchCorr)
library(cowplot)

setwd("C:/.../")

# raw data
data <- read.csv("Area_5_2026_03_10_18_08_43.csv") %>% 
        slice(-c(1:3))
colnames(data) <- make.unique(as.character(data[1,])) 
data <- data %>% slice(-1)
mzrt_v <- paste0(data$`Average Mz`, "@", data$`Average Rt(min)`)
rownames(data) <- mzrt_v
df <- data[,which(str_detect(colnames(data), "_"))] %>% t() %>% as.data.frame()

metadata <- do.call(rbind, str_split(rownames(df), "_")) %>% as.data.frame()
metadata$Label <- paste0(metadata[,1], "+", metadata[,2])
metadata$Label <- str_replace(metadata$Label, "Bacto\\+Peptone", "Media")
ds <- cbind(Label=metadata$Label, df) %>% as.data.frame()
ds[,-1] <- sapply(ds[,-1], as.numeric)
raw <- ds

# processed via App
data <- read.csv("Area_5_2026_03_10_18_08_43_filtered.csv") %>% 
        slice(-c(1:3))
colnames(data) <- make.unique(as.character(data[1,])) 
data <- data %>% slice(-1)
mzrt_v <- paste0(data$`Average Mz`, "@", data$`Average Rt(min)`)
rownames(data) <- mzrt_v
df <- data[,which(str_detect(colnames(data), "_"))] %>% t() %>% as.data.frame()

metadata <- do.call(rbind, str_split(rownames(df), "_")) %>% as.data.frame()
metadata$Label <- paste0(metadata[,1], "+", metadata[,2])
metadata$Label <- str_replace(metadata$Label, "Bacto\\+Peptone", "Media")
ds <- cbind(Label=metadata$Label, df) %>% as.data.frame()
ds[,-1] <- sapply(ds[,-1], as.numeric)
app <- ds

#................................................................
# Annotation Raw ----
#................................................................
target_df <- read.csv("annot table.csv")

mass <- calculateMass(target_df$Formula)
target_df$exactmass <- mass

df <- raw # dataset
cn <- colnames(df[,-1])
cn0 <- colnames(df[,-1])
cn2 <- t(data.frame(cn))
colnames(cn2) <- cn
peakIn <- peakInfo(PT = cn2, sep = "@", start = 1)
peakIn <- as.data.frame(cbind(peakIn, id = cn0)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

# set polarity, adduct, accuracy
parm <- Mass2MzParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+"), tolerance = 0, ppm = 5) 
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+"), tolerance = 0, ppm = 5, toleranceRt = 0.15) 
#MetaboCoreUtils::adducts() 

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_Compound)

#................................................................
# Volcano Raw ----
#................................................................
ds <- raw
df <- ds # dataset
ds <- subset(ds, Label != "Media")
df <- ds 

# MVI
ds_mvi <- df[,-1]
ds_mvi[ds_mvi == 0] <- NA

#noise <- as.numeric(quantile(1:min(ds_mvi, na.rm = T))[2]) # reduce noise by quantile: quantile(1:noise)[1] = 0% ... quantile(1:noise)[4] = 100%
noise <- 31.5 # from RAW data
sd <- as.numeric(30) # set sd value for random value generation as.numeric(quantile(1:noise)[2]) or noise*0.3
set.seed(1234)
ds_mvi[is.na(ds_mvi)] <- 0 # convert NA into 0
NAidx <- ds_mvi == 0 # NA index
imp.rand <- abs(rnorm(sum(NAidx), mean = noise, sd = sd)) # generate random values. other option: runif(sum(NAidx), noise-sd, noise+sd)
ds_mvi[NAidx] <- imp.rand
df <- as.data.frame(cbind(Label = df[,1], ds_mvi))

# Stat tests
gr <- as.factor(as.character(df$Label))
g <- t(combn(as.character(levels(gr)), 2))
g <- cbind(g[,2], g[,1])
g_sel <- g[which(g[,1] == "Pd+Pd"),]
#g_sel <- rbind(g_sel, c("Pd+Pd", "Pd+srf"))
# add manually 
#g_sel <- rbind(g_sel, c("Pd+Pd", "Pd+srf"))
#g_sel <- g
#idx <- which(rowSums(g == "Media") > 0 ) 
#g_sel <- g_sel[-idx, ] %>% t() %>% as.matrix()

df_l <- list()
names_comp <- c()

for (i in 1:nrow(g_sel)) {
  df_l[[i]] <- as.data.frame(subset(df, df$Label %in% g_sel[i,]))
  df_l[[i]][,-1] <- sapply(df_l[[i]][,-1], as.numeric)
  df_l[[i]][,1] <- as.factor(as.character(df_l[[i]][,1]))
  names_comp[i] <- paste0(g_sel[i,1]," / ",g_sel[i,2])
}

names(df_l) <- names_comp

stat_test_l <- lapply(1:length(df_l), function(x) sapply(2:ncol(df_l[[x]]), function(y) t.test(df_l[[x]][,y] ~ df_l[[x]][,1])$p.value)) # t.test # wilcox.test
stat_test_l <- lapply(1:length(stat_test_l), function(x) p.adjust(stat_test_l[[x]], "BH"))

stat_test_df <- do.call(rbind, stat_test_l)
colnames(stat_test_df) <- colnames(df_l[[1]][,-1])
stat_test_df <- as.data.frame(cbind(Groups = names(df_l), stat_test_df))
stat_test_df[,-1] <- sapply(stat_test_df[,-1], as.numeric)
stat_test_df <- data.frame(sapply(stat_test_df, function(x) ifelse(is.nan(x), 1, x)))

stat_test_df_t <- as.data.frame(t(stat_test_df))
colnames(stat_test_df_t) <- stat_test_df_t[1,]
stat_test_df_t <- as.data.frame(cbind(Compound = colnames(df), stat_test_df_t))
stat_test_df_t <- stat_test_df_t[-1,]

stat_test_df_v <- reshape2::melt(stat_test_df_t, value.name = "Adj.p-value", variable.name = "Groups", id=1)

# FC
FOLD.CHANGE.MG <- function(x, f, j, aggr_FUN = colMeans, combi_FUN = {function(x,y) "-"(x,y)})({
  f <- as.factor(f)
  i <- split(1:nrow(x), f)
  x <- sapply(i, function(i)({ aggr_FUN(x[i,])}))
  x <- t(x)
  x <- log2(x)
  j <- j
  ret <- combi_FUN(x[j[1,],], x[j[2,],])
  rownames(ret) <- paste(j[1,], j[2,], sep = ' / ')
  t(ret)
})

fdr <- FOLD.CHANGE.MG(x = df[,-1], f = df[,1], j = t(g_sel)) %>% as.data.frame()
fdr0 <- fdr

fdr <- as.data.frame(sapply(fdr, function(x) ifelse(is.nan(x), 0, x)))
fdr <- as.data.frame(sapply(fdr, function(x) ifelse(is.infinite(x), 0, x)))
rownames(fdr) <- rownames(fdr0)

fdr <- as.data.frame(cbind(Compound = colnames(df[,-1]), fdr))
fdr_v <- reshape2::melt(fdr, value.name = "FC", variable.name = "Groups", id=1)

# Combine for VP
data_stat <- full_join(stat_test_df_v, fdr_v, by = c("Groups", "Compound"))
data_stat[,3:4] <- sapply(data_stat[,3:4], as.numeric)
data_stat$`Adj.p-value.log` <- -log10(data_stat$`Adj.p-value`)

#to_del_med <- which(str_detect(unique(as.character(data_stat$Groups)), "Media"))
#data_stat <- subset(data_stat, data_stat$Groups %in% c(unique(as.character(data_stat$Groups))[-to_del_med]))
data_stat$Target <- ifelse(data_stat$Compound %in% md$id, "Target", "No")
data_stat$Name <- ifelse(data_stat$Compound %in% md$id, md$target_Compound[match(data_stat$Compound, md$id)], "No")
data_stat$TargetNew <- ifelse(data_stat$Compound == "374.22937@4.870", "Target_New", "No")
md_addH <- subset(md, adduct %in% c("[M+H]+", "[M+2H]2+"))
data_stat$Adduct <- ifelse(data_stat$Compound %in% md_addH$id, "M+H", md$adduct[match(data_stat$Compound, md$id)])
data_stat_raw <- data_stat

p <- ggplot(data_stat, aes(x = FC, y = `Adj.p-value.log`)) +

  geom_point(
    data = dplyr::filter(data_stat, Target == "Target"),
    shape = 21, fill = NA, colour = "yellow",
    size = 5, stroke = 3
  ) +

  geom_point(
    data = dplyr::filter(data_stat, Adduct == "M+H"),
    shape = 21, fill = NA, colour = "green",
    size = 5, stroke = 3
  ) +
  
    geom_point(
    data = dplyr::filter(data_stat, TargetNew == "Target_New"),
    shape = 21, fill = NA, colour = "blue",
    size = 5, stroke = 3
  ) +
  
  geom_point(aes(fill = Groups),
             shape = 21, colour = "black",
             alpha = 0.85, size = 5, stroke = 0.4) +

  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_fill_jama() +
  labs(x = "log2(FC)", y = "-log10(FDR)") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) 

p <- p+facet_wrap(vars(Groups), scales = "fixed", nrow = 1, strip.position = "top", labeller = label_wrap_gen(width = 15))
p
p_raw <- p

#................................................................
# Annotation App ----
#................................................................
target_df <- read.csv("annot table.csv")

mass <- calculateMass(target_df$Formula)
target_df$exactmass <- mass

df <- app # dataset
cn <- colnames(df[,-1])
cn0 <- colnames(df[,-1])
cn2 <- t(data.frame(cn))
colnames(cn2) <- cn
peakIn <- peakInfo(PT = cn2, sep = "@", start = 1)
peakIn <- as.data.frame(cbind(peakIn, id = cn0)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

# set polarity, adduct, accuracy
parm <- Mass2MzParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+"), tolerance = 0, ppm = 5) 
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+"), tolerance = 0, ppm = 5, toleranceRt = 0.15) 
#MetaboCoreUtils::adducts() 

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_Compound)

#................................................................
# Volcano App ----
#................................................................
ds <- app
df <- ds # dataset
ds <- subset(ds, Label != "Media")
df <- ds 

# MVI
ds_mvi <- df[,-1]
ds_mvi[ds_mvi == 0] <- NA

#noise <- as.numeric(quantile(1:min(ds_mvi, na.rm = T))[2]) # reduce noise by quantile: quantile(1:noise)[1] = 0% ... quantile(1:noise)[4] = 100%
noise <- 31.5 # from RAW data
sd <- as.numeric(30) # set sd value for random value generation as.numeric(quantile(1:noise)[2]) or noise*0.3
set.seed(1234)
ds_mvi[is.na(ds_mvi)] <- 0 # convert NA into 0
NAidx <- ds_mvi == 0 # NA index
imp.rand <- abs(rnorm(sum(NAidx), mean = noise, sd = sd)) # generate random values. other option: runif(sum(NAidx), noise-sd, noise+sd)
ds_mvi[NAidx] <- imp.rand
df <- as.data.frame(cbind(Label = df[,1], ds_mvi))

# Stat tests
gr <- as.factor(as.character(df$Label))
g <- t(combn(as.character(levels(gr)), 2))
g <- cbind(g[,2], g[,1])
g_sel <- g[which(g[,1] == "Pd+Pd"),]
#g_sel <- rbind(g_sel, c("Pd+Pd", "Pd+srf"))
# add manually 
#g_sel <- rbind(g_sel, c("Pd+Pd", "Pd+srf"))
#g_sel <- g
#idx <- which(rowSums(g == "Media") > 0 ) 
#g_sel <- g_sel[-idx, ] %>% t() %>% as.matrix()

df_l <- list()
names_comp <- c()

for (i in 1:nrow(g_sel)) {
  df_l[[i]] <- as.data.frame(subset(df, df$Label %in% g_sel[i,]))
  df_l[[i]][,-1] <- sapply(df_l[[i]][,-1], as.numeric)
  df_l[[i]][,1] <- as.factor(as.character(df_l[[i]][,1]))
  names_comp[i] <- paste0(g_sel[i,1]," / ",g_sel[i,2])
}

names(df_l) <- names_comp

stat_test_l <- lapply(1:length(df_l), function(x) sapply(2:ncol(df_l[[x]]), function(y) t.test(df_l[[x]][,y] ~ df_l[[x]][,1])$p.value)) # t.test # wilcox.test
stat_test_l <- lapply(1:length(stat_test_l), function(x) p.adjust(stat_test_l[[x]], "BH"))

stat_test_df <- do.call(rbind, stat_test_l)
colnames(stat_test_df) <- colnames(df_l[[1]][,-1])
stat_test_df <- as.data.frame(cbind(Groups = names(df_l), stat_test_df))
stat_test_df[,-1] <- sapply(stat_test_df[,-1], as.numeric)
stat_test_df <- data.frame(sapply(stat_test_df, function(x) ifelse(is.nan(x), 1, x)))

stat_test_df_t <- as.data.frame(t(stat_test_df))
colnames(stat_test_df_t) <- stat_test_df_t[1,]
stat_test_df_t <- as.data.frame(cbind(Compound = colnames(df), stat_test_df_t))
stat_test_df_t <- stat_test_df_t[-1,]

stat_test_df_v <- reshape2::melt(stat_test_df_t, value.name = "Adj.p-value", variable.name = "Groups", id=1)

# FC
FOLD.CHANGE.MG <- function(x, f, j, aggr_FUN = colMeans, combi_FUN = {function(x,y) "-"(x,y)})({
  f <- as.factor(f)
  i <- split(1:nrow(x), f)
  x <- sapply(i, function(i)({ aggr_FUN(x[i,])}))
  x <- t(x)
  x <- log2(x)
  j <- j
  ret <- combi_FUN(x[j[1,],], x[j[2,],])
  rownames(ret) <- paste(j[1,], j[2,], sep = ' / ')
  t(ret)
})

fdr <- FOLD.CHANGE.MG(x = df[,-1], f = df[,1], j = t(g_sel)) %>% as.data.frame()
fdr0 <- fdr

fdr <- as.data.frame(sapply(fdr, function(x) ifelse(is.nan(x), 0, x)))
fdr <- as.data.frame(sapply(fdr, function(x) ifelse(is.infinite(x), 0, x)))
rownames(fdr) <- rownames(fdr0)

fdr <- as.data.frame(cbind(Compound = colnames(df[,-1]), fdr))
fdr_v <- reshape2::melt(fdr, value.name = "FC", variable.name = "Groups", id=1)

# Combine for VP
data_stat <- full_join(stat_test_df_v, fdr_v, by = c("Groups", "Compound"))
data_stat[,3:4] <- sapply(data_stat[,3:4], as.numeric)
data_stat$`Adj.p-value.log` <- -log10(data_stat$`Adj.p-value`)

#to_del_med <- which(str_detect(unique(as.character(data_stat$Groups)), "Media"))
#data_stat <- subset(data_stat, data_stat$Groups %in% c(unique(as.character(data_stat$Groups))[-to_del_med]))
data_stat$Target <- ifelse(data_stat$Compound %in% md$id, "Target", "No")
data_stat$Name <- ifelse(data_stat$Compound %in% md$id, md$target_Compound[match(data_stat$Compound, md$id)], "No")
data_stat$TargetNew <- ifelse(data_stat$Compound == "374.22937@4.87", "Target_New", "No")
md_addH <- subset(md, adduct %in% c("[M+H]+", "[M+2H]2+"))
data_stat$Adduct <- ifelse(data_stat$Compound %in% md_addH$id, "M+H", md$adduct[match(data_stat$Compound, md$id)])
data_stat_app <- data_stat

p <- ggplot(data_stat, aes(x = FC, y = `Adj.p-value.log`)) +

  geom_point(
    data = dplyr::filter(data_stat, Target == "Target"),
    shape = 21, fill = NA, colour = "yellow",
    size = 5, stroke = 3
  ) +

  geom_point(
    data = dplyr::filter(data_stat, Adduct == "M+H"),
    shape = 21, fill = NA, colour = "green",
    size = 5, stroke = 3
  ) +
  
    geom_point(
    data = dplyr::filter(data_stat, TargetNew == "Target_New"),
    shape = 21, fill = NA, colour = "blue",
    size = 5, stroke = 3
  ) +
  
  geom_point(aes(fill = Groups),
             shape = 21, colour = "black",
             alpha = 0.85, size = 5, stroke = 0.4) +

  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_fill_jama() +
  labs(x = "log2(FC)", y = "-log10(FDR)") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) 

p <- p+facet_wrap(vars(Groups), scales = "fixed", nrow = 1, strip.position = "top", labeller = label_wrap_gen(width = 15))
p
p_app <- p
#ggplotly(p)
#................................................................
# Cowplot ----
#................................................................

library(cowplot)
p_raw <- p_raw + coord_cartesian(xlim = c(-20,20), ylim = c(0,4.2))
p_app <- p_app + coord_cartesian(xlim = c(-20,20), ylim = c(0,4.2))
plot_grid(p_raw, p_app, labels = c('A', 'B'), label_size = 25, nrow = 2)

#................................................................
# lollipop plot ----
#................................................................

data_stat_lpp <- rbind(cbind(data_stat_app, Label = "App"), cbind(data_stat_raw, Label = "Raw")) %>% as.data.frame()
plot_data <- subset(data_stat_lpp, data_stat_lpp$Name != "No")
plot_data$Groups %>% unique()
plot_data <- subset(plot_data, plot_data$Groups %in% c("Pd+Pd / Pd+001srf", "Pd+Pd / Pd+Bs", "Pd+Pd / Pd+dPPS"))

plot_data <- plot_data %>%
  mutate(Y_Label = paste0(Name, " ", Adduct, "")) %>%
  
  mutate(Label = factor(Label, levels = c("App", "Raw")))

# 3. Generate the plot using your exact column names
ggplot(plot_data, aes(x = `Adj.p-value.log`, y = Y_Label, color = Label)) +
  
  # Add a vertical threshold line (FDR = 0.05 is ~ 1.3 on a log scale)
  geom_vline(xintercept = 1.301, linetype = "dashed", color = "gray50", size = 1) +
  
  # Draw the sticks, dodged side-by-side
  geom_linerange(aes(xmin = 0, xmax = `Adj.p-value.log`), 
                 position = position_dodge(width = 0.6), 
                 linewidth = 1, alpha = 0.8) +
  
  # Draw the dots at the end of the sticks
  geom_point(aes(x = `Adj.p-value.log`), 
             position = position_dodge(width = 0.6), 
             size = 3) +
  
  # Separate the plot by your 'Groups' column
  # scales = "free_y" ensures each facet only shows compounds present in that specific comparison
  facet_wrap(~ Groups, scales = "free_y", ncol = 1) +
  
  # Set specific colors for Raw and App
  scale_color_npg() +
  
  #scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    breaks = c(0, 1.301, 2, 3, 4, 5), 
    labels = c("1.0", "0.05", "0.01", "0.001", "0.0001", "0.00001")
  )+
  
  labs(
    title = "",
    x = "FDR",
    y = "",
    color = ""
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(margin = margin(t = 0)),
        axis.text.y = element_text(margin = margin(r = 0)),
        #axis.ticks.y = element_line(color = "gray50", linewidth = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold", size = 11)
  )

#................................................................
#### Compare Total Annotations ----
#................................................................
# raw data
data <- read.csv("Area_5_2026_03_10_18_08_43.csv") %>% 
        slice(-c(1:3))
colnames(data) <- make.unique(as.character(data[1,])) 
data <- data %>% slice(-1)
mzrt_v <- paste0(data$`Average Mz`, "@", data$`Average Rt(min)`)
rownames(data) <- mzrt_v
df <- data[,which(str_detect(colnames(data), "_"))] %>% t() %>% as.data.frame()

metadata <- do.call(rbind, str_split(rownames(df), "_")) %>% as.data.frame()
metadata$Label <- paste0(metadata[,1], "+", metadata[,2])
metadata$Label <- str_replace(metadata$Label, "Bacto\\+Peptone", "Media")
ds <- cbind(Label=metadata$Label, df) %>% as.data.frame()
ds[,-1] <- sapply(ds[,-1], as.numeric)
raw <- ds

target_df <- read.csv("annot table.csv")

mass <- calculateMass(target_df$Formula)
target_df$exactmass <- mass

df <- raw # dataset
cn <- colnames(df[,-1])
cn0 <- colnames(df[,-1])
cn2 <- t(data.frame(cn))
colnames(cn2) <- cn
peakIn <- peakInfo(PT = cn2, sep = "@", start = 1)
peakIn <- as.data.frame(cbind(peakIn, id = cn0)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

# set polarity, adduct, accuracy
parm <- Mass2MzParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+"), tolerance = 0, ppm = 5) 
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+"), tolerance = 0, ppm = 5, toleranceRt = 0.15) 
#MetaboCoreUtils::adducts() 

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_Compound)

# processed via App
data <- read.csv("Area_5_2026_03_10_18_08_43_filtered.csv") %>% 
        slice(-c(1:3))
colnames(data) <- make.unique(as.character(data[1,])) 
data <- data %>% slice(-1)
mzrt_v <- paste0(data$`Average Mz`, "@", data$`Average Rt(min)`)
rownames(data) <- mzrt_v
df <- data[,which(str_detect(colnames(data), "_"))] %>% t() %>% as.data.frame()

metadata <- do.call(rbind, str_split(rownames(df), "_")) %>% as.data.frame()
metadata$Label <- paste0(metadata[,1], "+", metadata[,2])
metadata$Label <- str_replace(metadata$Label, "Bacto\\+Peptone", "Media")
ds <- cbind(Label=metadata$Label, df) %>% as.data.frame()
ds[,-1] <- sapply(ds[,-1], as.numeric)
app <- ds

target_df <- read.csv("annot table.csv")

mass <- calculateMass(target_df$Formula)
target_df$exactmass <- mass

df <- app # dataset
cn <- colnames(df[,-1])
cn0 <- colnames(df[,-1])
cn2 <- t(data.frame(cn))
colnames(cn2) <- cn
peakIn <- peakInfo(PT = cn2, sep = "@", start = 1)
peakIn <- as.data.frame(cbind(peakIn, id = cn0)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

# set polarity, adduct, accuracy
parm <- Mass2MzParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+"), tolerance = 0, ppm = 5) 
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+"), tolerance = 0, ppm = 5, toleranceRt = 0.15) 
#MetaboCoreUtils::adducts() 

matched_features <- matchValues(peakIn, target_df, parm)
md2 <- matchedData(matched_features)
md2 <- as.data.frame(md2)
md2 <- na.omit(md2)
unique(md2$target_Compound)

# Compare
#unique(md2$target_compound) %>% as.data.frame() %>% View()
intersect(md$target_Compound, md2$target_Compound)
setdiff(md2$target_Compound, md$target_Compound)

setdiff(target_df$Compound, md$target_Compound)
