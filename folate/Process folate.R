#................................................................
#### Settings ----
#................................................................

library(tidyverse)
library(dplyr)
library(tidyr)
library(MetaboAnnotation)
library(forcats)
library(MetaboCoreUtils)
library(ggplot2)
library(data.table)
library(ggsci)
library(cowplot)

setwd("C:/.../")

#................................................................
#### Annotation App ----
#................................................................

df <- read_csv("xcms_table_filtered.csv") %>% as.data.frame()
peakIn <- as.data.frame(cbind(mz = df$mzmed, rt = df$rtmed, id = df$feature)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- round(as.numeric(peakIn$rt)/60, 4)

target_df <- read.csv("All_affected_cycles_NEG upd.csv")

mass <- calculateMass(target_df$formula)
target_df$exactmass <- mass
target_df <- target_df[,c(1,2,9,4)]

# set polarity, adduct, accuracy
parm <- Mass2MzParam(adducts = c("[M-H]-", "[M-2H]2-", "[M+Na-2H]-","[M+Cl]-","[M+K-2H]-","[M+C2H3N-H]-","[M+CHO2]-", "[M+C2H3O2]-","[M+Br]-", "[2M-H]-", "[M]-"), tolerance = 0, ppm = 5) 
parm <- Mass2MzRtParam(adducts = c("[M-H]-", "[M-2H]2-", "[M+Na-2H]-","[M+Cl]-","[M+K-2H]-","[M+C2H3N-H]-","[M+CHO2]-", "[M+C2H3O2]-","[M+Br]-", "[2M-H]-", "[M]-"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.12)
#MetaboCoreUtils::adducts(polarity = c("positive", "negative")) 
#MetaboCoreUtils::adducts(polarity = c("negative")) 

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_compound)
unique(md$target_compound) %>% length()

# Plot
library(ggsci)

adduct_summary <- tibble(adduct = md$adduct) %>%
  dplyr::count(adduct, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  arrange(desc(n))

total_n <- sum(adduct_summary$n)

compound_summary <- tibble(target_Compound = md$target_compound) %>%
  dplyr::count(target_Compound, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  arrange(desc(n))

ggplot(adduct_summary, aes(x = reorder(adduct, n), y = n)) +
  geom_col(color = "black", aes(fill = adduct), size = 1, width = 0.5) +
  #coord_flip() +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, percent)), hjust = 0.5, vjust = -0.5, size = 5.5) +
  labs(x = NULL, y = "Count")+
  #labs(x = NULL, y = "Count", title = paste("Adduct counts (Total =", total_n, ")")) +
  expand_limits(y = max(adduct_summary$n) * 1.15) + theme_classic(base_size = 16) + theme(legend.position = "none")+
  scale_fill_npg() 

adduct_summary_app <- adduct_summary

# plot occurance
occ_tbl <- md %>%
  mutate(target_compound = trimws(as.character(target_compound))) %>%
  filter(!is.na(target_compound), target_compound != "") %>%
  dplyr::count(target_compound, name = "occurrence") %>%
  arrange(desc(occurrence), target_compound)

occ_dist <- occ_tbl %>%
  dplyr::count(occurrence, name = "n_compounds") %>%
  arrange(occurrence) %>%
  mutate(
    occurrence_f = factor(occurrence, levels = occurrence),
    percent = 100 * n_compounds / sum(n_compounds),
    lab = sprintf("%d (%.1f%%)", n_compounds, percent)
  )

ggplot(occ_dist, aes(x = occurrence_f, y = n_compounds)) +
  geom_col(color = "black", aes(fill = occurrence_f), size = 1, width = 0.7) +
  geom_text(aes(label = lab), vjust = -0.3, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = "Adducts per compound",
    y = "Number of compounds"
  ) +
  theme_classic(base_size = 20) + theme(legend.position = "none") + scale_fill_aaas()

occ_dist_app <- occ_dist

#................................................................
#### Annotation Raw ----
#................................................................

df <- read_csv("xcms_table.csv") %>% as.data.frame()
peakIn <- as.data.frame(cbind(mz = df$mzmed, rt = df$rtmed, id = df$feature)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- round(as.numeric(peakIn$rt)/60, 4)

target_df <- read.csv("All_affected_cycles_NEG upd.csv")

mass <- calculateMass(target_df$formula)
target_df$exactmass <- mass
target_df <- target_df[,c(1,2,9,4)]

# set polarity, adduct, accuracy
parm <- Mass2MzParam(adducts = c("[M-H]-", "[M-2H]2-", "[M+Na-2H]-","[M+Cl]-","[M+K-2H]-","[M+C2H3N-H]-","[M+CHO2]-", "[M+C2H3O2]-","[M+Br]-", "[2M-H]-", "[M]-"), tolerance = 0, ppm = 5) 
parm <- Mass2MzRtParam(adducts = c("[M-H]-", "[M-2H]2-", "[M+Na-2H]-","[M+Cl]-","[M+K-2H]-","[M+C2H3N-H]-","[M+CHO2]-", "[M+C2H3O2]-","[M+Br]-", "[2M-H]-", "[M]-"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.12)
#MetaboCoreUtils::adducts(polarity = c("positive", "negative")) 
#MetaboCoreUtils::adducts(polarity = c("negative")) 

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_compound)
unique(md$target_compound) %>% length()

# Plot
library(ggsci)

adduct_summary <- tibble(adduct = md$adduct) %>%
  dplyr::count(adduct, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  arrange(desc(n))

total_n <- sum(adduct_summary$n)

compound_summary <- tibble(target_Compound = md$target_compound) %>%
  dplyr::count(target_Compound, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  arrange(desc(n))

ggplot(adduct_summary, aes(x = reorder(adduct, n), y = n)) +
  geom_col(color = "black", aes(fill = adduct), size = 1, width = 0.5) +
  #coord_flip() +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, percent)), hjust = 0.5, vjust = -0.5, size = 5.5) +
  labs(x = NULL, y = "Count")+
  #labs(x = NULL, y = "Count", title = paste("Adduct counts (Total =", total_n, ")")) +
  expand_limits(y = max(adduct_summary$n) * 1.15) + theme_classic(base_size = 16) + theme(legend.position = "none")+
  scale_fill_npg() 

adduct_summary_raw <- adduct_summary

# plot occurance
occ_tbl <- md %>%
  mutate(target_compound = trimws(as.character(target_compound))) %>%
  filter(!is.na(target_compound), target_compound != "") %>%
  dplyr::count(target_compound, name = "occurrence") %>%
  arrange(desc(occurrence), target_compound)

occ_dist <- occ_tbl %>%
  dplyr::count(occurrence, name = "n_compounds") %>%
  arrange(occurrence) %>%
  mutate(
    occurrence_f = factor(occurrence, levels = occurrence),
    percent = 100 * n_compounds / sum(n_compounds),
    lab = sprintf("%d (%.1f%%)", n_compounds, percent)
  )

ggplot(occ_dist, aes(x = occurrence_f, y = n_compounds)) +
  geom_col(color = "black", aes(fill = occurrence_f), size = 1, width = 0.7) +
  geom_text(aes(label = lab), vjust = -0.3, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = "Adducts per compound",
    y = "Number of compounds"
  ) +
  theme_classic(base_size = 20) + theme(legend.position = "none") + scale_fill_aaas()

occ_dist_raw <- occ_dist

# Combine App & Raw ----
library(cowplot)
occ_dist <- rbind(cbind(occ_dist_raw, Label = "Raw"), cbind(occ_dist_app, Label = "App"))

p2 <- ggplot(occ_dist, aes(x = occurrence_f, y = n_compounds, fill = Label)) +
      geom_col(color = "black", position = position_dodge(width = 0.8, preserve = "single"), size = 1, width = 0.7) +
      #geom_text(aes(label = paste0(round(percent,1), " %")), position = position_dodge(width = 0.8), vjust = -0.3, size = 5) +
      #geom_text(aes(label = n_compounds), position = position_dodge(width = 0.8), vjust = -0.3, size = 5) + 
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)), n.breaks = 7) +
      labs(
        x = "Adducts per Compound",
        y = "Count"
      ) +
      #+ coord_flip()
      theme_classic(base_size = 16) + scale_fill_npg() +
      theme(legend.position = c(0.85, 0.85), #legend.justification = c(1, 1), 
           # legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
           legend.text = element_text(size = 15),         # Makes the group names larger
        legend.key.size = unit(1.0, "cm"),
            legend.title = element_blank()) 
  
adduct_summary <- rbind(cbind(adduct_summary_raw, Label = "Raw"), cbind(adduct_summary_app, Label = "App"))

p1 <- ggplot(adduct_summary, aes(x = reorder(adduct, -n), y = n, fill = Label)) +
      geom_col(color = "black",size = 1, width = 0.5, position = position_dodge(width = 0.8, preserve = "single"),) +
      #coord_flip() +
      #geom_text(aes(label = paste0(round(percent,1), " %")), hjust = 0.5, vjust = -0.5, size = 5.5, position = position_dodge(width = 0.8)) +
      #geom_text(aes(label = n), hjust = 0.5, vjust = -0.5, size = 5.5, position = position_dodge(width = 0.8)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0)), n.breaks = 7) +
      #labs(x = "Adduct Type", y = "Count")+
      labs(x = NULL, y = "Count")+    
      #labs(x = NULL, y = "Count", title = paste("Adduct counts (Total =", total_n, ")")) +
      expand_limits(y = max(adduct_summary$n) * 1.15) + theme_classic(base_size = 16) + theme(legend.position = "none")+
      scale_fill_npg() +theme(legend.position = c(0.85, 0.85), #legend.justification = c(1, 1), 
           # legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
           legend.text = element_text(size = 18),         # Makes the group names larger
        legend.key.size = unit(1.0, "cm"), #axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_blank()) 

shared_legend <- get_legend(
  p1 + theme(legend.box.margin = margin(0, 0, 0, -150)) # Adds a little breathing room
)

p1 <- p1 + theme(legend.position = "none")+ coord_flip()
p2 <- p2 + theme(legend.position = "none")+ coord_flip()

plot_row <- plot_grid(
  p2, 
  p1, 
  labels = c('A', 'B'), 
  label_size = 25, 
  nrow = 1
)

plot_grid(
  plot_row, 
  shared_legend, 
  rel_widths = c(3, .01), 
  nrow = 1
)

#................................................................
#### Raw stats for MetaboAnalyst ----
#................................................................

df <- read_csv("xcms_table.csv") %>% as.data.frame()
mzrt <- cbind(mz = df$mzmed, rt = df$rtmed) %>% as.data.frame()
mzrt$cn <- paste0(round(mzrt$mz, 5), "@", round(mzrt$rt/60, 2))
df <- df[,-c(1:3)] %>% t() %>% as.data.frame()
colnames(df) <- mzrt$cn

tbl <- "xcms_table.csv"
bn <- tools::file_path_sans_ext(basename(tbl))
bn2 <- bn

Label <- do.call(rbind, str_split(rownames(df), "_"))[,3]
df <- cbind(Label, df) %>% as.data.frame()
df[,-1] <- sapply(df[,-1], as.numeric)
ds <- df

# prepare data
df <- ds
df <- subset(df, df$Label != "QC")
df$Label <- factor(df$Label, levels=c("KO-6", "KO-22", "PM-31", "PM-37", "WT"))
df$Label <- as.factor(df$Label)

# MVI
ds_mvi <- df[,-1]
ds_mvi[ds_mvi == 0] <- NA

noise <- as.numeric(quantile(1:min(ds_mvi, na.rm = T))[2]) # reduce noise by quantile: quantile(1:noise)[1] = 0% ... quantile(1:noise)[4] = 100%
sd <- as.numeric(30) # set sd value for random value generation as.numeric(quantile(1:noise)[2]) or noise*0.3
set.seed(1234)
ds_mvi[is.na(ds_mvi)] <- 0 # convert NA into 0
NAidx <- ds_mvi == 0 # NA index
imp.rand <- abs(rnorm(sum(NAidx), mean = noise, sd = sd)) # generate random values. other option: runif(sum(NAidx), noise-sd, noise+sd)
ds_mvi[NAidx] <- imp.rand
df <- as.data.frame(cbind(Label = df[,1], ds_mvi))

# Stat Test 
gr <- as.factor(df$Label)
g <- t(combn(as.character(levels(gr)), 2))
g <- as.data.frame(as.matrix(g))

g_sel <- which(g[,2] == "WT")
g_sel <- g[g_sel,]

df_l <- list()
names_comp <- c()

for (i in 1:nrow(g_sel)) {
  df_l[[i]] <- as.data.frame(subset(df, df$Label %in% g_sel[i,]))
  df_l[[i]][,-1] <- sapply(df_l[[i]][,-1], as.numeric)
  df_l[[i]][,1] <- as.factor(as.character(df_l[[i]][,1]))
  names_comp[i] <- paste0(g_sel[i,1]," / ",g_sel[i,2])
}

names(df_l) <- names_comp

stat_test_l <- lapply(1:length(df_l), function(x) sapply(2:ncol(df_l[[x]]), function(y) t.test(log2(df_l[[x]][,y]+1.1) ~ df_l[[x]][,1], var.equal = F)$p.value)) # t.test # wilcox.test
stat_test_l <- lapply(1:length(stat_test_l), function(x) p.adjust(stat_test_l[[x]], "BH"))

stat_test_l2 <- lapply(1:length(df_l), function(x) sapply(2:ncol(df_l[[x]]), function(y) t.test(log2(df_l[[x]][,y]+1.1) ~ df_l[[x]][,1], var.equal = F)$statistic)) 

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

stat_test_df2 <- do.call(rbind, stat_test_l2)
colnames(stat_test_df2) <- colnames(df_l[[1]][,-1])
stat_test_df2 <- as.data.frame(cbind(Groups = names(df_l), stat_test_df2))
stat_test_df2[,-1] <- sapply(stat_test_df2[,-1], as.numeric)
stat_test_df2 <- data.frame(sapply(stat_test_df2, function(x) ifelse(is.nan(x), 0, x)))

stat_test_df_t2 <- as.data.frame(t(stat_test_df2))
colnames(stat_test_df_t2) <- stat_test_df_t2[1,]
stat_test_df_t2 <- as.data.frame(cbind(Compound = colnames(df), stat_test_df_t2))
stat_test_df_t2 <- stat_test_df_t2[-1,]

stat_test_df_v2 <- reshape2::melt(stat_test_df_t2, value.name = "t-value", variable.name = "Groups", id=1)

stat_test_df_v_pv <- left_join(stat_test_df_v, stat_test_df_v2, by = c("Compound", "Groups"))

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

# Combine
data_stat <- full_join(stat_test_df_v_pv, fdr_v, by = c("Groups", "Compound"))
data_stat[,3:5] <- sapply(data_stat[,3:5], as.numeric)
length(which(data_stat$`Adj.p-value`<0.05))
data_stat$`Adj.p-value.log` <- -log10(data_stat$`Adj.p-value`)

# Create Files for MetaboAnalyst
out_dir <- "MetaboAnalyst_mummichog_inputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data_mummi <- data_stat %>%
  mutate(
    Compound = as.character(Compound),
    .nums = str_extract_all(Compound, "[0-9]+\\.?[0-9]*"),
    mz = as.numeric(vapply(.nums, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))),
    rt = as.numeric(vapply(.nums, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))),
    p.value = as.numeric(`Adj.p-value`),
    t.score = as.numeric(`t-value`)
  ) %>%
  filter(is.finite(mz), is.finite(p.value), is.finite(t.score)) 

groups <- sort(unique(data_mummi$Groups))

for (g in groups) {

  out <- data_mummi %>%
    filter(Groups == g) %>%
    arrange(`p.value`) %>%   
    dplyr::transmute(
      `m.z`    = mz,
      `p.value` = p.value,
      `t.score` = t.score,
      `rt` = round(rt*60, 5)
    )

  # make safe filename from group string
  fname <- g %>%
    str_replace_all("[/\\\\]+", "_vs_") %>%
    str_replace_all("[^A-Za-z0-9._-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")

  file_out <- file.path(out_dir, paste0("mummichog_", bn2, " ", fname, ".txt"))

  # tab-delimited .txt
  write_tsv(out, file_out)
}

message("Saved files to: ", normalizePath(out_dir))

#................................................................
#### App stats for MetaboAnalyst ----
#................................................................

df <- read_csv("xcms_table_filtered.csv") %>% as.data.frame()
mzrt <- cbind(mz = df$mzmed, rt = df$rtmed) %>% as.data.frame()
mzrt$cn <- paste0(round(mzrt$mz, 5), "@", round(mzrt$rt/60, 2))
df <- df[,-c(1:3)] %>% t() %>% as.data.frame()
colnames(df) <- mzrt$cn

tbl <- "xcms_table_filtered.csv"
bn <- tools::file_path_sans_ext(basename(tbl))
bn2 <- str_remove(bn, "xcms_table_")
bn2

Label <- do.call(rbind, str_split(rownames(df), "_"))[,3]
df <- cbind(Label, df) %>% as.data.frame()
df[,-1] <- sapply(df[,-1], as.numeric)
ds <- df

# prepare data
df <- ds
df <- subset(df, df$Label != "QC")
df$Label <- factor(df$Label, levels=c("KO-6", "KO-22", "PM-31", "PM-37", "WT"))
df$Label <- as.factor(df$Label)

# MVI
ds_mvi <- df[,-1]
ds_mvi[ds_mvi == 0] <- NA

noise <- as.numeric(quantile(1:min(ds_mvi, na.rm = T))[2]) # reduce noise by quantile: quantile(1:noise)[1] = 0% ... quantile(1:noise)[4] = 100%
sd <- as.numeric(30) # set sd value for random value generation as.numeric(quantile(1:noise)[2]) or noise*0.3
set.seed(1234)
ds_mvi[is.na(ds_mvi)] <- 0 # convert NA into 0
NAidx <- ds_mvi == 0 # NA index
imp.rand <- abs(rnorm(sum(NAidx), mean = noise, sd = sd)) # generate random values. other option: runif(sum(NAidx), noise-sd, noise+sd)
ds_mvi[NAidx] <- imp.rand
df <- as.data.frame(cbind(Label = df[,1], ds_mvi))

# Stat Test 
gr <- as.factor(df$Label)
g <- t(combn(as.character(levels(gr)), 2))
g <- as.data.frame(as.matrix(g))

g_sel <- which(g[,2] == "WT")
g_sel <- g[g_sel,]

df_l <- list()
names_comp <- c()

for (i in 1:nrow(g_sel)) {
  df_l[[i]] <- as.data.frame(subset(df, df$Label %in% g_sel[i,]))
  df_l[[i]][,-1] <- sapply(df_l[[i]][,-1], as.numeric)
  df_l[[i]][,1] <- as.factor(as.character(df_l[[i]][,1]))
  names_comp[i] <- paste0(g_sel[i,1]," / ",g_sel[i,2])
}

names(df_l) <- names_comp

stat_test_l <- lapply(1:length(df_l), function(x) sapply(2:ncol(df_l[[x]]), function(y) t.test(log2(df_l[[x]][,y]+1.1) ~ df_l[[x]][,1], var.equal = F)$p.value)) # t.test # wilcox.test
stat_test_l <- lapply(1:length(stat_test_l), function(x) p.adjust(stat_test_l[[x]], "BH"))

stat_test_l2 <- lapply(1:length(df_l), function(x) sapply(2:ncol(df_l[[x]]), function(y) t.test(log2(df_l[[x]][,y]+1.1) ~ df_l[[x]][,1], var.equal = F)$statistic)) 

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

stat_test_df2 <- do.call(rbind, stat_test_l2)
colnames(stat_test_df2) <- colnames(df_l[[1]][,-1])
stat_test_df2 <- as.data.frame(cbind(Groups = names(df_l), stat_test_df2))
stat_test_df2[,-1] <- sapply(stat_test_df2[,-1], as.numeric)
stat_test_df2 <- data.frame(sapply(stat_test_df2, function(x) ifelse(is.nan(x), 0, x)))

stat_test_df_t2 <- as.data.frame(t(stat_test_df2))
colnames(stat_test_df_t2) <- stat_test_df_t2[1,]
stat_test_df_t2 <- as.data.frame(cbind(Compound = colnames(df), stat_test_df_t2))
stat_test_df_t2 <- stat_test_df_t2[-1,]

stat_test_df_v2 <- reshape2::melt(stat_test_df_t2, value.name = "t-value", variable.name = "Groups", id=1)

stat_test_df_v_pv <- left_join(stat_test_df_v, stat_test_df_v2, by = c("Compound", "Groups"))

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

# Combine
data_stat <- full_join(stat_test_df_v_pv, fdr_v, by = c("Groups", "Compound"))
data_stat[,3:5] <- sapply(data_stat[,3:5], as.numeric)
length(which(data_stat$`Adj.p-value`<0.05))
data_stat$`Adj.p-value.log` <- -log10(data_stat$`Adj.p-value`)

# Create Files for MetaboAnalyst
out_dir <- "MetaboAnalyst_mummichog_inputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data_mummi <- data_stat %>%
  mutate(
    Compound = as.character(Compound),
    .nums = str_extract_all(Compound, "[0-9]+\\.?[0-9]*"),
    mz = as.numeric(vapply(.nums, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))),
    rt = as.numeric(vapply(.nums, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))),
    p.value = as.numeric(`Adj.p-value`),
    t.score = as.numeric(`t-value`)
  ) %>%
  filter(is.finite(mz), is.finite(p.value), is.finite(t.score)) 

groups <- sort(unique(data_mummi$Groups))

for (g in groups) {

  out <- data_mummi %>%
    filter(Groups == g) %>%
    arrange(`p.value`) %>%   
    dplyr::transmute(
      `m.z`    = mz,
      `p.value` = p.value,
      `t.score` = t.score,
      `rt` = round(rt*60, 5)
    )

  # make safe filename from group string
  fname <- g %>%
    str_replace_all("[/\\\\]+", "_vs_") %>%
    str_replace_all("[^A-Za-z0-9._-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")

  file_out <- file.path(out_dir, paste0("mummichog_", bn2, " ", fname, ".txt"))

  # tab-delimited .txt
  write_tsv(out, file_out)
}

message("Saved files to: ", normalizePath(out_dir))

#................................................................
#### Functional Analysis Output from MetaboAnalyst----
#................................................................

raw <- read_csv("mummichog_pathway_raw.csv")
colnames(raw)[1] <- "Pathway"
app <- read_csv("mummichog_pathway_filt.csv")
colnames(app)[1] <- "Pathway"
app$`P(Fisher)`

for_plot <- rbind(cbind(app, Label = "App"), cbind(raw, Label = "Raw")) %>% as.data.frame()
for_plot <- for_plot %>%
  group_by(Pathway) %>%
  # Keep the pathway ONLY if App < 0.1 AND Raw < 0.1
  filter(
    any(Label == "App" & `P(Fisher)` < 0.5) & 
    any(Label == "Raw" & `P(Fisher)` < 0.5)
  ) %>%
  ungroup()

# 1. Prepare the Data for Mirroring
plot_data_mirror <- for_plot %>%
  # Calculate standard log transformation
  mutate(neg_log_p = -log10(`P(Fisher)`)) %>%
  
  # THE MIRROR TRICK: Make 'App' values negative so they draw to the left
  mutate(plot_val = ifelse(Label == "App", -neg_log_p, neg_log_p)) %>%
  
  # Sort pathways purely by how well the App did, so the App's best is at the top
  group_by(Pathway) %>%
  mutate(App_Score = max(ifelse(Label == "App", neg_log_p, -999), na.rm = TRUE)) %>%
  ungroup() %>%
  
  # Filter to keep it clean (optional: only keep if at least one side is significant)
  filter(any(neg_log_p > -log10(0.05))) %>%
  
  # Reorder the factor based on the App's score
  mutate(Pathway = fct_reorder(Pathway, App_Score))

# 2. Generate the Mirror Plot
ggplot(plot_data_mirror, aes(x = plot_val, y = Pathway, fill = Label)) +
  
  # Add the significance threshold lines for BOTH sides (0.05 is ~1.301 on log scale)
  geom_vline(xintercept = -(-log10(0.05)), linetype = "dashed", color = "gray50", linewidth = 1) + # App threshold (Left)
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "gray50", linewidth = 1) +  # Raw threshold (Right)
  
  # Draw the mirrored bars
  geom_col(width = 0.7, alpha = 1, color = "black") +

  # Center line separating the two sides
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  
  scale_fill_npg() +
  
  # THE SYMMETRIC ABSOLUTE SCALE
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.05)),
    # We place breaks symmetrically around 0
    breaks = c(-3, -2, -1.301, 0, 1.301, 2, 3),
    # The labels are perfectly mirrored absolute p-values!
    labels = c("0.001", "0.01", "0.05", "1.0", "0.05", "0.01", "0.001")
  ) +
  
  labs(
    title = "",
    x = "← App Better (Absolute P-Value) Raw Better →",
    y = "",
    fill = ""
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(margin = margin(t = 0), face = "bold"),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom"
  )                    
                    
#................................................................
#### Compare Total Annotations ----
#................................................................

# App
df <- read_csv("xcms_table_filtered.csv") %>% as.data.frame()

target_df <- read.csv("All_affected_cycles_NEG upd.csv")
mass <- calculateMass(target_df$formula)
target_df$exactmass <- mass
target_df <- target_df[,c(1,2,9,4)]

peakIn <- as.data.frame(cbind(mz = df$mzmed, rt = df$rtmed, id = df$feature)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- round(as.numeric(peakIn$rt)/60, 4)

parm <- Mass2MzRtParam(adducts = c("[M-H]-", "[M-2H]2-", "[M+Na-2H]-","[M+Cl]-","[M+K-2H]-","[M+C2H3N-H]-","[M+CHO2]-", "[M+C2H3O2]-","[M+Br]-", "[2M-H]-", "[M]-"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.12)

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_compound)
unique(md$target_compound) %>% length()

# Raw
df <- read_csv("xcms_table.csv") %>% as.data.frame()

target_df <- read.csv("All_affected_cycles_NEG upd.csv")
mass <- calculateMass(target_df$formula)
target_df$exactmass <- mass
target_df <- target_df[,c(1,2,9,4)]

peakIn <- as.data.frame(cbind(mz = df$mzmed, rt = df$rtmed, id = df$feature)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- round(as.numeric(peakIn$rt)/60, 4)

parm <- Mass2MzRtParam(adducts = c("[M-H]-", "[M-2H]2-", "[M+Na-2H]-","[M+Cl]-","[M+K-2H]-","[M+C2H3N-H]-","[M+CHO2]-", "[M+C2H3O2]-","[M+Br]-", "[2M-H]-", "[M]-"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.12)

matched_features <- matchValues(peakIn, target_df, parm)
md2 <- matchedData(matched_features)
md2 <- as.data.frame(md2)
md2 <- na.omit(md2)
unique(md2$target_compound)
unique(md2$target_compound) %>% length()

# Compare
#unique(md2$target_compound) %>% as.data.frame() %>% View()
intersect(md$target_compound, md2$target_compound)
setdiff(md2$target_compound, md$target_compound)

setdiff(target_df$compound, md$target_compound)
