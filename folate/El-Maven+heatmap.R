#...........................................................
# El-Maven Output Stats ----
#...........................................................

library(stringr)
library(dplyr)
library(data.table)
library(reshape2)
library(plotly) 
library(ggsci)
library(ggplot2)
library(tidyr)

setwd("C:/..")

df <- as.data.frame(fread("Results POS-NEG all replicates.csv"))
rownames(df) <- df[,1]
df <- df[,-1]
df$Label <- factor(df$Label, levels=c("KO-6", "KO-22", "PM-31", "PM-37", "WT"))
df$Label <- as.factor(df$Label)

# MVI 
# Min Sample
ds_mvi <- as.data.frame(t(df[,-1]))
ds_mvi[ds_mvi == 0] <- NA

impute.MinSample <- function (dataSet.mvs) ({
  nSamples = dim(dataSet.mvs)[2]
  dataSet.imputed = dataSet.mvs
  minValue.samples = apply(dataSet.imputed, 2, min, na.rm = T)
  for (i in 1:(nSamples)) {
    dataSet.imputed[which(is.na(dataSet.mvs[, i])), i] = minValue.samples[i]
  }
  return(dataSet.imputed)
})

ds_mvi <- as.data.frame(t(impute.MinSample(ds_mvi)))
#df <- as.data.frame(cbind(Label = df$Label, ds_mvi))

# Noise Imputation based on Min Value
ds_mvi <- df[,-1]
ds_mvi[ds_mvi == 0] <- NA
noise <- as.numeric(quantile(1:min(ds_mvi, na.rm = T))[2]) # reduce noise by quantile: quantile(1:noise)[1] = 0% ... quantile(1:noise)[4] = 100%
sd <- as.numeric(quantile(1:noise)[2]) # set sd value for random value generation as.numeric(quantile(1:noise)[2]) or noise*0.3
set.seed(1234)
ds_mvi[is.na(ds_mvi)] <- 0 # convert NA into 0
NAidx <- ds_mvi == 0 # NA index
imp.rand <- abs(rnorm(sum(NAidx), mean = noise, sd = sd)) # generate random values. other option: runif(sum(NAidx), noise-sd, noise+sd)
ds_mvi[NAidx] <- imp.rand
df <- as.data.frame(cbind(Label = df$Label, ds_mvi))

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

# t.test # wilcox.test
stat_test_l <- lapply(1:length(df_l), function(x) sapply(2:ncol(df_l[[x]]), function(y) t.test(df_l[[x]][,y] ~ df_l[[x]][,1])$p.value)) # t.test # wilcox.test
# wilcox_test
#stat_test_l <- lapply(1:length(df_l), function(x) sapply(2:ncol(df_l[[x]]), function(y) pvalue(wilcox_test(df_l[[x]][,y] ~ df_l[[x]][,1])))) # wilcox_test

# adjust
stat_test_l <- lapply(1:length(stat_test_l), function(x) p.adjust(stat_test_l[[x]], "BH"))
#stat_test_l <- lapply(1:length(stat_test_l), function(x) p.adjust(stat_test_l[[x]], "bonf"))
#stat_test_l <- lapply(1:length(stat_test_l), function(x) p.adjust(stat_test_l[[x]], "BY"))

# data frame of results
stat_test_df <- do.call(rbind, stat_test_l)
colnames(stat_test_df) <- colnames(df_l[[1]][,-1])
#rownames(stat_test_df) <- names(df_l)
stat_test_df <- as.data.frame(cbind(Groups = names(df_l), stat_test_df))
stat_test_df[,-1] <- sapply(stat_test_df[,-1], as.numeric)
stat_test_df <- data.frame(sapply(stat_test_df, function(x) ifelse(is.nan(x), 1, x)))

stat_test_df_t <- as.data.frame(t(stat_test_df))
colnames(stat_test_df_t) <- stat_test_df_t[1,]
stat_test_df_t <- as.data.frame(cbind(Compound = colnames(df), stat_test_df_t))
stat_test_df_t <- stat_test_df_t[-1,]

stat_test_df_v <- reshape2::melt(stat_test_df_t, value.name = "Adj.p-value", variable.name = "Groups", id=1)

# FC
FOLD.CHANGE.MG <- function(x, f, aggr_FUN = colMeans, combi_FUN = {function(x,y) "-"(x,y)})({
  f <- as.factor(f)
  i <- split(1:nrow(x), f)
  x <- sapply(i, function(i)({ aggr_FUN(x[i,])}))
  x <- t(x)
  x <- log2(x)
  j <- combn(levels(f), 2)
  ret <- combi_FUN(x[j[1,],], x[j[2,],])
  rownames(ret) <- paste(j[1,], j[2,], sep = ' / ')
  t(ret)
})

fdr <- FOLD.CHANGE.MG(df[,-1], df[,1])
fdr0 <- fdr
fdr <- as.data.frame(fdr[,which(str_detect(colnames(fdr), "WT"))])
fdr <- as.data.frame(sapply(fdr, function(x) ifelse(is.nan(x), 0, x)))
fdr <- as.data.frame(sapply(fdr, function(x) ifelse(is.infinite(x), 0, x)))
rownames(fdr) <- rownames(fdr0)

fdr <- as.data.frame(cbind(Compound = colnames(df[,-1]), fdr))
fdr_v <- reshape2::melt(fdr, value.name = "FC", variable.name = "Groups", id=1)

# Combine
data_stat <- full_join(stat_test_df_v, fdr_v, by = c("Groups", "Compound"))
data_stat[,3:4] <- sapply(data_stat[,3:4], as.numeric)
data_stat$`Adj.p-value.log` <- -log10(data_stat$`Adj.p-value`)

data_stat$Compound_short <- data_stat$Compound
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "Pyr_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "Pur_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "Gly_TCA_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "SerGly_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "dTMP_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "SAM_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "His_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "Bet_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "PPP_")
data_stat$Compound_short <- do.call(rbind, str_split(data_stat$Compound_short, "_"))[,2]

data_stat <- full_join(stat_test_df_v, fdr_v, by = c("Groups", "Compound"))
data_stat[,3:4] <- sapply(data_stat[,3:4], as.numeric)
data_stat$`Adj.p-value.log` <- -log10(data_stat$`Adj.p-value`)

data_stat$Pathway <- data_stat$Compound
data_stat$Pathway[which(str_detect(data_stat$Pathway, "Add-Gly-TCA"))] <- "Energy"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "Add-Pyr"))] <- "Pyrimidine"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "Add-Pur"))] <- "Purine"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "Add-SerGly"))] <- "Serine_Glycine"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "Gly_TCA_"))] <- "Energy"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "Pur_"))] <- "Purine"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "Pyr_"))] <- "Pyrimidine"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "dTMP_"))] <- "Pyrimidine"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "dTMP_"))] <- "Pyrimidine"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "SAM_"))] <- "SAM"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "His_"))] <- "Histidine"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "Bet_"))] <- "Betaine"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "SerGly_"))] <- "Serine_Glycine"
data_stat$Pathway[which(str_detect(data_stat$Pathway, "PPP_"))] <- "PPP"

data_stat$Compound_short <- data_stat$Compound
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "Pyr_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "Pur_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "Gly_TCA_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "SerGly_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "dTMP_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "SAM_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "His_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "Bet_")
data_stat$Compound_short <- str_remove(data_stat$Compound_short, "PPP_")
data_stat$Compound_short <- do.call(rbind, str_split(data_stat$Compound_short, "_"))[,2]

data_stat$FC <- as.numeric(round(data_stat$FC, 1)) 
data_stat$`Adj.p-value.log` <- as.numeric(round(data_stat$`Adj.p-value.log`, 2)) 
data_stat$`Signigicant FC` <- ifelse(data_stat$`Adj.p-value` <= 0.05, data_stat$FC, 0)
data_stat$`Signigicant FC` <- ifelse(data_stat$`Signigicant FC` <= -1 | data_stat$`Signigicant FC` >= 1, data_stat$FC, 0)
data_stat$`Signigicant FC2` <- ifelse(data_stat$`Signigicant FC` != 0, T, F)
data_stat$FC_abs <- abs(data_stat$FC)
data_stat$`Adj.p-value2` <- as.numeric(round(data_stat$`Adj.p-value`, 3))

#...........................................................
# heatmap ----
#...........................................................
ids <- read.csv("../All_IDs_for_Affected_NEG_POS_csv.csv")
ids <- ids[c(1, 5,8, 6)]
colnames(ids)[1] <- "Compound"

data_stat_ids <- dplyr::left_join(data_stat, ids, by = "Compound")
which(is.na(data_stat_ids))

fwrite(data_stat_ids[,-c(9,10)], "data stats + ids.csv", row.names = T)

df_lf <- as.data.frame(fread("data stats + ids.csv"))[,-1]

df_all_media <- df_lf
df_all_media$Media <- "Low Folate"
df_all_media$Media <- as.factor(df_all_media$Media)
df_all_media$Pathway1 <- ifelse(str_detect(df_all_media$Compound, "Gly_TCA_"), "Energy", NA)
df_all_media$Pathway2 <- ifelse(str_detect(df_all_media$Compound, "Pur_"), "Purine", NA)
df_all_media$Pathway3 <- ifelse(str_detect(df_all_media$Compound, "Pyr_"), "Pyrimidine", NA)
df_all_media$Pathway4 <- ifelse(str_detect(df_all_media$Compound, "SerGly_"), "Ser-Gly", NA)
df_all_media$Pathway5 <- ifelse(str_detect(df_all_media$Compound, "SAM_"), "SAM", NA)
df_all_media$Pathway6 <- ifelse(str_detect(df_all_media$Compound, "PPP_"), "PPP", NA)

df_all_media <- df_all_media %>% pivot_longer(cols=c('Pathway1', 'Pathway2', 'Pathway3', 'Pathway4', 'Pathway5', 'Pathway6'),
                                              names_to='Pathways',
                                              values_to='Pathway_all') %>% as.data.frame()
df_all_media <- na.omit(df_all_media)
df_all_media$Groups <- factor(df_all_media$Groups, levels=c("KO-6 / WT", "KO-22 / WT", "PM-31 / WT", "PM-37 / WT"))


df_all_media_sub <- df_all_media
df_all_media_sub$Compound <- str_remove(df_all_media_sub$Compound, "Pur_")
df_all_media_sub$Compound <- str_remove(df_all_media_sub$Compound, "Pyr_")
df_all_media_sub$Compound <- str_remove(df_all_media_sub$Compound, "SerGly_")
df_all_media_sub$Compound <- str_remove(df_all_media_sub$Compound, "dTMP_")
df_all_media_sub$Compound <- str_remove(df_all_media_sub$Compound, "PPP_")
df_all_media_sub$Compound <- str_remove(df_all_media_sub$Compound, "Gly_TCA_")
df_all_media_sub$Compound <- str_remove(df_all_media_sub$Compound, "SAM_")
df_all_media_sub$Compound <- do.call(rbind, str_split(df_all_media_sub$Compound, "_"))[,1]
df_all_media_sub$Compound <- str_remove(df_all_media_sub$Compound, "a")
df_all_media_sub$Compound <- str_remove(df_all_media_sub$Compound, "b")
df_all_media_sub$Compound <- str_remove(df_all_media_sub$Compound, "A")
df_all_media_sub$Compound <- as.numeric(do.call(rbind, str_split(df_all_media_sub$Compound, "/"))[,1])
df_all_media_sub <- df_all_media_sub[order(df_all_media_sub$Compound, decreasing = F), ]

df_all_media_sub$Compound_short <- str_replace(df_all_media_sub$Compound_short, "bisphosphate", "2P")
df_all_media_sub$Compound_short <- str_replace(df_all_media_sub$Compound_short, "phosphate", "P")
df_all_media_sub$Compound_short <- str_replace(df_all_media_sub$Compound_short, "S-Adenosyl-L-methionine", "SAM")
df_all_media_sub$Compound_short <- str_replace(df_all_media_sub$Compound_short, "S-5-Adenosyl-L-homocysteine", "SAH")

# change order
ind_del_kg <- c(which(df_all_media_sub$Compound_short == "Ketoglutarate" & df_all_media_sub$Pathway_all == "Purine"),
                which(df_all_media_sub$Compound_short == "Ketoglutarate" & df_all_media_sub$Pathway_all == "Pyrimidine"))
df_all_media_sub <- df_all_media_sub[-ind_del_kg,]

df_all_media_sub[df_all_media_sub$Compound_short == "3-phospho-D-glycerate" & df_all_media_sub$Pathway_all == "Ser-Gly", ]$Compound <- 1
df_all_media_sub[df_all_media_sub$Compound_short == "Glutamic acid" & df_all_media_sub$Pathway_all == "Ser-Gly", ]$Compound <- 3
df_all_media_sub[df_all_media_sub$Compound_short == "Ketoglutarate" & df_all_media_sub$Pathway_all == "Ser-Gly", ]$Compound <- 4
df_all_media_sub[df_all_media_sub$Compound_short == "Serine" & df_all_media_sub$Pathway_all == "Ser-Gly", ]$Compound <- 6
df_all_media_sub[df_all_media_sub$Compound_short == "Glycine" & df_all_media_sub$Pathway_all == "Ser-Gly", ]$Compound <- 7

df_all_media_sub[df_all_media_sub$Compound_short == "5-PRPP" & df_all_media_sub$Pathway_all == "Purine", ]$Compound <- 1
df_all_media_sub[df_all_media_sub$Compound_short == "Glutamine" & df_all_media_sub$Pathway_all == "Purine", ]$Compound <- 2.1
df_all_media_sub[df_all_media_sub$Compound_short == "Glutamic acid" & df_all_media_sub$Pathway_all == "Purine", ]$Compound <- 2.5
df_all_media_sub[df_all_media_sub$Compound_short == "PRA" & df_all_media_sub$Pathway_all == "Purine", ]$Compound <- 2.9
df_all_media_sub[df_all_media_sub$Compound_short == "Glycine" & df_all_media_sub$Pathway_all == "Purine", ]$Compound <- 2.8
df_all_media_sub[df_all_media_sub$Compound_short == "Aspartate" & df_all_media_sub$Pathway_all == "Purine", ]$Compound <- 7.1
df_all_media_sub[df_all_media_sub$Compound_short == "ADP" & df_all_media_sub$Pathway_all == "Purine", ]$Compound <- 14
df_all_media_sub[df_all_media_sub$Compound_short == "ATP" & df_all_media_sub$Pathway_all == "Purine", ]$Compound <- 14.5
df_all_media_sub[df_all_media_sub$Compound_short == "Aspartate" & df_all_media_sub$Pathway_all == "Purine", ]$Compound <- 7.1
df_all_media_sub[df_all_media_sub$Compound_short == "Aspartate" & df_all_media_sub$Pathway_all == "Purine", ]$Compound <- 7.1

df_all_media_sub[df_all_media_sub$Compound_short == "Glutamic acid" & df_all_media_sub$Pathway_all == "Pyrimidine", ]$Compound <- 17.5
df_all_media_sub[df_all_media_sub$Compound_short == "5-PRPP" & df_all_media_sub$Pathway_all == "Pyrimidine", ]$Compound <- 7
df_all_media_sub[df_all_media_sub$Compound_short == "UTP" & df_all_media_sub$Pathway_all == "Pyrimidine", ]$Compound <- 17.4

df_all_media_sub[df_all_media_sub$Compound_short == "Methionine" & df_all_media_sub$Pathway_all == "SAM", ]$Compound <- 6

df_all_media_sub[df_all_media_sub$Compound_short == "D-glucose 6-P" & df_all_media_sub$Pathway_all == "PPP", ]$Compound <- 1
df_all_media_sub[df_all_media_sub$Compound_short == "NADP" & df_all_media_sub$Pathway_all == "PPP", ]$Compound <- 1.5
df_all_media_sub[df_all_media_sub$Compound_short == "NADPH" & df_all_media_sub$Pathway_all == "PPP", ]$Compound <- 1.8
df_all_media_sub[df_all_media_sub$Compound_short == "D-gluconate 6-P" & df_all_media_sub$Pathway_all == "PPP", ]$Compound <- 2

df_all_media_sub[df_all_media_sub$Compound_short == "Lactic acid" & df_all_media_sub$Pathway_all == "Energy", ]$Compound <- 10.8
df_all_media_sub[df_all_media_sub$Compound_short == "NADH" & df_all_media_sub$Pathway_all == "Energy", ]$Compound <- 13.9
df_all_media_sub[df_all_media_sub$Compound_short == "Succinyl-CoA" & df_all_media_sub$Pathway_all == "Energy", ]$Compound <- 14.1
df_all_media_sub[df_all_media_sub$Compound_short == "ADP" & df_all_media_sub$Pathway_all == "Energy", ]$Compound <- 14.5
df_all_media_sub[df_all_media_sub$Compound_short == "ATP" & df_all_media_sub$Pathway_all == "Energy", ]$Compound <- 14.6
df_all_media_sub[df_all_media_sub$Compound_short == "GDP" & df_all_media_sub$Pathway_all == "Energy", ]$Compound <- 14.7
df_all_media_sub[df_all_media_sub$Compound_short == "GTP" & df_all_media_sub$Pathway_all == "Energy", ]$Compound <- 14.8

df_all_media_sub2 <- df_all_media_sub
df_all_media_sub2$Significance <- ifelse(df_all_media_sub2$`Signigicant FC` != 0, "*", NA)

# Horizont
df_all_media_sub3 <- subset(df_all_media_sub2, df_all_media_sub2$Groups %in% c("KO-6 / WT", "KO-22 / WT", "PM-31 / WT", "PM-37 / WT"))
df_all_media_sub3$Media <- ifelse(df_all_media_sub3$Media == "Regular", "High Folate", "Low Folate")
df_all_media_sub3 <- subset(df_all_media_sub2, df_all_media_sub2$Media %in% c("Low Folate"))
df_all_media_sub3$Media <- factor(df_all_media_sub3$Media, levels=c("Low Folate", "High Folate"))

p<- ggplot(df_all_media_sub3, aes(x=tidytext::reorder_within(Compound_short, Compound, Pathway_all), group = Groups, y=Media, fill=FC)) + 
  geom_tile(color="black", linewidth=0.1, height = 1, width = 1) +
  geom_text(aes(label = Significance), color = "black", size = 4)+
  #scale_alpha_discrete(range = c(1, 0.3))+
  xlab("")+
  ylab("")+
  tidytext::scale_x_reordered()+
  ggh4x::facet_grid2(cols = vars(Pathway_all), rows = vars(Groups), scales = "free", space = "free")+
  theme_bw(base_size = 12)+ theme(text = element_text(family = "Verdana"),legend.position = 'bottom', axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #facet_wrap(vars(index, Analyte), scales = "free")+
  #geom_text(aes(label = as.factor(Method)), position=position_dodge(width=0.9), angle = 90, hjust = 'left', alpha = 1, color = "red")+
  #ggh4x::facet_grid2(cols = vars(degree), rows = vars(index), scales = "free_y", independent = "y")+
  #scale_fill_brewer("Paired")+
  scale_fill_gradient2(midpoint=0, low="#1874CD", mid="white", high="#CD2626")+
  theme(plot.margin = margin(3,.8,2,.8, "cm"),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    #axis.ticks.x=element_blank(),
    legend.title = element_text(size = 15),
    plot.title = element_text(hjust=0.5, size=15, face = "bold"), strip.text = element_text(size = 10),
    legend.text=element_text(size=12), axis.title=element_text(size=15)) + scale_y_discrete(labels = NULL, breaks = NULL)

p+  labs(fill = "FC:", alpha = "") + theme(panel.spacing=unit(1,"lines"))

p+  labs(fill = "FC:", alpha = "") +  theme(panel.spacing=unit(1,"lines"))+
  geom_rect(data = subset(df_all_media_sub3, Compound_short %in% c("IMP")), 
                                               fill = NA, colour = "black", size = 1.5, xmin = 13.5, xmax = 13.5, 
                                               ymin = -Inf,ymax = Inf, inherit.aes = F) +
  geom_rect(data = subset(df_all_media_sub3, Compound_short %in% c("FGAR")), 
            fill = NA, colour = "black", size = 1.5, xmin = 6.5, xmax = 6.5, 
            ymin = -Inf,ymax = Inf, inherit.aes = F)  + 
  geom_rect(data = subset(df_all_media_sub3, Compound_short %in% c("dTMP")), 
            fill = NA, colour = "black", size = 1.5, xmin = 11.5, xmax = 11.5, 
            ymin = -Inf,ymax = Inf, inherit.aes = F) +
  geom_rect(data = subset(df_all_media_sub3, Compound_short %in% c("Methionine")), 
            fill = NA, colour = "black", size = 1.5, xmin = 4.5, xmax = 4.5, 
            ymin = -Inf,ymax = Inf, inherit.aes = F) +
  geom_rect(data = subset(df_all_media_sub3, Compound_short == "Glycine" & Pathway_all == "Ser-Gly"), 
            fill = NA, colour = "black", size = 1.5, xmin = 4.5, xmax = 4.5, 
            ymin = -Inf,ymax = Inf, inherit.aes = F) 
