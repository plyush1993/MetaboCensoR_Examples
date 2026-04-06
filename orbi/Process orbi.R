#................................................................
#### Settings ----
#................................................................

library(tidyverse)
library(MetaboAnnotation)
library(MetaboCoreUtils)
library(ggplot2)
library(data.table)
library(ggsci)
library(RCy3)
library(cowplot)

setwd("C:/.../")

#................................................................
#### Compare Annotations Total----
#................................................................

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

#target_df <- read.csv("annot table - rt.csv")
#mass <- calculateMass(target_df$Formula)
#target_df$exactmass <- mass
#target_df <- target_df[,c(1,2,3,4)]


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

#target_df <- read.csv("annot table - rt.csv")
#mass <- calculateMass(target_df$Formula)
#target_df$exactmass <- mass
#target_df <- target_df[,c(1,2,3,4)]

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

#................................................................
#### Annotation App ----
#................................................................

df <- read_csv("orbi_iimn_gnps_quant_filtered (2x H2O adducts).csv") %>% as.data.frame()
peakIn <- as.data.frame(cbind(mz = df$`row m/z`, rt = df$`row retention time`, id = df$`row ID`)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

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

# set polarity, adduct, accuracy
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+", "[M+H-Hexose-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.01)
#MetaboCoreUtils::adducts(polarity = c("positive", "negative")) 
#MetaboCoreUtils::adducts(polarity = c("negative")) 

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_Compound)
unique(md$target_Compound) %>% length()

md <- subset(md, !md$target_Compound %in% c("Wulignan A1", "p-coumaraldehyde (Q27103652)|Phenylacrylic acid"))
unique(md$target_Compound) %>% length()

# Plot
library(ggsci)

adduct_summary <- tibble(adduct = md$adduct) %>%
  dplyr::count(adduct, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  arrange(desc(n))

total_n <- sum(adduct_summary$n)

compound_summary <- tibble(target_Compound = md$target_Compound) %>%
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
  mutate(target_compound = trimws(as.character(target_Compound))) %>%
  filter(!is.na(target_Compound), target_Compound != "") %>%
  dplyr::count(target_Compound, name = "occurrence") %>%
  arrange(desc(occurrence), target_Compound)

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

df <- read_csv("orbi_iimn_gnps_quant.csv") %>% as.data.frame()
peakIn <- as.data.frame(cbind(mz = df$`row m/z`, rt = df$`row retention time`, id = df$`row ID`)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

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

# set polarity, adduct, accuracy
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+", "[M+H-Hexose-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.01)
#MetaboCoreUtils::adducts(polarity = c("positive", "negative")) 
#MetaboCoreUtils::adducts(polarity = c("negative")) 

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_Compound)
unique(md$target_Compound) %>% length()

md <- subset(md, !md$target_Compound %in% c("Wulignan A1", "p-coumaraldehyde (Q27103652)|Phenylacrylic acid"))
unique(md$target_Compound) %>% length()

# Plot
library(ggsci)

adduct_summary <- tibble(adduct = md$adduct) %>%
  dplyr::count(adduct, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  arrange(desc(n))

total_n <- sum(adduct_summary$n)

compound_summary <- tibble(target_Compound = md$target_Compound) %>%
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
  mutate(target_compound = trimws(as.character(target_Compound))) %>%
  filter(!is.na(target_Compound), target_Compound != "") %>%
  dplyr::count(target_Compound, name = "occurrence") %>%
  arrange(desc(occurrence), target_Compound)

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
#................................................................
# Combine App & Raw ----
#................................................................
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
#### Cytoscape ----
#................................................................

library(RCy3)
cytoscapePing()
cytoscapeVersionInfo()

# Steroids
openSession("app 2x H2O adds single.cys")
openSession("raw single.cys")

#openSession("raw.cys")
#openSession("App Output (1).cys")

imported_styles <- importVisualStyles(filename = "styles.xml")
setVisualStyle(imported_styles[1])

data <- read.csv("canopus_formula_summary + MetaboAnnotation RT.csv", check.names = FALSE)

loadTableData(
  data = data, 
  table.key.column = "shared name", 
  data.key.column = "mappingFeatureId"
)

#setNodeColorDefault('#FFFFFF', style.name = "ClassDefault") 
#setEdgeColorDefault('#999999', style.name = "ClassDefault") 
#setNodeBorderColorDefault('#111111', style.name = "ClassDefault")

selectNodes(nodes = "Steroids", by.col = "NPC.superclass")

# color
sel_nodes <- getSelectedNodes()
setNodeColorBypass(node.names = sel_nodes, new.colors = "#FF0000")
clearSelection()

# delete
#clearNodePropertyBypass(node.names = sel_nodes, 'NODE_FILL_COLOR')

# more advanced with border
setNodeBorderColorMapping(
  table.column = 'MetaboAnnotation_status', 
  table.column.values = c(TRUE, FALSE), 
  colors = c('#0000FF', '#000000'), 
  mapping.type = 'd',
  style.name = "ClassDefault_0"
)

###########################################31
# Flavonoids
openSession("app 2x H2O adds single.cys")
openSession("raw single.cys")

#openSession("raw.cys")
#openSession("App Output (1).cys")

imported_styles <- importVisualStyles(filename = "styles.xml")
setVisualStyle(imported_styles[1])

data <- read.csv("canopus_formula_summary + MetaboAnnotation RT.csv", check.names = FALSE)

loadTableData(
  data = data, 
  table.key.column = "shared name", 
  data.key.column = "mappingFeatureId"
)

#setNodeColorDefault('#FFFFFF', style.name = "ClassDefault") 
#setEdgeColorDefault('#999999', style.name = "ClassDefault") 
#setNodeBorderColorDefault('#111111', style.name = "ClassDefault")

selectNodes(nodes = "Flavonoids", by.col = "NPC.superclass")

# color
sel_nodes <- getSelectedNodes()
setNodeColorBypass(node.names = sel_nodes, new.colors = "#FF0000")
clearSelection()

# delete
#clearNodePropertyBypass(node.names = sel_nodes, 'NODE_FILL_COLOR')

# more advanced with border
setNodeBorderColorMapping(
  table.column = 'MetaboAnnotation_status', 
  table.column.values = c(TRUE, FALSE), 
  colors = c('#0000FF', '#000000'), 
  mapping.type = 'd',
  style.name = "ClassDefault_0"
)

#................................................................
#### Merge SIRIUS by MetaboAnnotation with Annotation ----
#................................................................
sirius <- read_csv("canopus_formula_summary.csv")

target_df <- read.csv("annotation.csv")

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

peakIn <- cbind(mz = as.numeric(sirius$ionMass), rt = as.numeric(sirius$retentionTimeInMinutes), map_id = sirius$mappingFeatureId) %>% as.data.frame()

# set polarity, adduct, accuracy
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+", "[M+H-Hexose-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.01)
#MetaboCoreUtils::adducts() 

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md$MetaboAnnotation_status <- ifelse(is.na(md$target_Compound), F, T)
md <- md[,c(3:4,10, 14)]

sirius_ma <- left_join(sirius, md, by = c("mappingFeatureId" = "map_id"))
write_csv(sirius_ma, "canopus_formula_summary + MetaboAnnotation RT.csv")

#................................................................
#### Filter SIRIUS Output by probability ----
#................................................................

sirius <- read_csv("canopus_formula_summary.csv")
colnames(sirius)

sirius$`NPC#pathway` <- ifelse(sirius$`NPC#pathway Probability` > 0.9, sirius$`NPC#pathway`, "")
sirius$`NPC#superclass` <- ifelse(sirius$`NPC#superclass Probability` > 0.9, sirius$`NPC#superclass`, "")
sirius$`NPC#class` <- ifelse(sirius$`NPC#class Probability` > 0.9, sirius$`NPC#class`, "")

sirius$`ClassyFire#superclass` <- ifelse(sirius$`ClassyFire#superclass probability` > 0.9, sirius$`ClassyFire#superclass`, "")
sirius$`ClassyFire#class` <- ifelse(sirius$`ClassyFire#class Probability` > 0.9, sirius$`ClassyFire#class`, "")
sirius$`ClassyFire#subclass` <- ifelse(sirius$`ClassyFire#subclass Probability` > 0.9, sirius$`ClassyFire#subclass`, "")
sirius$`ClassyFire#level 5` <- ifelse(sirius$`ClassyFire#level 5 Probability` > 0.9, sirius$`ClassyFire#level 5`, "")
sirius$`ClassyFire#most specific class` <- ifelse(sirius$`ClassyFire#most specific class Probability` > 0.9, sirius$`ClassyFire#most specific class`, "")

sirius <- sirius %>% mutate(across(everything(), ~ tidyr::replace_na(as.character(.x), "")))
write_csv(sirius, "canopus_formula_summary FILT.csv")

#................................................................
#### Percentage of Predicted classes ----
#................................................................

# raw
df <- read_csv("orbi_iimn_gnps_quant.csv")
sirius <- read_csv("canopus_formula_summary FILT.csv") 
na_count_per_row <- rowSums(is.na(sirius))
to_del <- which(na_count_per_row >= 10)
tot_an <- intersect(sirius$mappingFeatureId[-to_del], df$`row ID`) %>% length()
round(tot_an/nrow(df)*100, 1)

# app
df <- read_csv("orbi_iimn_gnps_quant_filtered (2x H2O adducts).csv")
sirius <- read_csv("canopus_formula_summary FILT.csv") 
na_count_per_row <- rowSums(is.na(sirius))
to_del <- which(na_count_per_row >= 10)
tot_an <- intersect(sirius$mappingFeatureId[-to_del], df$`row ID`) %>% length()
round(tot_an/nrow(df)*100, 1)

###################################1
