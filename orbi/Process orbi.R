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

setwd("C:/Users/plyush/OneDrive - University of Haifa/Desktop/Data/Metabolomics Evaluation/Orbi_POS_DDA_MZML/mzMine_sirius/App/")

#................................................................
#### Compare Annotations Total----
#................................................................

target_df <- read.csv("annotation.csv")
mass <- calculateMass(target_df$Formula)
target_df$exactmass <- mass
target_df <- target_df[,c(1,2,3,4,6)]

df <- read_csv("orbi_iimn_gnps_quant.csv") %>% as.data.frame()
peakIn <- as.data.frame(cbind(mz = df$`row m/z`, rt = df$`row retention time`, id = df$`row ID`)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

#target_df <- read.csv("annot table - rt.csv")
#mass <- calculateMass(target_df$Formula)
#target_df$exactmass <- mass
#target_df <- target_df[,c(1,2,3,4)]


parm <- Mass2MzParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+"), tolerance = 0, ppm = 5) 
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.1)

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md <- na.omit(md)
unique(md$target_Compound)
unique(md$target_Compound) %>% length()
md <- subset(md, !md$target_Compound %in% c("Wulignan A1", "p-coumaraldehyde (Q27103652)|Phenylacrylic acid"))
cat("Wulignan A1 & p-coumaraldehyde (Q27103652)|Phenylacrylic acid are detected in blank")
unique(md$target_Compound) %>% length()

df <- read_csv("orbi_iimn_gnps_quant_filtered.csv") %>% as.data.frame()
peakIn <- as.data.frame(cbind(mz = df$`row m/z`, rt = df$`row retention time`, id = df$`row ID`)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

#target_df <- read.csv("annot table - rt.csv")
#mass <- calculateMass(target_df$Formula)
#target_df$exactmass <- mass
#target_df <- target_df[,c(1,2,3,4)]

#parm <- Mass2MzParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+"), tolerance = 0, ppm = 5) 
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.1)

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

df <- read_csv("orbi_iimn_gnps_quant_filtered.csv") %>% as.data.frame()
peakIn <- as.data.frame(cbind(mz = df$`row m/z`, rt = df$`row retention time`, id = df$`row ID`)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

target_df <- read.csv("annotation.csv")

mass <- calculateMass(target_df$Formula)
target_df$exactmass <- mass
target_df <- target_df[,c(1,2,3,4,6)]

# set polarity, adduct, accuracy
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.1)
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

# Plot2
assign_isomer_id <- function(rt, tol = 0.1) ({
  if (length(rt) == 0 || all(is.na(rt))) return(rep(NA_integer_, length(rt)))
  
  # Sort RTs to find the gaps, but keep track of original row order
  o <- order(rt)
  rt_sorted <- rt[o]
  
  # Create cluster IDs based on gaps larger than the tolerance
  cluster <- cumsum(c(1, diff(rt_sorted) > tol))
  
  # Return the cluster IDs mapped back to their original row order
  return(cluster[order(o)])
})

rt_tol <- 0.1

# 2. Add the isomer_id to the main dataframe
md2 <- md %>%
  mutate(rt = as.numeric(rt)) %>%
  filter(!is.na(target_Compound), !is.na(rt)) %>%
  group_by(target_Compound) %>%
  mutate(isomer_id = assign_isomer_id(rt, tol = rt_tol)) %>%
  ungroup()

# 3. Calculate Isomers: Grouped ONLY by Compound
# (Counts how many unique isomers exist per compound formula)
sum_rt <- md2 %>%
  group_by(target_Compound) %>%
  summarise(detect = n_distinct(isomer_id), .groups = "drop") %>%
  mutate(metric = "Isomers")

# 4. Calculate Adducts: Grouped by Compound AND Isomer
# (Counts distinct adducts for EACH specific isomer separately)
sum_adduct <- md2 %>%
  filter(!is.na(adduct)) %>%
  group_by(target_Compound, isomer_id) %>%
  summarise(detect = n_distinct(adduct), .groups = "drop") %>%
  # Drop the isomer_id column now so it binds cleanly with sum_rt later
  select(-isomer_id) %>% 
  mutate(metric = "Adducts")

sum_adduct_tot_app <- tibble(adduct = md2$adduct) %>%
  dplyr::count(adduct, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  arrange(desc(n))

# Combine for plotting exactly like you had it
plot_df <- bind_rows(sum_adduct, sum_rt) %>%
  mutate(one = 1)

det_levels <- sort(unique(plot_df$detect))
det_levels <- det_levels[is.finite(det_levels)]
det_levels_chr <- as.character(det_levels)

#pal_det <- setNames(ggsci::pal_npg()(length(det_levels_chr)), det_levels_chr)
pal_det <- setNames(viridisLite::viridis(length(det_levels_chr)), det_levels_chr)

plot_df <- plot_df %>%
  mutate(detect_f = factor(as.character(detect), levels = det_levels_chr))

plot_df2 <- plot_df %>%
  mutate(detect_num = as.numeric(as.character(detect_f))) %>%
  arrange(metric, detect_num, target_Compound, ) %>%   # sort within each metric
  select(-detect_num)

max_y <- plot_df2 %>%
  count(metric, name = "n_compounds") %>%
  pull(n_compounds) %>%
  max()
step <- max(1, ceiling(max_y / 6))

plot_df2_app <- plot_df2

ggplot(plot_df2, aes(x = metric, y = 1, fill = detect_f)) +
  geom_col(color = "black", width = 0.65, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = pal_det, drop = FALSE, name = "Occurrence: ") +
  labs(x = NULL, y = "Compounds") +
  scale_y_continuous(
    limits = c(0, max_y),
    breaks = seq(0, max_y, by = step),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom") + labs(caption = paste0("Total number of detected compounds including possible isomers: ", nrow(distinct(md[,c(2,4)]))))

#................................................................
#### Annotation Raw ----
#................................................................

df <- read_csv("orbi_iimn_gnps_quant.csv") %>% as.data.frame()
peakIn <- as.data.frame(cbind(mz = df$`row m/z`, rt = df$`row retention time`, id = df$`row ID`)) # "mz" column name necessary
peakIn$mz <- as.numeric(peakIn$mz)
peakIn$rt <- as.numeric(peakIn$rt)

target_df <- read.csv("annotation.csv")

mass <- calculateMass(target_df$Formula)
target_df$exactmass <- mass
target_df <- target_df[,c(1,2,3,4,6)]

# set polarity, adduct, accuracy
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.1)
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

# Plot2
assign_isomer_id <- function(rt, tol = 0.1) ({
  if (length(rt) == 0 || all(is.na(rt))) return(rep(NA_integer_, length(rt)))
  
  # Sort RTs to find the gaps, but keep track of original row order
  o <- order(rt)
  rt_sorted <- rt[o]
  
  # Create cluster IDs based on gaps larger than the tolerance
  cluster <- cumsum(c(1, diff(rt_sorted) > tol))
  
  # Return the cluster IDs mapped back to their original row order
  return(cluster[order(o)])
})

rt_tol <- 0.1

# 2. Add the isomer_id to the main dataframe
md2 <- md %>%
  mutate(rt = as.numeric(rt)) %>%
  filter(!is.na(target_Compound), !is.na(rt)) %>%
  group_by(target_Compound) %>%
  mutate(isomer_id = assign_isomer_id(rt, tol = rt_tol)) %>%
  ungroup()

# 3. Calculate Isomers: Grouped ONLY by Compound
# (Counts how many unique isomers exist per compound formula)
sum_rt <- md2 %>%
  group_by(target_Compound) %>%
  summarise(detect = n_distinct(isomer_id), .groups = "drop") %>%
  mutate(metric = "Isomers")

# 4. Calculate Adducts: Grouped by Compound AND Isomer
# (Counts distinct adducts for EACH specific isomer separately)
sum_adduct <- md2 %>%
  filter(!is.na(adduct)) %>%
  group_by(target_Compound, isomer_id) %>%
  summarise(detect = n_distinct(adduct), .groups = "drop") %>%
  # Drop the isomer_id column now so it binds cleanly with sum_rt later
  select(-isomer_id) %>% 
  mutate(metric = "Adducts")

sum_adduct_tot_raw <- tibble(adduct = md2$adduct) %>%
  dplyr::count(adduct, name = "n") %>%
  mutate(percent = 100 * n / sum(n)) %>%
  arrange(desc(n))

# Combine for plotting exactly like you had it
plot_df <- bind_rows(sum_adduct, sum_rt) %>%
  mutate(one = 1)

det_levels <- sort(unique(plot_df$detect))
det_levels <- det_levels[is.finite(det_levels)]
det_levels_chr <- as.character(det_levels)

#pal_det <- setNames(ggsci::pal_npg()(length(det_levels_chr)), det_levels_chr)
pal_det <- setNames(viridisLite::viridis(length(det_levels_chr)), det_levels_chr)

plot_df <- plot_df %>%
  mutate(detect_f = factor(as.character(detect), levels = det_levels_chr))

plot_df2 <- plot_df %>%
  mutate(detect_num = as.numeric(as.character(detect_f))) %>%
  arrange(metric, detect_num, target_Compound, ) %>%   # sort within each metric
  select(-detect_num)

max_y <- plot_df2 %>%
  count(metric, name = "n_compounds") %>%
  pull(n_compounds) %>%
  max()
step <- max(1, ceiling(max_y / 6))

plot_df2_raw <- plot_df2

ggplot(plot_df2, aes(x = metric, y = 1, fill = detect_f)) +
  geom_col(color = "black", width = 0.65, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = pal_det, drop = FALSE, name = "Occurrence: ") +
  labs(x = NULL, y = "Compounds") +
  scale_y_continuous(
    limits = c(0, max_y),
    breaks = seq(0, max_y, by = step),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom") + labs(caption = paste0("Total number of detected compounds including possible isomers: ", nrow(distinct(md[,c(2,4)]))))

#................................................................
# Combine App & Raw ----
#................................................................

add_tot <- rbind(cbind(sum_adduct_tot_app, Label = "App"), cbind(sum_adduct_tot_raw, Label = "Raw")) %>% as.data.frame()

library(cowplot)
plot_df2_cp <- rbind(cbind(plot_df2_raw, Label = "Raw"), cbind(plot_df2_app, Label = "App"))

p2<- ggplot(subset(plot_df2_cp, plot_df2_cp$metric == "Isomers"), aes(x = detect_f, fill = Label)) +
  geom_bar(color = "black", size = 1, width = 0.8, position = position_dodge(width = 1, preserve = "single"), stat="count") +
  scale_fill_npg() +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), n.breaks = 5) +
  labs(x = "Isomers per compound", y = "Count") +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom")  

p1 <- ggplot(subset(plot_df2_cp, plot_df2_cp$metric == "Adducts"), aes(x = detect_f, fill = Label)) +
  geom_bar(color = "black", size = 1, width = 0.8, position = position_dodge(width = 1, preserve = "single"), stat="count") +
  scale_fill_npg() +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), n.breaks = 7) +
  labs(x = "Adducts per compound", y = "Count") +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom") +theme(legend.position = c(0.85, 0.85), #legend.justification = c(1, 1), 
           # legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
           legend.text = element_text(size = 18),         # Makes the group names larger
        legend.key.size = unit(1.0, "cm"), #axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_blank()) 

shared_legend <- get_legend(
  p1 + theme(legend.box.margin = margin(0, 0, 0, -150)) # Adds a little breathing room
)

p3 <- ggplot(add_tot, aes(x = reorder(adduct, -n), y = n, fill = Label)) +
  geom_col(color = "black", position = position_dodge(width = 1), width = 0.8, linewidth = 0.8, size = 1) +
  labs(x = "", y = "Count", fill = "Data Source") +
  expand_limits(y = max(add_tot$n) * 1.05) + 
  theme_classic(base_size = 16) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), n.breaks = 7) +
  theme(
    legend.position = "bottom") +
  scale_fill_npg()

p1 <- p1 + theme(legend.position = "none")+ coord_flip()
p2 <- p2 + theme(legend.position = "none")+ coord_flip()
p3 <- p3 + theme(legend.position = "none")+ coord_flip()

plot_row <- plot_grid(
  p1, 
  p3, 
  p2,
  rel_widths = c(2.5, 2.5, 1.5), 
  labels = c('A', 'B', 'C'), 
  label_size = 25, 
  nrow = 1
)

plot_grid(
  plot_row, 
  shared_legend, 
  rel_widths = c(4, .01), 
  nrow = 1
)

plot_row <- plot_grid(
  p1, 
  p3, 
  rel_widths = c(1,1), 
  labels = c('A', 'B'), 
  label_size = 25, 
  nrow = 1
)

plot_grid(
  plot_row, 
  shared_legend, 
  rel_widths = c(4, .01), 
  nrow = 1
)

#................................................................
# Final detected annotated compounds ----
#................................................................
# 92 unique + (82 by 1 isomer, 9 - 2, 1-3)
cmp <- subset(plot_df2_cp, plot_df2_cp$metric == "Isomers")
total_features_by_label <- cmp %>%
  group_by(Label) %>%
  summarise(total_features = sum(detect, na.rm = TRUE))

print(total_features_by_label)

#................................................................
#### Cytoscape ----
#................................................................

library(RCy3)
cytoscapePing()
cytoscapeVersionInfo()

# Steroids
openSession("app single.cys")
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
openSession("app single.cys")
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

peakIn <- cbind(mz = as.numeric(sirius$ionMass), rt = as.numeric(sirius$retentionTimeInMinutes), map_id = sirius$mappingFeatureId) %>% as.data.frame()

# set polarity, adduct, accuracy
parm <- Mass2MzParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+"), tolerance = 0, ppm = 5) 
parm <- Mass2MzRtParam(adducts = c("[M+H]+", "[M+2H]2+", "[M+K]+", "[M+Na]+", "[M+NH4]+", "[M+2Na]2+", "[M+H+K]2+", "[M+H+Na]2+", "[M+2Na-H]+", "[M+H2O+H]+", "[M+2K-H]+", "[M+H-H2O]+"), 
                       tolerance = 0, ppm = 5, toleranceRt = 0.1)
#MetaboCoreUtils::adducts() 

matched_features <- matchValues(peakIn, target_df, parm)
md <- matchedData(matched_features)
md <- as.data.frame(md)
md$MetaboAnnotation_status <- ifelse(is.na(md$target_Compound), F, T)
md <- md[,c(3:4,10, 14)]

sirius_ma <- left_join(sirius, md, by = c("mappingFeatureId" = "map_id"))
write_csv(sirius_ma, "canopus_formula_summary + MetaboAnnotation RT.csv")

###################################1