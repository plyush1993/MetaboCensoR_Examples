# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# ==============================================================================
# ALGORITHM VISUALIZATION ----
# ==============================================================================

# ==============================================================================
# 1. SETUP: Scenario Definition (Neutral Mass = 300 Da, Base RT = 5.0 min) 
# ==============================================================================
mz_tol_da <- 0.01
rt_tol_min <- 0.05
ref_rt <- 5.0

# Define features belonging to the COMPOUND GRAPH
compound_features <- tibble::tibble(
  label = c("M+H", "M+Na", "M+NH4", "13C"),
  # Calculate theoretical masses
  mz_theo = c(301.0073, 322.9892, 318.0338, 302.0106),
  # Add slight random error to mz and rt to make it realistic
  mz = mz_theo + runif(4, -0.002, 0.002),
  rt = ref_rt  + c(0.0, 0.008, 0.018, -0.008),
  # Relative intensities (for point size)
  intensity = c(100, 40, 15, 10),
  group = c("Cluster", "Target Adduct (M+Na)", "Cluster", "Cluster")
)

# Define scattered BACKGROUND NOISE features (Singletons)
noise_features <- tibble::tibble(
  label = paste0("S", 1:5),
  mz = c(305.1, 322.9892, 330.0, 310.2, 298.8),
  rt = c(4.965, 5.015, 4.99, 4.975, 5.015),
  intensity = c(25, 68, 5, 12, 8),
  group = c("Noise", "Target False Adduct", "Noise", "Noise", "Noise")
)

# Combine datasets
all_features <- bind_rows(compound_features, noise_features)

# ==============================================================================
# 2. DEFINE EDGES: The connections forming the graph
# ==============================================================================
# We define edges between the precursor (M+H) and other valid adducts
edges <- tibble::tibble(
  # Find coordinates for M+H
  x     = all_features$rt[all_features$label == "M+H"],
  y     = all_features$mz[all_features$label == "M+H"],
  # Ends connect to other compound labels
  xend  = compound_features$rt[2:4],
  yend  = compound_features$mz[2:4],
  type  = c("Adduct Shift", "Adduct Shift", "Isotope Shift")
)

# ==============================================================================
# 3. DEFINE TOLERANCE WINDOW: Visualize the algorithm's decision space
# ==============================================================================
# Let's show the search window around the expected M+Na location
target_add_mz <- 322.9892 # Theoretical M+Na
win_rt_min <- ref_rt - rt_tol_min
win_rt_max <- ref_rt + rt_tol_min
win_mz_min <- target_add_mz - mz_tol_da
win_mz_max <- target_add_mz + mz_tol_da

# ==============================================================================
# 4. PLOTTING
# ==============================================================================

final_plot <- ggplot() +
  
  # --- Step 1: Tolerance Window ---
  # Draw acceptable search region around expected M+Na
  annotate("rect", xmin = 4.975, xmax = 5.025, 
           ymin = win_mz_min, ymax = win_mz_min, 
           fill = "#3498db", alpha = 0.15, color = "#2980b9", linetype = "dashed", size = 1) +
  # Label the tolerance constraints
  annotate("text", x = 5.025, y = win_mz_max - 0.005, 
           label = "m/z Tol: ±0.005 Da\nRT Tol: ±0.05 min", color = "#2980b9", 
           size = 3.5, hjust = 0, fontface = "bold") +

  # --- Step 2: Connections (Edges) ---
  # Only drawn between validated clustered features
  geom_segment(data = edges, aes(x=x, y=y, xend=xend, yend=yend, linetype = type), 
               color = "gray40", size = 0.8) +

  # --- Step 3: Shift Arrows ---
  # Visual explanation of mass shifts
  geom_curve(aes(x = 5, xend = 4.997, y = 301.1, yend = 322.8), linewidth = 0.8,
             curvature = -0.1, arrow = arrow(length = unit(0.2, "cm")), color = "#16a085") +
  annotate("text", x = 4.980, y = 315, label = "Δ m/z = 21.9819\n(H → Na shift)", 
           color = "#16a085", size = 3.5, hjust = 0, fontface = "bold") +
  
  # --- Step 4: Features (Nodes) ---
  # Size by intensity, color/shape by cluster membership
  geom_point(data = all_features, aes(x = rt, y = mz, fill = group, size = intensity), 
             shape = 21, color = "black", stroke = 1) +
  
  # --- Step 5: Labels & Aesthetics ---
  # Annotate adduct types using ggrepel to prevent overlap
  geom_text_repel(data = compound_features, aes(x = rt, y = mz, label = label), 
                  size = 4.5, fontface = "bold",
                  box.padding = 0.8,       # Space around the text box
                  point.padding = 0.6,     # Space between point and line
                  min.segment.length = 0,  # Always draw the line segment connecting text to point
                  nudge_y = -1.5,             # Give them a slight upward nudge initially
                  seed = 1610) +             # Set a seed so the random placement is reproducible for your manuscript!
  
  # Define scales
  scale_size_continuous(range = c(2, 9), guide = "none") + # Map intensity to size
  scale_fill_manual(values = c(
    "Cluster" = "#A3E4D7",                # Faded light teal
    "Target Adduct (M+Na)" = "#1abc9c",   # VIBRANT TEAL
    "Noise" = "#F5B7B1",                  # Faded light red
    "Target False Adduct" = "#e74c3c"     # VIBRANT RED
  )) +
  scale_linetype_manual(values = c("Adduct Shift" = "solid", "Isotope Shift" = "dotted")) +
  
  # Final touches (Manuscript Ready)
  theme_classic(base_size = 14) +
  coord_cartesian(xlim = c(4.96, 5.04), ylim = c(295, 330)) +
  labs(
    title = "mz-rt match",
    #subtitle = "Grouping related features by characteristic shifts within m/z-RT tolerance windows",
    x = "Retention Time (min)",
    y = "m/z"
  ) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold.italic", size = 20),
    plot.subtitle = element_text(size = 11, color = "gray30")
  )

# ==============================================================================
# 5. VIEW AND SAVE
# ==============================================================================
print(final_plot)
a <- final_plot




library(ggplot2)
library(dplyr)

# ==============================================================================
# Helper Function: Generate a variable with an EXACT correlation to another 
# ==============================================================================
force_correlation <- function(x, target_r, target_mean, target_sd) {
  # 1. Standardize the input variable
  x_std <- scale(x)[, 1]
  
  # 2. Generate random noise and make it perfectly orthogonal (uncorrelated) to x
  noise <- rnorm(length(x))
  noise_resid <- residuals(lm(noise ~ x_std))
  noise_std <- scale(noise_resid)[, 1]
  
  # 3. Combine them mathematically to hit the exact target_r
  y_std <- target_r * x_std + sqrt(1 - target_r^2) * noise_std
  
  # 4. Scale to your desired mean and standard deviation
  y_final <- (y_std * target_sd) + target_mean
  
  # 5. LC-MS intensities can't be negative, so we floor it at 0 just in case
  y_final <- pmax(y_final, 0)
  
  return(y_final)
}

# ==============================================================================
# 1. Simulate Feature Intensities Across 40 Samples
# ==============================================================================
set.seed(1234) # Keeps the random noise consistent every time you run it
n_samples <- 40

# 1. Base Precursor (M+H)
mh_intensity <- runif(n_samples, min = 10000, max = 100000)

# 2. True Adduct (M+Na): 
# Target r = 0.98. Mean is ~40% of M+H mean.
mna_intensity <- force_correlation(
  x = mh_intensity, 
  target_r = 0.9, 
  target_mean = mean(mh_intensity) * 0.40, 
  target_sd = 4000
)

# 3. Co-eluting Noise: 
# Target r = 0.10 (Basically no correlation). Random mean and spread.
noise_intensity <- force_correlation(
  x = mh_intensity, 
  target_r = 0.10, 
  target_mean = 35000, 
  target_sd = 4000
)

cor(mh_intensity, mna_intensity)   
cor(mh_intensity, noise_intensity) 

# Combine into a long-format dataframe for ggplot
df_scatter <- data.frame(
  MH_Intensity = rep(mh_intensity, 2),
  Target_Intensity = c(mna_intensity, noise_intensity),
  Feature = factor(rep(c("Target Adduct (M+Na)", "Target False Adduct"), each = n_samples),
                   levels = c("Target Adduct (M+Na)", "Target False Adduct"))
)

# ==============================================================================
# 2. Plotting the Scatterplot
# ==============================================================================
p_scatter <- ggplot(df_scatter, aes(x = MH_Intensity, y = Target_Intensity, 
                                    fill = Feature, color = Feature)) +
  
  # 1. Add linear regression trendlines with shaded confidence intervals (se = TRUE)
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linetype = "dashed", linewidth = 1) +
  
  # 2. Add the scatter points (shape = 21 allows both fill and a black border)
  geom_point(shape = 21, size = 4.5, color = "black", stroke = 1, alpha = 0.8) +
  
  # 3. Add statistical annotations directly onto the plot
  annotate("text", x = 70000, y = 16000, hjust = 0, size = 3.5, fontface = "bold", color = "#16a085",
           label = "r > 0.80\n(Valid M+Na)") +
           
  annotate("text", x = 14000, y = 43000, hjust = 0, size = 3.5, fontface = "bold", color = "#c0392b",
           label = "r < 0.8\n(Rejected)") +

  # 4. Apply your specific color theme
  scale_fill_manual(values = c("Target Adduct (M+Na)" = "#1abc9c", "Target False Adduct" = "#e74c3c")) +
  scale_color_manual(values = c("Target Adduct (M+Na)" = "black", "Target False Adduct" = "black")) +
  
  # 5. Manuscript formatting
  theme_classic(base_size = 14) +
  labs(
    title = "Correlation match",
    subtitle = "",
    x = "Precursor (M+H) Intensity",
    y = "Feature Intensity"
  ) +
  theme(
    plot.title = element_text(face = "bold.italic", size = 20),
    legend.position = "none",
    legend.title = element_blank(),
    plot.subtitle = element_text(color = "gray30")
  )

# ==============================================================================
# 3. View the Plot
# ==============================================================================
print(p_scatter)
b <- p_scatter



# combine ----
library(cowplot)
plot_row <- plot_grid(
  a, 
  b, 
  labels = c('', ''), 
  #rel_widths = c(2.5, 1.5), 
  label_size = 25, 
  nrow = 1
)
plot_row

plot_row <- plot_grid(a, b, ncol = 2, align = 'h', axis = 'b')

# Define the coordinates for the central zoom polygon
zoom_coords <- data.frame(
  # X coordinates: [TopLeft(PlotA), TopRight(PlotB), BottomRight(PlotB), BottomLeft(PlotA)]
  x = c(0.49, 0.53, 0.53, 0.49), 
  # Y coordinates: [TopEdge(PlotA), TopEdge(PlotB), BottomEdge(PlotB), BottomEdge(PlotA)]
  y = c(0.72, 0.78, 0.16, 0.67)
)

final_horizontal_zoom <- ggdraw(plot_row) +
  
  # The shaded polygon bridging the middle gap
  geom_polygon(data = zoom_coords, aes(x = x, y = y), 
               fill = "#2980b9", alpha = 0.15) 
  
  # Top border line: Connects TopLeft (0.72) to TopRight (0.80)
#  draw_line(x = c(0.50, 0.53), y = c(0.72, 0.80), 
 #           color = "#2980b9", linewidth = 0.8, linetype = "dashed") +
  
  # Bottom border line: Connects BottomLeft (0.65) to BottomRight (0.16)
  # THE FIX: Changed the first y coordinate from 0.72 to 0.65
 # draw_line(x = c(0.50, 0.53), y = c(0.65, 0.16), 
  #          color = "#2980b9", linewidth = 0.8, linetype = "dashed")

# IMPORTANT: Make sure your RStudio plot window is pulled WIDE before running this line!
print(final_horizontal_zoom)

# ==============================================================================
# SPECTRA PLOT NEW PEPTIDIC PRODUCT ----
# ==============================================================================

# Annotation of new peptidic product ----
library(tidyverse)

# 1) Paste spectrum here
spec_txt <- "
70.066078 31.339817
84.043587 350.908539
84.046219 75.183083
86.095894 6808.047363
87.097610 45.304710
87.099564 84.926041
87.101440 63.816311
102.054718 859.977539
103.059036 10.554861
113.071083 61.380577
116.069366 13.964893
120.080208 47.090916
129.065399 214.020111
129.068527 26.143578
130.067459 11.853921
132.101151 5351.151855
132.116333 38.971794
133.097702 9.255801
133.104141 199.892822
133.107300 22.084015
141.065033 243.086563
152.109558 8.606271
155.081116 9.418183
166.084991 12.341068
169.131180 10.392478
170.115723 156.699097
197.127335 5551.207520
197.143463 53.099068
197.148865 29.391228
198.112915 18.024454
198.125870 13.477745
198.131546 433.398834
199.180984 32.476494
200.173401 8.443889
201.121674 9.742949
215.138702 10000.000000
215.156113 46.603771
215.163437 19.973043
215.170654 11.204391
216.143356 783.008301
225.122345 1313.349487
225.145889 11.691538
226.126938 53.748596
227.084213 10.717243
227.172012 8.119123
243.132263 1333.972046
243.153198 11.366774
244.135269 25.169285
245.184052 1904.259277
245.198608 65.927284
245.204636 11.691538
246.187988 103.275253
248.645264 14.289658
255.169434 9.093418
261.088684 13.315363
272.131317 8.281507
310.213226 14.289658
314.652557 8.768654
338.208557 10.554861
356.221680 25.169285
374.144531 7.956741
375.187164 21.109722
"

spec <- readr::read_table2(I(spec_txt), col_names = c("mz", "int")) %>%
  mutate(rel_int = 100 * int / max(int, na.rm = TRUE)) %>%
  arrange(mz)

# 2) Define expected ELL fragments
# --- monoisotopic constants ---
H    <- 1.00727646688
H2O  <- 18.010564684

# residue masses (monoisotopic)
E <- 129.042593
L <- 113.084064

# --- add b1 and y3 here ---
frags <- tribble(
  ~frag,              ~mz_theor,
  "b1 (E + H)",        E + H,                         # 130.0498695
  "b2 (EL + H)",      E + L + H,                     # 243.1339335
  "b3 (ELL + H)",    E + 2*L + H,                   # 356.2179975
  "y1 (L + H2O + H)",  L + H2O + H,                   # 132.1019052
  "y2 (LL + H2O + H)", 2*L + H2O + H,                 # 245.1859692
  "y3 (ELL + H2O + H)",E + 2*L + H2O + H              # 374.2285622  (full-length y-ion)
)

# 3) Match expected fragments to observed peaks
# Choose ONE tolerance mode:
tol_mode <- "ppm"     # "ppm" or "da"
tol_ppm  <- 20        # e.g., 5–20 ppm
tol_da   <- 0.02      # e.g., 0.01–0.05 Da

match_one <- function(mz0, spec, tol_mode = "ppm", tol_ppm = 10, tol_da = 0.02) {

  idx <- which.min(abs(spec$mz - mz0))

  mz_obs   <- spec$mz[idx]
  diff_da  <- mz_obs - mz0
  diff_ppm <- diff_da / mz0 * 1e6

  hit <- if (tolower(tol_mode) == "ppm") {
    abs(diff_ppm) <= tol_ppm
  } else {
    abs(diff_da) <= tol_da
  }

  tibble(
    mz_obs   = mz_obs,
    diff_da  = diff_da,
    diff_ppm = diff_ppm,
    hit      = hit,
    int      = spec$int[idx],
    rel_int  = spec$rel_int[idx]
  )
}

hits <- frags %>%
  mutate(
    match = purrr::map(
      mz_theor,
      match_one,
      spec = spec,
      tol_mode = tol_mode,
      tol_ppm = tol_ppm,
      tol_da = tol_da
    )
  ) %>%
  tidyr::unnest(match)

library(dplyr)
library(ggplot2)
library(ggrepel)
library(plotly)

# assumes you already have:
# spec: tibble/data.frame with columns mz, int
# hits: tibble/data.frame with columns frag, mz_theor, mz_obs, int, hit, diff_ppm (at least)

top_n_labels <- 10  # how many most abundant peaks to label with accurate mass

# --- label top abundant peaks with accurate m/z
top_peaks <- spec %>%
  slice_max(order_by = int, n = top_n_labels, with_ties = FALSE) %>%
  mutate(mz_lab = sprintf("%.4f", mz))

# --- label matched fragment hits (b/y)
hits_lab <- hits %>%
  filter(hit) %>%
  mutate(
    mz_obs = as.numeric(mz_obs),
    int    = as.numeric(int),
    #frag_lab = paste0(frag, "\n", sprintf("%.5f", mz_obs), "\n", sprintf("%.1f ppm", diff_ppm))
    frag_lab = paste0(frag, "\n", sprintf("%.4f", mz_obs), "\n", sprintf("%.1f ppm", diff_ppm))
  )

# --- base ggplot
gg <- ggplot(spec, aes(x = mz, y = int)) +
  geom_segment(aes(xend = mz, y = 0, yend = int,
                   text = sprintf("m/z: %.5f<br>Int: %.2f", mz, int))) +
  geom_point(aes(text = sprintf("m/z: %.5f<br>Int: %.2f", mz, int)), size = 2.5, shape = 21, fill = "black", color = "black") +

  # label top peaks (accurate mass)
  geom_text_repel(
    data = top_peaks,
    aes(label = mz_lab),
    size = 4.3,
    min.segment.length = 0,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15
  ) +

  # optional: mark hits with an extra point layer
  geom_point(
    data = hits_lab,
    aes(x = mz_obs, y = int,
        text = paste0(frag, "<br>m/z obs: ", sprintf("%.5f", mz_obs),
                      "<br>ppm: ", sprintf("%.1f", diff_ppm))),
    size = 5, fill = "firebrick3", shape = 21, color = "black"
  ) +

  # label hits (b/y)
  geom_label_repel(
    data = hits_lab,
    aes(x = mz_obs, y = int, label = frag_lab),
    size = 4.1,
    min.segment.length = 0,
    max.overlaps = Inf,
    box.padding = 0.8,       # Increased so labels don't overlap each other
    point.padding = 1.8,     # Decreased so labels aren't thrown too far away
    nudge_y = 500,         # (Optional) Try uncommenting this to push ALL labels up by a fixed amount
    label.size = 0.5,
    color = "firebrick3",        
    fill = "#ffffff"
  ) +

  theme_classic(base_size = 14) +
  labs(x = "m/z", y = "Intensity")

# convert to plotly and keep print(p)
#p <- ggplotly(gg, tooltip = "text")
#print(p)
gg

# Annotation of new peptidic product 2 ----

# --- base ggplot
gg <- ggplot(spec, aes(x = mz, y = int)) +
  
  # 1. Base peaks (gray and slightly thinner so the hits stand out more)
  geom_segment(aes(xend = mz, y = 0, yend = int,
                   text = sprintf("m/z: %.4f<br>Int: %.2f", mz, int)),
               color = "gray40", linewidth = 0.6) +
  
  # 2. Highlighted HIT peaks (FIXED: Added x = mz_obs and inherit.aes = FALSE)
  geom_segment(data = hits_lab,
               aes(x = mz_obs, xend = mz_obs, y = 0, yend = int),
               color = "firebrick3", linewidth = 1.2, inherit.aes = FALSE) +
  
  # 3. Points for hits
  geom_point(data = hits_lab,
             aes(x = mz_obs, y = int,
                 text = paste0(frag, "<br>m/z obs: ", sprintf("%.4f", mz_obs),
                               "<br>ppm: ", sprintf("%.1f", diff_ppm))),
             size = 3.5, fill = "firebrick3", shape = 21, color = "black", inherit.aes = FALSE) +

  # 4. Points for base peaks (only show dots for the top ones to reduce clutter)
  geom_point(data = top_peaks %>% filter(!mz %in% hits_lab$mz_obs), 
             aes(x = mz, y = int, text = sprintf("m/z: %.4f<br>Int: %.2f", mz, int)), 
             size = 2.5, shape = 21, fill = "black", color = "black") +

  # 5. Label hits (b/y ions) - Forced straight UP!
  geom_label_repel(
    data = hits_lab,
    aes(x = mz_obs, y = int, label = frag_lab),
    size = 3.8,
    direction = "y",                         
    nudge_y = max(spec$int) * 0.12,          
    min.segment.length = 0,
    segment.color = "firebrick3",            
    segment.linetype = 2,                    
    box.padding = 0.5,
    label.padding = 0.3,
    label.size = 0.5,
    color = "firebrick3",        
    fill = "#ffffff",
    fontface = "bold",
    inherit.aes = FALSE                      # Protects against missing global columns
  ) +
  
  # 6. Label top peaks (accurate mass) - Excludes hits to prevent overlap!
  geom_text_repel(
    data = top_peaks %>% filter(!mz %in% hits_lab$mz_obs), 
    aes(x = mz, y = int, label = mz_lab),
    size = 3.5,
    direction = "y",                         
    nudge_y = max(spec$int) * 0.05,
    min.segment.length = 0,
    segment.color = "gray50",
    color = "black"
  ) +

  # 7. Clean theme & expand Y axis so labels don't get cut off at the top
  theme_classic(base_size = 14) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + 
  labs(x = "m/z", y = "Intensity")

# convert to plotly and keep print(p)
# p <- ggplotly(gg, tooltip = "text")
# print(p)
gg

####################################################################1
# Plot Annotation Results ----
####################################################################1

# 1. read data
df <- data.frame(
  Software = c("nontarget", "nontarget", "mzMine", "mzMine", "MS1FA", "MS1FA", "CAMERA", "CAMERA", "CAMERA", "CAMERA"),
  Tool = c("nontarget", "MetaboCensoR", "mzMine", "MetaboCensoR", "MS1FA", "MetaboCensoR", "CAMERA", "MetaboCensoR", "CAMERA", "MetaboCensoR"),
  Label = c("Tool", "App", "Tool", "App", "Tool", "App", "Tool", "App", "Tool", "App"),
  Type = c("Adducts", "Adducts", "Adducts", "Adducts", "NL", "NL", "Adducts", "Adducts", "Isotopes", "Isotopes"),
  Dataset = c("inter", "inter", "orbi", "orbi", "orbi", "orbi", "folate", "folate", "folate", "folate"),
  total = c(523, 450, 2420, 2479, 33, 33, 957, 1119, 796, 924),
  relative = c(1.00, 0.86, 0.98, 1.00, 1.00, 1.00, 0.86, 1.00, 0.86, 1.00)
)

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggsci)

# 2. Ensure factors are ordered nicely for the plot
df$Label <- factor(df$Label, levels = c("Tool", "App"))
df$Dataset <- factor(df$Dataset, levels = c("orbi", "folate", "inter"))
df$Tool <- factor(df$Tool, levels = c("MS1FA", "mzMine", "CAMERA", "nontarget", "MetaboCensoR"))

# The '\n' adds a line break so the header labels look clean
df$Type <- paste0("(", df$Type, ")")
df$Comparison <- factor(df$Type, levels = c("(Adducts)", "(NL)", "(Isotopes)"))

# 3. Create the corrected plot
ggplot(df, aes(x = Tool, y = relative, fill = Label)) +
  
  geom_bar(stat = "identity", position = position_dodge(width = 1.2), size = 1, width = 0.75, color = "black") +
  
  geom_text(aes(label = total), 
            position = position_dodge(width = 0.95), 
            vjust = -0.5, size = 5, alpha = 1, show.legend = FALSE) +
  
  # FIX: Facet by Dataset AND the new Comparison column
  facet_grid(~ Dataset + Comparison, scales = "free_x", space = "free_x") +
  
  scale_fill_manual(values = rev(ggsci::pal_npg()(2))) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +

  theme_classic(base_size = 16) +
  labs(
    title = "",
    x = "",
    y = "Relative Scope",
    fill = "Feature Type",
    alpha = "Pipeline"
  ) +
  theme(legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    panel.spacing.x = unit(0.2, "lines") 
  )
#...................................................