# Modern plant n-alkane distributions
# Concentration-weighted mean n-alkane distributions by vegetation endmember
# Based on modern plant reference data (e.g., Santos et al., 2022)
#
# Input: Excel file with columns:
#   - Sample: plant species name
#   - Ecological form: vegetation type
#   - C17, C19, C21, C23, C25, C27, C29, C31, C33, C35: n-alkane concentrations (µg/g dw)
#
# Usage: set your file path in the `data_file` variable below, then run the script.

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

#  1. USER INPUT

data_file <- "your_data.xlsx"   # <-- path to your Excel file

# Define which samples belong to each endmember group
# Replace species names with those in your own dataset
group1_samples <- c("Agrostis", "Nardus", "Erica")       # EM1
group2_samples <- c("Juniperus", "Cytisus")               # EM2
group3_samples <- c("Antinoria", "Juncos")                # EM3

# 2. LOAD & PREPARE DATA 

chains <- c("C17", "C19", "C21", "C23", "C25", "C27", "C29", "C31", "C33", "C35")

df <- read_excel(data_file)

# Normalise each sample to relative abundance (sum = 1)
df_norm <- df %>%
  select(Sample, `Ecological form`, all_of(chains)) %>%
  rowwise() %>%
  mutate(total = sum(c_across(all_of(chains))),
         across(all_of(chains), ~ . / total)) %>%
  ungroup() %>%
  select(-total)

# Add total concentration per sample (used as weighting factor)
df_norm <- df_norm %>%
  left_join(
    df %>%
      select(Sample, all_of(chains)) %>%
      rowwise() %>%
      mutate(conc_total = sum(c_across(all_of(chains)))) %>%
      ungroup() %>%
      select(Sample, conc_total),
    by = "Sample"
  )

# 3. ASSIGN ENDMEMBER GROUPS 

df_norm <- df_norm %>%
  mutate(Group = case_when(
    Sample %in% group1_samples ~ "EM1 — Terrestrial grasses & heathland\n(Agrostis, Nardus, Erica)",
    Sample %in% group2_samples ~ "EM2 — Shrubs & trees\n(Juniperus, Cytisus)",
    Sample %in% group3_samples ~ "EM3 — Aquatic grasses\n(Antinoria, Juncos)",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group))

# 4. COMPUTE CONCENTRATION-WEIGHTED MEANS

df_mean <- df_norm %>%
  group_by(Group) %>%
  summarise(across(all_of(chains), ~ weighted.mean(., w = conc_total))) %>%
  pivot_longer(cols = all_of(chains),
               names_to = "Chain",
               values_to = "Abundance") %>%
  mutate(Chain = as.numeric(gsub("C", "", Chain)))

# 5. PLOT 

em_colors <- c(
  "EM1 — Terrestrial grasses & heathland\n(Agrostis, Nardus, Erica)" = "#E8A020",
  "EM2 — Shrubs & trees\n(Juniperus, Cytisus)"                        = "#2E7D32",
  "EM3 — Aquatic grasses\n(Antinoria, Juncos)"                        = "#5B9BD5"
)

df_mean$Group <- factor(df_mean$Group, levels = names(em_colors))

ggplot(df_mean, aes(x = factor(Chain), y = Abundance, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  facet_wrap(~ Group, nrow = 1) +
  scale_fill_manual(values = em_colors) +
  scale_y_continuous(limits = c(0, 0.6),
                     breaks = seq(0, 0.6, 0.1),
                     expand = c(0, 0)) +
  labs(x = "n-alkane carbon #",
       y = "Relative abundance",
       title = "Modern plant n-alkane distributions (Santos et al., 2022)") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 9, face = "bold"),
    strip.background = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    plot.title = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )


# Precipitation δD curve (δD_prc) reconstructed from NMF-derived terrestrial n-alkane δD
#
# Input: Excel file with columns:
#   - Age (col 1): sample age in cal yr BP
#   - dD_prc (‰ VSMOW): reconstructed precipitation δD (EM1 + EM2 weighted mean)
#   - dD_prc lo (ε-12): lower uncertainty bound (systematic ε bias −12‰)
#   - dD_prc hi (ε+12): upper uncertainty bound (systematic ε bias +12‰)
#
# Note: data are under review (manuscript and PANGAEA submission pending).
# This script is shared for reproducibility; input data will be made available
# upon acceptance.
#
# Usage: set your file path in `data_file` below, then run the script.

library(readxl)
library(ggplot2)

# 1. USER INPUT 

data_file <- "your_data.xlsx"   # <-- path to your Excel file

# Modern mean annual precipitation δD (reference line)
# Replace with the value appropriate for your study site
dD_MAP <- -56   # ‰ VSMOW

# 2. LOAD & PREPARE DATA 

df <- read_excel(data_file)

names(df)[1] <- "age"
df$dD_prc    <- df[["dD_prc (‰ VSMOW)"]]
df$dD_prc_lo <- df[["dD_prc lo (ε-12)"]]
df$dD_prc_hi <- df[["dD_prc hi (ε+12)"]]

# 3. PLOT 

ggplot(df, aes(x = age, y = dD_prc)) +
  # Systematic uncertainty band (ε ± 12‰)
  geom_ribbon(aes(ymin = dD_prc_lo, ymax = dD_prc_hi),
              fill = "#4C4CCC", alpha = 0.18) +
  # Modern MAP reference line
  geom_hline(yintercept = dD_MAP, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  annotate("text", x = 11500, y = dD_MAP + 1.2,
           label = "Modern \u03b4D MAP",
           hjust = 1, size = 3, colour = "grey40") +
  # Reconstruction curve
  geom_line(colour = "#3030CC", linewidth = 0.9) +
  # Sample points
  geom_point(shape = 16, size = 2.2, colour = "#3030CC") +
  scale_x_reverse(limits = c(12000, -1000),
                  breaks = seq(-1000, 12000, by = 1000),
                  expand = c(0, 0)) +
  labs(x = "Age (cal yr BP)",
       y = expression(delta * D[prc] ~ "(\u2030 VSMOW)"),
       title = expression("Reconstructed precipitation " * delta * D ~ "(EM1 + EM2)")) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12))


# NMF-derived terrestrial isotope curves (δ²H_terr and δ¹³C_terr)
# Concentration-weighted isotopic composition of terrestrial endmembers
# (EM1 + EM2) from NMF source unmixing of sediment n-alkane records.
#
# Input: Excel file with columns:
#   - Edad_calBP : sample age in cal yr BP
#   - dD         : δ²H of terrestrial endmembers (‰ VSMOW)
#   - d13C       : δ¹³C of terrestrial endmembers (‰ VPDB)
#
# Note: data are associated with a manuscript under review and a pending
# PANGAEA submission. This script is shared for reproducibility; input
# data will be made available upon acceptance.
#
# Usage: set your file path in `data_file` below, then run the script.

library(readxl)
library(ggplot2)

# 1. USER INPUT 

data_file <- "your_data.xlsx"   # <-- path to your Excel file

# 2. LOAD DATA 

df <- read_excel(data_file)

# 3. SHARED THEME 

x_scale <- scale_x_reverse(
  limits = c(12000, -1000),
  breaks = seq(-1000, 12000, by = 1000),
  expand = c(0, 0)
)

base_theme <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

# 4. δ²H PLOT (purple) 

p_dD <- ggplot(df, aes(x = Edad_calBP, y = dD)) +
  geom_line(colour = "#7B2D8B", linewidth = 0.9) +
  geom_point(colour = "#7B2D8B", shape = 16, size = 2) +
  x_scale +
  labs(
    x = "Age (cal yr BP)",
    y = expression(delta^2 * H[terr] ~ "(\u2030 VSMOW)")
  ) +
  base_theme

# 5. δ¹³C PLOT (green) 

p_d13C <- ggplot(df, aes(x = Edad_calBP, y = d13C)) +
  geom_line(colour = "#2E7D32", linewidth = 0.9) +
  geom_point(colour = "#2E7D32", shape = 16, size = 2) +
  x_scale +
  labs(
    x = "Age (cal yr BP)",
    y = expression(delta^13 * C[terr] ~ "(\u2030 VPDB)")
  ) +
  base_theme
















