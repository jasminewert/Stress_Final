---
title: "HMZ_final"
author: "Jasmin Ewert"
date: "`r format(Sys.Date(), '%B %d, %Y')`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Categorical covariates
## Change NW package (mgp)
## Choose colors (heatmap, violin plots)
## Repeat sensitivity analysis


## 1. Load Required Libraries

```{r load-libraries, echo=FALSE, message=FALSE, warning=FALSE}
# Data Organization
library(tidyverse)
library(dplyr)
library(readr)
library(tidyr)

# Plotting
library(ggplot2)
library(ggpubr)
library(cowplot)
library(purrr)

# Network analysis
library(qgraph)
library(bootnet)
```

## 2. Import & Data Organization

### 2.1 Import

```{r load-data, echo=FALSE, message=FALSE, warning=FALSE}
Stress_Masterfile <- read_delim("Stress_Masterfile.csv", 
                                delim = ";", escape_double = FALSE, trim_ws = TRUE)
```

### 2.2 Remove Empty Columns

```{r remove-na-rows, results='hide', message=FALSE, warning=FALSE}
Stress_Masterfile <- Stress_Masterfile %>%
  filter(rowSums(is.na(.)) != ncol(.))
```

### 2.3 Identify and Clean Analyte Columns

```{r clean-analyte-columns, results='hide', message=FALSE, warning=FALSE}
# Identify analyte columns (adjust indices if needed)
analyte_columns <- colnames(Stress_Masterfile)[62:161]

# Replace spaces in column names with underscores
colnames(Stress_Masterfile) <- gsub(" ", "_", colnames(Stress_Masterfile))
analyte_columns <- gsub(" ", "_", analyte_columns)
```

### 2.4 Keep Only Complete Bio_IDs

TOBIAS: Code could be made simpler - We don't need to show the filtered ones // Also this is only Analytes, next is Psychometrics - either combine or specifiy the title -> If we need a flow chart, we should combine all of the IDs we exclude into a list (so we 
can use them for a flow chart)

#### 2.4.1 Analyte Columns

```{r complete-bio-ids, echo=TRUE, message=FALSE, warning=FALSE, results='hide'}
Stress_Masterfile_cleaned <- Stress_Masterfile %>%
  group_by(Bio_ID) %>%
  filter(across(all_of(analyte_columns), ~!is.na(.))) %>%
  ungroup()
```

### 2.4.2 Psychometric Ratings (GHQ and PSS)

```{r filter-t0-t1, echo=TRUE, message=FALSE, warning=FALSE, results='markup'}
Stress_Masterfile_filtered <- Stress_Masterfile_cleaned %>%
  group_by(Bio_ID) %>%
  filter(
    any(session == "T0" & !is.na(GHQ_total) & !is.na(PSS_Global_total)) &
    any(session == "T1" & !is.na(GHQ_total) & !is.na(PSS_Global_total))
  ) %>%
  ungroup()

# Print number of participants with both T0 and T1 GHQ + PSS
num_participants <- Stress_Masterfile_filtered %>% 
  distinct(Bio_ID) %>% 
  nrow()

cat("Number of participants with both T0 and T1 GHQ + PSS:", num_participants, "\n")
```

### 2.5 Recalculate BMI

```{r recalculate-bmi, include=FALSE}
Stress_Masterfile_filtered <- Stress_Masterfile_filtered %>%
  mutate(BMI_T0_new = as.numeric(Weight_T0) / (as.numeric(Height_T0) / 100)^2)
```

### 2.6 Calculate Delta Scores

```{r calculate-deltas, echo=TRUE, results='hide', message=FALSE, warning=FALSE}
delta_data <- Stress_Masterfile_filtered %>%
  group_by(Bio_ID) %>%
  summarize(across(c(GHQ_total, PSS_Global_total, STADI_Trait_total, WHO_total, 
                     GHQ_Somatic_symptoms, GHQ_Anxiety_insomnia, GHQ_Severe_depression, 
                     PSS_Helplessness, PSS_Self_efficacy, STADI_Trait_agitation, 
                     STADI_Trait_euthymia, STADI_Trait_dysthymia, LEC_Exposure, MIMIS_Exposure, 
                     all_of(analyte_columns)), 
                   ~ .[session == "T1"] - .[session == "T0"], 
                   .names = "Delta_{col}")) %>%
  ungroup()
```

## 3. Outlier Removal

TOBIAS: Here we def. wanna show the results in the output to make sure things are correct and include them in the supplement if needed

```{r outlier-removal, echo=TRUE, message=FALSE, warning=FALSE, results='markup'}
# Start with original data
delta_data <- delta_data

# Replace values > 2 SD from mean with NA
remove_outliers_sd <- function(x) {
  if (all(is.na(x))) return(x)
  mu <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  outliers <- abs(x - mu) > 2 * sd_x
  x[outliers] <- NA
  return(x)
}

# Identify delta variables
delta_vars <- delta_data %>%
  select(starts_with("Delta_")) %>%
  colnames()

# Keep a copy before cleaning
delta_data_before <- delta_data

# Remove outliers → store result in delta_data_final
delta_data_final <- delta_data %>%
  mutate(across(all_of(delta_vars), remove_outliers_sd))

# Count how many participants remain (should be the same unless entire rows became NA)
cat("Number of participants after outlier handling:\n")
print(nrow(delta_data_final))

# Summarize how many outliers were removed per variable
outlier_counts <- map_int(delta_vars, function(var) {
  sum(is.na(delta_data_final[[var]]) & !is.na(delta_data_before[[var]]))
})

outlier_summary <- data.frame(
  Variable = delta_vars,
  Outliers_Replaced = outlier_counts
)

# Print summary
print(outlier_summary)
```

## 4. Extract and Prepare Covariates

### 4.1 Covariate Extraction

TOBIAS: These variables should be ordinal - numeric handeling only for regression?

```{r prepare-covariates, echo=TRUE, message=FALSE, warning=FALSE}
# Create covariate dataset with sex and smoking as factors
covariates <- Stress_Masterfile_filtered %>%
  filter(session == "T0") %>%
  dplyr::select(Bio_ID, age_T0, sex_T0, BMI_T0_new, smoke) %>%
  mutate(
    age_T0 = as.numeric(age_T0),
    sex_factor = factor(sex_T0, levels = c("male", "female")),
    smoke_factor = factor(smoke, levels = c("no", "yes"))
  ) %>%
  dplyr::select(Bio_ID, age_T0, BMI_T0_new, sex_factor, smoke_factor)

# Remove old versions (if present) before merging to avoid .x/.y suffixes
delta_data_final <- delta_data_final %>%
  select(-any_of(c("age_T0", "BMI_T0_new", "sex_factor", "smoke_factor"))) %>%
  left_join(covariates, by = "Bio_ID")
```

### 4.2 Export

```{r export-csvs, echo=FALSE, message=FALSE, warning=FALSE}
#write.csv(Stress_Masterfile_cleaned, "Stress_Masterfile_cleaned.csv", row.names = FALSE)
#write.csv(delta_data, "delta_data.csv", row.names = FALSE)
#write.csv(delta_data_filtered, "delta_data_filtered.csv", row.names = FALSE)
#write.csv(delta_data_final, "delta_data_final.csv", row.names = FALSE)
```

### 5. Abbreviate Variable Names

```{r abbreviations, echo=TRUE, results='hide', message=FALSE, warning=FALSE}
# Create abbreviation map for analytes
abbrev_map <- c(
  "Delta_Neurofilament_light_polypeptide" = "NfL",
  "Delta_Interleukin-12" = "IL-12",
  "Delta_Interleukin-4" = "IL-4",
  "Delta_Neutrophil_elastase" = "ELA2",
  "Delta_Interleukin-10" = "IL-10",
  "Delta_Protein_S100-A11" = "S100A11",
  "Delta_Oncostatin-M" = "OSM",
  "Delta_Transforming_growth_factor_beta-1" = "TGF-\u03b21",
  "Delta_Thrombospondin-2" = "THBS2",
  "Delta_Interferon_beta" = "IFN-\u03b2",
  "Delta_Interferon_gamma" = "IFN-\u03b3",
  "Delta_Serum_amyloid_A-1_protein" = "SAA1",
  "Delta_Interferon_alpha-1/13" = "IFN-\u03b1",
  "Delta_Serum_amyloid_A-2_protein" = "SAA2",
  "Delta_Vascular_cell_adhesion_protein_1" = "VCAM-1",
  "Delta_Glial_fibrillary_acidic_protein" = "GFAP",
  "Delta_Neuropeptide_Y" = "NPY",
  "Delta_FK506-binding_protein_5" = "FKBP5",
  "Delta_Lipoprotein_lipase" = "LPL",
  "Delta_NOS" = "NOS",
  "Delta_Metalloproteinase_inhibitor_1" = "TIMP-1",
  "Delta_Metalloproteinase_inhibitor_2" = "TIMP-2",
  "Delta_C-C_motif_chemokine_7" = "CCL7",
  "Delta_Gro-gamma" = "CXCL3",
  "Delta_Brain-derived_neurotrophic_factor" = "BDNF",
  "Delta_Metalloproteinase_inhibitor_3" = "TIMP-3",
  "Delta_High_mobility_group_protein_B1" = "HMGB1",
  "Delta_Interleukin-6" = "IL-6",
  "Delta_Leptin" = "LEP",
  "Delta_C-C_motif_chemokine_2" = "CCL2",
  "Delta_Matrix_metalloproteinase-9" = "MMP-9",
  "Delta_Myeloperoxidase" = "MPO",
  "Delta_Vascular_endothelial_growth_factor_A" = "VEGF-A",
  "Delta_Platelet_endothelial_cell_adhesion_molecule" = "PECAM-1",
  "Delta_Platelet_factor_4" = "CXCL4",
  "Delta_Complement_C3" = "C3",
  "Delta_Stromelysin-1" = "MMP-3",
  "Delta_Matrilysin" = "MMP-7",
  "Delta_Neutrophil-activating_peptide_2" = "CXCL7",
  "Delta_A_disintegrin_and_metalloproteinase_with_thrombospondin_motifs_4" = "ADAMTS4",
  "Delta_C5a_anaphylatoxin" = "C5a",
  "Delta_C-C_motif_chemokine_4" = "CCL4",
  "Delta_Fibroblast_growth_factor_21" = "FGF21",
  "Delta_Neutrophil_collagenase" = "MMP-8",
  "Delta_C-X-C_motif_chemokine_5" = "CXCL5",
  "Delta_Growth-regulated_alpha_protein" = "CXCL1",
  "Delta_Fibroblast_growth_factor_2" = "FGF2",
  "Delta_Interleukin-1_beta" = "IL-1\u03b2",
  "Delta_C-C_motif_chemokine_3" = "CCL3",
  "Delta_von_Willebrand_factor" = "vWF",
  "Delta_Lipopolysaccharide-binding_protein" = "LBP",
  "Delta_iNOS" = "iNOS",
  "Delta_Gro-beta" = "CXCL2",
  "Delta_Thyroid_Stimulating_Hormone" = "TSH",
  "Delta_A_disintegrin_and_metalloproteinase_with_thrombospondin_motifs_5" = "ADAMTS5",
  "Delta_A_disintegrin_and_metalloproteinase_with_thrombospondin_motifs_1" = "ADAMTS1",
  "Delta_A_disintegrin_and_metalloproteinase_with_thrombospondin_motifs_13" = "ADAMTS13",
  "Delta_COX-1" = "COX-1",
  "Delta_NACHT,_LRR_and_PYD_domains-containing_protein_3" = "NLRP3",
  "Delta_A_disintegrin_and_metalloproteinase_with_thrombospondin_motifs_9" = "ADAMTS9",
  "Delta_Interleukin-8" = "IL-8",
  "Delta_E-selectin" = "E-Selectin",
  "Delta_C-X-C_motif_chemokine_6" = "CXCL6",
  "Delta_Interleukin-17A" = "IL-17A",
  "Delta_Transforming_growth_factor_beta-3" = "TGF-\u03b23",
  "Delta_CD40_ligand" = "CD40L",
  "Delta_A_disintegrin_and_metalloproteinase_with_thrombospondin_motifs_12" = "ADAMTS12",
  "Delta_Apolipoprotein(a)" = "Apo(a)",
  "Delta_Prostaglandin_G/H_synthase_2" = "PTGS2",
  "Delta_Fibronectin" = "FN1",
  "Delta_C-X-C_motif_chemokine_10" = "CXCL10",
  "Delta_P-selectin" = "P-Selectin",
  "Delta_Transforming_growth_factor_beta-2" = "TGF-\u03b22",
  "Delta_72_kDa_type_IV_collagenase" = "MMP-2",
  "Delta_C-reactive_protein" = "CRP",
  "Delta_Complement_C4" = "C4",
  "Delta_Macrophage_metalloelastase" = "MMP-12",
  "Delta_Interleukin-1_alpha" = "IL-1\u03b1",
  "Delta_Insulin" = "INS",
  "Delta_Corticotropin" = "ACTH",
  "Delta_Collagenase_3" = "MMP-13",
  "Delta_Annexin_A1" = "ANXA1",
  "Delta_Annexin_A2" = "ANXA2",
  "Delta_Eotaxin" = "CCL11",
  "Delta_C-C_motif_chemokine_5" = "CCL5",
  "Delta_Interleukin-18" = "IL-18",
  "Delta_Tumor_necrosis_factor" = "TNF",
  "Delta_beta-nerve_growth_factor" = "\u03b2-NGF",
  "Delta_Glial_cell_line-derived_neurotrophic_factor" = "GDNF",
  "Delta_Ferritin" = "FTH1",
  "Delta_Thrombopoietin" = "THPO",
  "Delta_Metalloproteinase_inhibitor_4" = "TIMP-4",
  "Delta_Appetite-regulating_hormone" = "Ghrelin",
  "Delta_High_mobility_group_protein_B2" = "HMGB2",
  "Delta_Insulin-like_growth_factor_I" = "IGF-I",
  "Delta_Granulocyte_colony-stimulating_factor" = "G-CSF",
  "Delta_Fibroblast_growth_factor_22" = "FGF22",
  "Delta_Protein_S100-A4" = "S100A4",
  "Delta_Indoleamine_2,3-dioxygenase_1" = "IDO-1",
  "Delta_Tryptophan_2,3-dioxygenase" = "TDO"
)
```

## 5. GHQ & PSS Violin Plots 

```{r ghq-pss-violin-plots, message=FALSE, warning=FALSE, results='markup'}
# Get Bio_IDs with valid GHQ and PSS delta 
valid_ghq_ids <- delta_data_final %>%
  filter(!is.na(Delta_GHQ_total)) %>%
  pull(Bio_ID)

valid_pss_ids <- delta_data_final %>%
  filter(!is.na(Delta_PSS_Global_total)) %>%
  pull(Bio_ID)

# Filter original masterfile for those participants
ghq_data <- Stress_Masterfile_cleaned %>%
  filter(Bio_ID %in% valid_ghq_ids & session %in% c("T0", "T1")) %>%
  dplyr::select(Bio_ID, session, GHQ_total) %>%
  mutate(Measure = "GHQ Score Total", Value = GHQ_total) %>%
  dplyr::select(-GHQ_total)

pss_data <- Stress_Masterfile_cleaned %>%
  filter(Bio_ID %in% valid_pss_ids & session %in% c("T0", "T1")) %>%
  dplyr::select(Bio_ID, session, PSS_Global_total) %>%
  mutate(Measure = "PSS Score Total", Value = PSS_Global_total) %>%
  dplyr::select(-PSS_Global_total)

# Combine into one dataframe
violin_data <- bind_rows(ghq_data, pss_data) %>%
  rename(`Time Point` = session) %>%
  mutate(`Time Point` = factor(`Time Point`, levels = c("T0", "T1")))

# Violin plot function (with paired t-test)
plot_violin <- function(df, measure) {
  df_filtered <- df %>% filter(Measure == measure)
  y_max <- max(df_filtered$Value, na.rm = TRUE)
  y_lim_upper <- y_max + 10

  ggplot(df_filtered, 
         aes(x = `Time Point`, y = Value, fill = `Time Point`)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "black") +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.2, color = "black") +
    stat_compare_means(method = "t.test", 
                       comparisons = list(c("T0", "T1")), 
                       label = "p.signif", 
                       paired = TRUE,         # <--- Paired t-test added here
                       tip.length = 0.01, 
                       size = 7,
                       label.y = y_max + 6,
                       bracket.size = 1) +
    labs(x = "Time Point", y = measure) +
    ylim(NA, y_lim_upper) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 14),
      legend.position = "none",
      plot.title = element_blank()
    ) +
    scale_fill_manual(values = c("T0" = "#afea7b", "T1" = "#eab67b"), drop = FALSE)
}

# Plot
p1 <- plot_violin(violin_data, "GHQ Score Total")
p2 <- plot_violin(violin_data, "PSS Score Total")

print(p1)
print(p2)
```

## 6. Linear Regression Analysis: GHQ and PSS vs. Analytes

### 6.1 Run Regression

```{r linear-regression, echo=TRUE, message=FALSE, warning=FALSE, results='markup'}
# Identify analyte columns (starting from Neurofilament)
start_col <- which(colnames(delta_data_final) == "Delta_Neurofilament_light_polypeptide")
analyte_columns <- colnames(delta_data_final)[start_col:(start_col + 99)]

# Function to run regression and return results
run_regression <- function(dep_var, analyte_name) {
  regression_data <- delta_data_final %>%
    dplyr::select(Bio_ID, all_of(dep_var), all_of(analyte_name), age_T0, BMI_T0_new, smoke_factor, sex_factor) %>%
    rename(Dependent_Var = all_of(dep_var), 
           Analyte = all_of(analyte_name), 
           Age = age_T0, 
           BMI = BMI_T0_new, 
           Smoke = smoke_factor, 
           Sex = sex_factor) %>%
    drop_na()

  # Fit linear model with categorical covariates
  model <- lm(Dependent_Var ~ Analyte + Age + BMI + Smoke + Sex, data = regression_data)
  
  coef_value <- summary(model)$coefficients["Analyte", "Estimate"]
  p_value <- summary(model)$coefficients["Analyte", "Pr(>|t|)"]
  n_subjects <- nrow(regression_data)
  
  data.frame(
    Analyte = analyte_name,
    Dependent_Variable = dep_var,
    Coefficient = coef_value,
    P_Value = p_value,
    N_Subjects = n_subjects
  )
}

# Run regressions separately for GHQ and PSS
regression_summary_GHQ <- data.frame()
regression_summary_PSS <- data.frame()

for (analyte in analyte_columns) {
  regression_summary_GHQ <- rbind(regression_summary_GHQ, run_regression("Delta_GHQ_total", analyte))
  regression_summary_PSS <- rbind(regression_summary_PSS, run_regression("Delta_PSS_Global_total", analyte))
}

# Multiple correction
regression_summary_GHQ$FDR_BH <- p.adjust(regression_summary_GHQ$P_Value, method = "BH")
regression_summary_GHQ$FDR_Bonferroni <- p.adjust(regression_summary_GHQ$P_Value, method = "bonferroni")

regression_summary_PSS$FDR_BH <- p.adjust(regression_summary_PSS$P_Value, method = "BH")
regression_summary_PSS$FDR_Bonferroni <- p.adjust(regression_summary_PSS$P_Value, method = "bonferroni")
```

### 6.2 Combine GHQ and PSS Regression Results

```{r combine-regression-tables, echo=FALSE, message=FALSE, warning=FALSE}
regression_summary_table <- bind_rows(
  regression_summary_GHQ %>% 
    dplyr::select(Analyte, Dependent_Variable, Coefficient, P_Value, N_Subjects, FDR_BH),
  regression_summary_PSS %>% 
    dplyr::select(Analyte, Dependent_Variable, Coefficient, P_Value, N_Subjects, FDR_BH)
)

# Rename FDR column for consistency
colnames(regression_summary_table)[colnames(regression_summary_table) == "FDR_BH"] <- "FDR_P_Value"
```

## 7. Heatmap

### 7.1 Data Preparation Heat Map

```{r regression-heatmap, echo=TRUE, message=FALSE, warning=FALSE, results='markup'}
# Apply abbreviations and prepare heatmap data
heatmap_data <- regression_summary_table %>%
  dplyr::select(Analyte, Dependent_Variable, Coefficient, P_Value) %>%
  mutate(
    # Apply significance labels
    Significance = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # Use abbreviations for Analytes
    Analyte = ifelse(Analyte %in% names(abbrev_map),
                     abbrev_map[Analyte],
                     Analyte),
    # Factor Dependent Variable
    Dependent_Variable = factor(Dependent_Variable,
                                levels = c("Delta_GHQ_total", "Delta_PSS_Global_total"))
  )

# Sort analytes based on GHQ effect size (after abbreviation)
ghq_sorted <- heatmap_data %>%
  filter(Dependent_Variable == "Delta_GHQ_total") %>%
  arrange(desc(Coefficient)) %>%
  pull(Analyte) %>%
  unique()

# Final formatting for plot
heatmap_data <- heatmap_data %>%
  mutate(
    Dependent_Variable = recode(Dependent_Variable,
                                "Delta_GHQ_total" = "GHQ",
                                "Delta_PSS_Global_total" = "PSS"),
    Analyte = factor(Analyte, levels = ghq_sorted)
  )

```

### 7.2 Plotting Heatmap

```{r regression-heatmap-plot, echo=TRUE, message=FALSE, warning=FALSE, results='markup', fig.height=14, fig.width=9}
heatmap_plot <- ggplot(heatmap_data, aes(x = Dependent_Variable, y = Analyte, fill = Coefficient)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#0014a8", mid = "#ebf5ff", high = "#e32636",
    midpoint = 0, name = "Regression Coefficient"
  ) +  
  geom_text(aes(label = Significance), color = "black", size = 4.5, vjust = 0.65, fontface = "bold", family = "sans") +
  labs(title = "", x = "", y = "Analytes") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = "sans", face = "bold", color = "black"),
    axis.text.y = element_text(size = 7.5, vjust = 0.5, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 18, color = "black"),
    legend.title = element_text(size = 8, face = "bold", color = "black"),
    legend.text = element_text(size = 9, color = "black")
  )

# Show plot
print(heatmap_plot)
```

## 8. Significant Results and Regression Plots

```{r significant-results-inline, echo=TRUE, message=FALSE, warning=FALSE, results='markup', fig.width=6, fig.height=4}
# Extract significant correlations 
significant_results <- regression_summary_table %>%
  filter(P_Value < 0.05) %>%
  mutate(
    Dependent_Variable = recode(
      Dependent_Variable,
      "Delta_GHQ_total" = "GHQ",
      "Delta_PSS_Global_total" = "PSS"
    ),
    Significance = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Analyte_Clean = Analyte,  # ✅ This line is now correct
    Analyte_Label = ifelse(Analyte_Clean %in% names(abbrev_map),
                           abbrev_map[Analyte_Clean],
                           Analyte)
  )


# Plotting
for (i in 1:nrow(significant_results)) {

  analyte_name <- significant_results$Analyte_Clean[i]
  dep_var <- significant_results$Dependent_Variable[i]
  abbrev_label <- significant_results$Analyte_Label[i]

  # Check if the analyte exists in the data
  if (!(analyte_name %in% names(delta_data_final))) {
    cat("\n Skipping:", analyte_name, "→ Not found in `delta_data_final`.\n")
    next
  }

  dep_var_col <- ifelse(dep_var == "GHQ", "Delta_GHQ_total", "Delta_PSS_Global_total")

  plot_data <- delta_data_final %>%
    dplyr::select(Bio_ID, all_of(dep_var_col), all_of(analyte_name)) %>%
    rename(Dependent_Var = all_of(dep_var_col), Analyte = all_of(analyte_name))

  if (nrow(plot_data) == 0) {
    cat("\n Skipping:", analyte_name, "→ No valid data points!\n")
    next
  }

  p_value <- significant_results$P_Value[i]
  p_value_text <- ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", round(p_value, 3)))

  annotation_x <- max(plot_data$Analyte, na.rm = TRUE) * 0.9
  annotation_y <- max(plot_data$Dependent_Var, na.rm = TRUE) * 0.9

  p <- ggplot(plot_data, aes(x = Analyte, y = Dependent_Var)) +
    geom_point(alpha = 0.7, size = 1.5, color = "black") +
    geom_smooth(method = "lm", color = "red", linetype = "solid") +
    geom_segment(aes(x = min(Analyte, na.rm = TRUE), xend = max(Analyte, na.rm = TRUE),
                     y = min(Dependent_Var, na.rm = TRUE), yend = min(Dependent_Var, na.rm = TRUE)), color = "black") +
    geom_segment(aes(x = min(Analyte, na.rm = TRUE), xend = min(Analyte, na.rm = TRUE),
                     y = min(Dependent_Var, na.rm = TRUE), yend = max(Dependent_Var, na.rm = TRUE)), color = "black") +
    annotate("text", x = annotation_x, y = annotation_y, label = p_value_text,
             hjust = 1, vjust = 1, size = 6, fontface = "bold") +
    labs(
      x = bquote(.(abbrev_label) ~ "[" ~ Delta ~ T[1]-T[0] ~ "]"),
      y = bquote(.(dep_var) ~ "Total Score [" ~ Delta ~ T[1]-T[0] ~ "]"),
      title = NULL
    ) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      panel.grid = element_blank()
    )

  print(p)
}
```

## 9. Network Analysis

### 9.1 Network Analysis

###### Oridnal Data -> Change sex & Smoking to categorical and then use mgm instead of ggModSelect for nw estimation

```{r network-analysis-estimate, echo=TRUE, message=FALSE, warning=FALSE, results='markup'}
# Select analytes that showed significant correlation
selected_analytes <- significant_results$Analyte_Clean
network_vars <- c(selected_analytes, "Delta_GHQ_total", "Delta_PSS_Global_total")
covariate_vars <- c("age_T0", "BMI_T0_new", "smoke_factor", "sex_factor")
all_vars <- unique(c(network_vars, covariate_vars))

# Convert factors to integers before estimation (required for mgm)
network_input_full <- delta_data_final[, names(delta_data_final) %in% all_vars] %>%
  mutate(
    sex_factor = as.integer(sex_factor),
    smoke_factor = as.integer(smoke_factor)
  )

# Estimate full network
nw_full <- estimateNetwork(
  network_input_full,
  default = "mgm",
  lambda = 0.1,
  threshold = "none"
)

# Fix: Add dimnames so we can subset by name
dimnames(nw_full$graph) <- list(colnames(network_input_full), colnames(network_input_full))

# Now subset the network to analytes + outcomes only
nw <- nw_full
nw$graph <- nw$graph[network_vars, network_vars]
nw$labels <- nw$labels[network_vars]
```

### 9.2 Network Plotting

```{r network-plotting, echo=TRUE, message=FALSE, warning=FALSE, results='markup', fig.height=4, fig.width=6}
# Raw labels from graph (should match column names like "Delta_IL-6", "Delta_GHQ_total", etc.)
raw_labels <- colnames(nw$graph)

# Manually handle GHQ and PSS renaming first
raw_labels <- recode(raw_labels,
                     "Delta_GHQ_total" = "GHQ",
                     "Delta_PSS_Global_total" = "PSS")

# Then apply abbreviation map
formatted_labels <- ifelse(raw_labels %in% names(abbrev_map),
                           abbrev_map[raw_labels],
                           raw_labels)
# Aesthetics
node_colors <- rep("lightblue", length(formatted_labels))
node_size <- 11
node_border_width <- 2.5
node_border_color <- "black"
edge_colors <- ifelse(nw$graph > 0, "red", ifelse(nw$graph < 0, "blue", "gray"))
label_font <- ifelse(formatted_labels %in% c("GHQ", "PSS"), 1, 2)

# Plot network
plot(nw,
     layout = "spring",
     title = "",
     labels = formatted_labels,
     label.cex = 1.2,
     label.color = "black",
     label.font = label_font,
     edge.labels = FALSE,
     edge.color = edge_colors,
     node.color = node_colors,
     vsize = node_size,
     border.width = node_border_width,
     border.color = node_border_color
)
```

### 9.2 Stability Analysis

TOBIAS: Title and documentation - show only output

```{r sensitivity-analysis, echo=TRUE, message=FALSE, warning=FALSE, results='hold'}
boot1 <- bootnet(nw, default = "ggmModSelect", nBoots = 1000, type = "nonparametric")

plot(boot1, order = "sample", labels = FALSE)
plot(boot1, "edge", plot = "difference")

boot2 <- bootnet(nw, default = "ggmModSelect", nBoots = 1000, type = "case")
corStability(boot2)
```
