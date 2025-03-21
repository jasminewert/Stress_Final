---
title: "HMZ_final"
author: "Jasmin Ewert"
date: "2025-03-21"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
## Load Required Libraries

```{r load-libraries}
library(tidyverse)
library(ggplot2)
library(readr)
```

## Import

```{r load-data}
Stress_Masterfile <- read_delim("Stress_Masterfile.csv", 
                                delim = ";", escape_double = FALSE, trim_ws = TRUE)

# View(Stress_Masterfile) 
```

## Remove Rows with All NAs

```{r remove-na-rows}
Stress_Masterfile <- Stress_Masterfile %>%
  filter(rowSums(is.na(.)) != ncol(.))
```

## Identify and Clean Analyte Columns

```{r clean-analyte-columns}
# Identify analyte columns (adjust indices if needed)
analyte_columns <- colnames(Stress_Masterfile)[62:161]

# Replace spaces in column names with underscores
colnames(Stress_Masterfile) <- gsub(" ", "_", colnames(Stress_Masterfile))
analyte_columns <- gsub(" ", "_", analyte_columns)
```

## Keep Only Complete Bio_IDs

```{r complete-bio-ids}
bio_ids_to_keep <- Stress_Masterfile %>%
  group_by(Bio_ID) %>%
  filter(if_all(all_of(analyte_columns), ~ !is.na(.))) %>%
  pull(Bio_ID) %>%
  unique()

Stress_Masterfile_cleaned <- Stress_Masterfile %>%
  filter(Bio_ID %in% bio_ids_to_keep)

removed_bio_ids <- setdiff(unique(Stress_Masterfile$Bio_ID), bio_ids_to_keep)
cat("Removed Bio_IDs due to missing analyte values:\n")
print(removed_bio_ids)
```

## Filter for Participants with Both T0 and T1

```{r filter-t0-t1}
Stress_Masterfile_filtered <- Stress_Masterfile_cleaned %>%
  group_by(Bio_ID) %>%
  filter(any(session == "T0" & !is.na(GHQ_total)) &
         any(session == "T1" & !is.na(GHQ_total))) %>%
  ungroup()

# Print number of participants with both T0 and T1
num_participants <- Stress_Masterfile_filtered %>% 
  distinct(Bio_ID) %>% 
  nrow()
cat("Number of participants with both T0 and T1:", num_participants, "\n")
```

## Recalculate BMI

```{r recalculate-bmi}
Stress_Masterfile_filtered <- Stress_Masterfile_filtered %>%
  mutate(BMI_T0_new = as.numeric(Weight_T0) / (as.numeric(Height_T0) / 100)^2)
```

## Calculate Delta Scores

```{r calculate-deltas}
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

## Filter Participants Missing Analytes (Neurofilament Light Polypeptide)

```{r filter-nfl}
delta_data_filtered <- delta_data %>%
  filter(!is.na(Delta_Neurofilament_light_polypeptide))

removed_bio_ids <- setdiff(delta_data$Bio_ID, delta_data_filtered$Bio_ID)
cat("Removed Bio_IDs due to missing Delta_Neurofilament_light_polypeptide values:\n")
print(removed_bio_ids)

delta_data_final <- delta_data_filtered
```

## Extract and Prepare Covariates

```{r prepare-covariates}
covariates <- Stress_Masterfile_filtered %>%
  filter(session == "T0") %>%
  select(Bio_ID, age_T0, sex_T0, BMI_T0_new, smoke, meds) %>%
  mutate(
    age_T0 = as.numeric(age_T0),
    sex_numeric = ifelse(sex_T0 == "male", 0, ifelse(sex_T0 == "female", 1, NA)),
    smoke_numeric = ifelse(smoke == "yes", 1, ifelse(smoke == "no", 0, NA)),
    meds_numeric = ifelse(meds == "yes", 1, ifelse(meds == "no", 0, NA))
  ) %>%
  select(Bio_ID, age_T0, sex_numeric, BMI_T0_new, smoke_numeric, meds_numeric)

delta_data_final <- delta_data_final %>%
  left_join(covariates, by = "Bio_ID")
#View(delta_data_final)
```

## Export Cleaned Datasets

```{r export-csvs}
#write.csv(Stress_Masterfile_cleaned, "Stress_Masterfile_cleaned.csv", row.names = FALSE)
#write.csv(delta_data, "delta_data.csv", row.names = FALSE)
#write.csv(delta_data_filtered, "delta_data_filtered.csv", row.names = FALSE)
#write.csv(delta_data_final, "delta_data_final.csv", row.names = FALSE)
```

## GHQ & PSS Violin Plots

```{r ghq-pss-violin-plots, message=FALSE, warning=FALSE}
library(cowplot)

# Extract relevant data
violin_data <- Stress_Masterfile_cleaned %>%
  filter(session %in% c("T0", "T1")) %>%
  select(Bio_ID, session, GHQ_total, PSS_Global_total) %>%
  pivot_longer(cols = c(GHQ_total, PSS_Global_total), 
               names_to = "Measure", values_to = "Value") %>%
  rename(`Time Point` = session)

# Function to remove outliers 
remove_outliers <- function(df) {
  df %>%
    group_by(Measure, `Time Point`) %>%
    mutate(Mean = mean(Value, na.rm = TRUE),
           SD = sd(Value, na.rm = TRUE)) %>%
    filter(Value >= (Mean - 2 * SD) & Value <= (Mean + 2 * SD)) %>%
    ungroup()
}

# Apply outlier removal
violin_data <- remove_outliers(violin_data)

# Clean factor names
violin_data$Measure <- recode(violin_data$Measure, 
                              "GHQ_total" = "GHQ Score Total",
                              "PSS_Global_total" = "PSS Score Total")

violin_data <- violin_data %>%
  mutate(`Time Point` = factor(`Time Point`, levels = c("T0", "T1")))

# Custom violin plot function
plot_violin <- function(df, measure) {
  ggplot(df %>% filter(Measure == measure), 
         aes(x = `Time Point`, y = Value, fill = `Time Point`)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "black") +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.2, color = "black") +
    labs(x = "Time Point", y = measure) +
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
    scale_fill_manual(values = c("T0" = "#FF7070", "T1" = "#709BFF"), drop = FALSE)
}

# Generate and print both plots separately
p1 <- plot_violin(violin_data, "GHQ Score Total")
p2 <- plot_violin(violin_data, "PSS Score Total")

print(p1)
print(p2)
```

## Linear Regression Analysis: GHQ and PSS vs. Analytes

```{r linear-regression-ghq-pss, message=FALSE, warning=FALSE}
library(dplyr)

# Identify analyte columns (starting from Neurofilament)
start_col <- which(colnames(delta_data_final) == "Delta_Neurofilament_light_polypeptide")
analyte_columns <- colnames(delta_data_final)[start_col:(start_col + 99)]

# Function to run regression and return results
run_regression <- function(dep_var, analyte_name) {
  regression_data <- delta_data_final %>%
    select(Bio_ID, all_of(dep_var), all_of(analyte_name), age_T0, BMI_T0_new, smoke_numeric, sex_numeric) %>%
    rename(Dependent_Var = all_of(dep_var), Analyte = all_of(analyte_name), 
           Age = age_T0, BMI = BMI_T0_new, Smoke = smoke_numeric, Sex = sex_numeric)
  
  # Outlier removal: 2 SD
  regression_data_clean <- regression_data %>%
    filter(
      Dependent_Var >= (mean(Dependent_Var, na.rm = TRUE) - 2 * sd(Dependent_Var, na.rm = TRUE)) &
      Dependent_Var <= (mean(Dependent_Var, na.rm = TRUE) + 2 * sd(Dependent_Var, na.rm = TRUE)) &
      Analyte >= (mean(Analyte, na.rm = TRUE) - 2 * sd(Analyte, na.rm = TRUE)) &
      Analyte <= (mean(Analyte, na.rm = TRUE) + 2 * sd(Analyte, na.rm = TRUE))
    )
  
  model <- lm(Dependent_Var ~ Analyte + Age + BMI + Smoke + Sex, data = regression_data_clean)
  
  coef_value <- summary(model)$coefficients["Analyte", "Estimate"]
  p_value <- summary(model)$coefficients["Analyte", "Pr(>|t|)"]
  n_subjects <- nrow(regression_data_clean)
  
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

# view(regression_summary_GHQ)
# view(regression_summary_PSS)

# Export
# write.csv(regression_summary_GHQ, "HMZ_Publication/Regression_Summary_Table_GHQ_FDR.csv", row.names = FALSE)
# write.csv(regression_summary_PSS, "HMZ_Publication/Regression_Summary_Table_PSS_FDR.csv", row.names = FALSE)
```

## Combine GHQ and PSS Regression Results

```{r combine-regression-tables}
regression_summary_table <- bind_rows(
  regression_summary_GHQ %>% 
    select(Analyte, Dependent_Variable, Coefficient, P_Value, N_Subjects, FDR_BH),
  regression_summary_PSS %>% 
    select(Analyte, Dependent_Variable, Coefficient, P_Value, N_Subjects, FDR_BH)
)

# Rename FDR column for consistency
colnames(regression_summary_table)[colnames(regression_summary_table) == "FDR_BH"] <- "FDR_P_Value"
```

## Heatmap of Regression Coefficients

```{r regression-heatmap-cleaned, message=FALSE, warning=FALSE}
library(ggplot2)
library(dplyr)

# Prepare and clean data for heatmap
heatmap_data <- regression_summary_table %>%
  select(Analyte, Dependent_Variable, Coefficient, P_Value) %>%
  mutate(
    Significance = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Dependent_Variable = factor(Dependent_Variable, 
                                levels = c("Delta_GHQ_total", "Delta_PSS_Global_total")),
    # Clean analyte names
    Analyte = gsub("Delta_", "", Analyte),
    Analyte = gsub("_", " ", Analyte),
    Analyte = gsub("\\.", " ", Analyte)
  )

# Sort analytes based on GHQ coefficient (after cleaning names)
ghq_sorted <- heatmap_data %>%
  filter(Dependent_Variable == "Delta_GHQ_total") %>%
  arrange(desc(Coefficient)) %>%
  pull(Analyte) %>%
  unique()

# Apply final formatting
heatmap_data <- heatmap_data %>%
  mutate(
    Dependent_Variable = recode(Dependent_Variable, 
                                "Delta_GHQ_total" = "GHQ", 
                                "Delta_PSS_Global_total" = "PSS"),
    Analyte = factor(Analyte, levels = ghq_sorted)
  )

# Plot the heatmap
heatmap_plot <- ggplot(heatmap_data, aes(x = Dependent_Variable, y = Analyte, fill = Coefficient)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#0014a8", mid = "#ebf5ff", high = "#e32636", midpoint = 0, name = "Regression Coefficient") +  
  geom_text(aes(label = Significance), color = "black", size = 5, vjust = 0.75) +
  labs(title = "Regression Coefficients Heatmap", x = "", y = "Analytes") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 7, face = "bold", vjust = 0.5),
        axis.text.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"))

# Save and display
ggsave("HMZ_Publication/Regression_Heatmap.png", plot = heatmap_plot, width = 8, height = 12, dpi = 300)
print(heatmap_plot)
```
## Significant Results and Inline Regression Plots

```{r significant-results-inline, echo=FALSE, message=FALSE, warning=FALSE, fig.width=6, fig.height=4}
library(ggplot2)
library(dplyr)

#### Significant Results ####

# Extract significant correlations (p < 0.05)
significant_results <- regression_summary_table %>%
  filter(P_Value < 0.05) %>%
  mutate(Analyte = gsub("Delta_", "", Analyte),
         Analyte = gsub("_", " ", Analyte),
         Analyte = gsub("\\.", " ", Analyte),
         Dependent_Variable = recode(Dependent_Variable,
                                     "Delta_GHQ_total" = "GHQ",
                                     "Delta_PSS_Global_total" = "PSS"),
         Significance = case_when(
           P_Value < 0.001 ~ "***",
           P_Value < 0.01 ~ "**",
           P_Value < 0.05 ~ "*",
           TRUE ~ ""
         )) %>%
  select(Analyte, Dependent_Variable, Coefficient, P_Value, Significance)

# Clean analyte names to match delta_data_final
significant_results <- significant_results %>%
  mutate(Analyte_Clean = paste0("Delta_", gsub(" ", "_", Analyte)))

# Identify available analytes
available_analytes <- intersect(significant_results$Analyte_Clean, colnames(delta_data_final))
missing_analytes <- setdiff(significant_results$Analyte_Clean, colnames(delta_data_final))

cat("\n✅ Available analytes for plotting:", length(available_analytes), "\n")
print(available_analytes)

cat("\n❌ Missing analytes (not found in delta_data_final):", length(missing_analytes), "\n")
print(missing_analytes)

#### Inline Plots ####

for (i in 1:nrow(significant_results)) {

  analyte_name <- significant_results$Analyte_Clean[i]
  dep_var <- significant_results$Dependent_Variable[i]

  if (!(analyte_name %in% available_analytes)) {
    cat("\n⚠ Skipping:", analyte_name, "→ Not found in `delta_data_final`.\n")
    next
  }

  dep_var_col <- ifelse(dep_var == "GHQ", "Delta_GHQ_total", "Delta_PSS_Global_total")

  plot_data <- delta_data_final %>%
    select(Bio_ID, all_of(dep_var_col), all_of(analyte_name)) %>%
    rename(Dependent_Var = all_of(dep_var_col), Analyte = all_of(analyte_name)) %>%
    filter(
      Dependent_Var >= (mean(Dependent_Var, na.rm = TRUE) - 2 * sd(Dependent_Var, na.rm = TRUE)) &
      Dependent_Var <= (mean(Dependent_Var, na.rm = TRUE) + 2 * sd(Dependent_Var, na.rm = TRUE)) &
      Analyte >= (mean(Analyte, na.rm = TRUE) - 2 * sd(Analyte, na.rm = TRUE)) &
      Analyte <= (mean(Analyte, na.rm = TRUE) + 2 * sd(Analyte, na.rm = TRUE))
    )

  if (nrow(plot_data) == 0) {
    cat("\n⚠ Skipping:", analyte_name, "→ No valid data points after filtering!\n")
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
      title = "",
      x = bquote(.(significant_results$Analyte[i]) ~ "[" ~ Delta ~ T[1]-T[0] ~ "]"),
      y = bquote(.(dep_var) ~ "Total Score [" ~ Delta ~ T[1]-T[0] ~ "]")
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


## Network Analysis

```{r network-analysis, echo=TRUE, message=FALSE, warning=FALSE}
# Load required libraries
library(qgraph)
library(bootnet)
library(ppcor)

# Define Variables
selected_analytes <- significant_results$Analyte_Clean  # List of 15 significant analytes
network_vars <- c(selected_analytes, "Delta_GHQ_total", "Delta_PSS_Global_total")  # Excluding covariates
covariate_vars <- c("age_T0", "BMI_T0_new", "smoke_numeric", "sex_numeric")  # Covariates

# Extract Data
network_data <- delta_data_final[, names(delta_data_final) %in% network_vars]
covariate_data <- delta_data_final[, names(delta_data_final) %in% covariate_vars]

# Residualize the data to control for covariates
residualized_data <- network_data  # Copy data for modification

for (var in network_vars) {
  if (var != "Bio_ID") {  # Avoid modifying the Bio_ID column
    formula <- as.formula(paste0("`", var, "` ~ ", paste(covariate_vars, collapse = " + ")))
    residualized_data[[var]] <- resid(lm(formula, data = delta_data_final))
  }
}

# Function to remove outliers (values > ±2 SD)
remove_outliers <- function(data) {
  for (col in names(data)) {
    if (is.numeric(data[[col]])) {  # Only apply to numeric columns
      mean_val <- mean(data[[col]], na.rm = TRUE)
      sd_val <- sd(data[[col]], na.rm = TRUE)
      # Set values outside mean ± 2*SD to NA
      data[[col]][data[[col]] < (mean_val - 2 * sd_val) | 
                  data[[col]] > (mean_val + 2 * sd_val)] <- NA
    }
  }
  return(data)
}

# Apply outlier removal
residualized_data <- remove_outliers(residualized_data)

# Remove rows where any column has NA
residualized_data <- na.omit(residualized_data)

# Compute Partial Pearson Correlation Matrix
cor_matrix <- pcor(residualized_data, method = "pearson")$estimate

# Format labels
formatted_labels <- gsub("Delta_", "", colnames(cor_matrix))
formatted_labels <- gsub("_", " ", formatted_labels)
formatted_labels[formatted_labels == "GHQ total"] <- "GHQ"
formatted_labels[formatted_labels == "PSS Global total"] <- "PSS"
formatted_labels[formatted_labels == "Interleukin-4"] <- "IL-4"
formatted_labels[formatted_labels == "Vascular cell adhesion protein 1"] <- "VCAM-1"
formatted_labels[formatted_labels == "Platelet factor 4"] <- "CXCL4"
formatted_labels[formatted_labels == "Neutrophil collagenase"] <- "MMP-8"
formatted_labels[formatted_labels == "Growth-regulated alpha protein"] <- "CXCL1"
formatted_labels[formatted_labels == "Lipopolysaccharide-binding protein"] <- "LBP"
formatted_labels[formatted_labels == "Gro-beta"] <- "CXCL2"
formatted_labels[formatted_labels == "A disintegrin and metalloproteinase with thrombospondin motifs 9"] <- "ADAMTS9"
formatted_labels[formatted_labels == "Interleukin-17A"] <- "IL-17A"
formatted_labels[formatted_labels == "CD40 ligand"] <- "CD40L"
formatted_labels[formatted_labels == "Prostaglandin G/H synthase 2"] <- "PTGS2"
formatted_labels[formatted_labels == "C-X-C motif chemokine 10"] <- "CXCL10"
formatted_labels[formatted_labels == "Annexin A2"] <- "ANXA2"
formatted_labels[formatted_labels == "Indoleamine 2,3-dioxygenase 1"] <- "IDO-1"

# Set custom edge colors: red for positive, blue for negative
edge_colors <- ifelse(cor_matrix > 0, "red", 
                      ifelse(cor_matrix < 0, "blue", "grey"))

# Plot the network
qgraph(cor_matrix,
       layout = "spring",
       theme = "colorblind",
       labels = formatted_labels,
       label.cex = 2,
       color = "skyblue",
       minimum = 0.1,
       cut = 0.3,
       vsize = 8,
       esize = 20,
       edge.labels = TRUE,
       edge.label.cex = 1.2,
       title = "Partial Pearson Network with Coefficients",
       legend.cex = 1.2,
       edge.color = edge_colors)  # Correct color mapping
```
### We removed outliers more than 2xstd from mean. We entirely removed the Bio_ID. This reduces the sample size from 82 to 48. Adapt the outlier removal by using different outlier removal strategies? 

### Outlier Removal only for the respective Analyte

```{r network-analysis-variablewise-outlier, echo=TRUE, message=FALSE, warning=FALSE}
# Load required libraries
library(qgraph)
library(bootnet)
library(ppcor)

# Define Variables
selected_analytes <- significant_results$Analyte_Clean  # List of 15 significant analytes
network_vars <- c(selected_analytes, "Delta_GHQ_total", "Delta_PSS_Global_total")
covariate_vars <- c("age_T0", "BMI_T0_new", "smoke_numeric", "sex_numeric")

# Extract Data
network_data <- delta_data_final[, names(delta_data_final) %in% network_vars]
covariate_data <- delta_data_final[, names(delta_data_final) %in% covariate_vars]

# Residualize the data to control for covariates
residualized_data <- network_data  # Copy for modification

for (var in network_vars) {
  formula <- as.formula(paste0("`", var, "` ~ ", paste(covariate_vars, collapse = " + ")))
  residualized_data[[var]] <- resid(lm(formula, data = delta_data_final))
}

# Function to remove outliers only for the specific variable (±2 SD)
remove_outliers_variablewise <- function(data) {
  for (col in names(data)) {
    if (is.numeric(data[[col]])) {
      m <- mean(data[[col]], na.rm = TRUE)
      s <- sd(data[[col]], na.rm = TRUE)
      data[[col]][data[[col]] < (m - 2 * s) | data[[col]] > (m + 2 * s)] <- NA
    }
  }
  return(data)
}

# Apply outlier removal
residualized_data <- remove_outliers_variablewise(residualized_data)

# Compute partial correlations using pairwise complete observations
cor_matrix <- cor(residualized_data, method = "pearson", use = "pairwise.complete.obs")

# Format labels
formatted_labels <- gsub("Delta_", "", colnames(cor_matrix))
formatted_labels <- gsub("_", " ", formatted_labels)
formatted_labels[formatted_labels == "GHQ total"] <- "GHQ"
formatted_labels[formatted_labels == "PSS Global total"] <- "PSS"
formatted_labels[formatted_labels == "Interleukin-4"] <- "IL-4"
formatted_labels[formatted_labels == "Vascular cell adhesion protein 1"] <- "VCAM-1"
formatted_labels[formatted_labels == "Platelet factor 4"] <- "CXCL4"
formatted_labels[formatted_labels == "Neutrophil collagenase"] <- "MMP-8"
formatted_labels[formatted_labels == "Growth-regulated alpha protein"] <- "CXCL1"
formatted_labels[formatted_labels == "Lipopolysaccharide-binding protein"] <- "LBP"
formatted_labels[formatted_labels == "Gro-beta"] <- "CXCL2"
formatted_labels[formatted_labels == "A disintegrin and metalloproteinase with thrombospondin motifs 9"] <- "ADAMTS9"
formatted_labels[formatted_labels == "Interleukin-17A"] <- "IL-17A"
formatted_labels[formatted_labels == "CD40 ligand"] <- "CD40L"
formatted_labels[formatted_labels == "Prostaglandin G/H synthase 2"] <- "PTGS2"
formatted_labels[formatted_labels == "C-X-C motif chemokine 10"] <- "CXCL10"
formatted_labels[formatted_labels == "Annexin A2"] <- "ANXA2"
formatted_labels[formatted_labels == "Indoleamine 2,3-dioxygenase 1"] <- "IDO-1"

# Define edge colors
edge_colors <- ifelse(cor_matrix > 0, "red",
                      ifelse(cor_matrix < 0, "blue", "grey"))

# Plot the network
qgraph(cor_matrix,
       layout = "spring",
       theme = "colorblind",
       labels = formatted_labels,
       label.cex = 2,
       color = "skyblue",
       minimum = 0.1,
       cut = 0.3,
       vsize = 8,
       esize = 20,
       edge.labels = TRUE,
       edge.label.cex = 1.2,
       title = "Partial Pearson Network (Variable-wise Outlier Removal)",
       legend.cex = 1.2,
       edge.color = edge_colors)
```


```{r spearman-network-plot, echo=TRUE, message=FALSE, warning=FALSE}
# Load libraries
library(qgraph)
library(ppcor)

# Define variables
selected_analytes <- significant_results$Analyte_Clean
network_vars <- c(selected_analytes, "Delta_GHQ_total", "Delta_PSS_Global_total")
covariate_vars <- c("age_T0", "BMI_T0_new", "smoke_numeric", "sex_numeric")

# Extract data
network_data <- delta_data_final[, names(delta_data_final) %in% network_vars]
covariate_data <- delta_data_final[, names(delta_data_final) %in% covariate_vars]

# Residualize data to adjust for covariates
residualized_data <- network_data
for (var in colnames(network_data)) {
  formula <- as.formula(paste0("`", var, "` ~ ", paste(covariate_vars, collapse = " + ")))
  residualized_data[[var]] <- resid(lm(formula, data = delta_data_final))
}

# Compute Spearman partial correlation
cor_result <- pcor(residualized_data, method = "spearman")

# Extract correlation matrix and ensure it's square and numeric
cor_matrix <- as.matrix(cor_result$estimate)

# Set custom edge colors
edge_colors <- ifelse(cor_matrix > 0, "firebrick", "steelblue")

# Clean labels for plot
formatted_labels <- gsub("Delta_", "", colnames(cor_matrix))
formatted_labels <- gsub("_", " ", formatted_labels)
formatted_labels[formatted_labels == "GHQ total"] <- "GHQ"
formatted_labels[formatted_labels == "PSS Global total"] <- "PSS"
formatted_labels[formatted_labels == "Interleukin-4"] <- "IL-4"
formatted_labels[formatted_labels == "Vascular cell adhesion protein 1"] <- "VCAM-1"
formatted_labels[formatted_labels == "Platelet factor 4"] <- "CXCL4"
formatted_labels[formatted_labels == "Neutrophil collagenase"] <- "MMP-8"
formatted_labels[formatted_labels == "Growth-regulated alpha protein"] <- "CXCL1"
formatted_labels[formatted_labels == "Lipopolysaccharide-binding protein"] <- "LBP"
formatted_labels[formatted_labels == "Gro-beta"] <- "CXCL2"
formatted_labels[formatted_labels == "A disintegrin and metalloproteinase with thrombospondin motifs 9"] <- "ADAMTS9"
formatted_labels[formatted_labels == "Interleukin-17A"] <- "IL-17A"
formatted_labels[formatted_labels == "CD40 ligand"] <- "CD40L"
formatted_labels[formatted_labels == "Prostaglandin G/H synthase 2"] <- "PTGS2"
formatted_labels[formatted_labels == "C-X-C motif chemokine 10"] <- "CXCL10"
formatted_labels[formatted_labels == "Annexin A2"] <- "ANXA2"
formatted_labels[formatted_labels == "Indoleamine 2,3-dioxygenase 1"] <- "IDO-1"

# Plot the network with increased node spacing
qgraph(cor_matrix,
       layout = "spring",
       layout.par = list(repulse.rad = 20), # directly sets radius of repulsion; default ~10
       # repulsion = 1.4,
       theme = "colorblind",
       labels = formatted_labels,
       label.cex = 1.2,
       color = "skyblue",
       border.width = 1.2,
       border.color = "grey40",
       minimum = 0.15,
       cut = 0.3,
       vsize = 7,                             
       esize = 18,
       # edge.labels = TRUE,
       edge.label.cex = 1.0,
       edge.label.bg = "white",
       edge.label.position = 0.5,             
       edge.label.margin = 0.01,              
       edge.color = edge_colors,
       title = "",
       legend.cex = 1.2
)

```
