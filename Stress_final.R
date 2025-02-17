library(tidyverse)
library(ggplot2)
library(readr)

Stress_Masterfile <- read_delim("Stress_Masterfile.csv", 
                                delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(Stress_Masterfile)

# Remove rows where all values are NA
Stress_Masterfile <- Stress_Masterfile %>%
  filter(rowSums(is.na(.)) != ncol(.))

# Ensure analyte columns are correctly identified
analyte_columns <- colnames(Stress_Masterfile)[62:161]  # Adjust column indices if needed


# Replace spaces in column names with underscores
colnames(Stress_Masterfile) <- gsub(" ", "_", colnames(Stress_Masterfile))
analyte_columns <- gsub(" ", "_", analyte_columns)


# Remove Bio_IDs with NA in any analyte column
bio_ids_to_keep <- Stress_Masterfile %>%
  group_by(Bio_ID) %>%
  filter(if_all(all_of(analyte_columns), ~ !is.na(.))) %>%
  pull(Bio_ID) %>%
  unique()

# Keep only complete Bio_IDs
Stress_Masterfile_cleaned <- Stress_Masterfile %>%
  filter(Bio_ID %in% bio_ids_to_keep)

removed_bio_ids <- setdiff(unique(Stress_Masterfile$Bio_ID), bio_ids_to_keep)
cat("Removed Bio_IDs due to missing analyte values:\n")
print(removed_bio_ids)

# Filtering for T0 and T1
Stress_Masterfile_filtered <- Stress_Masterfile_cleaned %>%
  group_by(Bio_ID) %>%
  filter(any(session == "T0" & !is.na(GHQ_total)) &
           any(session == "T1" & !is.na(GHQ_total))) %>%
  ungroup()

# Recalculate BMI 
Stress_Masterfile_filtered <- Stress_Masterfile_filtered %>%
  mutate(BMI_T0_new = as.numeric(Weight_T0) / (as.numeric(Height_T0) / 100)^2)

# Calculate delta values
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


# write.csv(Stress_Masterfile_cleaned, "Stress_Masterfile_cleaned.csv", row.names = FALSE)
# write.csv(delta_data, "delta_data.csv", row.names = FALSE)


# Remove rows with NA in Delta_Neurofilament_light_polypeptide
delta_data_filtered <- delta_data %>%
  filter(!is.na(Delta_Neurofilament_light_polypeptide))

# Print removed Bio_IDs
removed_bio_ids <- setdiff(delta_data$Bio_ID, delta_data_filtered$Bio_ID)
cat("Removed Bio_IDs due to missing Delta_Neurofilament_light_polypeptide values:\n")
print(removed_bio_ids)

# write.csv(delta_data_filtered, "delta_data_filtered.csv", row.names = FALSE)

delta_data_final <- delta_data_filtered

# Extract covariates
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

# Merge covariates into delta_data_final
delta_data_final <- delta_data_final %>%
  left_join(covariates, by = "Bio_ID")

# write.csv(delta_data_final, "delta_data_final.csv", row.names = FALSE)


### GHQ & PSS Plots
library(cowplot)

# Extract relevant data from the cleaned dataset
violin_data <- Stress_Masterfile_cleaned %>%
  filter(session %in% c("T0", "T1")) %>%
  select(Bio_ID, session, GHQ_total, PSS_Global_total) %>%
  pivot_longer(cols = c(GHQ_total, PSS_Global_total), 
               names_to = "Measure", values_to = "Value")

# Rename "session" to "Time Point"
violin_data <- violin_data %>% rename(`Time Point` = session)

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

# Rename Measure values for clean titles
violin_data$Measure <- recode(violin_data$Measure, 
                              "GHQ_total" = "GHQ Score Across Timepoints",
                              "PSS_Global_total" = "PSS Score Across Timepoints")

violin_data <- violin_data %>%
  mutate(`Time Point` = factor(`Time Point`, levels = c("T0", "T1")))

# Violin Plot
plot_violin <- function(df, measure) {
  ggplot(df %>% filter(Measure == measure), 
         aes(x = `Time Point`, y = Value, fill = `Time Point`)) +
    geom_violin(trim = FALSE, alpha = 0.7, color = "black") +  # Add border color
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.2, color = "black") +  # Adjust jitter
    labs(title = measure, x = "Time Point", y = "Score") +
    theme_minimal(base_size = 16) +  
    theme(
      panel.grid = element_blank(),  # Remove grid
      plot.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.position = "top"
    ) +
    scale_fill_manual(values = c("T0" = "#FF7070", "T1" = "#709BFF"), drop = FALSE)  # Ensure correct mapping
}

# Create the violin plots
p1 <- plot_violin(violin_data, "GHQ Score Across Timepoints")
p2 <- plot_violin(violin_data, "PSS Score Across Timepoints")

# Combine plots using cowplot
combined_plot <- plot_grid(p1, p2, labels = c("A", "B"), ncol = 2, label_size = 20)

# Print the combined plot
print(combined_plot)



##### Linear Regression ####

library(dplyr)

start_col <- which(colnames(delta_data_final) == "Delta_Neurofilament_light_polypeptide")
analyte_columns <- colnames(delta_data_final)[start_col:(start_col + 99)]

# Empty dataframe to store results
regression_summary_table <- data.frame(
  Analyte = character(),
  Dependent_Variable = character(),
  Coefficient = numeric(),
  P_Value = numeric(),
  N_Subjects = integer(),
  stringsAsFactors = FALSE
)

# Function to run regression and store results
run_regression <- function(dep_var, analyte_name) {
  # Prepare the regression dataset
  regression_data <- delta_data_final %>%
    select(Bio_ID, all_of(dep_var), all_of(analyte_name), age_T0, BMI_T0_new, smoke_numeric, sex_numeric) %>%
    rename(Dependent_Var = all_of(dep_var), Analyte = all_of(analyte_name), 
           Age = age_T0, BMI = BMI_T0_new, Smoke = smoke_numeric, Sex = sex_numeric)
  
  # Remove outliers based on 2x standard deviation
  regression_data_clean <- regression_data %>%
    filter(
      Dependent_Var >= (mean(Dependent_Var, na.rm = TRUE) - 2 * sd(Dependent_Var, na.rm = TRUE)) &
        Dependent_Var <= (mean(Dependent_Var, na.rm = TRUE) + 2 * sd(Dependent_Var, na.rm = TRUE)) &
        Analyte >= (mean(Analyte, na.rm = TRUE) - 2 * sd(Analyte, na.rm = TRUE)) &
        Analyte <= (mean(Analyte, na.rm = TRUE) + 2 * sd(Analyte, na.rm = TRUE))
    )
  
  # Run linear regression
  model <- lm(Dependent_Var ~ Analyte + Age + BMI + Smoke + Sex, data = regression_data_clean)
  
  # Extract values
  coef_value <- summary(model)$coefficients["Analyte", "Estimate"]
  p_value <- summary(model)$coefficients["Analyte", "Pr(>|t|)"]
  n_subjects <- nrow(regression_data_clean)
  
  # Append results to dataframe
  return(data.frame(
    Analyte = analyte_name,
    Dependent_Variable = dep_var,
    Coefficient = coef_value,
    P_Value = p_value,
    N_Subjects = n_subjects
  ))
}

# Loop over each analyte for GHQ and PSS
for (analyte in analyte_columns) {
  regression_summary_table <- rbind(regression_summary_table, run_regression("Delta_GHQ_total", analyte))
  regression_summary_table <- rbind(regression_summary_table, run_regression("Delta_PSS_Global_total", analyte))
}

# Apply FDR correction 
regression_summary_table$FDR_P_Value <- p.adjust(regression_summary_table$P_Value, method = "BH")
view(regression_summary_table)

# Save the table as a CSV file
# write.csv(regression_summary_table, "HMZ_Publication/Regression_Summary_Table_FDR.csv", row.names = FALSE)


### Multiple Correction per Questionnaire ###

##### Linear Regression with BH & Bonferroni FDR Correction ####

library(dplyr)
start_col <- which(colnames(delta_data_final) == "Delta_Neurofilament_light_polypeptide")

analyte_columns <- colnames(delta_data_final)[start_col:(start_col + 99)]

regression_summary_GHQ <- data.frame(
  Analyte = character(),
  Dependent_Variable = character(),
  Coefficient = numeric(),
  P_Value = numeric(),
  N_Subjects = integer(),
  FDR_BH = numeric(),
  FDR_Bonferroni = numeric(),
  stringsAsFactors = FALSE
)

regression_summary_PSS <- data.frame(
  Analyte = character(),
  Dependent_Variable = character(),
  Coefficient = numeric(),
  P_Value = numeric(),
  N_Subjects = integer(),
  FDR_BH = numeric(),
  FDR_Bonferroni = numeric(),
  stringsAsFactors = FALSE
)

# Function to run regression and store results
run_regression <- function(dep_var, analyte_name) {
  # Prepare the regression dataset
  regression_data <- delta_data_final %>%
    select(Bio_ID, all_of(dep_var), all_of(analyte_name), age_T0, BMI_T0_new, smoke_numeric, sex_numeric) %>%
    rename(Dependent_Var = all_of(dep_var), Analyte = all_of(analyte_name), 
           Age = age_T0, BMI = BMI_T0_new, Smoke = smoke_numeric, Sex = sex_numeric)
  
  # Remove outliers based on 2x standard deviation
  regression_data_clean <- regression_data %>%
    filter(
      Dependent_Var >= (mean(Dependent_Var, na.rm = TRUE) - 2 * sd(Dependent_Var, na.rm = TRUE)) &
        Dependent_Var <= (mean(Dependent_Var, na.rm = TRUE) + 2 * sd(Dependent_Var, na.rm = TRUE)) &
        Analyte >= (mean(Analyte, na.rm = TRUE) - 2 * sd(Analyte, na.rm = TRUE)) &
        Analyte <= (mean(Analyte, na.rm = TRUE) + 2 * sd(Analyte, na.rm = TRUE))
    )
  
  # Run linear regression
  model <- lm(Dependent_Var ~ Analyte + Age + BMI + Smoke + Sex, data = regression_data_clean)
  
  # Extract values
  coef_value <- summary(model)$coefficients["Analyte", "Estimate"]
  p_value <- summary(model)$coefficients["Analyte", "Pr(>|t|)"]
  n_subjects <- nrow(regression_data_clean)
  
  # Return results as a data frame
  return(data.frame(
    Analyte = analyte_name,
    Dependent_Variable = dep_var,
    Coefficient = coef_value,
    P_Value = p_value,
    N_Subjects = n_subjects
  ))
}

# Loop over each analyte for GHQ and store results
for (analyte in analyte_columns) {
  regression_summary_GHQ <- rbind(regression_summary_GHQ, run_regression("Delta_GHQ_total", analyte))
}

# Loop over each analyte for PSS and store results
for (analyte in analyte_columns) {
  regression_summary_PSS <- rbind(regression_summary_PSS, run_regression("Delta_PSS_Global_total", analyte))
}

# Apply both Benjamini-Hochberg and Bonferroni corrections separately for GHQ and PSS
regression_summary_GHQ$FDR_BH <- p.adjust(regression_summary_GHQ$P_Value, method = "BH")
regression_summary_GHQ$FDR_Bonferroni <- p.adjust(regression_summary_GHQ$P_Value, method = "bonferroni")

regression_summary_PSS$FDR_BH <- p.adjust(regression_summary_PSS$P_Value, method = "BH")
regression_summary_PSS$FDR_Bonferroni <- p.adjust(regression_summary_PSS$P_Value, method = "bonferroni")


view(regression_summary_GHQ)
view(regression_summary_PSS)

# write.csv(regression_summary_GHQ, "HMZ_Publication/Regression_Summary_Table_GHQ_FDR.csv", row.names = FALSE)
# write.csv(regression_summary_PSS, "HMZ_Publication/Regression_Summary_Table_PSS_FDR.csv", row.names = FALSE)



#### Heatmap ####

library(ggplot2)
library(dplyr)

# Prepare data for heatmap
heatmap_data <- regression_summary_table %>%
  select(Analyte, Dependent_Variable, Coefficient, P_Value) %>%
  mutate(Significance = case_when(
    P_Value < 0.001 ~ "***",
    P_Value < 0.01 ~ "**",
    P_Value < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Convert Dependent_Variable to factor 
heatmap_data$Dependent_Variable <- factor(heatmap_data$Dependent_Variable, 
                                          levels = c("Delta_GHQ_total", "Delta_PSS_Global_total"))

# Sort analytes based on GHQ coefficient
ghq_sorted <- heatmap_data %>%
  filter(Dependent_Variable == "Delta_GHQ_total") %>%
  arrange(desc(Coefficient)) %>%
  mutate(Analyte = gsub("Delta_", "", Analyte),  # Remove "Delta" prefix
         Analyte = gsub("_", " ", Analyte),  # Replace underscores with spaces
         Analyte = gsub("\\.", " ", Analyte)) %>%  # Replace "." with spaces
  pull(Analyte)  

# Apply sorted order & clean names
heatmap_data <- heatmap_data %>%
  mutate(Analyte = gsub("Delta_", "", Analyte),  # Remove "Delta" prefix
         Analyte = gsub("_", " ", Analyte),  # Replace underscores with spaces
         Analyte = gsub("\\.", " ", Analyte),  # Replace "." with spaces
         Dependent_Variable = recode(Dependent_Variable, 
                                     "Delta_GHQ_total" = "GHQ", 
                                     "Delta_PSS_Global_total" = "PSS"),
         Analyte = factor(Analyte, levels = ghq_sorted))  # Preserve sorted order

# Heatmap plot 
heatmap_plot <- ggplot(heatmap_data, aes(x = Dependent_Variable, y = Analyte, fill = Coefficient)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#0014a8", mid = "#ebf5ff", high = "#e32636", midpoint = 0, name = "Regression Coefficient") +  
  geom_text(aes(label = Significance), color = "black", size = 3, vjust = 0.75) +  # Lowered symbols
  labs(title = "Regression Coefficients Heatmap", x = "Dependent Variable", y = "Analytes") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 7, face = "bold", vjust = 0.5),
        axis.text.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"))

# Save & display
# ggsave("HMZ_Publication/Regression_Heatmap.png", plot = heatmap_plot, width = 8, height = 12, dpi = 300)
print(heatmap_plot)



#### Significant Results ####


# Extract significant correlations (p < 0.05)
significant_results <- regression_summary_table %>%
  filter(P_Value < 0.05) %>%
  mutate(Analyte = gsub("Delta_", "", Analyte),  # Remove "Delta" prefix
         Analyte = gsub("_", " ", Analyte),  # Replace underscores with spaces
         Analyte = gsub("\\.", " ", Analyte),  # Remove dots
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

 print(significant_results)


#### Regression Plots ####

library(ggplot2)
library(dplyr)

# Define correct column names
dep_var_col <- "Delta_GHQ_total"
analyte_col <- "Delta_Neutrophil_collagenase"

# Extract relevant data
plot_data <- delta_data_final %>%
  select(Bio_ID, all_of(dep_var_col), all_of(analyte_col)) %>%
  rename(GHQ = all_of(dep_var_col), Neutrophil_Collagenase = all_of(analyte_col))

# Outlier removal: Keep values within 2 standard deviations
plot_data <- plot_data %>%
  filter(
    GHQ >= (mean(GHQ, na.rm = TRUE) - 2 * sd(GHQ, na.rm = TRUE)) &
      GHQ <= (mean(GHQ, na.rm = TRUE) + 2 * sd(GHQ, na.rm = TRUE)) &
      Neutrophil_Collagenase >= (mean(Neutrophil_Collagenase, na.rm = TRUE) - 2 * sd(Neutrophil_Collagenase, na.rm = TRUE)) &
      Neutrophil_Collagenase <= (mean(Neutrophil_Collagenase, na.rm = TRUE) + 2 * sd(Neutrophil_Collagenase, na.rm = TRUE))
  )

# Clean analyte name to match significant_results
clean_analyte_name <- gsub("Delta_", "", analyte_col)  # Remove "Delta"
clean_analyte_name <- gsub("_", " ", clean_analyte_name)  # Replace underscores with spaces
clean_analyte_name <- gsub("\\.", " ", clean_analyte_name)  # Replace dots with spaces

# Extract regression values from significant_results
signif_result <- significant_results %>%
  filter(Analyte == clean_analyte_name, Dependent_Variable == "GHQ")

# Check if a match was found
if (nrow(signif_result) == 0) {
  stop("Error: No matching entry found in significant_results for ", clean_analyte_name, " vs. GHQ")
}

# Get values from the table
n_subjects <- nrow(plot_data)  # Number of subjects after outlier removal
beta_value <- round(signif_result$Coefficient, 3)  # Round beta to 3 decimals
p_value <- signif_result$P_Value

# Format p-value for display (avoid scientific notation)
if (!is.na(p_value) && p_value < 0.001) {
  p_value_text <- "p < 0.001"
} else {
  p_value_text <- paste0("p = ", round(p_value, 3))
}

# Create annotation text
annotation_text <- paste0("N = ", n_subjects, 
                          "\nβ = ", beta_value, 
                          "\n", p_value_text)

# Create scatter plot with updated aesthetics
p <- ggplot(plot_data, aes(x = Neutrophil_Collagenase, y = GHQ)) +
  geom_point(alpha = 0.9, size = 2, color = "#a00053") +  # Smaller black data points
  geom_smooth(method = "lm", color = "#00a09d", linetype = "solid") +  # Blue regression line
  annotate("text", x = min(plot_data$Neutrophil_Collagenase), 
           y = max(plot_data$GHQ), label = annotation_text, 
           hjust = 0, vjust = 1, size = 5, fontface = "bold") +  # Annotation box
  labs(title = "Neutrophil Collagenase vs. GHQ",
       x = "Neutrophil Collagenase",
       y = "GHQ") +
  theme_minimal(base_size = 16) +  # Larger text for readability
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

# Save & display plot
#ggsave("HMZ_Publication/Neutrophil_Collagenase_vs_GHQ.png", plot = p, width = 6, height = 4, dpi = 300)
print(p)


#### Loop LR Plots ####

library(ggplot2)
library(dplyr)


# Rename the columns to replace periods with underscores
names(delta_data_final) <- gsub("\\.", "_", names(delta_data_final))

# Check the new column names
colnames(delta_data_final)


# Print all available column names in delta_data_final 
cat("\n🔍 Available columns in delta_data_final:\n", paste(names(delta_data_final), collapse = ", "), "\n\n")

# Loop through all significant correlations
for (i in 1:nrow(significant_results)) {
  
  # Extract analyte and dependent variable
  analyte_name <- significant_results$Analyte_Clean[i]  # Already cleaned
  dep_var <- significant_results$Dependent_Variable[i]
  
  # Map GHQ and PSS back to their original column names
  dep_var_col <- ifelse(dep_var == "GHQ", "Delta_GHQ_total", "Delta_PSS_Global_total")
  
  # Create multiple possible name formats
  possible_names <- c(
    gsub(" ", "_", analyte_name),    # Replace spaces with underscores
    gsub(" ", ".", analyte_name),    # Replace spaces with dots
    paste0("Delta_", analyte_name),  # Add "Delta_" prefix
    paste0("Delta_", gsub(" ", "_", analyte_name)),  # Delta_ with underscores
    paste0("Delta_", gsub(" ", ".", analyte_name))   # Delta_ with dots
  )
  
  # Find matching column name in delta_data_final
  analyte_col_match <- names(delta_data_final)[names(delta_data_final) %in% possible_names]
  
  # If no match found, print the problem and continue
  if (length(analyte_col_match) != 1) {
    cat("\n⚠ Skipping:", analyte_name, "- No unique match found in delta_data_final.\n")
    cat("🔍 Searched for:", paste(possible_names, collapse = ", "), "\n")
    cat("📝 Available analyte columns:\n", paste(names(delta_data_final), collapse = ", "), "\n")
    next
  }
  
  # Extract relevant data
  plot_data <- delta_data_final %>%
    select(Bio_ID, all_of(dep_var_col), all_of(analyte_col_match)) %>%
    rename(Dependent_Var = all_of(dep_var_col), Analyte = all_of(analyte_col_match))
  
  # Outlier removal: Keep values within 2 standard deviations
  plot_data <- plot_data %>%
    filter(
      Dependent_Var >= (mean(Dependent_Var, na.rm = TRUE) - 2 * sd(Dependent_Var, na.rm = TRUE)) &
        Dependent_Var <= (mean(Dependent_Var, na.rm = TRUE) + 2 * sd(Dependent_Var, na.rm = TRUE)) &
        Analyte >= (mean(Analyte, na.rm = TRUE) - 2 * sd(Analyte, na.rm = TRUE)) &  # Fixed typo here
        Analyte <= (mean(Analyte, na.rm = TRUE) + 2 * sd(Analyte, na.rm = TRUE))   # Fixed typo here
    )
  
  # Extract values from significant_results
  beta_value <- round(significant_results$Coefficient[i], 3)
  p_value <- significant_results$P_Value[i]
  n_subjects <- nrow(plot_data)
  
  # Format p-value for display
  if (!is.na(p_value) && p_value < 0.001) {
    p_value_text <- "p < 0.001"
  } else {
    p_value_text <- paste0("p = ", round(p_value, 3))
  }
  
  # Create annotation text
  annotation_text <- paste0("N = ", n_subjects, 
                            "\nβ = ", beta_value, 
                            "\n", p_value_text)
  
  # Create scatter plot
  p <- ggplot(plot_data, aes(x = Analyte, y = Dependent_Var)) +
    geom_point(alpha = 0.7, size = 2, color = "black") +  # Smaller black data points
    geom_smooth(method = "lm", color = "#1f78b4", linetype = "dashed") +  # Blue regression line
    annotate("text", x = min(plot_data$Analyte, na.rm = TRUE), 
             y = max(plot_data$Dependent_Var, na.rm = TRUE), label = annotation_text, 
             hjust = 0, vjust = 1, size = 5, fontface = "bold") +  # Annotation box
    labs(title = paste(analyte_name, "vs.", dep_var),
         x = analyte_name,
         y = dep_var) +
    theme_minimal(base_size = 16) +  # Larger text for readability
    theme(plot.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16))
  
  # Print the plot in RStudio
  print(p)
  
  # Define filename for saving
  filename <- paste0("HMZ_Publication/", gsub(" ", "_", analyte_name), "_vs_", dep_var, ".png")
  
  # Save the plot
  ggsave(filename, plot = p, width = 6, height = 4, dpi = 300)
  
  # Print status message
  cat("\n✅ Saved plot:", filename, "\n")
}




####### Significant Analytes <0.1 #######

library(dplyr)

# Get analytes from the regression table (p < 0.1)
significant_analytes_p0_1 <- c(
  "Delta_Neutrophil_collagenase", "Delta_Growth_regulated_alpha_protein", 
  "Delta_Gro_beta", "Delta_Indoleamine_2_3_dioxygenase_1", 
  "Delta_Lipopolysaccharide_binding_protein", "Delta_Vascular_cell_adhesion_protein_1", 
  "Delta_Annexin_A2", "Delta_C_X_C_motif_chemokine_10", 
  "Delta_A_disintegrin_and_metalloproteinase_with_thrombospondin_motifs_9", 
  "Delta_Interleukin_17A", "Delta_Prostaglandin_G_H_synthase_2", 
  "Delta_Interleukin_4", "Delta_CD40_ligand", "Delta_Eotaxin", 
  "Delta_Platelet_factor_4", "Delta_Fibroblast_growth_factor_22", 
  "Delta_C_X_C_motif_chemokine_6", "Delta_Brain_derived_neurotrophic_factor", 
  "Delta_A_disintegrin_and_metalloproteinase_with_thrombospondin_motifs_5", 
  "Delta_Interferon_gamma", "Delta_Metalloproteinase_inhibitor_1", 
  "Delta_iNOS", "Delta_Appetite_regulating_hormone", 
  "Delta_Neutrophil_activating_peptide_2", "Delta_Gro_gamma", 
  "Delta_P_selectin", "Delta_Indoleamine_2_3_dioxygenase_1", 
  "Delta_von_Willebrand_factor"
)

# Find the matching column names in delta_data_final
matching_columns_p0_1 <- names(delta_data_final)[names(delta_data_final) %in% significant_analytes_p0_1]

# Select the relevant data (GHQ, PSS, covariates, and analyte delta values)
significant_data_p0_1 <- delta_data_final %>%
  select(Bio_ID, Delta_GHQ_total, Delta_PSS_Global_total, 
         all_of(matching_columns_p0_1),  # Inserted analytes
         age_T0, sex_numeric, BMI_T0_new, smoke_numeric, meds_numeric)

# Save as CSV
# write.csv(significant_data_p0_1, "HMZ_Publication/Significant_Correlates_p0_1.csv", row.names = FALSE)






### Circular Plot ####

### CIRCULAR NETWORK PLOT ###

# Load necessary libraries
library(ggplot2)
library(scales)
library(igraph)
library(ggraph)

# Define the proteins and their respective functional groups
protein_groups <- c(
  "Neutrophil collagenase" = "inflammatory",
  "Growth regulated alpha protein" = "inflammatory",
  "Gro beta" = "inflammatory",
  "Indoleamine 2 3 dioxygenase 1" = "inflammatory",
  "Lipopolysaccharide binding protein" = "inflammatory",
  "Vascular cell adhesion protein 1" = "vascular",
  "Annexin A2" = "inflammatory",
  "C X C motif chemokine 10" = "inflammatory",
  "A disintegrin and metalloproteinase with thrombospondin motifs 9" = "metabolic",
  "Interleukin 17A" = "inflammatory",
  "Prostaglandin G H synthase 2" = "inflammatory",
  "Interleukin 4" = "inflammatory",
  "CD40 ligand" = "inflammatory",
  "Eotaxin" = "inflammatory",
  "Platelet factor 4" = "inflammatory",
  "Fibroblast growth factor 22" = "metabolic",
  "C X C motif chemokine 6" = "inflammatory",
  "Brain derived neurotrophic factor" = "metabolic",
  "A disintegrin and metalloproteinase with thrombospondin motifs 5" = "metabolic",
  "Interferon gamma" = "inflammatory",
  "Metalloproteinase inhibitor 1" = "inflammatory",
  "iNOS" = "inflammatory",
  "Appetite regulating hormone" = "metabolic",
  "Neutrophil activating peptide 2" = "inflammatory",
  "Gro gamma" = "inflammatory",
  "P selectin" = "vascular",
  "Indoleamine 2 3 dioxygenase 1" = "inflammatory",
  "von Willebrand factor" = "vascular"
)

# Define functional groups for analytes (Modify this based on your dataset)
functional_groups <- c(
  "inflammatory" = "skyblue1",
  "metabolic" = "orange1",
  "vascular" = "palegreen2"
)

# Create a data frame representing the relationships (edges) between analytes
edges <- data.frame(
  from = sample(names(protein_groups), length(protein_groups), replace = TRUE),  # Random connections (Modify if needed)
  to = names(protein_groups),  # Target nodes
  weight = runif(length(protein_groups), -0.5, 0.5)  # Random correlation values for example
)

# Create the graph object from the edges data frame
network_graph <- graph_from_data_frame(edges, directed = FALSE)

# Ensure all proteins in the network are assigned to a group
V(network_graph)$category <- protein_groups[V(network_graph)$name]

# Convert node categories to colors
node_colors <- setNames(functional_groups[V(network_graph)$category], V(network_graph)$name)

# Refined Circular Layout Visualization
ggraph(network_graph, layout = "circle") +  
  geom_edge_link(aes(edge_alpha = 0.7 * abs(weight), edge_width = weight, color = weight)) +
  scale_edge_width_continuous(range = c(0.3, 1.2)) +  # Adjust edge thickness
  scale_edge_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +  # Edge colors
  geom_node_point(aes(color = category), size = 10, alpha = 0.9) +  # Node colors
  scale_color_manual(values = functional_groups) +  
  geom_node_text(aes(label = name), size = 6, fontface = "bold", color = "black") +  # Labels
  theme_void() +
  theme(legend.position = "right", 
        legend.key.size = unit(0.1, "cm"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

# Save the plot
ggsave("Circular_Network_Plot.png", width = 10, height = 8, dpi = 300)

