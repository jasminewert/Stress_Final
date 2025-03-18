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

#write.csv(regression_summary_GHQ, "HMZ_Publication/Regression_Summary_Table_GHQ_FDR.csv", row.names = FALSE)
#write.csv(regression_summary_PSS, "HMZ_Publication/Regression_Summary_Table_PSS_FDR.csv", row.names = FALSE)



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
  geom_text(aes(label = Significance), color = "black", size = 5, vjust = 0.75) +  # Lowered symbols
  labs(title = "Regression Coefficients Heatmap", x = "", y = "Analytes") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 7, face = "bold", vjust = 0.5),
        axis.text.x = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"))

# Save & display
ggsave("HMZ_Publication/Regression_Heatmap.png", plot = heatmap_plot, width = 8, height = 12, dpi = 300)
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
# Ensure analyte names match the format in delta_data_final
significant_results <- significant_results %>%
  mutate(Analyte_Clean = paste0("Delta_", gsub(" ", "_", Analyte)))  # Add "Delta_" prefix and replace spaces with "_"

# Find which analytes exist in delta_data_final
available_analytes <- intersect(significant_results$Analyte_Clean, colnames(delta_data_final))
missing_analytes <- setdiff(significant_results$Analyte_Clean, colnames(delta_data_final))

# Print results
cat("\n✅ Available analytes for plotting:", length(available_analytes), "\n")
print(available_analytes)

cat("\n❌ Missing analytes (not found in delta_data_final):", length(missing_analytes), "\n")
print(missing_analytes)

 
##Loop

library(ggplot2)
library(dplyr)

# Loop through all significant analytes and generate plots
for (i in 1:nrow(significant_results)) {
  
  analyte_name <- significant_results$Analyte_Clean[i]
  dep_var <- significant_results$Dependent_Variable[i]
  
  # Check if analyte exists in delta_data_final
  if (!(analyte_name %in% available_analytes)) {
    cat("\n⚠ Skipping:", analyte_name, "→ Not found in `delta_data_final`.\n")
    next
  }
  
  # Select correct dependent variable column
  dep_var_col <- ifelse(dep_var == "GHQ", "Delta_GHQ_total", "Delta_PSS_Global_total")
  
  # Extract relevant data
  plot_data <- delta_data_final %>%
    select(Bio_ID, all_of(dep_var_col), all_of(analyte_name)) %>%
    rename(Dependent_Var = all_of(dep_var_col), Analyte = all_of(analyte_name))
  
  # Remove outliers (values within 2 SD)
  plot_data <- plot_data %>%
    filter(
      Dependent_Var >= (mean(Dependent_Var, na.rm = TRUE) - 2 * sd(Dependent_Var, na.rm = TRUE)) &
        Dependent_Var <= (mean(Dependent_Var, na.rm = TRUE) + 2 * sd(Dependent_Var, na.rm = TRUE)) &
        Analyte >= (mean(Analyte, na.rm = TRUE) - 2 * sd(Analyte, na.rm = TRUE)) &
        Analyte <= (mean(Analyte, na.rm = TRUE) + 2 * sd(Analyte, na.rm = TRUE))
    )
  
  # Skip if no valid data points remain
  if (nrow(plot_data) == 0) {
    cat("\n⚠ Skipping:", analyte_name, "→ No valid data points after filtering!\n")
    next
  }
  
  # Extract p-value
  p_value <- significant_results$P_Value[i]
  
  # Format p-value
  p_value_text <- ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", round(p_value, 3)))
  
  # Position p-value at upper-right
  annotation_x <- max(plot_data$Analyte, na.rm = TRUE) * 0.9  
  annotation_y <- max(plot_data$Dependent_Var, na.rm = TRUE) * 0.9  
  
  # Create scatter plot
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
      title = " ",  
      x = bquote(.(significant_results$Analyte[i]) ~ "[" ~ Delta ~ T[1]-T[0] ~ "]"),  
      y = bquote(.(dep_var) ~ "Total Score [" ~ Delta ~ T[1]-T[0] ~ "]")  # ✅ Updated y-axis title
    ) +
    theme_minimal(base_size = 16) +  
    theme(
      plot.title = element_blank(),  
      axis.text = element_text(size = 14),  
      axis.title = element_text(size = 16),  
      panel.grid = element_blank()  
    )
  
  # Print the plot
  print(p)
  
  # Save plot
  filename <- paste0("HMZ_Publication/", gsub(" ", "_", significant_results$Analyte[i]), "_vs_", dep_var, ".png")
  ggsave(filename, plot = p, width = 6, height = 4, dpi = 300)
  
  cat("\n✅ Saved plot:", filename, "\n")
} 

### Circular Plot ####

### CIRCULAR NETWORK PLOT ###

# Load necessary libraries
library(ggplot2)
library(scales)
library(igraph)
library(ggraph)

# Define the 15 significant proteins with their respective functional groups, matching delta_data_final
protein_groups <- c(
  "Delta_Neutrophil_collagenase" = "inflammatory",
  "Delta_Growth-regulated_alpha_protein" = "inflammatory",
  "Delta_Gro-beta" = "inflammatory",
  "Delta_Indoleamine_2,3-dioxygenase_1" = "inflammatory",
  "Delta_Lipopolysaccharide-binding_protein" = "inflammatory",
  "Delta_Vascular_cell_adhesion_protein_1" = "vascular",
  "Delta_Annexin_A2" = "inflammatory",
  "Delta_C-X-C_motif_chemokine_10" = "inflammatory",
  "Delta_A_disintegrin_and_metalloproteinase_with_thrombospondin_motifs_9" = "metabolic",
  "Delta_Interleukin-17A" = "inflammatory",
  "Delta_Prostaglandin_G/H_synthase_2" = "inflammatory",
  "Delta_Interleukin-4" = "inflammatory",
  "Delta_CD40_ligand" = "inflammatory",
  "Delta_Eotaxin" = "inflammatory",
  "Delta_Platelet_factor_4" = "inflammatory"
)

# Print to verify
print(protein_groups)

# Define functional groups for analytes (Modify this based on your dataset)
functional_groups <- c(
  "inflammatory" = "skyblue1",
  "metabolic" = "orange1",
  "vascular" = "palegreen2"
)


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggcorrplot)

# Step 1: Extract relevant protein data from delta_data_final
selected_proteins <- names(protein_groups)  # Get the 15 significant proteins

# Create a new dataset with only the selected proteins
correlation_data <- delta_data_final %>%
  select(all_of(selected_proteins))

# Step 2: Rename proteins (remove "Delta" & "_")
clean_names <- gsub("Delta_", "", selected_proteins)  # Remove "Delta_"
clean_names <- gsub("_", " ", clean_names)  # Replace "_" with spaces
colnames(correlation_data) <- clean_names  # Apply cleaned names

# Step 3: Remove outliers (values outside ±2 SD)
remove_outliers <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  x[x < (mean_x - 2 * sd_x) | x > (mean_x + 2 * sd_x)] <- NA
  return(x)
}

correlation_data <- correlation_data %>%
  mutate(across(everything(), remove_outliers))

# Step 4: Compute the correlation matrix
correlation_matrix <- cor(correlation_data, use = "pairwise.complete.obs", method = "pearson")

# Step 5: Generate the heatmap
heatmap_plot <- ggcorrplot(correlation_matrix, 
                           type = "lower", 
                           lab = TRUE, 
                           lab_size = 2, 
                           colors = c("#0014a8", "#ebf5ff", "#e32636"),
                           outline.color = "black",
                           show.legend = TRUE) +
  labs(title = "Correlation Matrix of Significant Proteins") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the heatmap
print(heatmap_plot)

# Save the heatmap
#ggsave("HMZ_Publication/Protein_Correlation_Heatmap.png", plot = heatmap_plot, width = 10, height = 8, dpi = 300)


### Network Plot


# Load necessary libraries
library(igraph)
library(ggraph)
library(tidyverse)

# Convert correlation matrix to a tidy format
correlation_long <- as.data.frame(as.table(correlation_matrix)) %>%
  rename(Protein1 = Var1, Protein2 = Var2, Correlation = Freq) %>%
  filter(Protein1 != Protein2)  # Remove self-loops

# Set a correlation threshold (e.g., abs(r) > 0.3)
threshold <- 0.25
edges <- correlation_long %>%
  filter(abs(Correlation) > threshold)

# Create nodes data frame
nodes <- data.frame(name = unique(c(edges$Protein1, edges$Protein2))) %>%
  mutate(category = protein_groups[name])  # Assign functional groups

# Define node colors based on category
functional_colors <- c("inflammatory" = "skyblue1", "metabolic" = "orange1", "vascular" = "palegreen2")

# Load stringr for text formatting
library(stringr)

# Modify the names manually with explicit line breaks (as done previously)
nodes <- nodes %>%
  mutate(label = case_when(
    name == "Prostaglandin G/H synthase 2" ~ "Prostaglandin\nG/H synthase 2",
    name == "A disintegrin and metalloproteinase with thrombospondin motifs 9" ~ "A disintegrin and metalloproteinase\nwith thrombospondin motifs 9",  # Explicit line breaks for three lines
    name == "Vascular cell adhesion protein 1" ~ "Vascular cell\nadhesion protein 1",
    name == "Lipopolysaccharide-binding protein" ~ "Lipopolysaccharide\nbinding protein",  # Now split into 2 lines
    name == "Neutrophil collagenase" ~ "Neutrophil\ncollagenase",
    name == "Growth-regulated alpha protein" ~ "Growth-regulated\nalpha protein",
    name == "Platelet factor 4" ~ "Platelet\nfactor 4",
    name == "C-X-C motif chemokine 10" ~ "C-X-C motif\nchemokine 10",
    name == "Indoleamine 2,3-dioxygenase 1" ~ "Indoleamine 2,3-\ndioxygenase 1",
    TRUE ~ name  # Keep other names unchanged
  ))

# Create igraph object
network_graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# Load stringr for text wrapping
library(stringr)

# Define the file path
file_path <- "C:/Users/jasmi/OneDrive/Dokumente/SomaLogic/HMZ_Publication/Protein_Network_Plot.png"

# Define the node sizes for each protein (Width, Height) with matching labels and adjustments
node_sizes <- data.frame(
  label = c("Interleukin-4", "Neutrophil collagenase", "Prostaglandin\nG/H synthase", 
            "Vascular cell\nadhesion protein 1", "CD40 ligand", "Eotaxin", 
            "Gro-beta", "Platelet\nfactor 4", "Indoleamine 2,3-\ndioxygenase 1", 
            "C-X-C motif\nchemokine 10", "Lipopolysaccharide\nbinding protein", "Interleukin-17",
            "Annexin A2", "Growth-regulated\nalpha protein", 
            "A disintegrin and metalloproteinase\nwith thrombospondin motifs 9"), 
  width = c(0.5, 0.6, 0.7, 0.8, 0.4, 0.3, 0.3, 0.4, 0.7, 0.6, 0.7, 0.8,
            0.4, 0.6, 1.1),  # Increased width significantly for "A disintegran...", reduced for "Growth-regulated alpha protein" and increased for "Vascular cell adhesion protein"
  height = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
             0.2, 0.2, 0.2)  # Height for each protein (adjusted as requested)
)

# Generate the network plot with customized node sizes for each protein
network_plot <- ggraph(network_graph, layout = "fr") +  # Using Fruchterman-Reingold layout for more space
  geom_edge_link(aes(edge_alpha = abs(Correlation), edge_width = abs(Correlation), edge_color = Correlation)) +
  scale_edge_width_continuous(range = c(1, 5), guide = "none") +  
  scale_edge_alpha_continuous(guide = "none") +  
  scale_edge_color_gradient2(
    low = "#0014a8", mid = "#ebf5ff", high = "#e32636", midpoint = 0,
    guide = "none"
  ) +  
  # Create rectangular nodes with width and height adjusted for each protein and black borders
  geom_tile(aes(x = x, y = y, 
                width = ifelse(label %in% node_sizes$label, 
                               node_sizes$width[match(label, node_sizes$label)], 0.5), 
                height = ifelse(label %in% node_sizes$label, 
                                node_sizes$height[match(label, node_sizes$label)], 0.2), 
                fill = category), 
            color = "black", # Add black border
            show.legend = FALSE) +  
  scale_fill_manual(values = functional_colors, guide = "none") +  
  # Add protein names inside the nodes
  geom_node_text(aes(x = x, y = y, label = label), size = 6, fontface = "bold", color = "black") +  
  theme_void() +
  theme(
    plot.margin = margin(400, 400, 400, 400),  # Increase margins significantly for more space
    legend.position = "none"  
  )

# Print the adjusted plot with customized node sizes and black borders
print(network_plot)

# Save the plot with customized node sizes and black borders (larger figure)
ggsave(file_path, plot = network_plot, width = 25, height = 25, dpi = 300, limitsize = FALSE)  # Save with larger figure and nodes


###

# Load necessary libraries
library(igraph)
library(ggraph)
library(tidyverse)

# Step 1: Add GHQ and PSS to the dataset
selected_variables <- c(names(protein_groups), "Delta_GHQ_total", "Delta_PSS_Global_total")

# Create dataset including GHQ and PSS
correlation_data <- delta_data_final %>%
  select(all_of(selected_variables))

# Clean variable names
clean_names <- gsub("Delta_", "", selected_variables)  
clean_names <- gsub("_", " ", clean_names)  
colnames(correlation_data) <- clean_names  

# Remove outliers
remove_outliers <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  x[x < (mean_x - 2 * sd_x) | x > (mean_x + 2 * sd_x)] <- NA
  return(x)
}

correlation_data <- correlation_data %>%
  mutate(across(everything(), remove_outliers))

# Step 2: Compute the correlation matrix and print it
correlation_matrix <- cor(correlation_data, use = "pairwise.complete.obs", method = "pearson")

# Print the correlation matrix
print(correlation_matrix)

# Step 3: Convert correlation matrix to a long format for network analysis
correlation_long <- as.data.frame(as.table(correlation_matrix)) %>%
  rename(Protein1 = Var1, Protein2 = Var2, Correlation = Freq) %>%
  filter(Protein1 != Protein2)  

# Step 4: Apply correlation threshold
threshold <- 0.25
edges <- correlation_long %>%
  filter(abs(Correlation) > threshold)

# Step 9: Generate the **final** corrected network plot with **extremely tiny** nodes
network_plot <- ggraph(network_graph, layout = layout_matrix) +
  geom_edge_link(aes(edge_alpha = abs(Correlation), edge_width = abs(Correlation), edge_color = Correlation)) +
  scale_edge_width_continuous(range = c(0.5, 2), guide = "none") +  
  scale_edge_alpha_continuous(guide = "none") +
  scale_edge_color_gradient2(
    low = "#0014a8", mid = "#ebf5ff", high = "#e32636", midpoint = 0,
    guide = "none"
  ) +
  # **Reduce node size to be as small as possible**
  geom_tile(aes(x = x, y = y, 
                width = ifelse(name %in% c("GHQ", "PSS"), 0.5, 0.3),  # **Super tiny nodes**
                height = ifelse(name %in% c("GHQ", "PSS"), 0.25, 0.2),  
                fill = category),  
            color = "black", show.legend = FALSE) +
  scale_fill_manual(values = functional_colors, guide = "none") +
  # **Keep text readable**
  geom_node_text(aes(x = x, y = y, label = label), size = 3.5, fontface = "bold", color = "black") +  
  theme_void() +
  theme(
    plot.margin = margin(50, 50, 50, 50),  
    legend.position = "none"
  )

# Step 10: Save the **final** corrected plot
ggsave("Protein_Network_Plot_Fixed.png", plot = network_plot, 
       width = 18, height = 18, dpi = 500, limitsize = FALSE)



