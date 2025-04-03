library(ranger)
library(tidyverse)
library(data.table)
library(tidymodels)
library(clusterProfiler)
library(org.Mm.eg.db)
library(future)
library(rstudioapi)  # For setting script directory

# Function to load and preprocess data
load_and_preprocess_data <- function(file_path, treatment_info) {
  counts <- fread(file_path)
  
  # Transpose and convert to a data frame
  transposed_counts <- as.data.frame(t(counts))
  colnames(transposed_counts) <- transposed_counts[1, ]
  transposed_counts <- transposed_counts[-1, ]
  transposed_counts$Sample <- rownames(transposed_counts)
  rownames(transposed_counts) <- NULL  # Clear row names to avoid conflicts
  
  nt <- as.data.frame(sapply(transposed_counts, as.numeric))
  
  # Filter genes based on coefficient of variation
  goodbob <- which(apply(nt, 2, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) >= 0.5)
  nt$Sample <- transposed_counts$Sample
  transposed_counts <- nt[, c(which(colnames(nt) == 'Sample'), goodbob)]
  
  # User-defined treatment groups
  treatment_data <- fread(treatment_info)
  
  # Ensure that the treatment file contains a 'Sample' column that matches the counts data
  if (!"Sample" %in% colnames(treatment_data)) {
    stop("The treatment file must contain a 'Sample' column that matches the sample names in the counts file.")
  }
  
  # Merge treatment information with the transposed counts data
  transposed_counts <- left_join(transposed_counts, treatment_data, by = "Sample")
  
  return(transposed_counts)
}

# Function to run Random Forest model and compute importance scores
run_random_forest <- function(data, outcome_var, num_trees = 500, num_seeds = 10) {
  plan(multisession, workers = 8)  # Enable parallel processing
  
  seeds <- sample(1:1000, num_seeds)  # Generate random seeds
  all_importance <- list()
  all_errors <- numeric(num_seeds)
  
  for (i in seq_along(seeds)) {
    set.seed(seeds[i])
    data_split <- initial_split(data, prop = 2/3, strata = all_of(outcome_var))
    train_data <- training(data_split)
    test_data <- testing(data_split)
    
    gene_recipe <- recipe(as.formula(paste(outcome_var, "~ .")), data = train_data) %>%
      update_role(Sample, new_role = "ID") %>%
      step_normalize(all_numeric_predictors())
    
    pgrz <- prep(gene_recipe, training = train_data)
    train_data_prepped <- bake(pgrz, new_data = train_data)
    test_data_prepped <- bake(pgrz, new_data = test_data)
    
    x_train <- train_data_prepped[, !names(train_data_prepped) %in% c(outcome_var, "Sample")]
    y_train <- train_data_prepped[[outcome_var]]
    
    mod <- ranger(y = y_train, x = x_train, num.trees = num_trees, importance = 'permutation')
    all_errors[i] <- mod$prediction.error
    
    importance_df <- data.frame(
      Gene = names(mod$variable.importance),
      Importance = mod$variable.importance
    )
    
    all_importance[[i]] <- importance_df
  }
  
  importance_combined <- bind_rows(all_importance, .id = "Run")
  importance_summary <- importance_combined %>%
    group_by(Gene) %>%
    summarize(Mean_Importance = mean(Importance, na.rm = TRUE), 
              SD_Importance = sd(Importance, na.rm = TRUE)) %>%
    arrange(desc(Mean_Importance))
  
  list(
    top_genes = head(importance_summary, 200),
    mean_oob_error = mean(all_errors),
    sd_oob_error = sd(all_errors)
  )
}

# Function to perform KEGG pathway enrichment
perform_kegg_analysis <- function(top_genes) {
  converted_ids <- bitr(top_genes$Gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  merged_genes <- top_genes %>% left_join(converted_ids, by = c("Gene" = "ENSEMBL"))
  enrich_results <- enrichKEGG(gene = merged_genes$ENTREZID, organism = "mmu")
  return(as.data.frame(enrich_results))
}

# Function to visualize results
generate_plots <- function(kegg_results) {
  ggplot(kegg_results[1:10, ], aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust), fill = p.adjust)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(x = "Gene Set", y = "-log10(Adjusted P-value)", title = "Top 10 Enriched KEGG Pathways") +
    coord_flip() +
    theme_minimal()
}

# Main execution
message("Select your counts table file")
file_path <- file.choose()

message("Select your treatment group file")
treatment_file <- file.choose()

# Specify the column name in the treatment file representing your outcome variable
message("Enter the outcome variable name (the column name in the treatment file for the experimental groups):")
outcome_var <- readline()

# Load and preprocess the data
transposed_counts <- load_and_preprocess_data(file_path, treatment_file)

# Run the Random Forest model
rf_results <- run_random_forest(transposed_counts, outcome_var)

# Perform KEGG pathway analysis on top genes
kegg_results <- perform_kegg_analysis(rf_results$top_genes)

# Generate plots for KEGG results
generate_plots(kegg_results)

# Print the top genes from Random Forest
print(rf_results$top_genes)

