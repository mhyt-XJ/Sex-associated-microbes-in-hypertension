# Modified Mediation Analysis Script with 95% Confidence Intervals

# Load required packages
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("mediation")) install.packages("mediation")
if (!require("broom")) install.packages("broom")
library(tidyverse)
library(mediation)
library(broom)

# Read data
data <- read.table("123.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure column names don't contain special characters
colnames(data) <- make.names(colnames(data))

# Define variable positions
outcome_cols <- 4  # Columns 4-5 are outcomes
cause_col <- 6      # Column 6 is the independent variable
feature_cols <- 53   # Columns 7+ are features/mediators

# Create empty dataframe to store results with CI columns
results <- data.frame(
  Outcome = character(),
  Cause = character(),
  Feature = character(),
  ACME = numeric(),
  ACME_95CI_lower = numeric(),
  ACME_95CI_upper = numeric(),
  ACME_p = numeric(),
  ADE = numeric(),
  ADE_95CI_lower = numeric(),
  ADE_95CI_upper = numeric(),
  ADE_p = numeric(),
  Prop_Mediated = numeric(),
  Total_Effect = numeric(),
  Total_Effect_95CI_lower = numeric(),
  Total_Effect_95CI_upper = numeric(),
  stringsAsFactors = FALSE
)

# Function to extract bootstrap CIs
extract_ci <- function(med_result, effect_type = "d0") {
  if (effect_type == "d0") {
    # ACME
    ci_lower <- med_result$d0.ci[1]
    ci_upper <- med_result$d0.ci[2]
  } else if (effect_type == "z0") {
    # ADE
    ci_lower <- med_result$z0.ci[1]
    ci_upper <- med_result$z0.ci[2]
  } else if (effect_type == "tau") {
    # Total Effect
    ci_lower <- med_result$tau.ci[1]
    ci_upper <- med_result$tau.ci[2]
  }
  return(c(ci_lower, ci_upper))
}

# Loop for mediation analysis
for (outcome_col in outcome_cols) {
  outcome_name <- colnames(data)[outcome_col]
  cause_name <- colnames(data)[cause_col]
  
  for (feature_col in feature_cols) {
    feature_name <- colnames(data)[feature_col]
    
    # Prepare data
    df <- data.frame(
      Y = data[, outcome_col],
      X = data[, cause_col],
      M = data[, feature_col]
    ) %>% na.omit()  # Remove missing values
    
    # Skip if insufficient data
    if (nrow(df) < 10) {
      warning(paste("Insufficient data for", outcome_name, "~", cause_name, "~", feature_name))
      next
    }
    
    tryCatch({
      # Build models
      # Step 1: Effect of independent variable on mediator
      model_m <- lm(M ~ X, data = df)
      
      # Step 2: Effect of independent variable and mediator on outcome
      model_y <- lm(Y ~ X + M, data = df)
      
      # Perform mediation analysis with increased bootstrap samples
      med <- mediate(model_m, model_y, 
                     treat = "X", mediator = "M",
                     boot = TRUE, sims = 1000)  # Increased to 1000 for better CI estimation
      
      # Extract confidence intervals
      acme_ci <- extract_ci(med, "d0")
      ade_ci <- extract_ci(med, "z0")
      total_ci <- extract_ci(med, "tau")
      
      # Store results
      temp_result <- data.frame(
        Outcome = outcome_name,
        Cause = cause_name,
        Feature = feature_name,
        ACME = med$d0,
        ACME_95CI_lower = acme_ci[1],
        ACME_95CI_upper = acme_ci[2],
        ACME_p = med$d0.p,
        ADE = med$z0,
        ADE_95CI_lower = ade_ci[1],
        ADE_95CI_upper = ade_ci[2],
        ADE_p = med$z0.p,
        Prop_Mediated = med$n0,
        Total_Effect = med$tau.coef,
        Total_Effect_95CI_lower = total_ci[1],
        Total_Effect_95CI_upper = total_ci[2],
        stringsAsFactors = FALSE
      )
      
      # Combine results
      results <- rbind(results, temp_result)
      
    }, error = function(e) {
      message(paste("Error in analysis for", outcome_name, "~", cause_name, "~", feature_name))
      message(paste("Error message:", e$message))
    })
  }
}

# Multiple testing correction
results$ACME_p_adj <- p.adjust(results$ACME_p, method = "fdr")
results$ADE_p_adj <- p.adjust(results$ADE_p, method = "fdr")

# Filter significant mediation relationships (ACME p < 0.05)
significant_results <- results %>% 
  filter(ACME_p_adj < 0.05) %>% 
  arrange(ACME_p_adj)

# Print summary of significant results
cat("Significant mediation results (FDR corrected p < 0.05):\n")
print(paste("Number of significant mediation effects:", nrow(significant_results)))

if (nrow(significant_results) > 0) {
  cat("\nTop significant mediation effects:\n")
  print(significant_results %>% 
          select(Outcome, Cause, Feature, ACME, ACME_95CI_lower, ACME_95CI_upper, ACME_p_adj) %>% 
          head(10))
}

# Save results
write.csv(results, "mediation_results_all_with_CI.csv", row.names = FALSE)

if (nrow(significant_results) > 0) {
  write.csv(significant_results, "mediation_results_significant_with_CI.csv", row.names = FALSE)
}

# Generate summary statistics
cat("\nSummary Statistics:\n")
cat("Total mediation tests performed:", nrow(results), "\n")
cat("Significant mediation effects (ACME FDR < 0.05):", nrow(significant_results), "\n")

# Save session info for reproducibility
sink("mediation_analysis_session_info.txt")
print(sessionInfo())
sink()