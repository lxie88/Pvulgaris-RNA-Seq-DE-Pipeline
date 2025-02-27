#!/usr/bin/env Rscript

###############################################################################
#                    SLEUTH-BASED DIFFERENTIAL EXPRESSION ANALYSIS            #
# This script uses Sleuth to analyze gene/transcript expression data generated
# by Kallisto. We run both likelihood ratio tests (LRT) and Wald tests at the  
# gene and transcript levels, and optionally assess genotype x treatment       
# interaction.                                                                 
#                                                                             
# Outputs:                                                                     
#   - CSV files of differentially expressed genes/transcripts                 
#   - Volcano plots (Wald test) for visualizing expression differences         
###############################################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LIBRARIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# We suppress messages to keep the console output clean.
suppressMessages({
  library(cowplot)
  library(sleuth)
  library(ggplot2)
  library(dplyr)     # Ensure we have dplyr for data wrangling
})

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PATH SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Define paths to your input (Kallisto outputs) and output (results) folders.
# Adjust these paths to match your environment.
# For example:
#   datapath <- "/home/username/output"
#   resultdir <- "/home/username/results"
# Currently set for HPC environment:
datapath <- "the Kallisto outputs directory"
resultdir <- "the result directory"

# Set the working directory to store output files in `resultdir`.
setwd(resultdir)

# ~~~~~~~~~~~~~~~~~~~~~~~~ SAMPLE METADATA LOADING ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Load a CSV that links sample IDs to conditions/genotypes.
# This CSV can contain columns:
#   - Run: sample ID
#   - condition: experimental condition (e.g., P_control, P_restriction)
#   - Cultivar: genotype (e.g., DOR 364, IAC Imperador)
sample_to_condition <- read.table("./sraRun.csv", header = TRUE, sep = ',')
head(sample_to_condition)

# Rename relevant columns to align with Sleuth's sample-metadata expectations:
sample_to_condition <- dplyr::select(
  sample_to_condition,
  sample    = Run,
  condition = condition,
  genotype  = Cultivar
)
head(sample_to_condition)

# Match sample IDs to folder paths where Kallisto outputs are located.
sample_id     <- dir(file.path(datapath))
kallisto_dirs <- file.path(datapath, sample_id)

# Add the full path to Kallisto outputs in our sample metadata.
sample_to_condition <- dplyr::mutate(sample_to_condition, path = kallisto_dirs)
print(sample_to_condition)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GENE MAPPING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# This table maps transcript IDs to gene IDs and gene descriptions.
# Because Phytozome's biomaRt is tricky, we pull a manually exported file 
# (Pvulgaris_gene.txt). Adjust the code if you're using a different annotation.
gene <- read.csv("../Pvulgaris_gene.txt")
head(gene)

# Rename the columns to what Sleuth expects:
#   - 'target_id': Kallisto output transcript IDs
#   - 'gene_id':   Gene identifiers
#   - 'gene_description': textual description of each gene
t2g <- dplyr::rename(
  gene,
  target_id        = Transcript.Name,
  gene_id          = Gene.Name,
  gene_description = Description
)

# ~~~~~~~~~~~~~~~~~~~~~~~~ GENOTYPE-SPECIFIC TESTS ~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Example: DOR 364 genotype under two phosphorus (P) treatments.

# 1) Filter the sample metadata to only include the genotype DOR 364:
DOR364 <- sample_to_condition[sample_to_condition$genotype == 'DOR 364', ]

# 2) Build a Sleuth object (so_DOR364) and run LRT at gene level:
#   - 'aggregation_column = gene_id' means we aggregate transcript-level 
#     abundances up to the gene level.
so_DOR364 <- sleuth_prep(
  DOR364,
  target_mapping          = t2g,
  aggregation_column      = 'gene_id',
  extra_bootstrap_summary = TRUE,
  read_bootstrap_tpm      = TRUE
)

# 3) Fit two models:
#   - Full model: includes an effect for condition (~ condition)
#   - Reduced model: intercept-only (~ 1)
so_DOR364 <- sleuth_fit(so_DOR364, ~ condition, 'full')
so_DOR364 <- sleuth_fit(so_DOR364, ~ 1,         'reduced')

# 4) Use the likelihood ratio test (LRT) to compare the two models.
so_DOR364 <- sleuth_lrt(so_DOR364, 'reduced', 'full')

# ~~~~~~ LRT results at the gene level (aggregated) ~~~~~~ #
DOR364_lrt_res_agg <- sleuth_results(
  so_DOR364, 
  'reduced:full', 
  'lrt', 
  show_all       = TRUE,
  pval_aggregate = TRUE
)
# Remove duplicates to avoid complications if multiple transcripts map to one gene.
DOR364_lrt_res_agg <- DOR364_lrt_res_agg[!duplicated(DOR364_lrt_res_agg$target_id), ]

cat("Total genes tested:\n")
print(nrow(DOR364_lrt_res_agg))

significant_DOR364_res_agg <- dplyr::filter(DOR364_lrt_res_agg, qval <= 0.05)

cat("DE genes (qval <= 0.05):\n")
print(nrow(significant_DOR364_res_agg))

# Write the results to CSV
write.csv(significant_DOR364_res_agg, "DE_genes_DOR364.csv")

# ~~~~~~ LRT results at the transcript level ~~~~~~ #
DOR364_lrt_res <- sleuth_results(
  so_DOR364, 
  'reduced:full', 
  'lrt', 
  show_all       = TRUE,
  pval_aggregate = FALSE
)

cat("Total transcripts tested:\n")
print(nrow(DOR364_lrt_res))

significant_DOR364_lrt <- dplyr::filter(DOR364_lrt_res, qval <= 0.05)

cat("DE transcripts (qval <= 0.05):\n")
print(nrow(significant_DOR364_lrt))

# Write the transcript-level results to CSV
write.csv(significant_DOR364_lrt, "DE_transcripts_DOR364.csv")

# ~~~~~~~~~~~~~~~~~~~ WALD TEST + VOLCANO PLOT (DOR364) ~~~~~~~~~~~~~~~~~~~~~ #
# Fit a Wald test to get per-condition coefficients (betas).
# Typically used to generate effect sizes/log fold changes for each condition.
models(so_DOR364)  # Check your available models/coefficients
so_DOR364_wt <- sleuth_wt(so_DOR364, 'conditionP_restriction')  # Compare P_restriction vs reference

DOR364_wt_res <- sleuth_results(
  so_DOR364_wt,
  'conditionP_restriction',
  test_type       = 'wt',
  pval_aggregate  = FALSE
)

# Merge the LRT results with beta (b), standard error (se_b), and mean_obs from Wald test
res <- merge(
  DOR364_lrt_res, 
  DOR364_wt_res[, c('target_id','b','se_b','mean_obs')], 
  by   = 'target_id', 
  sort = FALSE
)
res <- na.omit(res)

write.csv(res, "DE_combinedgenes_DOR364.csv")

res_sig <- dplyr::filter(res, qval <= 0.05)
write.csv(res_sig, "sig_DE_combinedgenes_DOR364.csv")

# Create a binary factor to mark significance for the volcano plot
sig <- ifelse(res$qval <= 0.05, "Significant", "Not Significant")

# Plot the volcano
p1 <- ggplot(res, aes(x = b, y = -log10(qval))) +
  geom_point(aes(col = sig), na.rm = FALSE) +
  ggtitle("DOR364 under P-control vs. P-restriction") +
  theme_bw() +
  scale_color_brewer(palette = "Set2") +
  xlab("Wald Beta") +
  ylab("-log10(qvalue)")

ggsave(p1, filename = './DOR364_volcano.png')

###############################################################################
#                   FUNCTION FOR REUSABLE SLEUTH ANALYSIS                     #
###############################################################################
# We often need to repeat the same steps for different subsets (e.g. genotypes,
# conditions). This function encapsulates the entire LRT+Wald pipeline. 
#
# Inputs:
#   data           : Data frame containing sample metadata (sample, condition, path)
#   var            : A label for naming output files (e.g. "DOR364", "P_control")
#   wt_coeffient   : The coefficient name for the Wald test (e.g., "conditionP_restriction")
#
# Outputs:
#   - CSV results for LRT gene-level, LRT transcript-level, combined LRT+Wald 
#   - Volcano plot of Wald test results 
###############################################################################

DE_sleuth <- function(data, var, wt_coeffient) {
  
  # 1) Prepare the Sleuth object
  so <- sleuth_prep(
    data, 
    target_mapping          = t2g,
    aggregation_column      = 'gene_id',
    extra_bootstrap_summary = TRUE,
    read_bootstrap_tpm      = TRUE
  )
  
  # 2) Fit full and reduced models
  so <- sleuth_fit(so, ~ condition, 'full')
  so <- sleuth_fit(so, ~ 1,         'reduced')
  
  # 3) Perform likelihood ratio test (LRT)
  so <- sleuth_lrt(so, 'reduced', 'full')
  
  # ~~~~~~ LRT: Gene level (aggregated) ~~~~~~ #
  lrt_res_agg <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE, pval_aggregate = TRUE)
  lrt_res_agg <- lrt_res_agg[!duplicated(lrt_res_agg$target_id), ]  # Remove duplicates
  
  cat("Total aggregated genes tested:\n")
  print(nrow(lrt_res_agg))
  
  significant_res_agg <- dplyr::filter(lrt_res_agg, qval <= 0.05)
  cat("Number of DE genes at aggregated level:\n")
  print(nrow(significant_res_agg))
  
  # Save the gene-level LRT results
  write.csv(significant_res_agg, paste0(var, "DE_genes.csv"))
  
  # ~~~~~~ LRT: Transcript level ~~~~~~ #
  lrt_res <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE, pval_aggregate = FALSE)
  
  cat("Total transcripts tested:\n")
  print(nrow(lrt_res))
  
  significant_lrt <- dplyr::filter(lrt_res, qval <= 0.05)
  cat("Number of DE transcripts:\n")
  print(nrow(significant_lrt))
  
  # Save the transcript-level LRT results
  write.csv(significant_lrt, paste0(var, "DE_transcript.csv"))
  
  # 4) Wald test for effect sizes (betas)
  models(so)  # Check available coefficients
  wt <- sleuth_wt(so, wt_coeffient)
  
  wt_res <- sleuth_results(wt, wt_coeffient, test_type = 'wt', pval_aggregate = FALSE)
  
  # Merge LRT table with the Wald test beta columns
  res <- merge(lrt_res, wt_res[, c('target_id','b','se_b','mean_obs')], by = 'target_id', sort = FALSE)
  res <- na.omit(res)
  
  write.csv(res, paste0(var, "DE_combinedgenes.csv"))
  
  res_sig <- dplyr::filter(res, qval <= 0.05)
  write.csv(res_sig, paste0(var, "sig_DE_combinedgenes.csv"))
  
  # Create significance indicator and plot
  sig <- ifelse(res$qval <= 0.05, "Significant", "Not Significant")
  
  p <- ggplot(res, aes(x = b, y = -log10(qval))) +
    geom_point(aes(col = sig), na.rm = FALSE) +
    theme_bw() +
    scale_color_brewer(palette = "Set2") +
    xlab("Wald Beta") +
    ylab("-log10(qvalue)")
  
  # Export the volcano plot
  ggsave(p, filename = paste0(var, "volcano.png"))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXAMPLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 1) Test for genotype "IAC Imperador" under two P treatments
data <- sample_to_condition[sample_to_condition$genotype=='IAC Imperador',]
DE_sleuth(data, var = 'IAC_Imperador_', wt_coeffient = 'conditionP_restriction')

# 2) Test for condition P_control across both genotypes
#    We remove the 'condition' column and rename genotype as 'condition' 
#    for a new comparison (basically flipping the factor).
data <- sample_to_condition[sample_to_condition$condition == 'P_control', ][, -2]
colnames(data)[2] <- 'condition'
DE_sleuth(data, 'P_control_', 'conditionIAC Imperador')

# 3) Test for condition P_restriction across both genotypes
data <- sample_to_condition[sample_to_condition$condition == 'P_restriction', ][, -2]
colnames(data)[2] <- 'condition'
DE_sleuth(data, 'P_restriction_', 'conditionIAC Imperador')

# ~~~~~~~~~~~~~~~~~~~~~~~ INTERACTION EFFECT TEST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Test for interaction of genotype x condition using a single model:
# Full model:   ~ genotype + condition + genotype:condition
# Reduced model:~ genotype + condition
so <- sleuth_prep(
  sample_to_condition,
  target_mapping          = t2g,
  aggregation_column      = 'gene_id',
  extra_bootstrap_summary = TRUE,
  read_bootstrap_tpm      = TRUE
)

so <- sleuth_fit(so, ~ genotype + condition + genotype:condition, 'full')
so <- sleuth_fit(so, ~ genotype + condition,                      'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

# LRT: Gene-level (aggregated)
lrt_res_agg <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE, pval_aggregate = TRUE)
cat("Total genes tested (interaction):\n")
print(nrow(lrt_res_agg))

significant_res_agg <- dplyr::filter(lrt_res_agg, qval <= 0.05)
cat("DE genes (interaction):\n")
print(nrow(significant_res_agg))
write.csv(significant_res_agg, "DE_interaction_genes.csv")

# LRT: Transcript-level
lrt_res <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE, pval_aggregate = FALSE)
cat("Total transcripts tested (interaction):\n")
print(nrow(lrt_res))

significant_lrt <- dplyr::filter(lrt_res, qval <= 0.05)
cat("DE transcripts (interaction):\n")
print(nrow(significant_lrt))
write.csv(significant_lrt, 'interact_DE_transcript.csv')

# ~~~~~~~~~~~~~~~~~~~~~ WALD TEST FOR INTERACTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
models(so)
wt <- sleuth_wt(so, 'genotypeIAC Imperador:conditionP_restriction')

wt_res <- sleuth_results(wt, 'genotypeIAC Imperador:conditionP_restriction', test_type = 'wt', pval_aggregate = FALSE)

# Combine LRT results with Wald betas
res <- merge(lrt_res, wt_res[, c('target_id','b','se_b','mean_obs')], by='target_id', sort=FALSE)
res <- na.omit(res)
write.csv(res, 'DE_intercombinedgenes.csv')

res_sig <- dplyr::filter(res, qval <= 0.05)
write.csv(res_sig, "sig_DE_intercombinedgenes.csv")

# Volcano plot for interaction effect
sig <- ifelse(res$qval <= 0.05, "Significant", "Not Significant")

p <- ggplot(res, aes(x = b, y = -log10(qval))) +
  geom_point(aes(col = sig), na.rm = FALSE) +
  theme_bw() +
  scale_color_brewer(palette = "Set2") +
  xlab("Wald Beta") +
  ylab("-log10(qvalue)")

ggsave(p, filename = 'inter_volcano.png')
