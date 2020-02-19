#!/usr/bin/env Rscript
suppressMessages({
    library('cowplot')
    library('sleuth')
    library(ggplot2)
})
########STEP 1 Prepare all data that is going to be used ######################
# set input and output dirs
# datapath is your kallisto output 
# datapath <- "/home/lx83659/Desktop/BINF8940/output"  
# datapath <- "C:/Users/meng/Desktop/BINF8940/output"
datapath <- "/work/gene8940/lx83659/FinalProject/output"
# remember to make results dir in the bash
# resultdir <- "/home/lx83659/Desktop/BINF8940/results"   
# resultdir <-  "C:/Users/meng/Desktop/BINF8940/results"
resultdir <- "/work/gene8940/lx83659/FinalProject/results"
setwd(resultdir)

#need to revise SRArun to combine sample information with kallisto output
sample_to_condition <- read.table("../sraRun.csv", header = TRUE,
                                  sep = ',')
head(sample_to_condition)
sample_to_condition <- dplyr::select(sample_to_condition,sample = Run, condition = condition, genotype = Cultivar)
head(sample_to_condition)
sample_id <- dir(file.path(datapath))
kallisto_dirs <- file.path(datapath, sample_id)
sample_to_condition <- dplyr::mutate(sample_to_condition, path = kallisto_dirs)
print(sample_to_condition) 
##############STEP2 Relate transcripts to geness #####################
# get gene ID and gene discription from phytozome. The biomaRt is offical designedfor ensemble database.The bioMart tool in phytozome has a bug which make it 
# difficult to write code for extracting those data. The port code should be 443, not 80 in their registry text.
# I have tried several times and failed. so I queried directly from website. 
# Personally, I think Phytozome is not user-friendly............
gene <- read.csv("../Pvulgaris_gene.txt")
head(gene)
t2g <- dplyr::rename(gene, target_id = Transcript.Name,
                     gene_id=Gene.Name, gene_description=Description) 

##############STEP3 Create slueth object and likelihood test for customerized models at gene level and transcripts level seperately#####################
# select DOR364 genotype under two p treatments
DOR364 <- sample_to_condition[sample_to_condition$genotype=='DOR 364',]
# prepare sleuth model 
# data is combined with the gene id and gene description together. 
so_DOR364<- sleuth_prep(DOR364, target_mapping = t2g,aggregation_column = 'gene_id',extra_bootstrap_summary = TRUE,read_bootstrap_tpm=TRUE)
# fit a full model that includes a parameter for condition 
so_DOR364 <- sleuth_fit(so_DOR364, ~condition, 'full')
# fit a reduced model that only includes an intercept
so_DOR364 <- sleuth_fit(so_DOR364, ~1, 'reduced')
so_DOR364 <- sleuth_lrt(so_DOR364, 'reduced', 'full')

# likelihood ratio test for models: basically you fit a full model with parameters and then you fit a reduced 
# model with parameter 

# (1) likelihood ratio test for gene level with setting pval_aggregate = T
DOR364_lrt_res_agg <- sleuth_results(so_DOR364, 'reduced:full', 'lrt', show_all = T,pval_aggregate = T)
# The ‘num_aggregated_transcripts’ column lists the number of transcripts used to make the gene determination. 
# ‘pval’ displays the p-value for the gene. 
# ‘qval’ displays the Benjamini-Hochberg-adjusted false discovery rate for the gene.
# remove duplicate genes 
DOR364_lrt_res_agg <- DOR364_lrt_res_agg[!duplicated(DOR364_lrt_res_agg$target_id), ]
print ("how many genes in total? ")
print(nrow(DOR364_lrt_res_agg))
significant_DOR364_res_agg<- dplyr::filter(DOR364_lrt_res_agg, qval <= 0.05)
print ("how many differentially expressed genes:")
print(nrow(significant_DOR364_res_agg))
significant_DOR364_res_agg[order(significant_DOR364_res_agg$qval,decreasing = F),]
# save DE genes as csv file 
write.csv(significant_DOR364_res_agg, "DE_genes_DOR364.csv")


# (2) likelihood ratio test for transcript level with setting pval_aggregate = F
DOR364_lrt_res <- sleuth_results(so_DOR364, 'reduced:full', 'lrt', show_all = T,pval_aggregate = FALSE)
head(DOR364_lrt_res)
print ("how many transcripts? ")
print(nrow(DOR364_lrt_res))
# significant diff expression transcripts 
significant_DOR364_lrt<- dplyr::filter(DOR364_lrt_res, qval <= 0.05)
head(significant_DOR364_lrt)
print ("how many differentially expressed transcripts:")
print(nrow(significant_DOR364_lrt))
write.csv(significant_DOR364_lrt, "DE_transcripts_DOR364.csv")
##############STEP4 Create slueth object and wald test for customerized models at gene level and transcripts level seperately#####################
# wald test:
# A Wald test is fitting the full model and then using a coefficient
#(or more) and its (or their) standard error to perform a test for difference from 0. 
# Genes with a coefficient (or more) significantly different from 0 are differentially expressed and their log2 fold-change is the coefficient.
# wald test for sleuth model
# Here we use walt test is to get beta value and used for ploting volcano plots for DE 
models(so_DOR364)
# make wt test model to get beta values
DOR364_wt <- sleuth_wt(so_DOR364,'conditionP_restriction')
DOR364_wt_res <- sleuth_results(DOR364_wt,'conditionP_restriction',test_type = 'wt',pval_aggregate = FALSE)
# merge beta vale with lrt results
res <- merge(DOR364_lrt_res, DOR364_wt_res[, c('target_id', 'b','se_b','mean_obs')], on='target_id', sort = FALSE)
res <- na.omit(res)
write.csv(res, "DE_combinedgenes_DOR364.csv")
res_sig <- dplyr::filter(res, qval <= 0.05)
write.csv(res_sig, "sig_DE_combinedgenes_DOR364.csv")
sig <- ifelse(res$qval <= 0.05, "Significant", "Not Significant")
#pdf(paste0(IAC,".pdf"))
p1 <- ggplot(res, aes(b, -log10(qval))) +
  geom_point(aes(col=sig), na.rm = F) +
  ggtitle(paste0("DOR364 under two P treatments volcano plot")) +
  theme_bw() + scale_color_brewer(palette = "Set2") +
  #geom_hline(yintercept = -log10(0.05), color = "#990000", linetype = "dashed", size = 0.1) +
  xlab("Wald Beta") + ylab("-log10(qvalue)");
ggsave(p1, filename = './DOR364_volcano.png' )

#################################################################################################
# since we will repeat the process for different model test, I wrote a function to simplify 
###############################################################################################

DE_sleuth <- function (data,var,wt_coeffient)
{ so<- sleuth_prep(data, target_mapping = t2g,aggregation_column = 'gene_id',extra_bootstrap_summary = TRUE,read_bootstrap_tpm=TRUE)
  so <- sleuth_fit(so, ~condition, 'full')
  so <- sleuth_fit(so, ~1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  # (1) likelihood ratio test for gene level 
  lrt_res_agg <- sleuth_results(so, 'reduced:full', 'lrt', show_all = T,pval_aggregate = T)
  lrt_res_agg <-lrt_res_agg[!duplicated(lrt_res_agg$target_id), ]
  print ("how many genes in total? ")
  print(nrow(lrt_res_agg))
  significant_res_agg<- dplyr::filter(lrt_res_agg, qval <= 0.05)
  print ("how many differentially expressed genes:")
  print(nrow(significant_res_agg))
  significant_res_agg[order(significant_res_agg$qval,decreasing = F),]
  # save DE genes as csv file 
  write.csv(significant_res_agg, paste0(var,"DE_genes.csv"))
  # (2) likelihood ratio test for transcript level with setting pval_aggregate = F
  lrt_res <- sleuth_results(so, 'reduced:full', 'lrt', show_all = T,pval_aggregate = FALSE)
  head(lrt_res)
  print ("how many transcripts? ")
  print(nrow(lrt_res))
  # significant diff expression transcripts 
  significant_lrt<- dplyr::filter(lrt_res, qval <= 0.05)
  head(significant_DOR364_lrt)
  print ("how many differentially expressed transcripts:")
  print(nrow(significant_lrt))
  # save DE genes as csv file 
  write.csv(significant_lrt, paste0(var,'DE_transcript.csv'))
  models(so)
  # make wt test model to get beta values
  wt <- sleuth_wt(so, paste0(wt_coeffient))
  wt_res <- sleuth_results(wt, paste0(wt_coeffient),test_type = 'wt',pval_aggregate = FALSE)
  # merge beta vale with lrt results
  res <- merge(lrt_res, wt_res[, c('target_id', 'b','se_b','mean_obs')], on='target_id', sort = FALSE)
  res <- na.omit(res)
  write.csv(res, paste0(var,'DE_combinedgenes.csv'))
  res_sig <- dplyr::filter(res, qval <= 0.05)
  write.csv(res_sig, paste0(var,"sig_DE_combinedgenes.csv"))
  sig <- ifelse(res$qval <= 0.05, "Significant", "Not Significant")
  p <- ggplot(res, aes(b, -log10(qval))) +
    geom_point(aes(col=sig), na.rm = F) +
    theme_bw() + scale_color_brewer(palette = "Set2") +
    xlab("Wald Beta") + ylab("-log10(qvalue)")
  ggsave(p, filename = paste0(var, 'volcano.png'))}

# 2. test for genotype IAC Imperor under two P treatments 
data <- sample_to_condition[sample_to_condition$genotype=='IAC Imperador',]
DE_sleuth(data, var = 'IAC Imperador',wt_coeffient = 'conditionP_restriction')


# 3. test for condition P_control for two genotypes 
data <- sample_to_condition[sample_to_condition$condition=='P_control',][, -2]
colnames(data)[2] <- 'condition'
DE_sleuth(data, 'P_control','conditionIAC Imperador')


# 4 test in P restriction condition for two genotypes
data <- sample_to_condition[sample_to_condition$condition=='P_restriction',][, -2]
colnames(data)[2] <- 'condition'
DE_sleuth(data, 'P_restriction','conditionIAC Imperador')


# 5 test for P x genotype interaction effect, need to write new one 
# test for interaction of genotype and control 
so<- sleuth_prep(sample_to_condition, target_mapping = t2g,aggregation_column = 'gene_id',extra_bootstrap_summary = TRUE,read_bootstrap_tpm=TRUE)
#so <- sleuth_prep(, extra_bootstrap_summary = TRUE,read_bootstrap_tpm=TRUE)
so <- sleuth_fit(so, ~genotype+condition+genotype:condition, 'full')
so <- sleuth_fit(so, ~genotype+condition, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
# likelihood test for models 
lrt_res_agg <- sleuth_results(so, 'reduced:full', 'lrt', show_all = T,pval_aggregate = T)
print ("how many genes in total? ")
print(nrow(lrt_res_agg))
significant_res_agg<- dplyr::filter(lrt_res_agg, qval <= 0.05)
print ("how many differentially expressed genes:")
print(nrow(significant_res_agg))
significant_res_agg[order(significant_res_agg$qval,decreasing = F),]
# save DE genes as csv file 
write.csv(significant_res_agg, paste0("DE_interaction_genes.csv"))
# likelihood ratio test for transcript level with setting pval_aggregate = F
lrt_res <- sleuth_results(so, 'reduced:full', 'lrt', show_all = T,pval_aggregate = FALSE)
head(lrt_res)
print ("how many transcripts? ")
print(nrow(lrt_res))
# significant diff expression transcripts 
significant_lrt<- dplyr::filter(lrt_res, qval <= 0.05)
head (significant_DOR364_lrt)
print ("how many differentially expressed transcripts:")
print(nrow(significant_lrt))
# save DE transcripts as csv file 
write.csv(significant_lrt,'interact_DE_transcript.csv')
models(so)
# make wt test model to get beta values
wt <- sleuth_wt(so,'genotypeIAC Imperador:conditionP_restriction')
wt_res_agg <- sleuth_results(wt,'genotypeIAC Imperador:conditionP_restriction',test_type = 'wt')
wt_res <- sleuth_results(wt,'genotypeIAC Imperador:conditionP_restriction',test_type = 'wt',pval_aggregate = FALSE)
res <- merge(lrt_res, wt_res[, c('target_id', 'b','se_b','mean_obs')], on='target_id', sort = FALSE)
res <- na.omit(res)
write.csv(res, 'DE_intercombinedgenes.csv')
res_sig <- dplyr::filter(res, qval <= 0.05)
write.csv(res_sig,"sig_DE_intercombinedgenes.csv")
sig <- ifelse(res$qval <= 0.05, "Significant", "Not Significant")
p <- ggplot(res, aes(b, -log10(qval))) +
  geom_point(aes(col=sig), na.rm = F) +
  theme_bw() + scale_color_brewer(palette = "Set2") +
  xlab("Wald Beta") + ylab("-log10(qvalue)")
ggsave(p, filename = 'inter_volcano.png')

