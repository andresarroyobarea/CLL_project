################################## CLL PROJECT #################################
library(dplyr)
library(DESeq2)
library(survival)
library(survminer)
library(cutpointr)
library(biomaRt)
library(stringr)
library(edgeR)
library(ggplot2)
library(Hmisc)
library(genefilter)
library(writexl)
library(psych)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# Functions
surv_plot <- function(fit, max_time, title){
  
  survival_plot <- ggsurvplot(
    fit,                       
    pval = TRUE,             
    conf.int = F,        
    xlim = c(0, max_time),       
    break.time.by = 500,
    ggtheme = theme_bw(),
    legend = "right",
    legend.labs = c("High", "Low"),
    title = title,
    legend.title = "Level",
    lynetype = "solid",
    risk.table = F) 
  
  survival_plot
}

add_ge_metadat <- function(gene_name){
  
  gene_ge <- cll_data_norm %>% filter(grepl(gene_name, rownames(cll_data_norm))) %>% as.vector() %>% unlist() %>% as.data.frame()
  gene_ge$icgc_donor_id <- rownames(gene_ge)
  rownames(gene_ge) <- NULL
  gene_ge <- gene_ge[, c(2, 1)]
  
  hcgc_gene_ge <- getBM(filters = "ensembl_gene_id", attributes= "hgnc_symbol", values = gene_name, mart =  mart)
  
  colnames(gene_ge)[2] <- paste0(hcgc_gene_ge, "_exp")
  cll_donor_data <- merge(cll_donor_data, gene_ge, by = "icgc_donor_id")
  
  cll_donor_data
}

preproc_surv <- function(data, gene_var){
  
  # Filter patients with NA values in SG
  SG_cll_donor_data <- data %>% filter(!is.na(survival_time_last_followup_SG))
  
  # Filter patients with NA values in SLP
  SLP_cll_donor_data <- data %>% filter(!is.na(relapse_interval_SLP))
  
  # Extract gene expression values.
  gene_exp_SLP <- SLP_cll_donor_data %>% dplyr::select(all_of(gene_var)) %>% pull()
  gene_exp_SG <- SG_cll_donor_data %>% dplyr::select(all_of(gene_var)) %>% pull()
  
  
  # Metrics calculation
  # ROC for SLP and SG
  cp_roc_gene_SLP <- cutpointr(SLP_cll_donor_data, !!sym(gene_var), relapse_type, direction = ">=", pos_class = 1,
                               neg_class = 0, method = maximize_metric, metric = sum_sens_spec, na.rm = T)
  
  cp_roc_gene_SG <- cutpointr(SG_cll_donor_data, !!sym(gene_var), vital_status, direction = ">=", pos_class = 1,
                              neg_class = 0, method = maximize_metric, metric = sum_sens_spec, na.rm = T)
  
  # Maxstat approach for SLP and SG
  cp_maxst_gene_SLP <- surv_cutpoint(SLP_cll_donor_data,
                                     time = "relapse_interval_SLP",
                                     event = "relapse_type",
                                     variables = gene_var,
                                     minprop = 0.1)
  
  
  
  cp_maxst_gene_SG <- surv_cutpoint(SG_cll_donor_data,
                                    time = "survival_time_last_followup_SG",
                                    event = "vital_status",
                                    variables = gene_var,
                                    minprop = 0.1)
  
  help("surv_cutpoint")
  # Gene name for new variables
  gene_name <- gsub("_.*", "", gene_var)
  
  # Quantiles, tertiles and quintiles
  quantiles <-unname(cll_donor_data %>% dplyr::select(all_of(gene_var)) %>% as.matrix() %>% quantile())
  tertiles <- unname(cll_donor_data %>% dplyr::select(all_of(gene_var)) %>% as.matrix() %>% quantile(., probs = seq(0, 1, 1/3)))
  quintiles <- unname(cll_donor_data %>% dplyr::select(all_of(gene_var)) %>% as.matrix() %>% quantile(., probs = seq(0, 1, 1/5)))
  
  
  # Cimplete table with media, mean, quantile, tertil and quartil divisions.
  cll_donor_data <- cll_donor_data %>% mutate(!!paste0(gene_name, "_level_ROC_SLP") := ifelse(!!sym(gene_var)  > cp_roc_gene_SLP$optimal_cutpoint, "high", "low"),
                                              
                                              !!paste0(gene_name, "_level_ROC_SG") := ifelse(!!sym(gene_var)  > cp_roc_gene_SG$optimal_cutpoint, "high", "low"),
                                              
                                              !!paste0(gene_name, "_level_maxstat_SLP") := ifelse(!!sym(gene_var)  > cp_maxst_gene_SLP$cutpoint$cutpoint, "high", "low"),
                                              
                                              !!paste0(gene_name, "_level_maxstat_SG") := ifelse(!!sym(gene_var)  > cp_maxst_gene_SG$cutpoint$cutpoint, "high", "low"), 
                                              
                                              !!paste0(gene_name, "_level_median") := ifelse(!!sym(gene_var)  > median(!!sym(gene_var) ) , "high", "low"),
                                              
                                              !!paste0(gene_name, "_level_mean") := ifelse(!!sym(gene_var)  > mean(!!sym(gene_var) ), "high", "low"),
                                              
                                              !!paste0(gene_name, "_level_quantiles") := ifelse(!!sym(gene_var)  > quantiles[4], "high" , ifelse(!!sym(gene_var) < quantiles[2], "low", NA)),
                                              
                                              !!paste0(gene_name, "_level_tertiles") := ifelse(!!sym(gene_var) > tertiles[3], "high" , ifelse(!!sym(gene_var)  < tertiles[2], "low", NA)),
                                              
                                              !!paste0(gene_name, "_level_quintiles") := ifelse(!!sym(gene_var) > quintiles[5], "high", ifelse(!!sym(gene_var)  <quintiles[2], "low", NA))
                                              
  )
  
  cll_donor_data
}


f.shapiro.stat <- function(x, n_diff_numbers = 3) {
  res <- ifelse(sum(!is.na(unique(x))) < n_diff_numbers, NA, as.numeric(shapiro.test(x)$p.value))
  return(res)
}



# 0. IMPORT RNA-SEQ DATA --------------------------------------------------------------------------------------------------------------------
cll_data <- read.table("./data/exp_seq_CLLE_ES.tsv", sep = "\t", header = T)

# 1. PREPROCESSING DATA TO RAW COUNTS MATRIX ------------------------------------------------------------------------------------------------
# Select only columns of interest (donor, specimen, gene ensembl id y raw read counts)
cll_data <- cll_data %>% dplyr::select(icgc_donor_id, icgc_specimen_id, gene_id, raw_read_count)

# 3 donor were evaluated 3 times (3 specimens - 3 samples from the same patient). We supposed that there were technical
# replicates. We kept only one for each patient (there with the lowest icgc_specimen_id).

# Specimens to remove
spec_remove <- c("SP16419", "SP16425", "SP16401", "SP16407", "SP16389", "SP16377")

cll_data <- subset(cll_data, !(cll_data$icgc_specimen_id %in% spec_remove))

rm(spec_remove)

# Saving donors and specimen info for all patients that had gene expression data (RNA-seq).
cll_rna_seq_donors_spec <- cll_data %>% dplyr::select(icgc_donor_id, icgc_specimen_id) %>% unique()
rownames(cll_rna_seq_donors_spec) <- NULL

# Gene expression was evaluated using ensembl id and ensembl id version genes so there are
# genes that were measured two times (at gene and transcript level). For further purposes, it is
# needed to retain all genes that were measured in common between all patients. In this case, 
# genes with id version was shared between all and they were retained.
cll_data_id_version <- subset(cll_data, grepl("\\.", cll_data$gene_id))

# Select the ensembl id without version to checking.
cll_data_id <- subset(cll_data, !grepl("\\.", cll_data$gene_id))

# Remove specfimen info.
cll_data_id <- cll_data_id %>% dplyr::select(icgc_donor_id, gene_id, raw_read_count)

cll_data_id_version <- cll_data_id_version %>% dplyr::select(icgc_donor_id, gene_id, raw_read_count)

# Split cll_data based on donor id to extract 304 dataframes.
cll_data_id_list <- split(cll_data_id, f = cll_data_id$icgc_donor_id)

cll_data_id_version_list <- split(cll_data_id_version, f = cll_data_id_version$icgc_donor_id)

# The 10 donors with not ensembl id version information.
ten_only_ensembl_id <- setdiff(names(cll_data_id_list), names(cll_data_id_version_list))
cll_data_id_list <- cll_data_id_list[names(cll_data_id_list) %in% ten_only_ensembl_id]

# For now, the 10 donors were excluded for further analysis and we maintain the 294 patients
# that shared commmon ensembl id version.
rm(cll_data_id, cll_data_id_list)
# Filtering the donors_specimen table to retain only the 294 patients
cll_rna_seq_donors_spec_id_version <- cll_data %>% dplyr::select(icgc_donor_id, icgc_specimen_id) %>% unique() %>% filter(!(icgc_donor_id %in% ten_only_ensembl_id))

# The raw_read_counts column name is changed to donor_id for each dataset.
for (i in seq_along(cll_data_id_version_list)){
  colnames(cll_data_id_version_list[[i]])[3] <- cll_data_id_version_list[[i]]$icgc_donor_id[1]
}

# Donor DO51973 has no info about the first gene in all donors (ENSG00000103066.8). To avoid other possible conflictive genes in this donor
# common genes shared between him/her and other random donor were obtained and they were considered the total number of genes (57819) 
common_genes <- intersect(cll_data_id_version_list$DO223366$gene_id, cll_data_id_version_list$DO51973$gene_id)

# All dataframes in the list were filtered using common genes.
cll_data_id_version_list <- lapply(cll_data_id_version_list, function(x){x %>% filter(x$gene_id %in% common_genes)})

# Donor column is not useful, so it was removed for all datasets.
cll_data_id_version_list <- lapply(cll_data_id_version_list, function(x){x <- x %>% dplyr::select(2, 3)})

# All data frames were merged by row.
cll_data_id_version <- do.call("cbind", cll_data_id_version_list)

# Ensembl id were fixed as rownames.
rownames(cll_data_id_version) <- cll_data_id_version$DO223355.gene_id

# ¿Are the ensembl id ordered in the same way for all donors? Like we found duplicated columns, we can affirm it.
table(duplicated(as.list(cll_data_id_version)) == TRUE)

# Duplicated ensembl id are in odd columns and ther were depleted.
cll_data_id_version <- cll_data_id_version[, -seq(1, 588, 2)]

# Properly colnames were fixed.
colnames(cll_data_id_version) <- gsub("\\..*", "", colnames(cll_data_id_version))

# Check if data was processed in a suitable way checking 3 randomly selected donors
table(cll_data %>% filter(icgc_donor_id == "DO52017", gene_id != "ENSG00000103066.8") %>% dplyr::select(raw_read_count) == cll_data_id_version$DO52017)
table(cll_data %>% filter(icgc_donor_id == "DO223355", gene_id != "ENSG00000103066.8") %>% dplyr::select(raw_read_count) == cll_data_id_version$DO223355)
table(cll_data %>% filter(icgc_donor_id == "DO52014", gene_id != "ENSG00000103066.8") %>% dplyr::select(raw_read_count) == cll_data_id_version$DO52014)

# There are not missing values in our data.
table(is.na(cll_data_id_version))

# Donors were saved for further analysis.
cll_donors <- colnames(cll_data_id_version)

# Specimens were saved 
cll_specimens_rna_seq <- cll_rna_seq_donors_spec_id_version$icgc_specimen_id

# Raw RNA-seq count matrix was saved.
write.table(cll_data_id_version, "./data/rna_seq_count_matrix_raw.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

rm(cll_data)

# 2. PREPROCESSING METADATA ---------------------------------------------------------------------------------------------------------

## 2.1 Donor data
cll_donor_data <- read.table("./metadata/clinical_raw_metadata/donor.tsv", sep = "\t", header = T)

# Last followup time is equal to survival time for this data.
table(cll_donor_data$donor_interval_of_last_followup == cll_donor_data$donor_survival_time)

# Select relevant information.
cll_donor_data <- cll_donor_data %>% dplyr::select(icgc_donor_id, donor_sex, donor_vital_status, disease_status_last_followup, donor_relapse_type, donor_relapse_interval, donor_tumour_stage_at_diagnosis,
                                            donor_survival_time)

# Rename columns
colnames(cll_donor_data)[2:8] <- gsub(pattern = "donor_", replacement = "", x = colnames(cll_donor_data)[2:8])
colnames(cll_donor_data)[8] <- "survival_time_last_followup_SG"

# Add specimen RNA-seq info.
cll_donor_data <- merge(cll_donor_data, cll_rna_seq_donors_spec, by = "icgc_donor_id")
cll_donor_data <- cll_donor_data %>% relocate(icgc_specimen_id, .after = icgc_donor_id)
colnames(cll_donor_data[2]) <- "icgc_rna_seq_specimen_id"

# Recode vital status: decesased (1) and alive (0)
cll_donor_data$vital_status <- as.numeric(ifelse(cll_donor_data$vital_status == "deceased", 1, 0))

# Recode sex : male (1) and female (0)
cll_donor_data$sex <- as.factor(ifelse(cll_donor_data$sex == "male", 1, 0))

# Unknown class in disease status
cll_donor_data$disease_status_last_followup[cll_donor_data$disease_status_last_followup == ""] <- "unknown"

# Unknown class in tumour stage at diagnosis.
cll_donor_data$tumour_stage_at_diagnosis[cll_donor_data$tumour_stage_at_diagnosis == ""] <- "unknown"

# Unknown class in "disease_status_last_followup".
cll_donor_data$disease_status_last_followup[cll_donor_data$disease_status_last_followup == ""] <- "unknown"

# There was a "relapse" individual that was merged with "progression" people.
cll_donor_data$disease_status_last_followup[cll_donor_data$disease_status_last_followup == "relapse"] <- "progression"

# Recode factor levels names.
cll_donor_data$disease_status_last_followup <- recode_factor(cll_donor_data$disease_status_last_followup,
                                                             `complete remission` = "complete_remission", 
                                                             `partial remission` = "partial_remission")

# Change relapse type codification. We suppose that individuals categorized as "progression" suffered a relapse while 
# individuals with no information, did not suffere any relaps episode. Relapse was codified as 1 and non-relapse as 0.
cll_donor_data$relapse_type <- as.numeric(ifelse(cll_donor_data$relapse_type == "progression (liquid tumours)", 1, 0))
colnames(cll_donor_data)[7] <- "relapse_interval_SLP"

# Add SLP
cll_donor_data <- cll_donor_data %>% mutate(relapse_interval_SLP = ifelse(relapse_type == 0, survival_time_last_followup_SG, relapse_interval_SLP))

# Library size
library_size <- as.data.frame(colSums(cll_data_id_version))
library_size$icgc_donor_id <- rownames(library_size)
rownames(library_size) <- NULL
library_size <- library_size[, c(2, 1)]
colnames(library_size)[2] <- "library_size"

# Filter donor data with samples used in RNA-seq.
cll_donor_data <- cll_donor_data[cll_donor_data$icgc_donor_id %in% colnames(cll_data_id_version), ]

# Add library size in donor data.
cll_donor_data <- merge(cll_donor_data, library_size, by = "icgc_donor_id")

# Factorise metadata columns.
cll_donor_data$disease_status_last_followup <- as.factor(cll_donor_data$disease_status_last_followup)
cll_donor_data$tumour_stage_at_diagnosis <- as.factor(cll_donor_data$tumour_stage_at_diagnosis)

# Save metadata
write.csv(cll_donor_data, "./metadata/clinical_prepr_metadata/cll_metadata.csv", sep = ",", col.names = TRUE,
          row.names = F)

# 3. NORMALIZATION WITH DESEQ2 -----------------------------------------------------------------------------------------------
dds_cll_data <- DESeqDataSetFromMatrix(countData = cll_data_id_version,
                                      colData = cll_donor_data,
                                      design = ~ disease_status_last_followup)


# Estimate size factors using dds object.
dds_cll_data <- estimateSizeFactors(dds_cll_data)

# Extract normalized count matrix. This matrix ONLY WILL BE USED TO EXTRACT INFO FOR SURVIVAL ANALYSIS. FOR DGE STUDIES, A PREVIOUS FILTERING IS RECOMMENDED.
cll_data_norm <- as.data.frame(counts(dds_cll_data, normalized = TRUE))

# Save RNA-seq count matrix normalized.
write.table(cll_data_norm, "./data/rna_seq_count_matrix_norm.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)


# 4. SURVIVAL ANALYSIS --------------------------------------------------------------------------------------------------------
# 4.1 Extract specific gene expression and add to metadata. ENSG00000136997.10 is the MYC ensembl identifier.
cll_donor_data <- cll_donor_data[, -11]
cll_donor_data <- add_ge_metadat("ENSG00000143384")
colnames(cll_donor_data)

# 4.2 Fix cutfoffs to dicide gene expression into high and low using the metrics used below.
cll_donor_data_surv <- cll_donor_data[, -c(11:20)]
cll_donor_data_surv <- preproc_surv(cll_donor_data, "MCL1_exp")

# Recode SLP and SG as numeric to plot results.
cll_donor_data_surv$vital_status <- as.numeric(cll_donor_data_surv$vital_status)
cll_donor_data_surv$vital_status <- ifelse(cll_donor_data_surv$vital_status == 2, 1, 0)
cll_donor_data_surv$relapse_type <- as.numeric(cll_donor_data_surv$relapse_type)
cll_donor_data_surv$relapse_type <- ifelse(cll_donor_data_surv$relapse_type == 2, 1, 0)

colnames(cll_donor_data_surv)

# Storing the variables used to divide gene expression into high and low for SLP
gene_level_slp <- colnames(cll_donor_data_surv)[c(12, 14, 16:20)]
gene_level_slp


# Storing the variables used to divide gene expression into high and low for SG.
gene_level_sg <- colnames(cll_donor_data_surv)[c(13, 15, 16:20)]
gene_level_sg


# Bucle to generate several plots in once time and save them into a file.
for(i in 1:length(gene_level_slp)){
  
  fit <- survfit(as.formula(paste0("Surv(relapse_interval_SLP, relapse_type) ~", gene_level_slp[i])), data = cll_donor_data_surv)
  
  png(paste0("SLP_", str_split(gene_level_slp[i], "_")[[1]][3], ".png"), width = 700, height = 600)
  
  p <- ggsurvplot(fit,
                  pval = TRUE, conf.int = FALSE, 
                  ggtheme = theme_bw(), 
                  palette = c("#FF0027", "#060606"),
                  xlim = c(0,10000),
                  break.x.by = 1000,
                  title = paste0("SLP_", str_split(gene_level_slp[i], "_")[[1]][3]),
                  legend = "right",
                  legend.labs = c("High", "Low"),
                  lynetype = "solid",
                  xlab = "Time (days)",
                  legend.title = "Gene expression level", 
                  pval.method = T)
  
  p <- p$plot+theme(plot.title = element_text(hjust = 0.5))
  
  print(
   p
  )
  dev.off()
}




# 5. GENE-GENE CORRELATION ANALYSIS --------------------------------------------------------------------------------------------------------

# This analysis was developed to evaluate possible correlations between the 6 genes and the other genes in the RNA-seq data.

# Genes of interest.
# MYC: ENSG00000136997
# HNRNPK: ENSG00000165119
# MCL1: ENSG00000143384
# NCL: ENSG00000115053
# PIEZO1: ENSG00000103335
# SAMHD1: ENSG00000101347

our_genes <- c("ENSG00000136997", "ENSG00000165119", "ENSG00000143384", "ENSG00000115053", "ENSG00000103335", "ENSG00000101347")

# Counts as column to evaluate correlation.
counts <- as.data.frame(t(as.matrix(cll_data_norm)))

# Ensembl ID code.
colnames(counts) <- gsub("\\..*", "", colnames(counts))

# Next piece of code was used to extract the correlation tables for each gene.

# Correlation analysis between a gene and the rest.
corr_an <- corr.test(counts[, "ENSG00000101347"], counts, method = "spearman", adjust = "fdr")

# Correlation values extraction.
corr <- as.data.frame(t(as.matrix(corr_an$r)))
corr$Gene_1 <- "SAMHD1"
colnames(corr)[1] <- "Correlation_Spearman"

# P-values extraction for each pairs of correlated genes.  
p_v <- as.data.frame(t(as.matrix(corr_an$p)))
colnames(p_v) <- "P_value"
  
# P-adj extraction for each pairs of correlated genes.  
p_adj <- as.data.frame(t(as.matrix(corr_an$p.adj)))
colnames(p_adj)[1] <- "P_adj"
p_adj$Gene_2 <- rownames(p_adj)

# Merging results
result <- merge(corr, p_v, by = 0)
colnames(result)[1] <- "Gene_2"
result <- merge(result, p_adj, by = "Gene_2")
result <- result[, c(3, 1, 2, 4, 5)]

# Selecting correlation pairs without NA values, with absolute R upper 0.2 and p-adjuted upper 0.05.
result <- result %>% filter(!(is.na(Correlation_Spearman)) & P_adj < 0.05 & abs(Correlation_Spearman) > 0.2 & Gene_2 != "ENSG00000101347") 

# Round correlation values
result$Correlation_Spearman <- round(result$Correlation_Spearman, 2)

# Substitution of Ensembl ID codes for Gene symbols.
symbols <- getBM(filters = "ensembl_gene_id", attributes = c("hgnc_symbol", "ensembl_gene_id"), values = result$Gene_2, mart =  mart)
colnames(symbols)[2] <- "Gene_2"                  
result <- merge(result, symbols, by = "Gene_2", all = TRUE)
result$hgnc_symbol[result$hgnc_symbol == ""] <- NA
result$hgnc_symbol <- ifelse(is.na(result$hgnc_symbol), result$Gene_2, result$hgnc_symbol)
result <- result %>% select(!Gene_2) %>% relocate(hgnc_symbol, .after = Gene_1) %>% rename(Gene_2 = hgnc_symbol)
  
# Order data by decreasing absolute R value.
result <- result %>% arrange(desc(abs(Correlation_Spearman)))

# Save results in an Excel file.
myc_results <- result
HNRNPK_results <- result
MCL1_results <- result
NCL_results <- result
PIEZO1_results <- result
SAMHD1_results <- result
  
rm(symbols, result, p_v, p_adj, corr_an, corr)

genes_corr_res <- list(MYC = myc_results, 
                       HNRNPK = HNRNPK_results,
                       MCL1 = MCL1_results,
                       NCL = NCL_results,
                       PIEZO1 = PIEZO1_results,
                       SAMHD1 = SAMHD1_results)

write_xlsx(genes_corr_res, "./results/CLL_genes_correlation_vs_all.xlsx")


