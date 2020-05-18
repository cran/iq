## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=TRUE, include=FALSE, echo=FALSE-------------------------------------
require("knitr")
local_file_exist <- file.exists("DIA-report-long-format.txt")

## ---- eval = local_file_exist, include = TRUE---------------------------------
library("iq") # if not already installed, run install.packages("iq") 

raw <- read.delim("DIA-Report-long-format.txt")

selected <- raw$F.ExcludedFromQuantification == "False" & 
            !is.na(raw$PG.Qvalue) & (raw$PG.Qvalue < 0.01) &
            !is.na(raw$EG.Qvalue) & (raw$EG.Qvalue < 0.01)

raw <- raw[selected,]

## ----eval = local_file_exist, echo=TRUE, message = TRUE, include = TRUE-------
sample_id  <- "R.FileName" 

secondary_id <- c("EG.Library", "FG.Id", "FG.Charge", "F.FrgIon", "F.Charge", "F.FrgLossType")

norm_data <- iq::preprocess(raw, 
                            sample_id  = sample_id, 
                            secondary_id = secondary_id)

protein_list <- iq::create_protein_list(norm_data)

result <- iq::create_protein_table(protein_list)

## ----eval = local_file_exist--------------------------------------------------
annotation_columns <- c("PG.Genes", "PG.ProteinNames")

extra_names <- iq::extract_annotation(rownames(result$estimate), 
                                      raw, 
                                      annotation_columns = annotation_columns)

write.table(cbind(Protein = rownames(result$estimate),
                  extra_names[, annotation_columns],
                  MaxLFQ_annotation = result$annotation,
                  result$estimate),
            "iq-MaxLFQ.txt", sep = "\t", row.names = FALSE)

## ----eval = local_file_exist--------------------------------------------------
#--------------------- Replacing ---------------------
# protein_list <- iq::create_protein_list(norm_data) #
# result <- iq::create_protein_table(protein_list)   #
#-----------------------------------------------------

result_faster <- iq::fast_MaxLFQ(norm_data)


## ----eval = local_file_exist--------------------------------------------------
cat("Max difference =", max(abs(result_faster$estimate - result$estimate), na.rm = TRUE), "\n")

cat("Identical NAs =", identical(is.na(result_faster$estimate), is.na(result$estimate)), "\n")

cat("Equal annotation =", identical(result_faster$annotation, result$annotation), "\n")

## ----eval=local_file_exist----------------------------------------------------
system.time({
    protein_list <- iq::create_protein_list(norm_data)
    result <- iq::create_protein_table(protein_list)
})

system.time({
    result_faster <- iq::fast_MaxLFQ(norm_data)
})

## ----eval=local_file_exist----------------------------------------------------
sample_id  <- "R.FileName" 

secondary_id <- c("EG.Library", "FG.Id", "FG.Charge", "F.FrgIon", "F.Charge", "F.FrgLossType")

annotation_columns <- c("PG.Genes", "PG.ProteinNames")

iq_dat <- iq::fast_read("DIA-report-long-format.txt",
                        sample_id  = sample_id, 
                        secondary_id = secondary_id,
                        filter_string_equal = c("F.ExcludedFromQuantification" = "False"),
                        annotation_col = annotation_columns)

iq_norm_data <- iq::fast_preprocess(iq_dat$quant_table)

result_fastest <- iq::fast_MaxLFQ(iq_norm_data, 
                                  row_names = iq_dat$protein[, 1], 
                                  col_names = iq_dat$sample)

## ----eval = local_file_exist--------------------------------------------------
cat("Max difference =", max(abs(result_fastest$estimate - result$estimate), na.rm = TRUE), "\n")

cat("Identical NAs =", identical(is.na(result_fastest$estimate), is.na(result$estimate)), "\n")

cat("Equal annotation =", identical(result_fastest$annotation, result$annotation), "\n")

## ----eval = local_file_exist--------------------------------------------------
iq_extra_names <- iq::extract_annotation(rownames(result_fastest$estimate), 
                                         iq_dat$protein, 
                                         annotation_columns = annotation_columns)

write.table(cbind(Protein = rownames(result_fastest$estimate),
                  iq_extra_names[, annotation_columns],
                  MaxLFQ_annotation = result_fastest$annotation,
                  result_fastest$estimate), 
            "iq-MaxLFQ-fast.txt", sep = "\t", row.names = FALSE)

## ----eval = local_file_exist--------------------------------------------------
sample_id  <- "R.FileName" 

secondary_id <- c("EG.Library", "FG.Id", "FG.Charge", "F.FrgIon", "F.Charge", "F.FrgLossType")

annotation_columns <- c("PG.Genes", "PG.ProteinNames")

system.time({
    
    # reading data
    raw <- read.delim("DIA-report-long-format.txt")

    # filtering
    selected <- raw$F.ExcludedFromQuantification == "False" & 
                !is.na(raw$PG.Qvalue) & raw$PG.Qvalue < 0.01 &
                !is.na(raw$EG.Qvalue) & raw$EG.Qvalue < 0.01

    raw <- raw[selected,]

    ## process

    norm_data <- iq::preprocess(raw, 
                                sample_id  = sample_id, 
                                secondary_id = secondary_id)

    protein_list <- iq::create_protein_list(norm_data)
    
    result <- iq::create_protein_table(protein_list)
    
})

system.time({
    iq_dat <- iq::fast_read("DIA-report-long-format.txt",
                            sample_id  = sample_id, 
                            secondary_id = secondary_id,
                            filter_string_equal = c("F.ExcludedFromQuantification" = "False"),
                            annotation_col = annotation_columns)

    iq_norm_data <- iq::fast_preprocess(iq_dat$quant_table)

    result_fastest <- iq::fast_MaxLFQ(iq_norm_data, 
                                      row_names = iq_dat$protein[, 1], 
                                      col_names = iq_dat$sample)
})

