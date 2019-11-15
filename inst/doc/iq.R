## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE, echo=FALSE------------------------------------
require("knitr")
local_file_exist <- file.exists("Bruderer15-DIA-longformat-compact.txt.gz")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("iq")

## ---- eval=TRUE, include = TRUE------------------------------------------
library("iq")

## ---- eval=FALSE, include = TRUE-----------------------------------------
#  raw <- read.delim("Bruderer15-DIA-longformat.txt")
#  
#  selected <- raw$F.ExcludedFromQuantification == "False" &
#              raw$F.FrgLossType == "noloss" &
#              (is.na(raw$PG.Qvalue) | raw$PG.Qvalue <= 0.01) &
#              (is.na(raw$EG.Qvalue) | raw$EG.Qvalue <= 0.01)
#  
#  raw <- raw[selected, c("R.Condition","PG.ProteinGroups", "EG.ModifiedSequence", "FG.Charge",
#                         "F.FrgIon", "F.Charge", "F.PeakArea", "PG.Genes", "PG.ProteinNames")]
#  
#  write.table(raw, "Bruderer15-DIA-longformat-compact.txt", sep = "\t", row.names = FALSE)

## ---- eval=local_file_exist, include=TRUE--------------------------------
raw <- read.delim(gzfile("Bruderer15-DIA-longformat-compact.txt.gz"))

## ----eval=local_file_exist, echo=TRUE, message = FALSE, include=TRUE-----
norm_data <- iq::preprocess(raw)
protein_list <- iq::create_protein_list(norm_data)
result <- iq::create_protein_table(protein_list)

## ---- eval=local_file_exist, include=TRUE, fig.width=6.5, fig.height=4.5----
hist(log2(raw[, "F.PeakArea"]), 100, main = "Histogram of log2 intensities", 
     col = "steelblue", border = "steelblue", freq = FALSE)

## ---- eval = FALSE, include=TRUE-----------------------------------------
#  write.table(cbind(Protein = rownames(result$estimate),
#                    result$estimate,
#                    annotation = result$annotation),
#              "output-maxLFQ.txt", sep = "\t", row.names = FALSE)

## ---- eval=local_file_exist, include=TRUE, fig.width=6.5, fig.height=4.5----
iq::plot_protein(protein_list$P00366, main = "Protein P00366", split = NULL)  

## ---- eval=local_file_exist, include=TRUE, fig.width=6.5, fig.height=4.5----
iq::plot_protein(rbind(protein_list$P00366, 
                       MaxLFQ = iq::maxLFQ(protein_list$P00366)$estimate), 
                 main = "MaxLFQ quantification of P00366", 
                 col = c(rep("gray", nrow(protein_list$P00366)), "green"), 
                 split = NULL)  

## ---- eval=local_file_exist, include=TRUE, fig.width=6.5, fig.height=4.5----
iq::plot_protein(protein_list$P00366, main = "Protein P00366", cex = 0.4)  

## ---- eval=local_file_exist, include=TRUE, fig.width=6.5, fig.height=4.5----
MaxLFQ_estimate <- iq::maxLFQ(protein_list$P12799)$estimate

ground_truth <-  log2(rep(c(200, 125.99, 79.37, 50, 4, 2.52, 1.59, 1), each = 3))
ground_truth <- ground_truth - mean(ground_truth) + mean(MaxLFQ_estimate)

iq::plot_protein(rbind(MaxLFQ = MaxLFQ_estimate,
                       Groundtruth = ground_truth), 
                 main = "P12799 - MaxLFQ versus groundtruth",  
                 split = 0.75, 
                 col = c("green", "gold"))  

## ---- eval=FALSE, include = TRUE-----------------------------------------
#  extra_names <- iq::extract_annotation(rownames(result$estimate),
#                                        raw,
#                                        annotation_columns = c("PG.Genes", "PG.ProteinNames"))
#  
#  write.table(cbind(Protein = rownames(result$estimate),
#                    extra_names[, c("PG.Genes", "PG.ProteinNames")],
#                    result$estimate,
#                    annotation = result$annotation),
#              "output-maxLFQ-annotation.txt", sep = "\t", row.names = FALSE)

## ---- eval=FALSE, include = TRUE-----------------------------------------
#  tab <- read.delim("./Mtb_feature_alignment_requant_filtered_max10_fixed_noUPS.tsv",
#                    stringsAsFactors = FALSE)
#  
#  tab$Condition[tab$Condition == "d20_6h"] <- "d20_06h"
#  
#  tab_list <- vector("list", nrow(tab))
#  for (i in 1:nrow(tab)) {
#      a <- unlist(strsplit(tab[i, "aggr_Fragment_Annotation"], ";"))
#      b <- unlist(strsplit(tab[i, "aggr_Peak_Area"], ";"))
#  
#      tab_list[[i]] <- NULL
#      for (j in 1:length(a)) {
#          tab[i, "aggr_Fragment_Annotation"] <- a[j]
#          tab[i, "aggr_Peak_Area"] <- b[j]
#          tab_list[[i]] <- rbind(tab_list[[i]], tab[i,])
#      }
#  }
#  
#  tab_extended <- do.call(rbind.data.frame, tab_list)
#  
#  quant <- as.double(tab_extended$aggr_Peak_Area)
#  
#  short_name <- paste(tab_extended$Condition, tab_extended$BioReplicate,
#                      tab_extended$Run, sep = "_")
#  
#  tab_small <- cbind(tab_extended[, c("ProteinName", "FullPeptideName", "Charge",
#                                      "aggr_Fragment_Annotation")], quant, short_name)
#  
#  write.table(tab_small, "Schubert15-OpenSWATH.txt", sep = "\t", row.names = FALSE)

## ---- eval=local_file_exist, include=TRUE--------------------------------
tab_small <- read.delim(gzfile("Schubert15-OpenSWATH.txt.gz"))

## ---- eval=local_file_exist, include=TRUE--------------------------------
head(tab_small)

## ---- eval=local_file_exist, include=TRUE, fig.width=6.5, fig.height=4.5----
hist(log2(tab_small[, "quant"]), 100, main = "Histogram of log2 intensities", 
     col = "steelblue", border = "steelblue", freq = FALSE)

## ---- eval=FALSE, include = TRUE-----------------------------------------
#  norm_data <- iq::preprocess(tab_small,
#                              primary_id = "ProteinName",
#                              secondary_id = c("FullPeptideName", "Charge",
#                                               "aggr_Fragment_Annotation"),
#                              sample_id = "short_name",
#                              intensity_col = "quant")
#  
#  protein_list <- iq::create_protein_list(norm_data)
#  
#  result <- iq::create_protein_table(protein_list)
#  
#  write.table(cbind(Protein = rownames(result$estimate),
#                    result$estimate,
#                    annotation = result$annotation),
#              "Schubert-output-maxLFQ.txt", sep = "\t", row.names = FALSE)

## ---- eval=local_file_exist, include = TRUE------------------------------
dda <- read.delim(gzfile("proteinGroups.txt.gz"))
dda <- subset(dda, Reverse == "")    # remove reversed entries
rownames(dda) <- dda[,"Protein.IDs"] # use protein group ids as rownames
lfq <- grep("^LFQ", colnames(dda))
dda_log2 <- log2(dda[, lfq])
dda_log2[dda_log2 == -Inf] <- NA

colnames(dda_log2) <- sprintf("C%02d", 1:24)

## ---- eval=local_file_exist, include = TRUE------------------------------
evidence <- read.delim(gzfile("evidence.txt.gz"), stringsAsFactors = FALSE)
rownames(evidence) <- evidence[, "id"]

# median normalization
ex <- paste0(paste0("sample", rep(1:8, each=3),"_R0"), rep(1:3, 8))
ex_median <- rep(NA, length(ex))
names(ex_median) <- ex
for (i in ex) {
    tmp <- subset(evidence, Experiment == i)
    ex_median[i] <- median(tmp[,"Intensity"], na.rm = TRUE)
}
f <- mean(ex_median)/ex_median
evidence[, "Intensity"] <- evidence[, "Intensity"] * f[evidence[, "Experiment"]]

# create a protein list
p_list <- list()

for (i in 1:nrow(dda)) {
    tmp <- unlist(strsplit(as.character(dda[i, "Evidence.IDs"]), ";"))
    a <- evidence[tmp,]
    b <- data.frame(cn = a[, "Experiment"], 
                    rn = paste(a[, "Modified.sequence"], a[, "Charge"], sep="_"), 
                    quant = a[, "Intensity"])
    b <- b[!is.na(b$quant),]
    m <- matrix(NA, nrow = nlevels(b$rn), ncol = length(ex), 
                dimnames = list(levels(b$rn), ex))
  
    if (nrow(b) > 0) {
        for (j in 1:nrow(b)) {
            rn <- as.character(b$rn[j])
            cn <- as.character(b$cn[j])
            if (is.na(m[rn, cn])) {
                m[rn, cn] <- b$quant[j]
            } else {
                m[rn, cn] <- m[rn, cn] + b$quant[j]
            }
        }
        p_list[[rownames(dda)[i]]] <- log2(m)
    } else {
        p_list[[rownames(dda)[i]]] <- NA
    }
}

## ---- eval=local_file_exist, include=TRUE, fig.width=6.5, fig.height=4.5----
w1 <- iq::maxLFQ(p_list$A1L0T0)$estimate
w2 <- as.numeric(dda_log2["A1L0T0", ])
w2 <- w2 - mean(w2, na.rm = TRUE) + mean(w1, na.rm = TRUE)
tmp <- rbind(p_list$A1L0T0,
             `MaxLFQ by iq` = w1,
             `MaxLFQ by MaxQuant` = w2)
colnames(tmp) <- sprintf("C%02d", 1:24)

iq::plot_protein(tmp,
                 main = "A1L0T0", cex = 0.5, split = 0.65, 
                 col = c(rep("gray", nrow(p_list$A1L0T0)), "green", "blue"))  

## ---- eval=FALSE, include = TRUE-----------------------------------------
#  output_mq <- iq::create_protein_table(p_list)
#  
#  write.table(cbind(Protein = rownames(output_mq$estimate),
#                    output_mq$estimate,
#                    annotation = output_mq$annotation),
#              "output_IQ_LFQ.txt", sep = "\t", row.names = FALSE)

## ---- eval=FALSE, include = TRUE-----------------------------------------
#  # default MaxLFQ
#  output <- iq::create_protein_table(protein_list)
#  
#  # median polish
#  output <- iq::create_protein_table(protein_list, method = "medpolish")
#  
#  # top 3
#  output <- iq::create_protein_table(protein_list, method = "topN", N = 3)
#  
#  # top 5
#  output <- iq::create_protein_table(protein_list, method = "topN", N = 5)
#  
#  # MeanInt in the original intensity space
#  output <- iq::create_protein_table(protein_list, method = "meanInt",
#                                     aggregation_in_log_space = FALSE)

