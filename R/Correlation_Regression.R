#' @import ggplot2
#' @import readxl
#' @import tidyr
#' @import dplyr
#' @import stringr
#' @import pheatmap
#' @import broom
#' @importFrom stats cor quantile
NULL

# un-upset lintr
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Genus_1", "Genus_2", "product_r2", "R", "Genus1", "Genus2", "X1", "X2", "rsID",
    "rho_xy", "Combined_Genus"
  ))
}

## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  coringTaxa creates correlation dataframe for taxa
#'
#' This function creates an output correlation data frame for all microbial taxa (or other organisms such as
#' viral or parasitic taxa)
#' @param microbeAbund the taxa abundance dataframe (rownames sample names and colnames taxa Genus/species/family)
#' @return A data frame of correlations between taxa
#' @export
#' @keywords taxa correlation
#' @examples
#' data(microbeAbund)
#' x <- coringTaxa(microbeAbund)
#'
# Function to run correlation on the taxa-taxa
coringTaxa <- function(microbeAbund) {
  # start parsing
  ppr <- rownames(microbeAbund)
  vapply_type <- array(data = 0.0, dim = nrow(microbeAbund))
  microbeAbundA <- as.data.frame(vapply(microbeAbund, as.numeric, vapply_type))
  rownames(microbeAbundA) <- ppr
  coringTaxaR <- cor(microbeAbundA)
  coringTaxaF <- as.data.frame(coringTaxaR)
  # Make a long R dataframe
  coringTaxaF$Genus_1 <- rownames(coringTaxaF)
  coringTaxalong <- coringTaxaF %>%
    pivot_longer(cols = 0:length(rownames(coringTaxaF)), names_to = "Genus_2", values_to = "R")
  # Remove genus correlations which are the same for Genus_1 and 2 values
  coringTaxalongParsed <- coringTaxalong %>% filter(Genus_1 != Genus_2)
  coringTaxalongParsed2 <- coringTaxalongParsed %>%
    rename(Genus1 = Genus_1, Genus2 = Genus_2)
  return(as.data.frame(coringTaxalongParsed2))
}

# correlationMicrobes<-coringTaxa(microbeAbund)


## SNPs

## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  RegSnp creates a dataframe of parsed long snp files
#'
#' This internal function takes the orginal snp dataframe and retruns a long parsed snp dataframe
#' @param SnpFile the snp file (rownames is sample number and colnames is the snps)
#' @param microbeAbund the microbial abundance file (rownames is sample number and colnames is
#'                     the microbe)
#' @return A long parsed datframe of snps
#' @keywords snp parsed long
#' @examples
#' \dontrun{
#' x <- RegSnp(SnpFile, microbeAbund)
#' }
#'
# function to make the long snp format
RegSnp <- function(SnpFile, microbeAbund) {
  SnpFilt <- rownames(microbeAbund)
  SnpAss_filt <- SnpFile[rownames(SnpFile) %in% SnpFilt, ]
  mycor <- cor(cbind(SnpAss_filt, microbeAbund))
  mycor[is.na(mycor)] <- 0
  # Remove initial coloumns for rs correlation
  n_snps <- length(SnpFile)
  n_mycor <- ncol(mycor)
  # Remove the first cor with snps and
  # keep the rest of cor os snps with taxa as it starts there
  mycor_parse <- as.data.frame(mycor[, seq.int(n_snps + 1, n_mycor)]) # Use seqint avoid drama
  # Remove all the lines after rs ends as you do not need the taxa for row
  mycor_parse2 <- as.data.frame(mycor_parse[seq.int(1, n_snps), ])
  mycor_parse2_R2 <- mycor_parse2^2
  mycor_parse2_R2$rsID <- rownames(mycor_parse2_R2)
  mycor_parse2_long <- mycor_parse2_R2 %>%
    pivot_longer(cols = !any_of("rsID"), names_to = "Genus1", values_to = "R2")
  return(mycor_parse2_long)
}

## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  allToAllProduct creates a dataframe of snp and taxa correlations
#'
#' This internal function takes the original snp dataframe and returns a long parsed snp dataframe
#' @param SnpFile the snp file (rownames is sample number and colnames is the snps)
#' @param microbeAbund he taxa abundance dataframe (rownames sample names and colnames taxa Genus/species/family)
#' @param rsID Default is NULL and will run across the who dataset unless specific rsID/SNP/chr_region is specified
#' @return A dataframe of correlations between snps and taxa
#' @keywords snptaxa Correlation
#' @export
#' @examples
#' data(microbeAbund)
#' data(SnpFile)
#' x <- allToAllProduct(SnpFile, microbeAbund, "chr1.171282963_T")
#'
# function to make R2
allToAllProduct <- function(SnpFile, microbeAbund, rsID = NULL) {
  r2_long <- RegSnp(SnpFile, microbeAbund)
  # Define function apply_product
  # Function to solve on one of the items
  apply_product <- function(rsID) {
    filtered <- filter(r2_long, rsID == !!rsID)
    vapply_type <- array(data = 0.0, dim = length(filtered$R2))
    mtx <- vapply(filtered$R2, `*`, vapply_type, filtered$R2)
    rownames(mtx) <- filtered$Genus1
    colnames(mtx) <- filtered$Genus1
    mtx <- as.data.frame(1.0 / mtx)
    mtx$rsID <- rsID
    mtx$Genus1 <- rownames(mtx)
    result <- pivot_longer(mtx,
      cols = !any_of(c("rsID", "Genus1")),
      names_to = "Genus2", values_to = "product_r2"
    )
    return(result)
  }
  # Get rows for rsID
  if (!is.null(rsID)) { # If value of rsID is not null do this:
    return(apply_product(rsID))
    # Now apply to rest
  } else {
    rsids <- unique(r2_long$rsID)
    result <- lapply(rsids, apply_product) %>%
      bind_rows() %>%
      filter(Genus1 != Genus2)
    return(result)
  }
}

# one_rs_id<-allToAllProduct(SnpFile,microbeAbund,"chr1.171282963_T")
# for_all_rsids<-allToAllProduct(SnpFile,microbeAbund)


## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  taxaSnpCor for estimation of the rho value between snp, taxa correlations across datasets
#'
#' This function produces a log heatmap +1 of the correlation rho values across snp, taxa dataframe.
#' @param for_all_rsids A dataframe result of correlation analysis between the snps and taxa dataframe,
#'  an output of allToAllProduct() function.
#' @param correlationMicrobes A dataframe of correlation between coringTaxa() function.
#' @param probs Default is NULL if other that all rho values are wanted the value can be subseted using c(x,y).
#' @return A data frame of correlations between taxa
#' @export
#' @keywords rho estimation
#' @examples
#' data(microbeAbund)
#' data(SnpFile)
#'
#' for_all_rsids <- allToAllProduct(SnpFile, microbeAbund)
#' correlationMicrobes <- coringTaxa(microbeAbund)
#' x <- taxaSnpCor(for_all_rsids, correlationMicrobes)
#'
taxaSnpCor <- function(for_all_rsids, correlationMicrobes, probs = NULL) {
  invprods_and_cors <- for_all_rsids %>%
    left_join(correlationMicrobes, by = c("Genus1", "Genus2")) %>%
    filter(is.finite(product_r2)) %>%
    mutate(rho_xy = product_r2 * R)
  # Deduplicate gene combinations (i.e., "A", "B" == "B", "A", per rsID)
  invprods_and_cors <- invprods_and_cors %>%
    mutate(
      X1 = ifelse(Genus1 < Genus2, Genus1, Genus2),
      X2 = ifelse(Genus1 > Genus2, Genus1, Genus2)
    )
  invprods_and_cors$duplicated <- duplicated(invprods_and_cors %>% select(rsID, X1, X2))
  invprods_and_cors <- invprods_and_cors %>%
    filter(!duplicated) %>%
    select(-X1, -X2, -duplicated)
  if (is.null(probs)) { # If value of rsID is not null do this:
    return(invprods_and_cors)
    # Now apply to rest
  } else {
    extremes <- quantile(invprods_and_cors$rho_xy, probs = probs) ## ask about this
    final_var <- invprods_and_cors %>% filter(rho_xy < extremes[1] | rho_xy > extremes[2])
    final_var$Genus1 <- gsub("X.", "", final_var$Genus1)
    final_vara <- final_var %>%
      mutate_if(is.character, str_trim)
    final_var_long <- final_vara[final_vara$Genus1 != final_vara$Genus2, ]
    return(final_var_long)
  }
}

# taxa_SNP_Cor<-taxaSnpCor(for_all_rsids,correlationMicrobes)
# taxa_SNP_Cor_lim<-taxaSnpCor(for_all_rsids,correlationMicrobes,probs = c(0.0001, 0.9999))


## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  mbQtlCorHeatmap for making heatmap for snp, taxa rho values
#'
#' This function produces a log heatmap +1 of the correlation rho values across snp, taxa datasets
#' @param final_var_long the long data frame of rho values created by the taxaSnpCor() function.
#' @param labels_col set to NULL ass default if TRUE, labels will appear on the heatmap.
#' @param ... all other parameters for pheatmap.
#' @return A data frame of correlations between taxa
#' @export
#' @keywords rho heatmap
#' @examples
#'
#' data(microbeAbund)
#' data(SnpFile)
#' for_all_rsids <- allToAllProduct(SnpFile, microbeAbund)
#' correlationMicrobes <- coringTaxa(microbeAbund)
#' taxaSnpCor(for_all_rsids, correlationMicrobes)
#' final_var_long <- taxaSnpCor(for_all_rsids, correlationMicrobes, probs = c(0.0001, 0.9999))
#' x <- mbQtlCorHeatmap(final_var_long)
#'
mbQtlCorHeatmap <- function(final_var_long, labels_col = NULL, ...) {
  final_var_long$Combined_Genus <- paste(final_var_long$Genus1, final_var_long$Genus2, sep = "_")
  final_var2 <- final_var_long %>% select(rsID, Combined_Genus, rho_xy)
  final_var_wide <- final_var2 %>%
    pivot_wider(
      names_from = Combined_Genus, values_from = rho_xy,
      values_fill = 0
    ) %>%
    as.data.frame()
  rownames(final_var_wide) <- final_var_wide$rsID
  final_var_wide2 <- final_var_wide %>% select(-rsID)
  if (is.null(labels_col)) {
    return(pheatmap(as.matrix(log(abs(final_var_wide2) + 1)), ...))
  } else {
    return(pheatmap(as.matrix(log(abs(final_var_wide2) + 1)), labels_col = FALSE, ...))
  }
}
# mbQtlCorHeatmap(taxa_SNP_Cor_lim,fontsize_col=4,fontsize_row=4)
