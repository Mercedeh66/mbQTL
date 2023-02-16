#' @import ggplot2
#' @import readxl
#' @import tidyr
#' @import dplyr
#' @import stringr
#' @importFrom stats as.formula glm p.adjust
#' @importFrom methods is
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".", "term", "p.value", "Pathogen", "genotype"
  ))
}

## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  binarizeMicrobe binarizes microbe abundace file based on user's cutoff
#'
#' This function creates a dataframe output produces a formatted dataframe prepared.
#' @param microbeAbund the taxa abundance dataframe (rownames sample names and colnames taxa Genus/species/family)
#' @param cutoff cutoff at which the user chose to call taxa positive or negative across samples (should be a numeric value
#' for normalized or count values).
#' @param selectmicrobe default is and all taxa are considered at the same time, if the user is interested in a
#' specific pathogen use name of the pathogen for example "Haemophilus".
#' @return A data frame of microbial abundance.
#' @keywords logitPlotDataframe logitParsing
#'
binarizeMicrobe <- function(microbeAbund, cutoff = NULL, selectmicrobe = NULL) {
  if (is.null(cutoff)) { # If value of rsID is not null do this:
    micrbial_abundance_2 <- microbeAbund %>% mutate_if(is.numeric, ~ 1 * (. != 0))
  } else {
    micrbial_abundance_2 <- microbeAbund %>% mutate_if(is.numeric, ~ 1 * (. > cutoff))
  }
  if (is.null(selectmicrobe)) { # If value of rsID is not null do this:
    micrbial_abundance_3 <- micrbial_abundance_2
    micrbial_abundance_3$ID <- rownames(micrbial_abundance_3)
    return(micrbial_abundance_3)
  } else {
    micrbial_abundance_3 <- micrbial_abundance_2 %>% select(all_of(selectmicrobe))
    micrbial_abundance_3$ID <- rownames(micrbial_abundance_3)
    return(micrbial_abundance_3)
  }
}

## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  prepareCorData prpares and joins snp-taxa and taxa-taxa correlation file.
#'
#' This function creates a dataframe output produces a formatted dataframe prepared.
#' @param microbeAbund the taxa abundance dataframe (rownames sample names and colnames taxa Genus/species/family)
#' @param SnpFile the snp dataframe (values 0,1,2 indicating zygosity), rownames sample names and colnames snp names.
#' @param cutoff default is NULL, hence anything above cutoff is considered positive, the cutoff at which the
#' specific or all taxa are considered positive for the pathogen (1 indicates positive and 0 negative).
#' @param selectmicrobe default is and all taxa are considered at the same time, if the user is interested in a
#' specific pathogen use name of the pathogen for example "Haemophilus".
#' @return A data frame which is a result of Logistic regression  products of individual snp, taxa relationships,
#' with P values and P value corrected values.
#' @keywords logitPlotDataframe
#'
## Function to produce dataframe for plot
prepareCorData <- function(microbeAbund, SnpFile, cutoff = NULL, selectmicrobe = NULL) {
  # Function binarizes based on user cutoff and either select microbe or look across all
  microbial_abundance_4 <- binarizeMicrobe(microbeAbund, cutoff = cutoff, selectmicrobe = selectmicrobe)
  n_snps <- ncol(SnpFile)
  SnpFile$ID <- rownames(SnpFile)
  SNPAss_filt <- SnpFile %>%
    pivot_longer(seq.int(1, n_snps), names_to = "rsID", values_to = "count")
  final_DF_Logit <- left_join(SNPAss_filt, microbial_abundance_4, by = "ID")
  return(final_DF_Logit)
}


## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  logRegSnpsTaxa Performs logistic regression analysis between taxa and SNPs and returns concordance statistics
#'
#' This function creates a dataframe output from the results of either a unique taxa and all snps or all taxa and
#' all snps in the dataset. The result is a dataframe with P values and FDRs of all regressions.
#' @param microbeAbund the taxa abundance dataframe (rownames sample names and colnames taxa Genus/species/family)
#' @param SnpFile the snp dataframe (values 0,1,2 indicating zygosity), rownames sample names and colnames snp names.
#' @param cutoff default is NULL, hence anything above cutoff is considered positive, the cutoff at which the
#' specific or all taxa are considered positive for the pathogen (1 indicates positive and 0 negative).
#' @param selectmicrobe default is and all taxa are considered at the same time, if the user is interested in a
#' specific pathogen use name of the pathogen for example "Haemophilus".
#' @return A data frame which is a result of Logistic regression  products of individual snp, taxa relationships,
#' with P values and P value corrected values (FDR, Bonferroni).
#' @export
#' @keywords taxa snp logistic_regression logit
#' @examples
#' data(microbeAbund)
#' data(SnpFile)
#' x <- logRegSnpsTaxa(microbeAbund, SnpFile, selectmicrobe = c("Haemophilus"))
#'
logRegSnpsTaxa <- function(microbeAbund, SnpFile, cutoff = NULL, selectmicrobe = NULL) {
  stopifnot("data frame is the expected input for SnpFile" = is(SnpFile, "data.frame"))
  stopifnot("data frame is the expected input for microbeAbund" = is(microbeAbund, "data.frame"))
  stopifnot("The value for 'cutoff'" = is.null(cutoff) || is.numeric(cutoff))
  stopifnot("The value for 'selectmicrobe'" = is.null(cutoff) || is.character(selectmicrobe))
  # Function binarizes based on user cutoff and either select microbe or look across all
  final_DF_Logit <- prepareCorData(microbeAbund, SnpFile,
    cutoff = cutoff,
    selectmicrobe = selectmicrobe
  )

  # if (is.null(selectmicrobe)){
  non_microb <- c("ID", "rsID", "count")
  microbes <- setdiff(colnames(final_DF_Logit), non_microb)
  chelo <- lapply(
    microbes,
    function(m) {
      x <- final_DF_Logit %>% group_by(rsID)
      ret <- do(x, tidy(glm(as.formula(paste(m, "~ count")), family = "binomial", data = .)))
      ret$microbe <- m
      return(ret)
    }
  ) %>%
    bind_rows()
  chelo_play <- filter(chelo, term == c("count"))
  # Add FDR coloumn or Bonf to the col
  chelo_play$FDR <- p.adjust(chelo_play$p.value, method = "fdr")
  chelo_play$BC <- p.adjust(chelo_play$p.value, method = "BH")
  # Sort by P values
  chelo_playF <- chelo_play %>% arrange(p.value)
  chelo_playFin <- chelo_playF %>% select(-term)
  return(chelo_playFin)
}


# log_link_res<-logRegSnpsTaxa(microbeAbund,SnpFile)
# log_link_res<-logRegSnpsTaxa(microbeAbund,SnpFile,selectmicrobe = c("Haemophilus"))


## Plot function for Logistic Regression
## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  logitPlotSnpTaxa produces bar plots for counts of ref vs alt vs het allells for
#'  particular rsID taxa combinations
#'
#' This function creates a dataframe output produces a formatted dataframe prepared.
#' @param SnpFile original snp file with (0,1,2 values for ref, het, alt genotypes), colnames
#'                SNP names, rownames, sample IDs.
#' @param microbeAbund original microbe abundance file (colnames microbe, rownames= sample IDs)
#' @param selectmicrobe name of the microbe of interest (for example individual significant microbes
#'                      associate with a snp).
#' @param rsID name of the snp of interest (for example individual significant snps associated with
#'             a microbe)
#' @param ref the name of reference genotype for example "GG"
#' @param alt the name of snp (variant) genotype for example "AA"
#' @param het the name of hetrozygote genotype for example "GA"
#' @param color the default is NULL and the color is set to c("#ffaa1e", "#87365b").
#' @param cutoff cutoff at which we call microbe present or absent
#' @return A bar ggplot comparing the counts of ref vs alt vs het genotype
#' @export
#' @keywords barplot logitplot
#' @examples
#' data(microbeAbund)
#' data(SnpFile)
#' x <- logitPlotSnpTaxa(microbeAbund, SnpFile,
#'   selectmicrobe = "Neisseria", rsID = "chr2.241072116_A",
#'   ref = NULL, alt = NULL, het = NULL, color = NULL, cutoff = NULL
#' )
logitPlotSnpTaxa <-
  function(microbeAbund, SnpFile, selectmicrobe = NULL,
           rsID, ref = NULL, alt = NULL, het = NULL, color = NULL, cutoff = NULL) {
    stopifnot("data frame is the expected input for SnpFile" = is(SnpFile, "data.frame"))
    stopifnot("data frame is the expected input for microbeAbund" = is(microbeAbund, "data.frame"))
    stopifnot("The value for 'rsID'" = is.null(cutoff) || is.character(rsID))
    stopifnot("The value for 'selectmicrobe'" = is.null(cutoff) || is.character(selectmicrobe))
    stopifnot("The value for 'ref'" = is.null(ref) || is.character(ref))
    stopifnot("The value for 'alt'" = is.null(alt) || is.character(alt))
    stopifnot("The value for 'het'" = is.null(het) || is.character(het))
    stopifnot("The value for 'color'" = is.null(color) || is.character(color))
    stopifnot("The value for 'cutoff'" = is.null(cutoff) || is.numeric(cutoff))
    final_DF_Logit <- prepareCorData(microbeAbund, SnpFile,
      cutoff = cutoff,
      selectmicrobe = selectmicrobe
    )

    final_DF_Logit_P <- as.data.frame(final_DF_Logit)
    final_DF_Logit_P_rs <- final_DF_Logit_P %>% filter(rsID == !!rsID)
    final_DF_Logit_plot <- final_DF_Logit_P_rs %>% select(all_of(c("ID", "count", "rsID", selectmicrobe)))

    # Ensure relevant columns are factors
    for (colname in c("count", "rsID", selectmicrobe)) {
      final_DF_Logit_plot[[colname]] <- as.factor(final_DF_Logit_plot[[colname]])
    }

    if (!is.null(ref)) {
      final_DF_Logit_plot_cased <- final_DF_Logit_plot %>%
        mutate(genotype = case_when(
          count == 2 ~ alt,
          count == 0 ~ ref,
          count == 1 ~ het
        ))
    } else {
      final_DF_Logit_plot_cased <- final_DF_Logit_plot %>%
        mutate(genotype = case_when(
          count == 2 ~ "2",
          count == 0 ~ "0",
          count == 1 ~ "1"
        ))
    }
    count_final_DF_Logit_plot_cased2 <-
      final_DF_Logit_plot_cased %>%
      mutate(Pathogen = case_when(
        final_DF_Logit_plot_cased[[selectmicrobe]] == 0 ~ "Absent",
        final_DF_Logit_plot_cased[[selectmicrobe]] == 1 ~ "Present"
      ))
    count_plot_file <- count_final_DF_Logit_plot_cased2 %>% count(Pathogen, genotype)
    if (is.null(color)) {
      # If value of rsID is not null do this:
      p <-
        ggplot(
          count_plot_file,
          aes(x = genotype, y = n, fill = Pathogen)
        ) +
        geom_bar(position = "dodge", stat = "identity") +
        scale_fill_manual(values = c("#ffaa1e", "#87365b")) +
        ylab("Genotype Count") +
        xlab(paste(rsID, "count"))
      return(p)
    } else {
      p <-
        ggplot(
          count_plot_file,
          aes(x = genotype, y = n, fill = Pathogen)
        ) +
        geom_bar(position = "dodge", stat = "identity") +
        scale_fill_manual(values = color) +
        ylab("Genotype Count") +
        xlab(paste(rsID, "count"))
      return(p)
    }
  }


# logit_plot <- logitPlotSnpTaxa(microbeAbund, SnpFile, selectmicrobe = "Neisseria",
#                                 rsID="chr2.241072116_A", ref="GG", alt="AA",het="AG")
