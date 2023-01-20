#' @import readxl
#' @import tidyr
#' @import dplyr
#' @import stringr
#' @import metagenomeSeq
NULL

#' ## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  metagenomeSeqToMbqtl Converts metagenomeSeq obj to compatible taxa dataframe
#'
#' This function takes and MRexperiement class object transforms it and makes the result dataframe compatible
#' with mbQTL taxa input file
#' @param meta_glom MRexperiement class obj from metagenomeSeq package.
#' @param norm A logical indicating whether or not to return normalized counts.
#' @param log TRUE/FALSE whether or not to log2 transform scale.
#' @param aggregate_taxa it is recommended that the normalization occurs at taxa level (default NULL)
#'                       however, if the user chooses to aggregate on the phyla/family/Genus or Species
#'                       level before normalization they have the option.
#' @return A data frame of normalized/not normalized counts compatible with mbQTL.
#' @keywords metagenomeSeq MRexperiment normalization
#' @export
#' @examples
#' data(metagenomeSeqObj)
#' x <- metagenomeSeqToMbqtl(metagenomeSeqObj, norm = TRUE, log = TRUE, aggregate_taxa = NULL)
metagenomeSeqToMbqtl <- function(meta_glom, norm, log, aggregate_taxa = NULL) {
  if (is.null(aggregate_taxa)) {
    obj <- meta_glom
    mat <- t(MRcounts(obj, norm = norm, log = log))
    return(as.data.frame(mat))
  } else {
    obj <- metagenomeSeq::aggTax(meta_glom, lvl = aggregate_taxa)
    mat <- t(MRcounts(obj, norm = norm, log = log))
    return(as.data.frame(mat))
  }
}

# compatible_metagenome_seq<-metagenomeSeqToMbqtl(metagenomeSeqObj, norm=TRUE, log=TRUE,
#                                                 aggregate_taxa= "Genus")
