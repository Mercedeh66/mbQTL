#' title `mbQTL` is a package for microbial QTL/GWAS Analysis
#' 
#' @description 
#' 
#' This package provides statistical methods for identifying significant relationships between
#' microbial/taxa and genetic SNP signatures. We use three models 1) linear regression between 
#' all taxa-snp. Main function is `linearTaxaSnp()`. 2) Correlation analysis between taxa-snp
#' across all taxa and snps. Main function is `taxaSnpCor()` and 3) Logistic regression analysis
#' between each taxa and each snp simultaneously or for a specific cases. Main function is
#' `logRegSnpsTaxa()`.
#' 
#' @seealso
#' The package vignette can be accessed with `vignette("mbQTL")`.
#' 
"_PACKAGE"