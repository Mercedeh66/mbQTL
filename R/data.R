#' `mbQTL` "CovFile"
#'
#' The "CovFile" is the covariate file for linear regression option of taxa and snp association.
#' The covariance file is generated randomly by assigning sex and site of collection to
#' the samples.) rownames are covariate and colnames samples.
#' @name CovFile
#' @docType data
#' @keywords data
NULL

#' `mbQTL` "SnpFile"
#' 
#' The "SnpFile" is randomly generated from GATK (Van der Auwera GA & O'Connor BD. (2020).
#' Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition).
#' O'Reilly Media) snp calls followed by plink (Purcell S, et al. (2007) PLINK: a
#' toolset for whole-genome association and population-based linkage analysis. American
#' Journal of Human Genetics) processing and it needs to be in (0,1,2) format representing
#' the zygosity of the snps.
#'
#' @name SnpFile
#' @docType data
#' @keywords data
NULL

#'  `mbQTL` "metagenomeSeqObj"
#'  
#' "MetagenomSeqObj" is an `MRexperiment` object format of the "microbeAbund" file.
#'
#' @name metagenomeSeqObj
#' @docType data
#' @keywords data
NULL

#' `mbQTL` "microbiomeAbund" File
#' 
#' This is the microbial Abudnance file generated from 16s it is either this file or
#' the "metaGenomeSeqObj".The "microbiomeAbund" file is a randomly generated file
#' format for total microbe presence (number of reads)(parasite/viral transcripts)
#' for specific species.This could be generated from select taxa results from 
#' `biom()` or `MRexperiment` objects as long as the samples are colnames and taxa
#'  are rownames.
#'
#' @name microbeAbund
#' @docType data
#' @keywords data
NULL
