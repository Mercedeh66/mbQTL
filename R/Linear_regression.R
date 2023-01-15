#' @import ggplot2
#' @import readxl
#' @import tidyr
#' @import dplyr
#' @import stringr
#' @import MatrixEQTL
#' @importFrom graphics hist

NULL

## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  linearTaxaSnp Performs linear regression analysis between taxa and SNPs and returns concordance statistics
#'
#' This function creates a dataframe output from the results all snps with all taxa linear regression analysis of all
#' snps in the dataset. The result is a dataframe with P values and FDRs of all regressions.MatrixeQRL core functions
#' are utilized to achieve this.
#' @param microbeAbund the taxa abundance dataframe (rownames sample names and colnames taxa Genus/species/family)
#' @param SnpFile the snp dataframe (values 0,1,2 indicating zygosity), rownames sample names and colnames snp names.
#' @param Covariate default is NULL, hence assumed non-existent. If covariates are available they need to be formatted
#' in the CovFile format, that is colnames are sample numbers matching samples in the microbe abundance and snp file
#'  and row names are the co-variates names (such as sex, disease etc).
#' @return A data frame which is a result of Linear Regression of all snp,taxa relationships,
#' with P values and P value corrected values.
#' @export
#' @keywords taxa snp linear regression LR
#' @examples
#' x <- linearTaxaSnp(microbeAbund, SnpFile, Covariate = CovFile)
#'
linearTaxaSnp <- function(microbeAbund, SnpFile, Covariate = NULL) {
  microbeAbund_i <- t(microbeAbund)
  SnpFile_i <- t(SnpFile)
  if (is.null(Covariate)) {
    useModel <- modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS
    pvOutputThreshold <- 5e-2
    errorCovariance <- numeric()
    snps <- SlicedData$new()
    SNPs_t_f2 <- as.matrix(SnpFile_i)
    snps$CreateFromMatrix(SNPs_t_f2)
    snps$fileOmitCharacters <- "NA" # denote missing values;
    snps$fileSliceSize <- 2000 # read file in pieces of 2,000 rows
    gene <- SlicedData$new()
    microbe_file_t_f2 <- as.matrix(microbeAbund_i)
    gene$CreateFromMatrix(microbe_file_t_f2)
    gene$fileDelimiter <- "\t" # the TAB character
    gene$fileOmitCharacters <- "NA" # denote missing values;
    gene$fileSliceSize <- 1000 # read file in pieces of 2,000 rows

    outPut_linear <- tempfile()
    # For P value dataframe
    me <- Matrix_eQTL_engine(
      snps = snps,
      gene = gene,
      output_file_name = outPut_linear,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      pvalue.hist = TRUE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE
    )

    kettle <- as.data.frame(cbind(
      me$all[["eqtls"]][["snps"]], me$all[["eqtls"]][["gene"]],
      me$all[["eqtls"]][["pvalue"]], me$all[["eqtls"]][["FDR"]]
    ))
    names(kettle) <- c("snps", "Genus", "pvalue", "FDR")
    return(kettle)
  } else {
    useModel <- modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS
    pvOutputThreshold <- 5e-2
    errorCovariance <- numeric()
    snps <- SlicedData$new()
    SNPs_t_f2 <- as.matrix(SnpFile_i)
    snps$CreateFromMatrix(SNPs_t_f2)
    snps$fileOmitCharacters <- "NA" # denote missing values;
    snps$fileSliceSize <- 2000 # read file in pieces of 2,000 rows

    gene <- SlicedData$new()
    microbe_file_t_f2 <- as.matrix(microbeAbund_i)
    gene$CreateFromMatrix(microbe_file_t_f2)
    gene$fileDelimiter <- "\t" # the TAB character
    gene$fileOmitCharacters <- "NA" # denote missing values;
    gene$fileSliceSize <- 1000 # read file in pieces of 1,000 rows

    cvrt <- SlicedData$new()
    covariate_f <- as.matrix(Covariate)
    cvrt$CreateFromMatrix(covariate_f)
    cvrt$fileDelimiter <- "\t" # the TAB character
    cvrt$fileOmitCharacters <- "NA" # denote missing values;
    cvrt$fileSliceSize <- 10 # read file in pieces of 1,0 rows

    outPut_linear <- tempfile()
    # For P value dataframe
    me <- Matrix_eQTL_engine(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name = outPut_linear,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      pvalue.hist = TRUE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE
    )
    kettle <- as.data.frame(cbind(
      me$all[["eqtls"]][["snps"]], me$all[["eqtls"]][["gene"]],
      me$all[["eqtls"]][["pvalue"]], me$all[["eqtls"]][["FDR"]]
    ))
    names(kettle) <- c("snps", "Genus", "pvalue", "FDR")
    return(kettle)
  }
}

# LinearAnalysisTaxaSNP<-linearTaxaSnp(microbeAbund,SnpFile,Covariate=CovFile)
# LinearAnalysisTaxaSNP2<-linearTaxaSnp(microbeAbund,SnpFile)


## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  histPvalueLm a histogram of Taxa and SNP linear regression analysis.
#'  This function creates a histogram object of all SNPs with all taxa Linear regression analysis p values.
#' @param LinearAnalysisTaxaSNP the data frame result created from the linearTaxaSnp() function.
#' @return A histogram object of p values observed from taxa and SNP Linear Regression analysis.
#' @export
#' @keywords taxa snp linear_regression plot histogram
#' @examples
#' x <- qqPlotLm(microbeAbund, SnpFile, Covariate = CovFile)
#'
histPvalueLm <- function(LinearAnalysisTaxaSNP) {
  return(hist(as.numeric(LinearAnalysisTaxaSNP$pvalue),
    col = "grey",
    main = "Histogram of Taxa-SNP LM P-Values", xlab = "P value",
    ylab = "Number of Samples"
  ))
}

# histPvalueLm(LinearAnalysisTaxaSNP)


## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, January 2023
#'  qqPlotLm creates QQ-Plot of all SNPs with all taxa Linear regression analysis
#'  This function creates QQ-Plot object of all SNPs with all taxa Linear regression analysis of
#'  expected versus observed P values
#' @param microbeAbund the taxa abundance dataframe (rownames sample names and colnames taxa Genus/species/family)
#' @param SnpFile the snp dataframe (values 0,1,2 indicating zygosity), rownames sample names and colnames snp names.
#' @param Covariate default is NULL, hence assumed non-existent. If covariates are available they need to be formatted
#' in the CovFile format, that is colnames are sample numbers matching samples in the microbe abundance and snp file
#'  and row names are the covariates names (such as sex, disease etc).
#' @return A QQplot object of expected versus obsesrved taxa and SNP Linear Regression analysis
#' @export
#' @keywords taxa snp linear_regression plot
#' @examples
#' x <- qqPlotLm(microbeAbund, SnpFile, Covariate = CovFile)
#'
qqPlotLm <- function(microbeAbund, SnpFile, Covariate = NULL) {
  microbeAbund_i <- t(microbeAbund)
  SnpFile_i <- t(SnpFile)
  if (is.null(Covariate)) {
    # useModel <- modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS
    # pvOutputThreshold <- 5e-2
    # errorCovariance <- numeric()
    snps <- SlicedData$new()
    SNPs_t_f2 <- as.matrix(SnpFile_i)
    snps$CreateFromMatrix(SNPs_t_f2)
    snps$fileOmitCharacters <- "NA" # denote missing values;
    snps$fileSliceSize <- 2000 # read file in pieces of 2,000 rows
    gene <- SlicedData$new()
    microbe_file_t_f2 <- as.matrix(microbeAbund_i)
    gene$CreateFromMatrix(microbe_file_t_f2)
    gene$fileDelimiter <- "\t" # the TAB character
    gene$fileOmitCharacters <- "NA" # denote missing values;
    gene$fileSliceSize <- 1000 # read file in pieces of 2,000 rows

    # outPut_linear <- tempfile()
    # For qqPlot
    filename <- tempfile()
    meq <- Matrix_eQTL_engine(
      snps = snps,
      gene = gene,
      # cvrt = cvrt1,
      output_file_name = filename,
      pvOutputThreshold = 1e-6,
      useModel = modelLINEAR,
      errorCovariance = numeric(),
      verbose = TRUE,
      pvalue.hist = "qqplot"
    )
    unlink(filename)
    return(plot(meq, pch = 16, cex = 0.7))
  } else {
    # useModel <- modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS
    # pvOutputThreshold <- 5e-2
    # errorCovariance <- numeric()
    snps <- SlicedData$new()
    SNPs_t_f2 <- as.matrix(SnpFile_i)
    snps$CreateFromMatrix(SNPs_t_f2)
    snps$fileOmitCharacters <- "NA" # denote missing values;
    snps$fileSliceSize <- 2000 # read file in pieces of 2,000 rows

    gene <- SlicedData$new()
    microbe_file_t_f2 <- as.matrix(microbeAbund_i)
    gene$CreateFromMatrix(microbe_file_t_f2)
    gene$fileDelimiter <- "\t" # the TAB character
    gene$fileOmitCharacters <- "NA" # denote missing values;
    gene$fileSliceSize <- 1000 # read file in pieces of 1,000 rows

    cvrt <- SlicedData$new()
    # covariate_f <- as.matrix(Covariate)
    cvrt$fileDelimiter <- "\t" # the TAB character
    cvrt$fileOmitCharacters <- "NA" # denote missing values;
    cvrt$fileSliceSize <- 10 # read file in pieces of 1,0 rows

    filename <- tempfile()
    meq <- Matrix_eQTL_engine(
      snps = snps,
      gene = gene,
      # cvrt = cvrt1,
      output_file_name = filename,
      pvOutputThreshold = 1e-6,
      useModel = modelLINEAR,
      errorCovariance = numeric(),
      verbose = TRUE,
      pvalue.hist = "qqplot"
    )
    unlink(filename)
    return(plot(meq, pch = 16, cex = 0.7))
  }
}


# qqPlotLm(microbeAbund,SnpFile,Covariate=CovFile)
