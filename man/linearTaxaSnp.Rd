% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Linear_regression.R
\name{linearTaxaSnp}
\alias{linearTaxaSnp}
\title{linearTaxaSnp Performs linear regression analysis between taxa and SNPs and returns concordance statistics}
\usage{
linearTaxaSnp(microbeAbund, SnpFile, Covariate = NULL)
}
\arguments{
\item{microbeAbund}{the taxa abundance dataframe (rownames sample names and colnames taxa Genus/species/family)}

\item{SnpFile}{the snp dataframe (values 0,1,2 indicating zygosity), rownames sample names and colnames snp names.}

\item{Covariate}{default is NULL, hence assumed non-existent. If covariates are available they need to be formatted
in the CovFile format, that is colnames are sample numbers matching samples in the microbe abundance and snp file
and row names are the co-variates names (such as sex, disease etc).}
}
\value{
A data frame which is a result of Linear Regression of all snp, taxa relationships,
with P values and P value corrected values.
}
\description{
This function creates a dataframe output from the results all snps with all taxa linear regression analysis of all
snps in the dataset. The result is a dataframe with P values and FDRs of all regressions. MatrixeQTL core functions
are utilized to achieve this. Note the main functions used are Matrix_eQTL_engine() assuming linear regression
with or without a covariate file.
}
\examples{
data(microbeAbund)
data(SnpFile)
data(CovFile)
x <- linearTaxaSnp(microbeAbund, SnpFile, Covariate = CovFile)

}
\keyword{LR}
\keyword{linear}
\keyword{regression}
\keyword{snp}
\keyword{taxa}
