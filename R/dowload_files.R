## Written by Mercedeh Movassagh <mercedeh@ds.dfci.harvard.edu>, Jan 2023

#' @importFrom utils read.table download.file
#' @importFrom data.table fread
#' @import R.utils
NULL

#' import Read internal Miranda file
#'
#' Reads internal Miranda file from extdata and returns it as a data.frame
#' @param fn filename
#' @return data.frame containing Miranda data
#' @examples
#' \donttest{
#' x <- importMbQTLFile("Mouse_miRanda.txt")
#' }
importMbQTLFile <- function(fn) {
  fn <- system.file("extdata", "mbQTLPrepFiles", fn, package = "mbQTL", mustWork = TRUE)
  ret1 <- read.table(fn, as.is = TRUE, sep = "\t")
  return(ret1)
}


#' downloadMbQTLFile Read internal Miranda file
#'
#' Reads internal mbQTL file from extdata and returns it as a data.frame
#' @param urlf URL of the specific chosen file
#' @return data.frame containing downloaded mbQTL file
#' @examples
#' \donttest{
#' x <- downloadMbQTLFile(
#'        "https://zenodo.org/record/4615670/files/mbQTL.txt.gz"
#'      )
#' }
downloadMbQTLFile <- function(urlf) {
  return(downloadMbQTLFile_(urlf))
}

# use internal function here so we can R.cache it without changing how it looks like for the user.
downloadMbQTLFile_ <- function(urlf) {
  # make effective statement error if it returns true or false for http with or without s
  assert_that(grepl("^https?://.*$", urlf))
  tmp <- tempfile(fileext = ".gz")
  op <- options(timeout = 99999)
  on.exit({
    unlink(tmp)
    options(op)
  })
  download.file(urlf, tmp)
  ret1 <- fread(tmp, sep = "\t", header = FALSE)
  return(as.data.frame(ret1))
}
