#The microbiomeAbund file is a randomly generated file format for total microbe presence (number of reads)
#(parasite/viral transcripts) for specific species. This could be generated from select taxa results from biom() or
#MRexperiment objects as long as the samples are colnames and taxa are rownames.
microbeAbund <- read.table("data-raw/Taxa_table_MVP_cor.txt", sep = "\t", row.names = 1, header = TRUE)
#The SNP file is randomly generated from GATK (Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud:
#Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media) snp calls followed by plink (Purcell S, et al. (2007)
#PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics)
# processing and it needs to be in (0,1,2) format representing the zygosity of the snps.
SnpFile <- read.table("data-raw/SNP_table_MVP_cor.txt", sep = "\t", row.names = 1, header = TRUE)
#The covariance file is generated randomly by assinging sex and site of collection to the samples.)
#rownames are covariate and colnames samples.
CovFile <- read.table("data-raw/covariates.txt", sep = "\t", row.names = 1, header = TRUE)
#MetagenomSeqObj is an MRexperiment obj format of the microbeAbund file.
metagenomeSeqObj <- readRDS("data-raw/metagenomeSeqReads.RDS")
usethis::use_data(microbeAbund, SnpFile, CovFile, metagenomeSeqObj, overwrite = TRUE)
