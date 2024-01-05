#'Read GSEA results
#'
#'date 06/01/2023
#'@param dir Directory of GSEA results
#'
#'@return Data Frame
#'@export

read_gsea <- function(dir) {
  files <- list.files(dir,pattern = 'gsea_report_[[:alnum:]]{1+}.tsv')
  neg <- read.table(paste0(dir,'/',files[1]),sep = '\t', header = TRUE)
  pos <- read.table(paste0(dir,'/',files[2]),sep = '\t', header = TRUE)
  gsea_rep <- as.data.frame(rbind(pos,neg))
}