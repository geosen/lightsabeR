#'Create a count table
#'
#'Creates a count table using multiple HTSeq output files as input
#'
#' 
#'
#'
#'date 08/03/2023
#'@param path Folder containing HTSeq output files 
#'@param pattern The pattern of file extensions to check. Defaults to '.txt'
#'@param remove_last_5 Whether to remove the last 5 rows containing ambiguous (etc.) counts or not
#'
#'
#'@return count table
#'@export

combine_htseq_files <- function(indir, pattern = '*.txt', remove_last_5 = TRUE) {
#Get list of all files in the directory

files<-list.files(path = indir, 
                  pattern = pattern)
#Merge the first to files and store
file1<-read.table(paste0(indir,files[1]),col.names = c("gene_id",files[1]))
file2<-read.table(paste0(indir,files[2]),col.names = c("gene_id",files[2]))
out.file<-merge(file1,file2,by="gene_id")

#For loop to merge content of remaining files
for(i in 3:length(files))
{
  file=read.table(paste0(indir,files[i]),col.names = c("gene_id",files[i]))
  out.file=merge(out.file,file,by="gene_id")
}
rownames(out.file) <- out.file$gene_id
out.file <- out.file[,-1]

if (remove_last_5){
  out.file <- out.file[6:dim(out.file)[1],]
}

out.file
}
