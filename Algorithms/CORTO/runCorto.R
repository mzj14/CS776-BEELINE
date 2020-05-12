library("corto")
library("optparse")

option_list <- list (
  make_option(c("-e","--expressionFile"), type = 'character',
              help= "Path to comma separated file containing gene-by-cell matrix with
              cell names as the first row and gene names as 
              the first column."),
  make_option(c("-s","--sourceFile"), type = 'character',
              help= "Path to line separated file containing gene names."),
  make_option(c("-o","--outFile"), , type = 'character',
              help= "outFile name to write the output ranked edges."),
  make_option(c("-b","--nBootstraps"), type = 'integer',
              help= "number of bootstraps."),
  make_option(c("-p","--pValue"), type = 'numeric',
              help= "p value."),
  make_option(c("-t","--nThreads"), type = 'integer',
              help= "number of threads.")
  )

parser <- OptionParser(option_list = option_list)

arguments <- parse_args(parser, positional_arguments = FALSE)

inmat = data.matrix(read.csv(arguments$expressionFile, header=TRUE, row.names=1, check.names=FALSE))

centroids = scan(arguments$sourceFile, character(), quote="")

regulon<-corto(inmat,centroids=centroids,nbootstraps=arguments$nBootstraps,p=arguments$pValue,nthreads=arguments$nThreads)

sink(arguments$outFile)

source_num <- length(regulon)

for (i in 1:source_num) {
  source_name <- names(regulon[i])
  source <- get(source_name, regulon[i])
  target_names <- names(source$tfmode)
  likelihoods <- source$likelihood
  target_num <- length(target_names)
  for (j in 1:target_num) {
    cat(source_name, target_names[j], likelihoods[j], sep = '\t', fill = TRUE)
  }
}
