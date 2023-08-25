suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  optparse::make_option(
    c("-i", "--counts"), type="character", default=NULL,
    help="path to counts table (rows: sgRNAs, columns: Samples)", metavar="character"
  ),
  optparse::make_option(
    c("-s", "--samplesheet"), type="character", default=NULL,
    help="path to sample sheet (columns: Index, Treat, Rep)", metavar="character"
  ),
  optparse::make_option(
    c("-t", "--treatment"), type="character", default=NULL,
    help="name of samples in treatment condition", metavar="character"
  ),
  optparse::make_option(
    c("-c", "--control"), type="character", default=NULL,
    help="name of samples in control condition", metavar="character"
  ),
  optparse::make_option(
    c("-z", "--T0"), type="character", default=FALSE,
    help="name of samples in T0 condition", metavar="character"
  )
);

runDeseq <- function(countsTable, colData) {
  colData$Treat <- as.factor(colData$Treat) %>% relevel(ref = 'T0')

  ### Create DESeq Object
  dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = countsTable[,colData %>% rownames],
      colData = colData,
      design = ~ 0 + Treat
  )

  ### Normalize counts
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds)

  ## Run DESeq tests
  dds <- DESeq2::DESeq(dds, test="Wald")

  ## print DESeq tests result names
  message(DESeq2::resultsNames(dds))

  return(dds)
}

writeResult <- function(res,resName){
  write.table(res, resName, sep="\t", quote=FALSE, col=TRUE, row=TRUE)
}

getResult <- function(
  dds, numerator, denominator, sgRNA2gene, cond='Treat', write=FALSE,name=NULL
){
  res <- DESeq2::results(
    dds,
    contrast=list(paste0(cond,numerator),paste0(cond,denominator)),
    listValues=c(1, -1)
  )

  res <- cbind(sgRNA2gene, res %>% data.frame)

  if(write) {
    writeResult(
      res, paste0(name, 'R', numerator, '_vs_', denominator, '.txt')
    )
  } else {
    return(res)
  }
}

getPheScores <- function(dds, Treatment, Control, sgRNA2gene, T0){
  ## rho  – treatment vs ctrl
  getResult(dds, Treatment, Control, sgRNA2gene, write=TRUE, name='rho')

  ## gamma  – ctrl vs T0
  if (T0){
    getResult(dds, Control, T0, sgRNA2gene, write=TRUE, name='gamma')
  }

  ## tau  – treatment vs T0
  if (T0){
    getResult(dds, Treatment, T0, sgRNA2gene, write=TRUE, name='tau')
  }

  ## kappa  – combination treatment vs. single treatment
  # not implemented
}


opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

message('Load count matrix and samplesheet!')
RawTable <- read.table(opt$counts, check.names=FALSE)
sgRNA2gene <- RawTable %>% dplyr::select(target)
countsTable <- RawTable %>% dplyr::select(!target)

colData <- read.table(
  opt$samplesheet, check.names=FALSE, sep=',', header=1
) %>% tibble::column_to_rownames('Index')

message('Run test!')
dds <- runDeseq(countsTable, colData)

message('Extract Normalized counts!')
normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
writeResult(normalized_counts, gsub(
  '.txt',"_DESeq2_normalized.txt", opt$counts
))

setwd(dirname(opt$counts))

message('Extract results!')
getPheScores(dds, opt$treatment, opt$control, sgRNA2gene,  opt$T0)

message('DONE!')
