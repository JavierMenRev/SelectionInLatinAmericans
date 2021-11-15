## libs
suppressMessages(library(readr))
suppressMessages(library(dplyr))

## file names
temp <- commandArgs(TRUE)
f1 <- as.character(temp[1])
f2 <- as.character(temp[2])

## read files
samples <- read.table(paste0(f1),stringsAsFactors=F,h=F,skip = 2)
ids <- read.table(paste0(f2),stringsAsFactors=F,h=F)
colnames(samples)[2] <- "IID"
colnames(ids) <- c("IID","FID")

## add ids
samples <- left_join(samples,ids,by="IID")
samples$INCLUDE <- 1
samples <- samples[,c("IID","FID","INCLUDE")]

## output
cat(format_tsv(samples,col_names = F))
