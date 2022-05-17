## Convert VCF (phased or unphased) into ChromoPainter format
## Note that this Chromopainter format might contain missing data and is therefore only suitable for AdaptMix 
## Do not use output from this R script to run other softwares that use ChromPainter format files as input

## libs
suppressMessages(library(data.table))
suppressMessages(library(stringr))

## args
args <- commandArgs(T)
file_in <- as.character(args[1])
file_out <- as.character(args[2])

## read vcf
vcf <- data.table::fread(paste0(file_in,".vcf"),skip="CHROM",stringsAsFactors = F, data.table = F)

## get number of haplotypes
pos_genos <- (which(colnames(vcf)=="FORMAT"))+1
num_haps <- length(pos_genos:ncol(vcf))*2

## get number of sites
num_sites <- nrow(vcf)

## get positions
pos_col <- which(colnames(vcf)=="POS")
positions <- vcf[,pos_col]
positions <- c("P",positions)

## get genotypes
vcf$Genotypes <- apply(vcf[,pos_genos:ncol(vcf)],1,paste,collapse="")
vcf$Genotypes <- stringr::str_replace_all(vcf$Genotypes, "[/]", "") ## for unphased data
vcf$Genotypes <- stringr::str_replace_all(vcf$Genotypes, "[|]", "") ## for phased data
vcf$Genotypes <- stringr::str_replace_all(vcf$Genotypes, "[.]", "?") ## for missing data

## save
write.table(num_haps,file=paste0(file_out,".chromopainter.haps"),row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(num_sites,file=paste0(file_out,".chromopainter.haps"),row.names=FALSE,col.names=FALSE,quote=FALSE,append=T)
write.table(t(positions),file=paste0(file_out,".chromopainter.haps"),row.names=FALSE,col.names=FALSE,quote=FALSE,append=T)
write.table(vcf$Genotypes,file=paste0(file_out,".chromopainter.haps"),row.names=FALSE,col.names=FALSE,quote=FALSE,append=T)