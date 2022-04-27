# Laura Tibbs-Cortes
# Code to create consensus sequence for accession genotyped multiple times (vcf format).
# Based on Python script at https://github.com/mdzievit/Genomic_Prediction/tree/master/ThePlantGenome_Version/Consensus_Script

# load libraries
# If any libraries are not installed, first run: install.packages("libraryname")
library(tidyverse)
library(data.table)

# Set your working directory to the "code" folder using setwd()

# read in input data
vcf.in <- fread("../data/example.dup.NAM.recode.vcf", skip="#CHROM") #  vcf file with original genotype data; skip metadata lines
dup.lines <- fread("../data/example.keep.dups.txt", header=F) # text file with names of all duplicate genotypes found in the vcf file

n.snps <- nrow(vcf.in) # get number of SNPs

# Get "short.id" for each genotype, which is the accession name.
# In this case, this is the part before the colon, but this may be different in different files!! 
# Make sure to use appropriate separator character for YOUR data here.
colnames(dup.lines)[1] <- "full.id" # full.id is the ID for each distinct sequence of the accession
dup.lines <- dup.lines %>%
  separate(col="full.id", into=c("short.id"), sep=":", remove=F, extra="drop")

out.vcf <- vcf.in[,1:9] # initialize output
for(short.name in unique(dup.lines$short.id)) { 
  print(paste("Starting", short.name))
  print(Sys.time())
  # pull the column names for the current taxa
  comb.cols <- dup.lines %>% filter(short.id==short.name) %>% select(full.id)
  
  # pull the relevant vcf columns
  current.vcf <- vcf.in %>%
    select(comb.cols$full.id)
  
  # go through and get consensus for each snp, following same method as consensus python script
  out <- tibble(a=rep(NA,n.snps),
                b=rep(NA,n.snps),
                het1=rep(NA,n.snps),
                het2=rep(NA,n.snps)) # make tibble to hold results for this taxa
  for(i in 1:nrow(current.vcf)) {

    current.snp <- unname(unlist(current.vcf[i,]))
    current.snp <- gsub(pattern="\\./\\.", replacement="n/n", current.snp) #./. gets treated as wildcard so replace with n/n
    
    out$a[i] <- sum(str_count(fixed("1/1"), current.snp))
    out$b[i] <- sum(str_count(fixed("0/0"), current.snp))
    out$het1[i] <- sum(str_count(fixed("1/0"), current.snp))
    out$het2[i] <- sum(str_count(fixed("0/1"), current.snp))
  }
  out <- out %>%
    mutate(total=a+b+het1+het2)
  out <- out %>%
    mutate(consensus=ifelse(total==0, "./.",
                            ifelse(a/total > 0.5, "1/1",
                                   ifelse(b/total > 0.5, "0/0",
                                          ifelse(het1/total >0.5, "1/0",
                                                 ifelse(het2/total >0.5, "0/1",
                                                        "./."))))))
  fwrite(out %>% select(consensus), paste0(short.name, "_prevcf.txt"), sep="\t")
  
  # join with output
  out.vcf <- cbind(out.vcf, out$consensus)
  colnames(out.vcf)[ncol(out.vcf)] <- short.name
}
fwrite(out.vcf, "consensus_vcf.txt", sep="\t")
