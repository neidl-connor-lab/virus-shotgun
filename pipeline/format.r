#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
options(stringsAsFactors=FALSE)

# helper functions
err <- function(message) {
  print(message)
  q(save="no", status=1)
}

# process arguments
args <- ArgumentParser()
args$add_argument("-o", "--ofile", required=TRUE, help="output CSV file")
args$add_argument("-v", "--vcf", required=TRUE, help="LoFreq VCF output")
args <- args$parse_args()

## check input -----------------------------------------------------------------
if(!file.exists(args$vcf)) {
  err(paste("VCF file not found:", args$vcf))
}

## load, format, and save ------------------------------------------------------
read.delim(args$vcf,
           comment.char="#",
           header=FALSE,
           col.names=c("Accession", "NT.position", "NT.ID", "NT.ref", "NT.alt", 
                       "QUAL", "FILTER", "INFO")) %>%
  # only keep passing SNVs
  filter(FILTER=="PASS") %>%
  # format SNV IDs and extract frequencies
  # pull and format depth data
  mutate(NT.ID=paste0(NT.position, "-", NT.ref, "-", NT.alt),
         Frequency=as.numeric(str_extract(INFO, "(?<=AF=)[0-9\\.]+")),
         InDel=str_detect(INFO, "INDEL")) %>%
  select(Accession, NT.ID, NT.position, Frequency, InDel, NT.ref, NT.alt) %>%
  # save
  write.csv(args$ofile, na="", row.names=FALSE)
