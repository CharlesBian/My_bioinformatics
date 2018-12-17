library(openxlsx)
library(readr)
library(survival)
library(survminer)
library(dplyr)

setwd(choose.dir())
getwd()

all_normal_countdata <- read.csv(choose.files())
all_tumor_countdata <- read.csv(choose.files())
pair_normal_countdata <- read.csv(choose.files())
pair_tumor_countdata <- read.csv(choose.files())

pheno <- read.csv(choose.files(),sep = '\t', header = TRUE)
name <- colnames(pheno[-1])


































