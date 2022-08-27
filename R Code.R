---
title: "Irene Hsueh's BS 858 Homework 9"
author: "Irene Hsueh"
date: "11/22/2021"
output:
  pdf_document: default
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
```

# Analyzing VCF File
```{r}
#Reading in Text File
vcf <- read.table("C:/Irene Hsueh's Documents/MS Applied Biostatistics/BS 858 - Statistical Genetics I/Class 10 - Sequencing Data/Homework/mendelian_2021.txt", header=TRUE, as.is=TRUE) %>% 
#Selecting and Renaming Variables
  dplyr::select(id = ID, 
                position = POS, 
                info = INFO, 
                ref = REF, 
                alt = ALT, 
                individual1 = IND01,
                individual2 = IND02, 
                individual3 = IND03,
                individual4 = IND04, 
                individual5 = IND05) 
head(vcf, 20)
```



```{r}
#Number of Alt Alleles
vcf$n_alt_alleles <- str_count(vcf$alt, ",") + 1
table(vcf$n_alt_alleles)
sum(vcf$n_alt_alleles > 1)


#Number of SNPs
vcf<- vcf %>%
  mutate(snp = ifelse(ref %in% c("A","C","G","T") & alt %in% c("A","C","G","T"), 
                      yes=TRUE, no=FALSE))
sum(vcf$snp==TRUE)


#Number of Alternate Alleles among The Five Individuals
alt_count <- function(x) {
	y <- ifelse(x == "1/1", yes=2, no=NA)
	y <- ifelse(x %in% c("0/1", "1/0"), yes=1, no=y)
	y <- ifelse(x %in% c("0/0"), yes=0, no=y)
	sum(y, na.rm=TRUE)
}

vcf$total_alt_alleles <- apply(vcf[,paste("individual",seq(1,5,1),sep="")], 1, alt_count)


#SNPs and Alternate Alleles 
sum(vcf$snp==TRUE & vcf$total_alt_alleles >= 1)


#Missing Rates
individual1_missing <- sum(vcf$individual1 =="./.")/nrow(vcf)*100
individual2_missing <- sum(vcf$individual2 =="./.")/nrow(vcf)*100
individual3_missing <- sum(vcf$individual3 =="./.")/nrow(vcf)*100
individual4_missing <- sum(vcf$individual4 =="./.")/nrow(vcf)*100
individual5_missing <- sum(vcf$individual5 =="./.")/nrow(vcf)*100


#Variant Sites
individual1_variant_sites <- sum(vcf$individual1 %in% c("1/1", "0/1", "1/0"))
individual2_variant_sites <- sum(vcf$individual2 %in% c("1/1", "0/1", "1/0"))
individual3_variant_sites <- sum(vcf$individual3 %in% c("1/1", "0/1", "1/0"))
individual4_variant_sites <- sum(vcf$individual4 %in% c("1/1", "0/1", "1/0"))
individual5_variant_sites <- sum(vcf$individual5 %in% c("1/1", "0/1", "1/0"))


#Singletons 
individual1_singletons <- sum(vcf$individual1 %in% c("0/1", "1/0") & vcf$total_alt_alleles==1)
individual2_singletons <- sum(vcf$individual2 %in% c("0/1", "1/0") & vcf$total_alt_alleles==1)
individual3_singletons <- sum(vcf$individual3 %in% c("0/1", "1/0") & vcf$total_alt_alleles==1)
individual4_singletons <- sum(vcf$individual4 %in% c("0/1", "1/0") & vcf$total_alt_alleles==1)
individual5_singletons <- sum(vcf$individual5 %in% c("0/1", "1/0") & vcf$total_alt_alleles==1)
```



### Transition-Tranversion Ratio
```{r}
#Individual 1 Ti/Tv Ratio
individual1_ti <- sum(vcf$individual1 %in% c("1/1", "0/1", "1/0") &
                          (vcf$ref=="G" & vcf$alt=="A" | 
                           vcf$ref=="A" & vcf$alt=="G" | 
                           vcf$ref=="T" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="T")) 
individual1_tv <- sum(vcf$individual1 %in% c("1/1", "0/1", "1/0") &
                          (vcf$ref=="A" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="A" | 
                           vcf$ref=="G" & vcf$alt=="T" | 
                           vcf$ref=="T" & vcf$alt=="G" |
                           vcf$ref=="A" & vcf$alt=="T" | 
                           vcf$ref=="T" & vcf$alt=="A" | 
                           vcf$ref=="G" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="G" )) 
individual1_ti_tv <- individual1_ti/individual1_tv


#Individual 2 Ti/Tv Ratio
individual2_ti <- sum(vcf$individual2 %in% c("1/1", "0/1", "1/0") &
                          (vcf$ref=="G" & vcf$alt=="A" | 
                           vcf$ref=="A" & vcf$alt=="G" | 
                           vcf$ref=="T" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="T")) 
individual2_tv <- sum(vcf$individual2 %in% c("1/1", "0/1", "1/0") &
                          (vcf$ref=="A" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="A" | 
                           vcf$ref=="G" & vcf$alt=="T" | 
                           vcf$ref=="T" & vcf$alt=="G" |
                           vcf$ref=="A" & vcf$alt=="T" | 
                           vcf$ref=="T" & vcf$alt=="A" | 
                           vcf$ref=="G" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="G" )) 
individual2_ti_tv <- individual2_ti/individual2_tv


#Individual 3 Ti/Tv Ratio
individual3_ti <- sum(vcf$individual3 %in% c("1/1", "0/1", "1/0") &
                          (vcf$ref=="G" & vcf$alt=="A" | 
                           vcf$ref=="A" & vcf$alt=="G" | 
                           vcf$ref=="T" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="T")) 
individual3_tv <- sum(vcf$individual3 %in% c("1/1", "0/1", "1/0") &
                          (vcf$ref=="A" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="A" | 
                           vcf$ref=="G" & vcf$alt=="T" | 
                           vcf$ref=="T" & vcf$alt=="G" |
                           vcf$ref=="A" & vcf$alt=="T" | 
                           vcf$ref=="T" & vcf$alt=="A" | 
                           vcf$ref=="G" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="G" )) 
individual3_ti_tv <- individual3_ti/individual3_tv


#Individual 4 Ti/Tv Ratio
individual4_ti <- sum(vcf$individual4 %in% c("1/1", "0/1", "1/0") &
                          (vcf$ref=="G" & vcf$alt=="A" | 
                           vcf$ref=="A" & vcf$alt=="G" | 
                           vcf$ref=="T" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="T")) 
individual4_tv <- sum(vcf$individual4 %in% c("1/1", "0/1", "1/0") &
                          (vcf$ref=="A" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="A" | 
                           vcf$ref=="G" & vcf$alt=="T" | 
                           vcf$ref=="T" & vcf$alt=="G" |
                           vcf$ref=="A" & vcf$alt=="T" | 
                           vcf$ref=="T" & vcf$alt=="A" | 
                           vcf$ref=="G" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="G" )) 
individual4_ti_tv <- individual4_ti/individual4_tv


#Individual 5 Ti/Tv Ratio
individual5_ti <- sum(vcf$individual5 %in% c("1/1", "0/1", "1/0") &
                          (vcf$ref=="G" & vcf$alt=="A" | 
                           vcf$ref=="A" & vcf$alt=="G" | 
                           vcf$ref=="T" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="T")) 
individual5_tv <- sum(vcf$individual5 %in% c("1/1", "0/1", "1/0") &
                          (vcf$ref=="A" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="A" | 
                           vcf$ref=="G" & vcf$alt=="T" | 
                           vcf$ref=="T" & vcf$alt=="G" |
                           vcf$ref=="A" & vcf$alt=="T" | 
                           vcf$ref=="T" & vcf$alt=="A" | 
                           vcf$ref=="G" & vcf$alt=="C" | 
                           vcf$ref=="C" & vcf$alt=="G" )) 
individual5_ti_tv <- individual5_ti/individual5_tv
```



### Associated Variants 
```{r}
associated_variants <- vcf %>% 
  filter(individual1=="0/0" & 
           individual5=="0/0" & 
           individual3=="0/1" &
           individual2=="0/1" &
           individual4=="0/1"
         )
```








