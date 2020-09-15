## ---------------------------
##
## Script name: MRMStaR/MRMSymR
##
## Purpose of script: Calculate symmetry of the peak 
##
## Author: Merve Alp
##
## Date Created: 2020-09-14
##
## Email: kezibanmerve.alp@mdc-berlin.de
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory 

setwd(choose.dir()) 

## Install and Import required packages

if   (!require("pacman")) install.packages("pacman")
pacman::p_load("dplyr",
               "tidyr",
               "reshape2",
               "stringr",
               "broom")

library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(broom)

## Import input file and modify column names

data              <- read.csv(file= "Input/Symmetry_Skyline_Export.csv", header=TRUE, sep=",", stringsAsFactors = F, as.is = T, dec = ".")
data$Protein.Name <- gsub(pattern = "^.*\\|",  replacement = "",  data$Protein.Name)###should be modified depending on how protein name look like,removes Uniprot_IDs


data[c("Retention.Time", "Start.Time", "End.Time")] <- sapply(data[c("Retention.Time", 
                                                                     "Start.Time", 
                                                                     "End.Time")], as.numeric)

## Calculate head and tail of a peak and divide sum of them by head*2. Return this value. If a peak is symmetric this value should be one

detect_sym  <- function(start, end, apex ) {
        head <- abs(apex - start)
        tail <- abs(apex - end)
        if(head != 0){
                sym <- (head + tail) / (2* head)  
                return(sym)
        }
        else 
                return(NA)
}

## Define sd of population as function, since I am not infering from sample to population here

sd.p = function(x){sd(x)*sqrt((length(x)-1)/length(x))} 

## Filter data

data       <- data %>%
              filter(Truncated == "False")   %>% 
              filter(Quantitative == "True") %>%
              filter(Standard.Type != "iRT") %>% 
              select(-c(Truncated, Quantitative, Standard.Type)) 


## Calculate symmetry values for all precursors

peakSym_table      <- data %>% 
                      group_by(Protein.Name, 
                               Peptide.Modified.Sequence, 
                               Fragment.Ion, 
                               Isotope.Label.Type, 
                               Replicate.Name) %>%
                      rowwise()%>%
                      summarise(peak.Sym = detect_sym(Start.Time, End.Time, Retention.Time))


## Calculate average symmetry score and std across transitions from same precursor


peakSym_table     <- peakSym_table %>% 
                     group_by(Protein.Name, 
                              Peptide.Modified.Sequence, 
                              Isotope.Label.Type, 
                              Replicate.Name) %>%
                     mutate(Sym.mean = mean(peak.Sym, na.rm = FALSE),
                            Sym.std  = sd.p(peak.Sym))


## Calculate CVs

peakSym_table    <- peakSym_table %>%
                    mutate(CV = Sym.std / Sym.mean)

## Export results in txt file

write.table( peakSym_table, file = 'peakSym_table.txt',
             quote = F, sep = '\t', dec = '.', row.names = F )


## Plot histograms

pdf("Output/Histograms.pdf")

hist(peakSym_table$peak.Sym, breaks=30, xlim=c(0,7), col=rgb(1,0,0,0.5), xlab="Symmetry Score", 
     ylab="Counts", main="Symmetry Value Distribution" )

hist(peakSym_table$CV, breaks=30, xlim=c(0,1), col=rgb(1,0,0,0.5), xlab="CV", 
     ylab="Counts", main="CV Value Distribution" )

dev.off()




