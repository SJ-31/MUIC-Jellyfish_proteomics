### This script can be used to get a msms.txt file from MaxQuant and transform it intoa msms.tab that can be used as 
# an input for Percolator
# Miguel Cosenza 
# v 0.1

## Packages ----

library(dplyr)
library(tidyverse)
library(here)

## Load msms.txt file ----
library("optparse")
parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input"))
parser <- add_option(parser, c("-o", "--output"))
args <- parse_args(parser)

msmsfile <- read_tsv(here(args$input),
                     na = c("NaN", NA, "")) %>% 
          janitor::clean_names()

# core format for percolator ----

pinfrommsms <- mutate(
          msmsfile, # original file 
          
          SpecId = paste(raw_file,
                         scan_number,
                         sequence,
                         charge,
                         scan_event_number,
                         sep = "-"),
          
          Label = ifelse(is.na(reverse),
                         yes = 1,
                         no = -1),
          
          ScanNr = scan_number,
          
          ExpMass = 1000,
          
          Mass = mass,
          
          m_z = m_z,
          
          rt = retention_time, 
          
          dm = ifelse(is.na(mass_error_ppm),
                    no = mass_error_ppm,
                    yes = 0),
          
          deltaM_da = ifelse(is.na(mass_error_da),
                             no = mass_error_da,
                             yes = 0),
          
          simpleDeltaM_ppm = simple_mass_error_ppm,
          
          missedCleavages = missed_cleavages,
          
          sequence_length = length,
          
          andromeda_score = score,
          
          delta_score = delta_score,
          
          precursor_intensity = precursor_intensity,
          
          number_of_matches = number_of_matches,
          
          intensity_coverage = intensity_coverage,
          
          peak_coverage = peak_coverage,
          
          Peptide = paste("_.",
                          sequence,
                          "._",
                          sep = ""),
          
          Protein = ifelse(is.na(proteins),
                           yes = paste("REV__",
                                       sequence,
                                       "_",
                                       sep = ""),
                           no = proteins)
          
) %>%
          
          dplyr::select(
                    
                    SpecId,
                    Label,
                    ScanNr, 
                    ExpMass,
                    Mass,
                    m_z,
                    rt,
                    dm,
                    deltaM_da,
                    simpleDeltaM_ppm,
                    missedCleavages,
                    sequence_length,
                    andromeda_score,
                    delta_score,
                    precursor_intensity,
                    number_of_matches,
                    intensity_coverage,
                    peak_coverage,
                    Peptide,
                    Protein,
                    
                    proteins,
                    sequence,
                    scan_number,
                    charge,
                    scan_event_number
                    
          ) 

## prep the charge variable for 'dummyfication' 

charge_var <- dplyr::select(msmsfile, 
                            proteins,
                            sequence,
                            scan_number,
                            scan_event_number,
                            charge) %>%
          mutate(Charge = factor(charge,
                                 levels = unique(charge) 
                                 %>% sort()))

charge_var2 <- model.matrix(~ 0 + Charge, data = charge_var)%>%
          as.data.frame() %>%
          bind_cols(charge_var) %>%
          select(charge, everything(), -Charge)

## final percolator-formated table  

pinfrommsms_tab <- left_join(pinfrommsms, charge_var2) %>% 
          
          dplyr::select(-c(
                    proteins,
                    sequence,
                    scan_number,
                    charge,
                    scan_event_number
          )) %>% 
          
          relocate(Peptide, Protein,
                   .after = last_col()) %>%
          na.omit()


# generate output file ----

write_delim(pinfrommsms_tab,
            args$output,
            delim = "\t")

#Note: PSMs with NAs in the 'precursor intensity' columns are excluded. NAs in this column make percolator crash

