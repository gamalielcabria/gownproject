# SingleM Summarise Visualization: Microeco Plot Alpha Diversity
# This is a script to analyse the results of SingleM summarise
# with options '--output-taxonomic-profile-with-extras'

########################### Libraries#############################
library(tidyverse)
library(microeco)

############################ Inputs ##############################
# Input files are from 'singleM_alphadiversity_01.R' 
# named as 'otu_table_<name>.csv' and 'taxa_table_<name>.csv'

# File location where the 'singleM_alphadiversity_01.R' outputs are
inputfolder <- '/home/glbcabria/Workbench/P2/singlem/summary'
setwd(inputfolder)

#This is <name> common to the otu_table and taxa_table
output <- 'singlem_GOWN2022' 

############################# Main ###############################
