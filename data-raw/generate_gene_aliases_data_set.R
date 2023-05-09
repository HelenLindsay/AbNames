library("tidyverse")
library("AbNames")
library("readxl")

existing <- ls()

source("data-raw/01_hgnc.R")
source("data-raw/02_ncbi.R")

library("biomaRt")
library("org.Hs.eg.db")

source("data-raw/03_entrez_ensembl.R")

unloadNamespace("biomaRt")

#source("data-raw/04_gene_aliases.R")
#source("data-raw/04_protein_databases.R")

unloadNamespace("org.Hs.eg.db")

#source("data-raw/07_manual_matches.R")
#source("data-raw/08_antibody_registry.R")
#source("data-raw/09_check_aliases.R")
