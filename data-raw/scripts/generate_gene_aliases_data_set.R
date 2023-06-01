library("tidyverse")
library("AbNames")
library("readxl")

existing <- ls()
downloads <- "../raw_data"

source("data-raw/01_hgnc.R")
source("data-raw/02_ncbi.R")

library("biomaRt")
library("org.Hs.eg.db")

source("data-raw/03_entrez_ensembl.R")

unloadNamespace("biomaRt")

source("data-raw/04_gene_aliases.R")
source("data-raw/05_protein_databases.R")

unloadNamespace("org.Hs.eg.db")

source("data-raw/06_manual_matches.R")
#source("data-raw/07_check_aliases.R")
