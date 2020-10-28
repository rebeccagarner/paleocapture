# 06_taxonomy.R
# Taxonomic analyses of free and captured metagenomes

setwd("paleocapture/")

# Load libraries
library(tidyverse)
library(vegan)
library(FactoMineR)
library(factoextra)
library(ggdendro)

#### Import IMG .map.txt files which relate MEGAHIT contig ID and IMG scaffold ID ####
map_06126_bottom <- read_tsv(file = "data/img/lp2017_06-126-bottom/3300032883.a.map.txt", col_names = c("contig_id", "scaffold_id"))
map_06126_top <- read_tsv(file = "data/img/lp2017_06-126-top/3300032795.a.map.txt", col_names = c("contig_id", "scaffold_id"))
map_06126_watercolumn <- read_tsv(file = "data/img/lp2017_06-126-watercolumn/3300032269.a.map.txt", col_names = c("contig_id", "scaffold_id"))

map_17054_bottom <- read_tsv(file = "data/img/lp2017_17-054-bottom/3300032797.a.map.txt", col_names = c("contig_id", "scaffold_id"))
map_17054_top <- read_tsv(file = "data/img/lp2017_17-054-top/3300032796.a.map.txt", col_names = c("contig_id", "scaffold_id"))
map_17054_watercolumn <- read_tsv(file = "data/img/lp2017_17-054-watercolumn/3300032328.a.map.txt", col_names = c("contig_id", "scaffold_id"))

map_17067_bottom <- read_tsv(file = "data/img/lp2017_17-067-bottom/3300032799.a.map.txt", col_names = c("contig_id", "scaffold_id"))
map_17067_top <- read_tsv(file = "data/img/lp2017_17-067-top/3300032798.a.map.txt", col_names = c("contig_id", "scaffold_id"))
map_17067_watercolumn <- read_tsv(file = "data/img/lp2017_17-067-watercolumn/3300032317.a.map.txt", col_names = c("contig_id", "scaffold_id"))


#### Extract IMG scaffold ID string lengths ####
map <- bind_rows(map_06126_bottom, map_06126_top, map_06126_watercolumn,
                 map_17054_bottom, map_17054_top, map_17054_watercolumn,
                 map_17067_bottom, map_17067_top, map_17067_watercolumn)

img_scaffold_ids <- map %>%
  select(-contig_id) %>%
  mutate(project_id = str_remove(scaffold_id, pattern = "_.*")) %>%
  distinct(project_id, .keep_all = TRUE) %>%
  mutate(scaffold_id_length = str_length(scaffold_id)) %>%
  select(-scaffold_id)


#### Import covstats of reads recruited at 90% identity and filter contigs with non-zero coverage ####
covstats_06126_swa_tsr <- read_tsv(file = "data/bbmap/covstats/lp2017_06-126-swa_tsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)
covstats_06126_swa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_06-126-swa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)
covstats_06126_tsa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_06-126-tsa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)

covstats_17054_swa_tsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-054-swa_tsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)
covstats_17054_swa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-054-swa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)
covstats_17054_tsa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-054-tsa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)

covstats_17067_swa_tsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-067-swa_tsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)
covstats_17067_swa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-067-swa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)
covstats_17067_tsa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-067-tsa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)

covstats_06126_swa <- read_tsv(file = "data/bbmap/covstats/lp2017_06-126-swa_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)
covstats_06126_tsa <- read_tsv(file = "data/bbmap/covstats/lp2017_06-126-tsa_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)
covstats_06126_bsa <- read_tsv(file = "data/bbmap/covstats/lp2017_06-126-bsa_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)

covstats_17054_swa <- read_tsv(file = "data/bbmap/covstats/lp2017_17-054-swa_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)
covstats_17054_tsa <- read_tsv(file = "data/bbmap/covstats/lp2017_17-054-tsa_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)
covstats_17054_bsa <- read_tsv(file = "data/bbmap/covstats/lp2017_17-054-bsa_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)

covstats_17067_swa <- read_tsv(file = "data/bbmap/covstats/lp2017_17-067-swa_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)
covstats_17067_tsa <- read_tsv(file = "data/bbmap/covstats/lp2017_17-067-tsa_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)
covstats_17067_bsa <- read_tsv(file = "data/bbmap/covstats/lp2017_17-067-bsa_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0)


#### Pull out scaffolds of interest with coverage information ####
scaffolds_06126_swa_tsr <- covstats_06126_swa_tsr %>% left_join(map_06126_watercolumn, by = "contig_id") %>% select(scaffold_id, Avg_fold)
scaffolds_06126_swa_bsr <- covstats_06126_swa_bsr %>% left_join(map_06126_watercolumn, by = "contig_id") %>% select(scaffold_id, Avg_fold)
scaffolds_06126_tsa_bsr <- covstats_06126_tsa_bsr %>% left_join(map_06126_top, by = "contig_id") %>% select(scaffold_id, Avg_fold)

scaffolds_17054_swa_tsr <- covstats_17054_swa_tsr %>% left_join(map_17054_watercolumn, by = "contig_id") %>% select(scaffold_id, Avg_fold)
scaffolds_17054_swa_bsr <- covstats_17054_swa_bsr %>% left_join(map_17054_watercolumn, by = "contig_id") %>% select(scaffold_id, Avg_fold)
scaffolds_17054_tsa_bsr <- covstats_17054_tsa_bsr %>% left_join(map_17054_top, by = "contig_id") %>% select(scaffold_id, Avg_fold)

scaffolds_17067_swa_tsr <- covstats_17067_swa_tsr %>% left_join(map_17067_watercolumn, by = "contig_id") %>% select(scaffold_id, Avg_fold)
scaffolds_17067_swa_bsr <- covstats_17067_swa_bsr %>% left_join(map_17067_watercolumn, by = "contig_id") %>% select(scaffold_id, Avg_fold)
scaffolds_17067_tsa_bsr <- covstats_17067_tsa_bsr %>% left_join(map_17067_top, by = "contig_id") %>% select(scaffold_id, Avg_fold)

scaffolds_06126_swa <- covstats_06126_swa %>% left_join(map_06126_watercolumn, by = "contig_id") %>% select(scaffold_id, Avg_fold)
scaffolds_06126_tsa <- covstats_06126_tsa %>% left_join(map_06126_top, by = "contig_id") %>% select(scaffold_id, Avg_fold)
scaffolds_06126_bsa <- covstats_06126_bsa %>% left_join(map_06126_bottom, by = "contig_id") %>% select(scaffold_id, Avg_fold)

scaffolds_17054_swa <- covstats_17054_swa %>% left_join(map_17054_watercolumn, by = "contig_id") %>% select(scaffold_id, Avg_fold)
scaffolds_17054_tsa <- covstats_17054_tsa %>% left_join(map_17054_top, by = "contig_id") %>% select(scaffold_id, Avg_fold)
scaffolds_17054_bsa <- covstats_17054_bsa %>% left_join(map_17054_bottom, by = "contig_id") %>% select(scaffold_id, Avg_fold)

scaffolds_17067_swa <- covstats_17067_swa %>% left_join(map_17067_watercolumn, by = "contig_id") %>% select(scaffold_id, Avg_fold)
scaffolds_17067_tsa <- covstats_17067_tsa %>% left_join(map_17067_top, by = "contig_id") %>% select(scaffold_id, Avg_fold)
scaffolds_17067_bsa <- covstats_17067_bsa %>% left_join(map_17067_bottom, by = "contig_id") %>% select(scaffold_id, Avg_fold)


#### Import IMG .phylodist.txt files with taxonomic annotations ####
phylodist_06126_bottom <- read_tsv(file = "data/img/lp2017_06-126-bottom/3300032883.a.phylodist.txt", col_names = c("gene_id", "homolog_gene_oid", "homolog_taxon_oid", "percent_identity", "lineage")) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length), gene_count = 1) %>% select(scaffold_id, gene_id, lineage, gene_count)
phylodist_06126_top <- read_tsv(file = "data/img/lp2017_06-126-top/3300032795.a.phylodist.txt", col_names = c("gene_id", "homolog_gene_oid", "homolog_taxon_oid", "percent_identity", "lineage")) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length), gene_count = 1) %>% select(scaffold_id, gene_id, lineage, gene_count)
phylodist_06126_watercolumn <- read_tsv(file = "data/img/lp2017_06-126-watercolumn/3300032269.a.phylodist.txt", col_names = c("gene_id", "homolog_gene_oid", "homolog_taxon_oid", "percent_identity", "lineage")) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length), gene_count = 1) %>% select(scaffold_id, gene_id, lineage, gene_count)

phylodist_17054_bottom <- read_tsv(file = "data/img/lp2017_17-054-bottom/3300032797.a.phylodist.txt", col_names = c("gene_id", "homolog_gene_oid", "homolog_taxon_oid", "percent_identity", "lineage")) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length), gene_count = 1) %>% select(scaffold_id, gene_id, lineage, gene_count)
phylodist_17054_top <- read_tsv(file = "data/img/lp2017_17-054-top/3300032796.a.phylodist.txt", col_names = c("gene_id", "homolog_gene_oid", "homolog_taxon_oid", "percent_identity", "lineage")) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length), gene_count = 1) %>% select(scaffold_id, gene_id, lineage, gene_count)
phylodist_17054_watercolumn <- read_tsv(file = "data/img/lp2017_17-054-watercolumn/3300032328.a.phylodist.txt", col_names = c("gene_id", "homolog_gene_oid", "homolog_taxon_oid", "percent_identity", "lineage")) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length), gene_count = 1) %>% select(scaffold_id, gene_id, lineage, gene_count)

phylodist_17067_bottom <- read_tsv(file = "data/img/lp2017_17-067-bottom/3300032799.a.phylodist.txt", col_names = c("gene_id", "homolog_gene_oid", "homolog_taxon_oid", "percent_identity", "lineage")) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length), gene_count = 1) %>% select(scaffold_id, gene_id, lineage, gene_count)
phylodist_17067_top <- read_tsv(file = "data/img/lp2017_17-067-top/3300032798.a.phylodist.txt", col_names = c("gene_id", "homolog_gene_oid", "homolog_taxon_oid", "percent_identity", "lineage")) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length), gene_count = 1) %>% select(scaffold_id, gene_id, lineage, gene_count)
phylodist_17067_watercolumn <- read_tsv(file = "data/img/lp2017_17-067-watercolumn/3300032317.a.phylodist.txt", col_names = c("gene_id", "homolog_gene_oid", "homolog_taxon_oid", "percent_identity", "lineage")) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length), gene_count = 1) %>% select(scaffold_id, gene_id, lineage, gene_count)


#### phylodist annotations for scaffolds with non-zero coverage ####
# Refine phylodist annotations by separating by taxonomic rank
# Divide Proteobacteria into phyla
phylodist_06126_swa_tsr <- scaffolds_06126_swa_tsr %>% left_join(phylodist_06126_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Lac Paula", metagenome = "swa_tsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_06126_swa_bsr <- scaffolds_06126_swa_bsr %>% left_join(phylodist_06126_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Lac Paula", metagenome = "swa_bsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_06126_tsa_bsr <- scaffolds_06126_tsa_bsr %>% left_join(phylodist_06126_top, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Lac Paula", metagenome = "tsa_bsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))

phylodist_17054_swa_tsr <- scaffolds_17054_swa_tsr %>% left_join(phylodist_17054_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Eightmile Lake", metagenome = "swa_tsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17054_swa_bsr <- scaffolds_17054_swa_bsr %>% left_join(phylodist_17054_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Eightmile Lake", metagenome = "swa_bsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17054_tsa_bsr <- scaffolds_17054_tsa_bsr %>% left_join(phylodist_17054_top, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Eightmile Lake", metagenome = "tsa_bsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))

phylodist_17067_swa_tsr <- scaffolds_17067_swa_tsr %>% left_join(phylodist_17067_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Grand lac Touradi", metagenome = "swa_tsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17067_swa_bsr <- scaffolds_17067_swa_bsr %>% left_join(phylodist_17067_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Grand lac Touradi", metagenome = "swa_bsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17067_tsa_bsr <- scaffolds_17067_tsa_bsr %>% left_join(phylodist_17067_top, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Grand lac Touradi", metagenome = "tsa_bsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))

phylodist_06126_swa <- scaffolds_06126_swa %>% left_join(phylodist_06126_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Lac Paula", metagenome = "swa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_06126_tsa <- scaffolds_06126_tsa %>% left_join(phylodist_06126_top, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Lac Paula", metagenome = "tsa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_06126_bsa <- scaffolds_06126_bsa %>% left_join(phylodist_06126_bottom, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Lac Paula", metagenome = "bsa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))

phylodist_17054_swa <- scaffolds_17054_swa %>% left_join(phylodist_17054_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Eightmile Lake", metagenome = "swa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17054_tsa <- scaffolds_17054_tsa %>% left_join(phylodist_17054_top, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Eightmile Lake", metagenome = "tsa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17054_bsa <- scaffolds_17054_bsa %>% left_join(phylodist_17054_bottom, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Eightmile Lake", metagenome = "bsa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))

phylodist_17067_swa <- scaffolds_17067_swa %>% left_join(phylodist_17067_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Grand lac Touradi", metagenome = "swa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17067_tsa <- scaffolds_17067_tsa %>% left_join(phylodist_17067_top, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Grand lac Touradi", metagenome = "tsa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17067_bsa <- scaffolds_17067_bsa %>% left_join(phylodist_17067_bottom, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Grand lac Touradi", metagenome = "bsa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))


#### Combine phylodist annotations for all metagenomes ####
# Remove scaffolds containing rRNA or tRNA genes
scaffolds_rna <- read_tsv("output/scaffolds_rna.tsv") %>%
  pull(scaffold_id)

phylodist_all <- bind_rows(phylodist_06126_swa_tsr,
                           phylodist_06126_swa_bsr,
                           phylodist_06126_tsa_bsr,
                           
                           phylodist_17054_swa_tsr,
                           phylodist_17054_swa_bsr,
                           phylodist_17054_tsa_bsr,
                           
                           phylodist_17067_swa_tsr,
                           phylodist_17067_swa_bsr,
                           phylodist_17067_tsa_bsr,
                           
                           phylodist_06126_swa,
                           phylodist_06126_tsa,
                           phylodist_06126_bsa,
                           
                           phylodist_17054_swa,
                           phylodist_17054_tsa,
                           phylodist_17054_bsa,
                           
                           phylodist_17067_swa,
                           phylodist_17067_tsa,
                           phylodist_17067_bsa) %>%
  filter(!scaffold_id %in% scaffolds_rna)

# Clear up space in R environment
rm(list = ls(pattern = "map"))
rm(list = ls(pattern = "covstats"))
rm(list = ls(pattern = "phylodist_\\d"))
rm(list = ls(pattern = "scaffold"))


#### Palettes ####
palette_shapes <- c(87, 84, 66, 25, 25, 8)
names(palette_shapes) <- c("swa", "tsa", "bsa", "swa_tsr", "swa_bsr", "tsa_bsr")

palette_lakes <- c("#ff9f1a", "#18dcff", "#7d5fff")
names(palette_lakes) <- c("Lac Paula", "Eightmile Lake", "Grand lac Touradi")


#### PCA of order-rank metagenome taxonomic compositions ####
# Site by order matrix in data frame format
phylodist_df_order <- phylodist_all %>%
  mutate(lake_metagenome = str_c(lake, metagenome, sep = ".")) %>%
  group_by(lake_metagenome, taxlevel_4) %>%
  summarize(sum_avgfold = sum(Avg_fold)) %>%
  filter(!grepl(pattern = "unclassified", x = taxlevel_4, ignore.case = TRUE)) %>%
  spread(taxlevel_4, sum_avgfold, fill = 0) %>%
  as.data.frame()
rownames(phylodist_df_order) <- phylodist_df_order$lake_metagenome
phylodist_df_order$lake_metagenome <- NULL

# Data transformation
phylodist_chord_order <- decostand(x = phylodist_df_order, method = "normalize")

# PCA on order-rank taxonomic profiles
phylodist_pca_order <- rda(phylodist_chord_order, scale = FALSE)
prcomp_sites_order <- scores(phylodist_pca_order)$sites
prcomp_species_order <- scores(phylodist_pca_order)$species
pc1_percent_var_order <- round(x = (eigenvals(x = phylodist_pca_order)[1])/(sum(eigenvals(x = phylodist_pca_order)))*100, digits = 1)
pc2_percent_var_order <- round(x = (eigenvals(x = phylodist_pca_order)[2])/(sum(eigenvals(x = phylodist_pca_order)))*100, digits = 1)

# Taxon variables contributing to PCA dimensions 1 and 2
phylodist_pca_res_order <- PCA(X = phylodist_chord_order, scale.unit = FALSE, ncp = 17, graph = FALSE)
fviz_contrib(phylodist_pca_res_order, choice = "var", axes = 1, top = 25)
fviz_contrib(phylodist_pca_res_order, choice = "var", axes = 2, top = 10)

tax_contributors <- phylodist_pca_res_order$var$contrib %>%
  as.data.frame() %>%
  rownames_to_column("taxon") %>%
  as_tibble() %>%
  filter(Dim.1 > 1.6 | Dim.2 > 0.39) %>%
  pull(taxon)

taxa_contribution_order <- prcomp_species_order %>%
  as_tibble(rownames = "tax_group") %>%
  filter(tax_group %in% tax_contributors)

phylodist_prcomp_order <- as.data.frame(prcomp_sites_order) %>%
  rownames_to_column("lake_metagenome") %>%
  as_tibble() %>%
  mutate(lake = str_replace(lake_metagenome, "\\..*", replacement = ""),
         metagenome = str_replace(lake_metagenome, ".*\\.", replacement = "")) %>%
  mutate(fill_shape = case_when(metagenome == "swa" ~ 1,
                                metagenome == "tsa" ~ 1,
                                metagenome == "bsa" ~ 1,
                                metagenome == "swa_tsr" ~ 1,
                                metagenome == "swa_bsr" ~ 0,
                                metagenome == "tsa_bsr" ~ 1))

(phylodist_pca_order_plot <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "lightgrey") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "lightgrey") +
    geom_point(aes(x = PC1, y = PC2, colour = lake, fill = ifelse(fill_shape == 1, lake, NA),
                   shape = factor(metagenome, levels = c("swa", "tsa", "bsa", "swa_tsr", "swa_bsr", "tsa_bsr"))),
               data = phylodist_prcomp_order, size = 5, alpha = 1, stroke = 1.5) +
    scale_shape_manual(values = palette_shapes, guide = "none") +
    scale_fill_manual(values = palette_lakes, guide = "none") +
    scale_colour_manual(values = palette_lakes, guide = "none") +
    geom_segment(aes(x = 0, y = 0, xend = 1*PC1, yend = 1*PC2),
                 data = taxa_contribution_order,
                 arrow = arrow(length = unit(x = 0.2, units = "cm")),
                 alpha = 0.7, colour = "gray60") +
    geom_text(aes(x = 1.1*PC1, y = 1.1*PC2), data = taxa_contribution_order,
              label = taxa_contribution_order$tax_group, colour = "gray60", size = 4) +
    labs(x = paste("PC1 (", pc1_percent_var_order, "%)", sep = ""),
         y = paste("PC2 (", pc2_percent_var_order, "%)", sep = "")) +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          legend.position = "none",
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5)))
#ggsave("figures/04b_paleocapture_pca_order_cov_rmrna.pdf", plot = phylodist_pca_order_plot, device = "pdf", width = 10, height = 10, units = "in")


#### Generate taxonomy heatmaps ####
phylodist_all_grouped_viruses <- phylodist_all %>%
  mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1,
                               taxlevel_2 == "Proteobacteria" ~ taxlevel_3,
                               TRUE ~ taxlevel_2))

# Sum coverage for each metagenome
phylodist_sum_cov <- phylodist_all_grouped_viruses %>%
  group_by(lake, metagenome) %>%
  summarize(sum_avgfold = sum(Avg_fold))

# Scale relative coverage to 100% metagenome profile
phylodist_all_sum_cov <- phylodist_all_grouped_viruses %>%
  left_join(phylodist_sum_cov, by = c("lake", "metagenome")) %>%
  mutate(relative_coverage = (Avg_fold/sum_avgfold * 100))

# Identify phyla with >= 10 % coverage
phyla_10pct_coverage <- phylodist_all_sum_cov %>%
  group_by(lake, metagenome, tax_group) %>%
  summarize(sum_relative_coverage = sum(relative_coverage)) %>%
  filter(sum_relative_coverage >= 10) %>%
  pull(tax_group) %>%
  unique()

# Draw dendrogram of phylum-rank taxa based on Bray-Curtis dissimilarity of gene composition
site_by_phylum_coverage <- phylodist_all_sum_cov %>%
  filter(tax_group %in% phyla_10pct_coverage) %>%
  group_by(lake, metagenome, tax_group) %>%
  summarize(sum_relative_coverage = sum(relative_coverage)) %>%
  ungroup() %>%
  spread(tax_group, sum_relative_coverage, fill = 0) %>%
  mutate(lake_metagenome = str_c(lake, metagenome, sep = ".")) %>%
  select(-lake, -metagenome) %>%
  as.data.frame()
rownames(site_by_phylum_coverage) <- site_by_phylum_coverage$lake_metagenome
site_by_phylum_coverage$lake_metagenome <- NULL

site_by_phylum_bray <- vegdist(t(site_by_phylum_coverage), method = "bray")
site_by_phylum_hclust <- hclust(site_by_phylum_bray, method = "ward.D2")
site_by_phylum_dendro <- dendro_data(site_by_phylum_hclust)
site_by_phylum_cut <- cutree(site_by_phylum_hclust, k = 3)

site_by_phylum_cut_df <- data.frame(cluster = site_by_phylum_cut) %>%
  rownames_to_column("label") %>%
  as_tibble() %>%
  left_join(label(site_by_phylum_dendro), by = "label")

(dendrogram_phylodist_phyla <- ggplot() + 
    geom_segment(data = segment(site_by_phylum_dendro), aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_text(data = site_by_phylum_cut_df, aes(x, y, label = label, hjust = 0)) +
    coord_flip() +
    scale_y_reverse(expand = c(1,0)) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(), axis.title = element_blank(),
          panel.background = element_blank(), panel.grid = element_blank(),
          legend.position = "none"))
#ggsave(filename = "figures/03_paleocapture_dendrogram_phylum_rmrna.pdf", plot = dendrogram_phylodist_phyla, device = "pdf", units = "in", width = 5, height = 10)

# Generate phylum-rank heatmap based on relative coverage
(heatmap_phylodist_phyla <- phylodist_all_sum_cov %>%
    filter(tax_group %in% phyla_10pct_coverage) %>%
    group_by(lake, metagenome, tax_group) %>%
    summarize(sum_relative_coverage = sum(relative_coverage)) %>%
    ggplot(aes(x = factor(metagenome, levels = c("swa", "tsa", "bsa", "swa_tsr", "swa_bsr", "tsa_bsr")),
               y = factor(tax_group, levels = c("Chloroflexi", "Firmicutes", "Deltaproteobacteria", "Euryarchaeota", "Cyanobacteria", "Planctomycetes", "Bacteroidetes", "Gammaproteobacteria", "Viruses", "Actinobacteria", "Alphaproteobacteria", "Betaproteobacteria")))) +
    facet_grid(~factor(lake, levels = c("Lac Paula", "Eightmile Lake", "Grand lac Touradi"))) +
    geom_tile(aes(fill = sum_relative_coverage), colour = "white") +
    #geom_text(aes(label = round(sum_relative_coverage, 1)), size = 2) +
    scale_x_discrete(expand = c(0,0)) +
    labs(fill = "Coverage (%)") +
    scale_fill_gradient(low = "white", high = "#505050") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "bottom", axis.title = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = -90, hjust = 0), axis.ticks = element_blank()))
#ggsave("figures/03_paleocapture_heatmap_phylum_rmrna.pdf", plot = heatmap_phylodist_phyla, device = "pdf", units = "in", width = 10, height = 8)


# Summarize range of phylum-rank relative coverage in a CSV table
summary_phylodist_phylum_cov <- phylodist_all_sum_cov %>%
  filter(tax_group %in% phyla_10pct_coverage) %>%
  group_by(lake, metagenome, tax_group) %>%
  summarize(sum_relative_coverage = round(sum(relative_coverage), 1)) %>%
  ungroup() %>%
  mutate(metagenome_type = case_when(!grepl("_", metagenome) ~ "free",
                                     grepl("_", metagenome) ~ "captured")) %>%
  spread(lake, sum_relative_coverage) %>%
  select(metagenome_type, metagenome, tax_group, `Lac Paula`, `Eightmile Lake`, `Grand lac Touradi`) %>%
  arrange(factor(metagenome, levels = c("swa", "tsa", "bsa", "swa_tsr", "swa_bsr", "tsa_bsr")))
# summary_phylodist_phylum_cov %>%
#   write_csv("output/summary_phyla_cov_rmrna.csv")


#### Identify orders with >= 5% coverage ####
orders_5pct_coverage <- phylodist_all_sum_cov %>%
  group_by(lake, metagenome, taxlevel_4) %>%
  summarize(sum_relative_coverage = sum(relative_coverage)) %>%
  filter(sum_relative_coverage >= 5) %>%
  filter(!grepl(pattern = "unclassified", x = taxlevel_4, ignore.case = TRUE)) %>%
  pull(taxlevel_4) %>%
  unique()

# Draw dendrogram of order-rank taxa based on Bray-Curtis dissimilarity of gene composition
site_by_order_coverage <- phylodist_all_sum_cov %>%
  filter(taxlevel_4 %in% orders_5pct_coverage) %>%
  group_by(lake, metagenome, taxlevel_4) %>%
  summarize(sum_relative_coverage = sum(relative_coverage)) %>%
  ungroup() %>%
  spread(taxlevel_4, sum_relative_coverage, fill = 0) %>%
  mutate(lake_metagenome = str_c(lake, metagenome, sep = ".")) %>%
  select(-lake, -metagenome) %>%
  as.data.frame()
rownames(site_by_order_coverage) <- site_by_order_coverage$lake_metagenome
site_by_order_coverage$lake_metagenome <- NULL

site_by_order_bray <- vegdist(t(site_by_order_coverage), method = "bray")
site_by_order_hclust <- hclust(site_by_order_bray, method = "ward.D2")
site_by_order_dendro <- dendro_data(site_by_order_hclust)
site_by_order_cut <- cutree(site_by_order_hclust, k = 3)

site_by_order_cut_df <- data.frame(cluster = site_by_order_cut) %>%
  rownames_to_column("label") %>%
  as_tibble() %>%
  left_join(label(site_by_order_dendro), by = "label")

(dendrogram_phylodist_order <- ggplot() + 
    geom_segment(data = segment(site_by_order_dendro), aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_text(data = site_by_order_cut_df, aes(x, y, label = label, hjust = 0)) +
    coord_flip() +
    scale_y_reverse(expand = c(1,0)) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(), axis.title = element_blank(),
          panel.background = element_blank(), panel.grid = element_blank(),
          legend.position = "none"))
#ggsave(filename = "figures/04a_paleocapture_dendrogram_order_rmrna.pdf", plot = dendrogram_phylodist_order, device = "pdf", units = "in", width = 5, height = 10)

(heatmap_phylodist_order <- phylodist_all_sum_cov %>%
    filter(taxlevel_4 %in% orders_5pct_coverage) %>%
    group_by(lake, metagenome, taxlevel_4) %>%
    summarize(sum_relative_coverage = sum(relative_coverage)) %>%
    ggplot(aes(x = factor(metagenome, levels = c("swa", "tsa", "bsa", "swa_tsr", "swa_bsr", "tsa_bsr")),
               y = factor(taxlevel_4, levels = c("Methanomicrobiales", "Nitrosomonadales", "Nitrospirales", "Desulfuromonadales", "Myxococcales", "Clostridiales", "Desulfobacterales", "Syntrophobacterales", "Candidatus Nanopelagicales", "Pelagibacterales", "Burkholderiales", "Corynebacteriales", "Planctomycetales", "Rhizobiales", "Synechococcales", "Caudovirales", "Sphingomonadales")))) +
    facet_grid(~factor(lake, levels = c("Lac Paula", "Eightmile Lake", "Grand lac Touradi"))) +
    geom_tile(aes(fill = sum_relative_coverage), colour = "white") +
    #geom_text(aes(label = round(sum_relative_coverage, 1)), size = 2) +
    scale_x_discrete(expand = c(0,0)) +
    labs(fill = "Coverage (%)") +
    scale_fill_gradient(low = "white", high = "#7F7F7F") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "bottom", axis.title = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = -90, hjust = 0), axis.ticks = element_blank()))
#ggsave("figures/04a_paleocapture_heatmap_order_rmrna.pdf", plot = heatmap_phylodist_order, device = "pdf", units = "in", width = 10, height = 10)


# Summarize range of order-rank relative coverage in a CSV table
summary_phylodist_order_cov <- phylodist_all_sum_cov %>%
  filter(taxlevel_4 %in% orders_5pct_coverage) %>%
  group_by(lake, metagenome, tax_group, taxlevel_4) %>%
  summarize(sum_relative_coverage = round(sum(relative_coverage), 1)) %>%
  ungroup() %>%
  mutate(metagenome_type = case_when(!grepl("_", metagenome) ~ "free",
                                     grepl("_", metagenome) ~ "captured")) %>%
  spread(lake, sum_relative_coverage) %>%
  select(metagenome_type, metagenome, tax_group, taxlevel_4, `Lac Paula`, `Eightmile Lake`, `Grand lac Touradi`) %>%
  arrange(factor(metagenome, levels = c("swa", "tsa", "bsa", "swa_tsr", "swa_bsr", "tsa_bsr")), tax_group)
# summary_phylodist_order_cov %>%
#   write_csv("output/summary_order_cov_rmrna.csv")
