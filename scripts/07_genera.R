# 07_genera.R
# Count genes by taxon

setwd("paleocapture/")

# Load libraries
library(tidyverse)
library(scales)

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
phylodist_06126_swa_tsr <- scaffolds_06126_swa_tsr %>% left_join(phylodist_06126_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Lac Paula", metagenome = "swa_tsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_06126_swa_bsr <- scaffolds_06126_swa_bsr %>% left_join(phylodist_06126_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Lac Paula", metagenome = "swa_bsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_06126_tsa_bsr <- scaffolds_06126_tsa_bsr %>% left_join(phylodist_06126_top, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Lac Paula", metagenome = "tsa_bsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))

phylodist_17054_swa_tsr <- scaffolds_17054_swa_tsr %>% left_join(phylodist_17054_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Eightmile Lake", metagenome = "swa_tsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17054_swa_bsr <- scaffolds_17054_swa_bsr %>% left_join(phylodist_17054_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Eightmile Lake", metagenome = "swa_bsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17054_tsa_bsr <- scaffolds_17054_tsa_bsr %>% left_join(phylodist_17054_top, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Eightmile Lake", metagenome = "tsa_bsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))

phylodist_17067_swa_tsr <- scaffolds_17067_swa_tsr %>% left_join(phylodist_17067_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Grand lac Touradi", metagenome = "swa_tsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17067_swa_bsr <- scaffolds_17067_swa_bsr %>% left_join(phylodist_17067_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Grand lac Touradi", metagenome = "swa_bsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17067_tsa_bsr <- scaffolds_17067_tsa_bsr %>% left_join(phylodist_17067_top, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Grand lac Touradi", metagenome = "tsa_bsr") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))

phylodist_06126_swa <- scaffolds_06126_swa %>% left_join(phylodist_06126_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Lac Paula", metagenome = "swa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_06126_tsa <- scaffolds_06126_tsa %>% left_join(phylodist_06126_top, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Lac Paula", metagenome = "tsa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_06126_bsa <- scaffolds_06126_bsa %>% left_join(phylodist_06126_bottom, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Lac Paula", metagenome = "bsa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))

phylodist_17054_swa <- scaffolds_17054_swa %>% left_join(phylodist_17054_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Eightmile Lake", metagenome = "swa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17054_tsa <- scaffolds_17054_tsa %>% left_join(phylodist_17054_top, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Eightmile Lake", metagenome = "tsa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17054_bsa <- scaffolds_17054_bsa %>% left_join(phylodist_17054_bottom, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Eightmile Lake", metagenome = "bsa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))

phylodist_17067_swa <- scaffolds_17067_swa %>% left_join(phylodist_17067_watercolumn, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Grand lac Touradi", metagenome = "swa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17067_tsa <- scaffolds_17067_tsa %>% left_join(phylodist_17067_top, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Grand lac Touradi", metagenome = "tsa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))
phylodist_17067_bsa <- scaffolds_17067_bsa %>% left_join(phylodist_17067_bottom, by = "scaffold_id") %>% filter(!is.na(gene_count)) %>% mutate(lake = "Grand lac Touradi", metagenome = "bsa") %>% separate(col = "lineage", into = c("taxlevel_1", "taxlevel_2", "taxlevel_3", "taxlevel_4", "taxlevel_5", "taxlevel_6", "taxlevel_7", "taxlevel_8"), sep = ";") %>% mutate(tax_group = case_when(taxlevel_1 == "Viruses" ~ taxlevel_1, taxlevel_2 == "Proteobacteria" ~ taxlevel_3, TRUE ~ taxlevel_2))


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


#### Calculate relative coverage ####
# Sum coverage for each metagenome
phylodist_sum_cov <- phylodist_all %>%
  group_by(lake, metagenome) %>%
  summarize(sum_avgfold = sum(Avg_fold))

# Scale relative coverage to 100% metagenome profile
phylodist_all_sum_cov <- phylodist_all %>%
  left_join(phylodist_sum_cov, by = c("lake", "metagenome")) %>%
  mutate(relative_coverage = (Avg_fold/sum_avgfold * 100))


#### Caudovirales ####
ngenes_order <- phylodist_all %>%
  group_by(lake, metagenome, taxlevel_1, taxlevel_2, taxlevel_3, taxlevel_4) %>%
  summarize(sum_gene_count = sum(gene_count))

ngenes_caudovirales <- ngenes_order %>%
  filter(taxlevel_4 == "Caudovirales") %>%
  ungroup() %>%
  group_by(lake, metagenome, taxlevel_4) %>%
  summarize(all_gene_count = sum(sum_gene_count)) %>%
  mutate(taxlevel_1 = "Viruses", tax_group = "Viruses") %>%
  rename(taxlevel_6 = taxlevel_4)


#### Identify genera with >= 1% coverage in free and captured SW metagenomes ####
genera_1pct_cov <- phylodist_all_sum_cov %>%
  filter(grepl("sw", metagenome)) %>%
  group_by(lake, metagenome, taxlevel_6) %>%
  summarize(sum_relative_coverage = sum(relative_coverage)) %>%
  filter(sum_relative_coverage >= 1) %>%
  filter(!grepl(pattern = "unclassified", x = taxlevel_6, ignore.case = TRUE)) %>%
  pull(taxlevel_6) %>%
  unique()

ngenes_genus <- phylodist_all %>%
  group_by(lake, metagenome, taxlevel_1, taxlevel_2, taxlevel_3, taxlevel_4, taxlevel_5, taxlevel_6) %>%
  summarize(sum_gene_count = sum(gene_count))

# Associate genera with coarser level taxonomy
genus_tax <- phylodist_all %>%
  filter(taxlevel_6 %in% genera_1pct_cov) %>%
  distinct(taxlevel_1, tax_group, taxlevel_6) %>%
  bind_rows(tibble(taxlevel_1 = "Viruses", tax_group = "Viruses", taxlevel_6 = "Caudovirales"))

# Count n genes by genera, metagenome, lake
ngenes_genus_sw <- ngenes_genus %>%
  filter(taxlevel_6 %in% genera_1pct_cov) %>%
  ungroup() %>%
  group_by(lake, metagenome, taxlevel_6) %>%
  summarize(all_gene_count = sum(sum_gene_count)) %>%
  left_join(genus_tax, by = "taxlevel_6")


#### Palettes ####
palette_shapes <- c(87, 84, 66, 25, 25, 8)
names(palette_shapes) <- c("swa", "tsa", "bsa", "swa_tsr", "swa_bsr", "tsa_bsr")

palette_lakes <- c("#ff9f1a", "#18dcff", "#7d5fff")
names(palette_lakes) <- c("Lac Paula", "Eightmile Lake", "Grand lac Touradi")


#### Plot n genes by genus as a scatter plot ####
# Plot ALL genera with coverage >= 1% in a SW metagenome
(ngenes_genus_all_plot <- ngenes_genus_sw %>%
    bind_rows(ngenes_caudovirales) %>%
    filter(!grepl("_", metagenome)) %>%
    ggplot(aes(x = factor(lake, levels = c("Lac Paula", "Eightmile Lake", "Grand lac Touradi")),
               y = all_gene_count,
               colour = factor(lake, levels = c("Lac Paula", "Eightmile Lake", "Grand lac Touradi")))) +
    geom_point(aes(shape = factor(metagenome, levels = c("swa", "tsa", "bsa"))), size = 3) +
    facet_wrap(taxlevel_1~tax_group + taxlevel_6, nrow = 2) +
    scale_shape_manual(values = palette_shapes, guide = "none") +
    scale_colour_manual(values = palette_lakes, guide = "none") +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    ylab("Number of genes") +
    theme(panel.grid.major.x = element_line(colour = "lightgrey", linetype = "dashed"),
          panel.background = element_rect(),
          axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_blank(), strip.placement = "out"))
#ggsave("figures/s03_paleocapture_ngenes_genera_all.pdf", plot = ngenes_genus_all_plot, device = "pdf", units = "in", width = 14, height = 16)

# Reorder genera by taxonomic groups and difference in n genes between SW and sediments
genera_in_order <- ngenes_genus_sw %>%
  bind_rows(ngenes_caudovirales) %>%
  filter(!grepl("_", metagenome)) %>%
  filter(taxlevel_6 %in% c("Candidatus Fonsibacter", "Candidatus Nanopelagicus", "Candidatus Planktophila", "Limnohabitans", "Polynucleobacter",
                           "Rhodoluna", "Sphingomonas", "Cyanobium", "Candidatus Methylopumilus", "Ilumatobacter",
                           "Synechococcus", "Aphanothece", "Mycobacterium", "Mycolicibacterium",
                           "Phaeodactylum", "Thalassiosira", "Rhodoferax", "Caudovirales")) %>%
  ungroup() %>%
  group_by(metagenome, taxlevel_6) %>%
  summarize(mean_gene_count = mean(all_gene_count)) %>%
  spread(metagenome, mean_gene_count) %>%
  mutate(sed = (tsa + bsa)/2) %>%
  mutate(swsed_difference = (swa - sed)/swa) %>%
  left_join(genus_tax, by = "taxlevel_6") %>%
  arrange(factor(taxlevel_1, levels = c("Viruses", "Bacteria", "Eukaryota")),
          factor(tax_group, levels = c("Viruses", "Cyanobacteria", "Betaproteobacteria", "Alphaproteobacteria", "Actinobacteria", "Bacillariophyta", "Ascomycota")),
          swsed_difference) %>%
  pull(taxlevel_6)

# Plot subset of genera mentionned in text
(ngenes_genus_subset_plot <- ngenes_genus_sw %>%
    bind_rows(ngenes_caudovirales) %>%
    filter(!grepl("_", metagenome)) %>%
    filter(taxlevel_6 %in% genera_in_order)  %>%
    ggplot(aes(x = factor(lake, levels = c("Lac Paula", "Eightmile Lake", "Grand lac Touradi")),
               y = all_gene_count,
               colour = factor(lake, levels = c("Lac Paula", "Eightmile Lake", "Grand lac Touradi")))) +
    geom_point(aes(shape = factor(metagenome, levels = c("swa", "tsa", "bsa"))), size = 3) +
    facet_grid(~factor(taxlevel_1, levels = c("Viruses", "Bacteria", "Eukaryota")) +
                 factor(tax_group, levels = c("Viruses", "Cyanobacteria", "Betaproteobacteria", "Alphaproteobacteria", "Actinobacteria", "Bacillariophyta", "Ascomycota")) +
                 factor(taxlevel_6, levels = genera_in_order)) +
    scale_shape_manual(values = palette_shapes, guide = "none") +
    scale_colour_manual(values = palette_lakes, guide = "none") +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    ylab("Number of genes") +
    theme(panel.grid.major.x = element_line(colour = "lightgrey", linetype = "dashed"),
          panel.background = element_rect(),
          axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_blank(), strip.placement = "out"))
#ggsave("figures/06_paleocapture_ngenes_genera.pdf", plot = ngenes_genus_subset_plot, device = "pdf", units = "in", width = 14, height = 8)

# Summarize range of order-rank relative coverage in a CSV table
summary_ngenes_genus <- ngenes_genus_sw %>%
  bind_rows(ngenes_caudovirales) %>%
  ungroup() %>%
  filter(!grepl("_", metagenome)) %>%
  spread(lake, all_gene_count) %>%
  select(metagenome, taxlevel_1, tax_group, taxlevel_6, `Lac Paula`, `Eightmile Lake`, `Grand lac Touradi`) %>%
  arrange(factor(metagenome, levels = c("swa", "tsa", "bsa")), taxlevel_1, tax_group)
# summary_ngenes_genus %>%
#   write_csv("output/05_paleocapture_summary_ngenes_genus.csv")
