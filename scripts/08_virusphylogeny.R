# 08_virusphylogeny.R
# Identify genes encoding Caudovirales terminase large subunit

setwd("paleocapture/")

# Load libraries
library(tidyverse)
library(seqinr)

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
rm(list = ls(pattern = "phylodist_\\d"))
rm(list = ls(pattern = "map"))
rm(list = ls(pattern = "covstats_\\d"))
rm(list = ls(pattern = "scaffolds_\\d"))


# Sum coverage for each metagenome
phylodist_sum_cov <- phylodist_all %>%
  group_by(lake, metagenome) %>%
  summarize(sum_avgfold = sum(Avg_fold))

# Scale relative coverage to 100% metagenome profile
phylodist_all_sum_cov <- phylodist_all %>%
  left_join(phylodist_sum_cov, by = c("lake", "metagenome")) %>%
  mutate(relative_coverage = (Avg_fold/sum_avgfold * 100))


#### Import IMG .ko.txt files ####
ko_06126_watercolumn <- read_tsv(file = "data/img/lp2017_06-126-watercolumn/3300032269.a.ko.txt", col_names = c("gene_id", "img_ko_flag", "ko_term", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length"), cols_only(gene_id = col_character(), ko_term = col_character())) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length)) %>% select(-project_id, -scaffold_id_length) %>% mutate(ko_term = str_remove(ko_term, "KO:"))
ko_06126_top <- read_tsv(file = "data/img/lp2017_06-126-top/3300032795.a.ko.txt", col_names = c("gene_id", "img_ko_flag", "ko_term", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length"), cols_only(gene_id = col_character(), ko_term = col_character())) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length)) %>% select(-project_id, -scaffold_id_length) %>% mutate(ko_term = str_remove(ko_term, "KO:"))
ko_06126_bottom <- read_tsv(file = "data/img/lp2017_06-126-bottom/3300032883.a.ko.txt", col_names = c("gene_id", "img_ko_flag", "ko_term", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length"), cols_only(gene_id = col_character(), ko_term = col_character())) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length)) %>% select(-project_id, -scaffold_id_length) %>% mutate(ko_term = str_remove(ko_term, "KO:"))

ko_17054_watercolumn <- read_tsv(file = "data/img/lp2017_17-054-watercolumn/3300032328.a.ko.txt", col_names = c("gene_id", "img_ko_flag", "ko_term", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length"), cols_only(gene_id = col_character(), ko_term = col_character())) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length)) %>% select(-project_id, -scaffold_id_length) %>% mutate(ko_term = str_remove(ko_term, "KO:"))
ko_17054_top <- read_tsv(file = "data/img/lp2017_17-054-top/3300032796.a.ko.txt", col_names = c("gene_id", "img_ko_flag", "ko_term", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length"), cols_only(gene_id = col_character(), ko_term = col_character())) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length)) %>% select(-project_id, -scaffold_id_length) %>% mutate(ko_term = str_remove(ko_term, "KO:"))
ko_17054_bottom <- read_tsv(file = "data/img/lp2017_17-054-bottom/3300032797.a.ko.txt", col_names = c("gene_id", "img_ko_flag", "ko_term", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length"), cols_only(gene_id = col_character(), ko_term = col_character())) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length)) %>% select(-project_id, -scaffold_id_length) %>% mutate(ko_term = str_remove(ko_term, "KO:"))

ko_17067_watercolumn <- read_tsv(file = "data/img/lp2017_17-067-watercolumn/3300032317.a.ko.txt", col_names = c("gene_id", "img_ko_flag", "ko_term", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length"), cols_only(gene_id = col_character(), ko_term = col_character())) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length)) %>% select(-project_id, -scaffold_id_length) %>% mutate(ko_term = str_remove(ko_term, "KO:"))
ko_17067_top <- read_tsv(file = "data/img/lp2017_17-067-top/3300032798.a.ko.txt", col_names = c("gene_id", "img_ko_flag", "ko_term", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length"), cols_only(gene_id = col_character(), ko_term = col_character())) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length)) %>% select(-project_id, -scaffold_id_length) %>% mutate(ko_term = str_remove(ko_term, "KO:"))
ko_17067_bottom <- read_tsv(file = "data/img/lp2017_17-067-bottom/3300032799.a.ko.txt", col_names = c("gene_id", "img_ko_flag", "ko_term", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length"), cols_only(gene_id = col_character(), ko_term = col_character())) %>% mutate(project_id = str_remove(gene_id, pattern = "_.*")) %>% left_join(img_scaffold_ids, by = "project_id") %>% mutate(scaffold_id = str_sub(gene_id, 1, scaffold_id_length)) %>% select(-project_id, -scaffold_id_length) %>% mutate(ko_term = str_remove(ko_term, "KO:"))


#### Join KO categories and protein products from table #####
ko_table <- read_csv("data/ko.csv") %>%
  distinct(ko_term, ko_protein)

ko_all <- bind_rows(ko_06126_watercolumn,
                    ko_06126_top,
                    ko_06126_bottom,
                    
                    ko_17054_watercolumn,
                    ko_17054_top,
                    ko_17054_bottom,
                    
                    ko_17067_watercolumn,
                    ko_17067_top,
                    ko_17067_bottom) %>%
  left_join(ko_table, by = "ko_term")

# Clear up space in R environment
rm(list = ls(pattern = "ko_\\d"))


#### Join phylodist and KO matrices ####
phylodist_ko <- phylodist_all_sum_cov %>%
  left_join(ko_all, by = c("gene_id", "scaffold_id")) %>%
  filter(!is.na(ko_term))


#### Identify Caudovirales terminase large subunit annotations
# Caudovirales terminase large subunit
terminase_caudovirales <- phylodist_ko %>%
  filter(taxlevel_4 == "Caudovirales") %>%
  filter(grepl("terminase", ko_protein, ignore.case = TRUE)) %>%
  filter(ko_term == "K06909")

# Write list of gene IDs with Caudovirales terminase large subunit annotations
# terminase_caudovirales %>%
#   distinct(gene_id) %>%
#   mutate(headers = str_c(">", gene_id)) %>%
#   select(headers) %>%
#   write_delim(path = "output/caudovirales_terminase_genes_ko.list", delim = "\n", col_names = FALSE)


#### Annotate phylogenetic tree nodes with isolation source and taxonomic assignment ####
# Import .fas file with genes of interest
terminase_fasta <- read.fasta("output/caudovirales_terminase_genes_ko.list.fas", seqtype = "AA", as.string = TRUE)

terminase_seqs <- tibble(header = names(terminase_fasta),
                         sequence = unlist(getSequence(terminase_fasta, as.string = TRUE)))

terminase_cat <- terminase_caudovirales %>%
  select(gene_id, lake, metagenome, taxlevel_8) %>%
  mutate(gene_count = 1) %>%
  spread(metagenome, gene_count, fill = 0) %>%
  mutate(isolation = case_when(swa == 1 & swa_tsr == 0 & swa_bsr == 0 ~ "SW",
                               swa == 1 & swa_tsr == 1 & swa_bsr == 0 ~ "SW SWa-TSr",
                               swa == 1 & swa_tsr == 0 & swa_bsr == 1 ~ "SW SWa-BSr",
                               swa == 1 & swa_tsr == 1 & swa_bsr == 1 ~ "SW SWa-TSr SWa-BSr",
                               
                               tsa == 1 & tsa_swr == 0 & tsa_bsr == 0 ~ "TS",
                               tsa == 1 & tsa_swr == 1 & tsa_bsr == 0 ~ "TS TSa-SWr",
                               tsa == 1 & tsa_swr == 0 & tsa_bsr == 1 ~ "TS TSa-BSr",
                               tsa == 1 & tsa_swr == 1 & tsa_bsr == 1 ~ "TS TSa-SWr TSa-BSr",
                               
                               bsa == 1 & bsa_swr == 0 & bsa_tsr == 0 ~ "BS",
                               bsa == 1 & bsa_swr == 1 & bsa_tsr == 0 ~ "BS BSa-SWr",
                               bsa == 1 & bsa_swr == 0 & bsa_tsr == 1 ~ "BS BSa-TSr",
                               bsa == 1 & bsa_swr == 1 & bsa_tsr == 1 ~ "BS BSa-SWr BSa-TSr")) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  mutate(header_cat = str_c(gene_id, taxlevel_8, lake, isolation, sep = " ")) %>%
  mutate(header_cat = str_replace_all(header_cat, " ", "_")) %>%
  select(gene_id, header_cat) %>%
  left_join(terminase_seqs, by = c("gene_id" = "header"))
#write.fasta(as.list(terminase_cat$sequence), terminase_cat$header_cat, "C:/Users/Gandalf/Desktop/terminase_cat_headers.fas")
