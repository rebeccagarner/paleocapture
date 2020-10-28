# 04_genecontent.R
# Calculate gene counts on scaffolds in free and captured metagenomes

setwd("paleocapture/")


# Load libraries
library(tidyverse)
library(gridExtra)
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
map <- rbind(map_06126_bottom, map_06126_top, map_06126_watercolumn,
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


#### Pull out scaffolds of interest ####
# Import list of scaffolds containing r/t RNA genes
scaffolds_rna <- read_tsv("output/scaffolds_rna.tsv") %>%
  pull(scaffold_id)

# Remove scaffolds containing r/t RNA genes in captured and free metagenomes
scaffolds_06126_swa_tsr <- covstats_06126_swa_tsr %>% left_join(map_06126_watercolumn, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_06126_swa_bsr <- covstats_06126_swa_bsr %>% left_join(map_06126_watercolumn, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_06126_tsa_bsr <- covstats_06126_tsa_bsr %>% left_join(map_06126_top, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)

scaffolds_17054_swa_tsr <- covstats_17054_swa_tsr %>% left_join(map_17054_watercolumn, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_17054_swa_bsr <- covstats_17054_swa_bsr %>% left_join(map_17054_watercolumn, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_17054_tsa_bsr <- covstats_17054_tsa_bsr %>% left_join(map_17054_top, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)

scaffolds_17067_swa_tsr <- covstats_17067_swa_tsr %>% left_join(map_17067_watercolumn, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_17067_swa_bsr <- covstats_17067_swa_bsr %>% left_join(map_17067_watercolumn, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_17067_tsa_bsr <- covstats_17067_tsa_bsr %>% left_join(map_17067_top, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)

scaffolds_06126_swa <- covstats_06126_swa %>% left_join(map_06126_watercolumn, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_06126_tsa <- covstats_06126_tsa %>% left_join(map_06126_top, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_06126_bsa <- covstats_06126_bsa %>% left_join(map_06126_bottom, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)

scaffolds_17054_swa <- covstats_17054_swa %>% left_join(map_17054_watercolumn, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_17054_tsa <- covstats_17054_tsa %>% left_join(map_17054_top, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_17054_bsa <- covstats_17054_bsa %>% left_join(map_17054_bottom, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)

scaffolds_17067_swa <- covstats_17067_swa %>% left_join(map_17067_watercolumn, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_17067_tsa <- covstats_17067_tsa %>% left_join(map_17067_top, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_17067_bsa <- covstats_17067_bsa %>% left_join(map_17067_bottom, by = "contig_id") %>% select(scaffold_id) %>% filter(!scaffold_id %in% scaffolds_rna)


#### Import files with genome feature information ####
gff_06126_bottom <- read_tsv(file = "data/img/lp2017_06-126-bottom/3300032883.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type)
gff_06126_top <- read_tsv(file = "data/img/lp2017_06-126-top/3300032795.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type)
gff_06126_watercolumn <- read_tsv(file = "data/img/lp2017_06-126-watercolumn/3300032269.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type)

gff_17054_bottom <- read_tsv(file = "data/img/lp2017_17-054-bottom/3300032797.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type)
gff_17054_top <- read_tsv(file = "data/img/lp2017_17-054-top/3300032796.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type)
gff_17054_watercolumn <- read_tsv(file = "data/img/lp2017_17-054-watercolumn/3300032328.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type)

gff_17067_bottom <- read_tsv(file = "data/img/lp2017_17-067-bottom/3300032799.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type)
gff_17067_top <- read_tsv(file = "data/img/lp2017_17-067-top/3300032798.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type)
gff_17067_watercolumn <- read_tsv(file = "data/img/lp2017_17-067-watercolumn/3300032317.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type)


#### Summarize gene counts for unique scaffolds ####
gene_freq_06126_bottom <- gff_06126_bottom %>% filter(type != "repeat_region" | type != "exon") %>% select(scaffold_id) %>% mutate(gene_count = 1) %>% group_by(scaffold_id) %>% summarize(sum_gene_count = sum(gene_count))
gene_freq_06126_top <- gff_06126_top %>% filter(type != "repeat_region" | type != "exon") %>% select(scaffold_id) %>% mutate(gene_count = 1) %>% group_by(scaffold_id) %>% summarize(sum_gene_count = sum(gene_count))
gene_freq_06126_watercolumn <- gff_06126_watercolumn %>% filter(type != "repeat_region" | type != "exon") %>% select(scaffold_id) %>% mutate(gene_count = 1) %>% group_by(scaffold_id) %>% summarize(sum_gene_count = sum(gene_count))

gene_freq_17054_bottom <- gff_17054_bottom %>% filter(type != "repeat_region" | type != "exon") %>% select(scaffold_id) %>% mutate(gene_count = 1) %>% group_by(scaffold_id) %>% summarize(sum_gene_count = sum(gene_count))
gene_freq_17054_top <- gff_17054_top %>% filter(type != "repeat_region" | type != "exon") %>% select(scaffold_id) %>% mutate(gene_count = 1) %>% group_by(scaffold_id) %>% summarize(sum_gene_count = sum(gene_count))
gene_freq_17054_watercolumn <- gff_17054_watercolumn %>% filter(type != "repeat_region" | type != "exon") %>% select(scaffold_id) %>% mutate(gene_count = 1) %>% group_by(scaffold_id) %>% summarize(sum_gene_count = sum(gene_count))

gene_freq_17067_bottom <- gff_17067_bottom %>% filter(type != "repeat_region" | type != "exon") %>% select(scaffold_id) %>% mutate(gene_count = 1) %>% group_by(scaffold_id) %>% summarize(sum_gene_count = sum(gene_count))
gene_freq_17067_top <- gff_17067_top %>% filter(type != "repeat_region" | type != "exon") %>% select(scaffold_id) %>% mutate(gene_count = 1) %>% group_by(scaffold_id) %>% summarize(sum_gene_count = sum(gene_count))
gene_freq_17067_watercolumn <- gff_17067_watercolumn %>% filter(type != "repeat_region" | type != "exon") %>% select(scaffold_id) %>% mutate(gene_count = 1) %>% group_by(scaffold_id) %>% summarize(sum_gene_count = sum(gene_count))


#### Count frequencies of genes per scaffold ####
scaffold_gene_count_06126_swa_tsr <- scaffolds_06126_swa_tsr %>% left_join(gene_freq_06126_watercolumn, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Lac Paula", metagenome = "swa_tsr")
scaffold_gene_count_06126_swa_bsr <- scaffolds_06126_swa_bsr %>% left_join(gene_freq_06126_watercolumn, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Lac Paula", metagenome = "swa_bsr")
scaffold_gene_count_06126_tsa_bsr <- scaffolds_06126_tsa_bsr %>% left_join(gene_freq_06126_top, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Lac Paula", metagenome = "tsa_bsr")

scaffold_gene_count_17054_swa_tsr <- scaffolds_17054_swa_tsr %>% left_join(gene_freq_17054_watercolumn, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Eightmile Lake", metagenome = "swa_tsr")
scaffold_gene_count_17054_swa_bsr <- scaffolds_17054_swa_bsr %>% left_join(gene_freq_17054_watercolumn, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Eightmile Lake", metagenome = "swa_bsr")
scaffold_gene_count_17054_tsa_bsr <- scaffolds_17054_tsa_bsr %>% left_join(gene_freq_17054_top, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Eightmile Lake", metagenome = "tsa_bsr")

scaffold_gene_count_17067_swa_tsr <- scaffolds_17067_swa_tsr %>% left_join(gene_freq_17067_watercolumn, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Grand lac Touradi", metagenome = "swa_tsr")
scaffold_gene_count_17067_swa_bsr <- scaffolds_17067_swa_bsr %>% left_join(gene_freq_17067_watercolumn, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Grand lac Touradi", metagenome = "swa_bsr")
scaffold_gene_count_17067_tsa_bsr <- scaffolds_17067_tsa_bsr %>% left_join(gene_freq_17067_top, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Grand lac Touradi", metagenome = "tsa_bsr")

scaffold_gene_count_06126_swa <- scaffolds_06126_swa %>% left_join(gene_freq_06126_watercolumn, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Lac Paula", metagenome = "swa")
scaffold_gene_count_06126_tsa <- scaffolds_06126_tsa %>% left_join(gene_freq_06126_top, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Lac Paula", metagenome = "tsa")
scaffold_gene_count_06126_bsa <- scaffolds_06126_bsa %>% left_join(gene_freq_06126_bottom, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Lac Paula", metagenome = "bsa")

scaffold_gene_count_17054_swa <- scaffolds_17054_swa %>% left_join(gene_freq_17054_watercolumn, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Eightmile Lake", metagenome = "swa")
scaffold_gene_count_17054_tsa <- scaffolds_17054_tsa %>% left_join(gene_freq_17054_top, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Eightmile Lake", metagenome = "tsa")
scaffold_gene_count_17054_bsa <- scaffolds_17054_bsa %>% left_join(gene_freq_17054_bottom, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Eightmile Lake", metagenome = "bsa")

scaffold_gene_count_17067_swa <- scaffolds_17067_swa %>% left_join(gene_freq_17067_watercolumn, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Grand lac Touradi", metagenome = "swa")
scaffold_gene_count_17067_tsa <- scaffolds_17067_tsa %>% left_join(gene_freq_17067_top, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Grand lac Touradi", metagenome = "tsa")
scaffold_gene_count_17067_bsa <- scaffolds_17067_bsa %>% left_join(gene_freq_17067_bottom, by = "scaffold_id") %>% mutate(sum_gene_count = coalesce(sum_gene_count, 0)) %>% count(sum_gene_count) %>% mutate(lake = "Grand lac Touradi", metagenome = "bsa")


#### Plot scaffold gene frequencies ####
palette_shape <- c(87, 84, 66, 25, 25, 8)
names(palette_shape) <- c("swa", "tsa", "bsa", "swa_tsr", "swa_bsr", "tsa_bsr")

palette_lakes <- c("#ff9f1a", "#18dcff", "#7d5fff")
names(palette_lakes) <- c("Lac Paula", "Eightmile Lake", "Grand lac Touradi")

(gene_freq_plot_06126_cap <- bind_rows(scaffold_gene_count_06126_swa_tsr,
                                       scaffold_gene_count_06126_swa_bsr,
                                       scaffold_gene_count_06126_tsa_bsr) %>%
    mutate(fill_shape = case_when(metagenome == "swa_tsr" ~ 1,
                                  metagenome == "swa_bsr" ~ 0,
                                  metagenome == "tsa_bsr" ~ 1)) %>%
    ggplot(aes(x = sum_gene_count, y = n, shape = metagenome, colour = lake, fill = ifelse(fill_shape == 1, lake, NA))) +
    geom_point(size = 3, alpha = 0.7) +
    scale_shape_manual(values = palette_shape, guide = "none") +
    scale_fill_manual(values = palette_lakes, guide = "none") +
    scale_color_manual(values = palette_lakes, guide = "none") +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(trans = pseudo_log_trans(1, 10),
                       breaks = c(0, 1, 10, 100)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()))

(gene_freq_plot_17054_cap <- bind_rows(scaffold_gene_count_17054_swa_tsr,
                                       scaffold_gene_count_17054_swa_bsr,
                                       scaffold_gene_count_17054_tsa_bsr) %>%
    mutate(fill_shape = case_when(metagenome == "swa_tsr" ~ 1,
                                  metagenome == "swa_bsr" ~ 0,
                                  metagenome == "tsa_bsr" ~ 1)) %>%
    ggplot(aes(x = sum_gene_count, y = n, shape = metagenome, colour = lake, fill = ifelse(fill_shape == 1, lake, NA))) +
    geom_point(size = 3, alpha = 0.7) +
    scale_shape_manual(values = palette_shape, guide = "none") +
    scale_fill_manual(values = palette_lakes, guide = "none") +
    scale_color_manual(values = palette_lakes, guide = "none") +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(trans = pseudo_log_trans(1, 10),
                       breaks = c(0, 1, 10, 100)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()))

(gene_freq_plot_17067_cap <- bind_rows(scaffold_gene_count_17067_swa_tsr,
                                       scaffold_gene_count_17067_swa_bsr,
                                       scaffold_gene_count_17067_tsa_bsr) %>%
    mutate(fill_shape = case_when(metagenome == "swa_tsr" ~ 1,
                                  metagenome == "swa_bsr" ~ 0,
                                  metagenome == "tsa_bsr" ~ 1)) %>%
    ggplot(aes(x = sum_gene_count, y = n, shape = metagenome, colour = lake, fill = ifelse(fill_shape == 1, lake, NA))) +
    geom_point(size = 3, alpha = 0.7) +
    scale_shape_manual(values = palette_shape, guide = "none") +
    scale_fill_manual(values = palette_lakes, guide = "none") +
    scale_color_manual(values = palette_lakes, guide = "none") +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(trans = pseudo_log_trans(1, 10),
                       breaks = c(0, 1, 10, 100)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()))

(gene_freq_plot_06126_free <- bind_rows(scaffold_gene_count_06126_swa,
                                       scaffold_gene_count_06126_tsa,
                                       scaffold_gene_count_06126_bsa) %>%
    ggplot(aes(x = sum_gene_count, y = n, shape = metagenome, colour = lake)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_shape_manual(values = palette_shape, guide = "none") +
    scale_fill_manual(values = palette_lakes, guide = "none") +
    scale_color_manual(values = palette_lakes, guide = "none") +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(trans = pseudo_log_trans(1, 10),
                       breaks = c(0, 1, 10, 100, 1000)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()))

(gene_freq_plot_17054_free <- bind_rows(scaffold_gene_count_17054_swa,
                                       scaffold_gene_count_17054_tsa,
                                       scaffold_gene_count_17054_bsa) %>%
    ggplot(aes(x = sum_gene_count, y = n, shape = metagenome, colour = lake)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_shape_manual(values = palette_shape, guide = "none") +
    scale_fill_manual(values = palette_lakes, guide = "none") +
    scale_color_manual(values = palette_lakes, guide = "none") +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(trans = pseudo_log_trans(1, 10),
                       breaks = c(0, 1, 10, 100, 1000)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()))

(gene_freq_plot_17067_free <- bind_rows(scaffold_gene_count_17067_swa,
                                       scaffold_gene_count_17067_tsa,
                                       scaffold_gene_count_17067_bsa) %>%
    ggplot(aes(x = sum_gene_count, y = n, shape = metagenome, colour = lake)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_shape_manual(values = palette_shape, guide = "none") +
    scale_fill_manual(values = palette_lakes, guide = "none") +
    scale_color_manual(values = palette_lakes, guide = "none") +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_x_continuous(trans = pseudo_log_trans(1, 10),
                       breaks = c(0, 1, 10, 100, 1000)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()))

#### Create composite plot of scaffold gene frequencies ####
(gene_freq_plot_all <- grid.arrange(gene_freq_plot_06126_free,
                                    gene_freq_plot_06126_cap,
                                    gene_freq_plot_17054_free,
                                    gene_freq_plot_17054_cap,
                                    gene_freq_plot_17067_free,
                                    gene_freq_plot_17067_cap,
                                    nrow = 2, ncol = 3,
                                    layout_matrix = rbind(c(1, 3, 5),
                                                          c(2, 4, 6)),
                                    left = "n scaffolds",
                                    bottom = "n genes on scaffolds"))
#ggsave(filename = "figures/s01_paleocapture_scaffolds_genes.pdf", plot = gene_freq_plot_all, device = "pdf", units = "in", width = 10, height = 8)


#### Calculate stats for  n genes on scaffolds among lakes and metagenomes ####
scaffold_gene_count_all <- bind_rows(scaffold_gene_count_06126_swa_tsr,
                                     scaffold_gene_count_06126_swa_bsr,
                                     scaffold_gene_count_06126_tsa_bsr,
                                     
                                     scaffold_gene_count_17054_swa_tsr,
                                     scaffold_gene_count_17054_swa_bsr,
                                     scaffold_gene_count_17054_tsa_bsr,
                                     
                                     scaffold_gene_count_17067_swa_tsr,
                                     scaffold_gene_count_17067_swa_bsr,
                                     scaffold_gene_count_17067_tsa_bsr,
                                     
                                     scaffold_gene_count_06126_swa,
                                     scaffold_gene_count_06126_tsa,
                                     scaffold_gene_count_06126_bsa,
                                     
                                     scaffold_gene_count_17054_swa,
                                     scaffold_gene_count_17054_tsa,
                                     scaffold_gene_count_17054_bsa,
                                     
                                     scaffold_gene_count_17067_swa,
                                     scaffold_gene_count_17067_tsa,
                                     scaffold_gene_count_17067_bsa)

#### Scaffold gene content and summary statistics ####
# Calculate n scaffolds with zero and above-zero n genes
scaffold_gene_content <- scaffold_gene_count_all %>%
  mutate(gene_content = case_when(sum_gene_count == 0 ~ "scaffolds_zero_genes",
                                  sum_gene_count > 0 ~ "scaffolds_with_genes")) %>%
  group_by(lake, metagenome, gene_content) %>%
  summarize(sum_scaffolds = sum(n)) %>%
  spread(gene_content, sum_scaffolds)

# Summarize n scaffolds, genes, and mean gene counts among lakes and metagenomes
(scaffold_gene_count_summary <- scaffold_gene_count_all %>%
    group_by(lake, metagenome) %>%
    summarize(nscaffolds_total = sum(n),
              ngenes_total = sum(n*sum_gene_count)) %>%
    left_join(scaffold_gene_content, by = c("lake", "metagenome")) %>%
    mutate(scaffolds_with_genes_percent = round(scaffolds_with_genes/nscaffolds_total * 100, 1),
           scaffolds_zero_genes_percent = round(scaffolds_zero_genes/nscaffolds_total * 100, 1)))

# Calculate percentages of scaffolds with and without gene content
scaffold_gene_count_summary %>%
  ungroup() %>%
  mutate(metagenome_type = case_when(!grepl("_", metagenome) ~ "free",
                                     grepl("_", metagenome) ~ "captured")) %>%
  group_by(metagenome_type) %>%
  summarize(mean_nscaffolds_with_genes = round(mean(scaffolds_with_genes_percent), 1),
            sd_nscaffolds_with_genes = round(sd(scaffolds_with_genes_percent), 1))

# Calculate scaffold gene content range (min and max)
scaffold_gene_count_range <- scaffold_gene_count_all %>%
  group_by(lake, metagenome) %>%
  summarize(min_gene_count = min(sum_gene_count),
            max_gene_count = max(sum_gene_count))

# Calculate percentages of scaffolds with single gene content
scaffold_gene_count_summarystats <- scaffold_gene_count_all %>%
  filter(sum_gene_count == 1) %>%
  left_join(scaffold_gene_count_summary, by = c("lake", "metagenome")) %>%
  mutate(scaffolds_single_gene_percent = round(n/nscaffolds_total * 100, 1)) %>%
  rename(scaffolds_single_gene = n) %>%
  select(-sum_gene_count) %>%
  left_join(scaffold_gene_count_range, by = c("lake", "metagenome")) %>%
  mutate(metagenome_type = case_when(!grepl("_", metagenome) ~ "free",
                                     grepl("_", metagenome) ~ "captured")) %>%
  select(lake, metagenome_type, metagenome, nscaffolds_total, ngenes_total,
         scaffolds_with_genes, scaffolds_with_genes_percent,
         scaffolds_zero_genes, scaffolds_zero_genes_percent,
         scaffolds_single_gene, scaffolds_single_gene_percent,
         min_gene_count, max_gene_count) %>%
  arrange(factor(lake, levels = c("Lac Paula", "Eightmile Lake", "Grand lac Touradi")),
          factor(metagenome_type, levels = c("free", "experimental")),
          factor(metagenome, levels = c("swa", "tsa", "bsa", "swa_tsr", "swa_bsr", "tsa_bsr")))
#write_csv(scaffold_gene_count_summarystats, "output/scaffolds_genes.csv", col_names = TRUE)
