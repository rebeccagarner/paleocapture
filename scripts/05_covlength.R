# 05_covlength.R
# Scaffold coverage vs. length

setwd("paleocapture/")


# Load libraries
library(tidyverse)
library(scales)
library(gridExtra)


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
covstats_06126_swa_tsr <- read_tsv(file = "data/bbmap/covstats/lp2017_06-126-swa_tsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold, Length) %>% filter(Avg_fold > 0)
covstats_06126_swa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_06-126-swa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold, Length) %>% filter(Avg_fold > 0)
covstats_06126_tsa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_06-126-tsa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold, Length) %>% filter(Avg_fold > 0)

covstats_17054_swa_tsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-054-swa_tsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold, Length) %>% filter(Avg_fold > 0)
covstats_17054_swa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-054-swa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold, Length) %>% filter(Avg_fold > 0)
covstats_17054_tsa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-054-tsa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold, Length) %>% filter(Avg_fold > 0)

covstats_17067_swa_tsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-067-swa_tsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold, Length) %>% filter(Avg_fold > 0)
covstats_17067_swa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-067-swa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold, Length) %>% filter(Avg_fold > 0)
covstats_17067_tsa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-067-tsa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold, Length) %>% filter(Avg_fold > 0)


#### Pull out scaffolds of interest ####
# Import list of scaffolds containing r/t RNA genes
scaffolds_rna <- read_tsv("output/scaffolds_rna.tsv") %>%
  pull(scaffold_id)

# Remove scaffolds containing r/t RNA genes in captured and free metagenomes
scaffolds_06126_swa_tsr <- covstats_06126_swa_tsr %>% left_join(map_06126_watercolumn, by = "contig_id") %>% select(-contig_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_06126_swa_bsr <- covstats_06126_swa_bsr %>% left_join(map_06126_watercolumn, by = "contig_id") %>% select(-contig_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_06126_tsa_bsr <- covstats_06126_tsa_bsr %>% left_join(map_06126_top, by = "contig_id") %>% select(-contig_id) %>% filter(!scaffold_id %in% scaffolds_rna)

scaffolds_17054_swa_tsr <- covstats_17054_swa_tsr %>% left_join(map_17054_watercolumn, by = "contig_id") %>% select(-contig_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_17054_swa_bsr <- covstats_17054_swa_bsr %>% left_join(map_17054_watercolumn, by = "contig_id") %>% select(-contig_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_17054_tsa_bsr <- covstats_17054_tsa_bsr %>% left_join(map_17054_top, by = "contig_id") %>% select(-contig_id) %>% filter(!scaffold_id %in% scaffolds_rna)

scaffolds_17067_swa_tsr <- covstats_17067_swa_tsr %>% left_join(map_17067_watercolumn, by = "contig_id") %>% select(-contig_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_17067_swa_bsr <- covstats_17067_swa_bsr %>% left_join(map_17067_watercolumn, by = "contig_id") %>% select(-contig_id) %>% filter(!scaffold_id %in% scaffolds_rna)
scaffolds_17067_tsa_bsr <- covstats_17067_tsa_bsr %>% left_join(map_17067_top, by = "contig_id") %>% select(-contig_id) %>% filter(!scaffold_id %in% scaffolds_rna)

# Clear up space in R environment
rm(list = ls(pattern = "map"))
rm(list = ls(pattern = "covstats"))
rm(scaffolds_rna)


#### Plot scaffold coverage vs. length for captured metagenomes ####
(covlength_06126_swa_tsr <- scaffolds_06126_swa_tsr %>%
   ggplot(aes(x = Length, y = Avg_fold)) +
   geom_hex() +
   scale_y_continuous(trans = "log1p") +
   scale_x_log10(labels = trans_format("log10", scales::math_format(10^.x))) +
   annotation_logticks(sides = "b") +
   theme_classic() %+replace%
   theme(axis.title = element_blank()))

(covlength_06126_swa_bsr <- scaffolds_06126_swa_bsr %>%
    ggplot(aes(x = Length, y = Avg_fold)) +
    geom_hex() +
    scale_y_continuous(trans = "log1p") +
    scale_x_log10(labels = trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    theme_classic() %+replace%
    theme(axis.title = element_blank()))

(covlength_06126_tsa_bsr <- scaffolds_06126_tsa_bsr %>%
    ggplot(aes(x = Length, y = Avg_fold)) +
    geom_hex() +
    scale_y_continuous(trans = "log1p") +
    scale_x_log10(labels = trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    theme_classic() %+replace%
    theme(axis.title = element_blank()))

(covlength_17054_swa_tsr <- scaffolds_17054_swa_tsr %>%
    ggplot(aes(x = Length, y = Avg_fold)) +
    geom_hex() +
    scale_y_continuous(trans = "log1p") +
    scale_x_log10(labels = trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    theme_classic() %+replace%
    theme(axis.title = element_blank()))

(covlength_17054_swa_bsr <- scaffolds_17054_swa_bsr %>%
    ggplot(aes(x = Length, y = Avg_fold)) +
    geom_hex() +
    scale_y_continuous(trans = "log1p") +
    scale_x_log10(labels = trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    theme_classic() %+replace%
    theme(axis.title = element_blank()))

(covlength_17054_tsa_bsr <- scaffolds_17054_tsa_bsr %>%
    ggplot(aes(x = Length, y = Avg_fold)) +
    geom_hex() +
    scale_y_continuous(trans = "log1p") +
    scale_x_log10(labels = trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    theme_classic() %+replace%
    theme(axis.title = element_blank()))

(covlength_17067_swa_tsr <- scaffolds_17067_swa_tsr %>%
    ggplot(aes(x = Length, y = Avg_fold)) +
    geom_hex() +
    scale_y_continuous(trans = "log1p") +
    scale_x_log10(labels = trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    theme_classic() %+replace%
    theme(axis.title = element_blank()))

(covlength_17067_swa_bsr <- scaffolds_17067_swa_bsr %>%
    ggplot(aes(x = Length, y = Avg_fold)) +
    geom_hex() +
    scale_y_continuous(trans = "log1p") +
    scale_x_log10(labels = trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    theme_classic() %+replace%
    theme(axis.title = element_blank()))

(covlength_17067_tsa_bsr <- scaffolds_17067_tsa_bsr %>%
    ggplot(aes(x = Length, y = Avg_fold)) +
    geom_hex() +
    scale_y_continuous(trans = "log1p") +
    scale_x_log10(labels = trans_format("log10", scales::math_format(10^.x))) +
    annotation_logticks(sides = "b") +
    theme_classic() %+replace%
    theme(axis.title = element_blank()))

# Combine separate plots
(covlength_all <- grid.arrange(covlength_06126_swa_tsr,
                               covlength_06126_swa_bsr,
                               covlength_06126_tsa_bsr,
                               
                               covlength_17054_swa_tsr,
                               covlength_17054_swa_bsr,
                               covlength_17054_tsa_bsr,
                               
                               covlength_17067_swa_tsr,
                               covlength_17067_swa_bsr,
                               covlength_17067_tsa_bsr,
                               
                               nrow = 3, ncol = 3,
                               layout_matrix = rbind(c(1, 4, 7),
                                                     c(2, 5, 8),
                                                     c(3, 6, 9)),
                               left = "Average fold",
                               bottom = "Scaffold length"))
#ggsave(filename = "figures/s02_paleocapture_scaffold_coverage_length.pdf", plot = covlength_all, device = "pdf", units = "in", width = 12, height = 10)
