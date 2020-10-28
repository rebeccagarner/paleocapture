# 02_rnagenes.R
# Identify scaffolds containing conserved genes (ribosomal or transfer rRNA genes)
# Count reads mapped (after removing scaffolds containing conserved genes)

# Set working directory
setwd("paleocapture/")

# Load libraries
library(tidyverse)
library(readxl)


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


#### Import .gff file containing scaffold feature information ####
# Join map data to gff
gff_06126_bottom <- read_tsv(file = "data/img/lp2017_06-126-bottom/3300032883.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type) %>% left_join(map_06126_bottom, by = "scaffold_id")
gff_06126_top <- read_tsv(file = "data/img/lp2017_06-126-top/3300032795.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type) %>% left_join(map_06126_top, by = "scaffold_id")
gff_06126_watercolumn <- read_tsv(file = "data/img/lp2017_06-126-watercolumn/3300032269.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type) %>% left_join(map_06126_watercolumn, by = "scaffold_id")

gff_17054_bottom <- read_tsv(file = "data/img/lp2017_17-054-bottom/3300032797.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type) %>% left_join(map_17054_bottom, by = "scaffold_id")
gff_17054_top <- read_tsv(file = "data/img/lp2017_17-054-top/3300032796.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type) %>% left_join(map_17054_top, by = "scaffold_id")
gff_17054_watercolumn <- read_tsv(file = "data/img/lp2017_17-054-watercolumn/3300032328.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type) %>% left_join(map_17054_watercolumn, by = "scaffold_id")

gff_17067_bottom <- read_tsv(file = "data/img/lp2017_17-067-bottom/3300032799.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type) %>% left_join(map_17067_bottom, by = "scaffold_id")
gff_17067_top <- read_tsv(file = "data/img/lp2017_17-067-top/3300032798.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type) %>% left_join(map_17067_top, by = "scaffold_id")
gff_17067_watercolumn <- read_tsv(file = "data/img/lp2017_17-067-watercolumn/3300032317.a.gff", col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes")) %>% select(scaffold_id, type) %>% left_join(map_17067_watercolumn, by = "scaffold_id")

# Combine all gff
gff_all <- bind_rows(gff_06126_bottom,
                     gff_06126_top,
                     gff_06126_watercolumn,
                     
                     gff_17054_bottom,
                     gff_17054_top,
                     gff_17054_watercolumn,
                     
                     gff_17067_bottom,
                     gff_17067_top,
                     gff_17067_watercolumn) %>%
  select(-contig_id)

# Clear up space in R environment
rm(list = ls(pattern = "gff_\\d"))


#### Filter scaffolds containing rRNA or tRNA genes ####
gff_rna <- gff_all %>%
  filter(type == "rRNA" | type == "tRNA") %>%
  distinct(scaffold_id)

# Write a tsv file containing the IDs of scaffolds containing r/t RNA genes
# gff_rna %>%
#   write_tsv("output/scaffolds_rna.tsv")


#### Import covstats of reads recruited at 90% identity and filter contigs with non-zero coverage ####
contigs_06126_swa_tsr <- read_tsv(file = "data/bbmap/covstats/lp2017_06-126-swa_tsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0) %>% left_join(map_06126_watercolumn, by = "contig_id") %>% mutate(lake = "Lac Paula", assembly = "watercolumn", experiment = "swa_tsr")
contigs_06126_swa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_06-126-swa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0) %>% left_join(map_06126_watercolumn, by = "contig_id") %>% mutate(lake = "Lac Paula", assembly = "watercolumn", experiment = "swa_bsr")
contigs_06126_tsa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_06-126-tsa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0) %>% left_join(map_06126_top, by = "contig_id") %>% mutate(lake = "Lac Paula", assembly = "top", experiment = "tsa_bsr")

contigs_17054_swa_tsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-054-swa_tsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0) %>% left_join(map_17054_watercolumn, by = "contig_id") %>% mutate(lake = "Eightmile Lake", assembly = "watercolumn", experiment = "swa_tsr")
contigs_17054_swa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-054-swa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0) %>% left_join(map_17054_watercolumn, by = "contig_id") %>% mutate(lake = "Eightmile Lake", assembly = "watercolumn", experiment = "swa_bsr")
contigs_17054_tsa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-054-tsa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0) %>% left_join(map_17054_top, by = "contig_id") %>% mutate(lake = "Eightmile Lake", assembly = "top", experiment = "tsa_bsr")

contigs_17067_swa_tsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-067-swa_tsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0) %>% left_join(map_17067_watercolumn, by = "contig_id") %>% mutate(lake = "Grand lac Touradi", assembly = "watercolumn", experiment = "swa_tsr")
contigs_17067_swa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-067-swa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0) %>% left_join(map_17067_watercolumn, by = "contig_id") %>% mutate(lake = "Grand lac Touradi", assembly = "watercolumn", experiment = "swa_bsr")
contigs_17067_tsa_bsr <- read_tsv(file = "data/bbmap/covstats/lp2017_17-067-tsa_bsr_90id_covstats.txt") %>% mutate(contig_id = str_replace(`#ID`, " flag.*", "")) %>% select(contig_id, Avg_fold) %>% filter(Avg_fold > 0) %>% left_join(map_17067_top, by = "contig_id") %>% mutate(lake = "Grand lac Touradi", assembly = "top", experiment = "tsa_bsr")

contigs_all <- bind_rows(contigs_06126_swa_tsr,
                         contigs_06126_swa_bsr,
                         contigs_06126_tsa_bsr,
                         
                         contigs_17054_swa_tsr,
                         contigs_17054_swa_bsr,
                         contigs_17054_tsa_bsr,
                         
                         contigs_17067_swa_tsr,
                         contigs_17067_swa_bsr,
                         contigs_17067_tsa_bsr)

# Clear up space in R environment
rm(list = ls(pattern = "map_\\d"))
rm(list = ls(pattern = "contigs_\\d"))


#### Import idxstats files summarizing n reads mapped to contigs ####
idxstats_06126_swa_tsr <- read_tsv("data/bbmap/idxstats/lp2017_06-126-swa_tsr_90id_mapped_idxstats.txt", col_names = c("reference_seq", "contig_length", "n_mapped", "n_unmapped")) %>% mutate(lake = "Lac Paula", assembly = "watercolumn", experiment = "swa_tsr") %>% separate(reference_seq, into = c("contig_id", "flag", "multi", "len"), sep = " ") %>% select(-flag, -multi, -len)
idxstats_06126_swa_bsr <- read_tsv("data/bbmap/idxstats/lp2017_06-126-swa_bsr_90id_mapped_idxstats.txt", col_names = c("reference_seq", "contig_length", "n_mapped", "n_unmapped")) %>% mutate(lake = "Lac Paula", assembly = "watercolumn", experiment = "swa_bsr") %>% separate(reference_seq, into = c("contig_id", "flag", "multi", "len"), sep = " ") %>% select(-flag, -multi, -len)
idxstats_06126_tsa_bsr <- read_tsv("data/bbmap/idxstats/lp2017_06-126-tsa_bsr_90id_mapped_idxstats.txt", col_names = c("reference_seq", "contig_length", "n_mapped", "n_unmapped")) %>% mutate(lake = "Lac Paula", assembly = "top", experiment = "tsa_bsr") %>% separate(reference_seq, into = c("contig_id", "flag", "multi", "len"), sep = " ") %>% select(-flag, -multi, -len)

idxstats_17054_swa_tsr <- read_tsv("data/bbmap/idxstats/lp2017_17-054-swa_tsr_90id_mapped_idxstats.txt", col_names = c("reference_seq", "contig_length", "n_mapped", "n_unmapped")) %>% mutate(lake = "Eightmile Lake", assembly = "watercolumn", experiment = "swa_tsr") %>% separate(reference_seq, into = c("contig_id", "flag", "multi", "len"), sep = " ") %>% select(-flag, -multi, -len)
idxstats_17054_swa_bsr <- read_tsv("data/bbmap/idxstats/lp2017_17-054-swa_bsr_90id_mapped_idxstats.txt", col_names = c("reference_seq", "contig_length", "n_mapped", "n_unmapped")) %>% mutate(lake = "Eightmile Lake", assembly = "watercolumn", experiment = "swa_bsr") %>% separate(reference_seq, into = c("contig_id", "flag", "multi", "len"), sep = " ") %>% select(-flag, -multi, -len)
idxstats_17054_tsa_bsr <- read_tsv("data/bbmap/idxstats/lp2017_17-054-tsa_bsr_90id_mapped_idxstats.txt", col_names = c("reference_seq", "contig_length", "n_mapped", "n_unmapped")) %>% mutate(lake = "Eightmile Lake", assembly = "top", experiment = "tsa_bsr") %>% separate(reference_seq, into = c("contig_id", "flag", "multi", "len"), sep = " ") %>% select(-flag, -multi, -len)

idxstats_17067_swa_tsr <- read_tsv("data/bbmap/idxstats/lp2017_17-067-swa_tsr_90id_mapped_idxstats.txt", col_names = c("reference_seq", "contig_length", "n_mapped", "n_unmapped")) %>% mutate(lake = "Grand lac Touradi", assembly = "watercolumn", experiment = "swa_tsr") %>% separate(reference_seq, into = c("contig_id", "flag", "multi", "len"), sep = " ") %>% select(-flag, -multi, -len)
idxstats_17067_swa_bsr <- read_tsv("data/bbmap/idxstats/lp2017_17-067-swa_bsr_90id_mapped_idxstats.txt", col_names = c("reference_seq", "contig_length", "n_mapped", "n_unmapped")) %>% mutate(lake = "Grand lac Touradi", assembly = "watercolumn", experiment = "swa_bsr") %>% separate(reference_seq, into = c("contig_id", "flag", "multi", "len"), sep = " ") %>% select(-flag, -multi, -len)
idxstats_17067_tsa_bsr <- read_tsv("data/bbmap/idxstats/lp2017_17-067-tsa_bsr_90id_mapped_idxstats.txt", col_names = c("reference_seq", "contig_length", "n_mapped", "n_unmapped")) %>% mutate(lake = "Grand lac Touradi", assembly = "top", experiment = "tsa_bsr") %>% separate(reference_seq, into = c("contig_id", "flag", "multi", "len"), sep = " ") %>% select(-flag, -multi, -len)

idxstats_all <- bind_rows(idxstats_06126_swa_tsr,
                          idxstats_06126_swa_bsr,
                          idxstats_06126_tsa_bsr,
                          
                          idxstats_17054_swa_tsr,
                          idxstats_17054_swa_bsr,
                          idxstats_17054_tsa_bsr,
                          
                          idxstats_17067_swa_tsr,
                          idxstats_17067_swa_bsr,
                          idxstats_17067_tsa_bsr)

# Clear up space in R environment
rm(list = ls(pattern = "idxstats_\\d"))

#### Join contigs and idxstats ####
contigs_idxstats <- contigs_all %>%
  left_join(idxstats_all, by = c("lake", "assembly", "experiment", "contig_id"))


#### Examine scaffolds containing versus not containing  r/t RNA genes ####
# Vectorize list of scaffolds containing r/t RNA genes
scaffolds_rna <- gff_rna %>%
  pull(scaffold_id)

contigs_idxstats_rmrna <- contigs_idxstats %>%
  filter(!scaffold_id %in% scaffolds_rna)


#### Summarize n mapped scaffolds ####
# n mapped scaffolds including scaffolds containing r/t RNA genes
(nscaffolds <- contigs_idxstats %>%
   mutate(scaffold_count = 1) %>%
   group_by(lake, experiment) %>%
   summarize(n_scaffolds = sum(scaffold_count)))

# n mapped scaffolds after removing scaffolds containing r/t RNA genes
(nscaffolds_rmrna <- contigs_idxstats_rmrna %>%
    mutate(scaffold_count = 1) %>%
    group_by(lake, experiment) %>%
    summarize(n_scaffolds_rmrna = sum(scaffold_count)))


#### Summarize n mapped reads ####
# n mapped reads INCLUDING reads mapped to scaffolds containing r/t RNA genes
(nmapped <- contigs_idxstats %>%
   group_by(lake, experiment) %>%
   summarize(sum_n_mapped = sum(n_mapped)))

# n mapped reads EXCLUDING reads mapped to scaffolds containing r/t RNA genes
(nmapped_rmrna <- contigs_idxstats_rmrna %>%
    group_by(lake, experiment) %>%
    summarize(sum_n_mapped_rmrna = sum(n_mapped)))


#### Import read tracking data ####
read_tracking <- read_xlsx("data/lp2017_metagenomes_qaqc_read_tracking.xlsx", col_names = TRUE) %>%
  mutate(unassembled_metagenome = case_when(env_fraction == "surface water" ~ "swr",
                                            env_fraction == "top sediment" ~ "tsr",
                                            env_fraction == "bottom sediment" ~ "bsr")) %>%
  select(lake, unassembled_metagenome, n_reads_R1_p_trimmed)


#### Summarize mapping stats ####
(scaffolds_mapped_summary <- nscaffolds %>%
   left_join(nscaffolds_rmrna, by = c("lake", "experiment")) %>%
   left_join(nmapped, by = c("lake", "experiment")) %>%
   left_join(nmapped_rmrna, by = c("lake", "experiment")) %>%
   arrange(factor(lake, levels = c("Lac Paula", "Eightmile Lake", "Grand lac Touradi")),
           factor(experiment, levels = c("swa_tsr", "swa_bsr", "tsa_bsr"))) %>%
   separate(experiment, into = c("assembled_metagenome", "unassembled_metagenome"), sep = "_", remove = FALSE) %>%
   left_join(read_tracking, by = c("lake", "unassembled_metagenome")) %>%
   mutate(percent_mapped_rmrna = round(sum_n_mapped_rmrna/n_reads_R1_p_trimmed * 100, 2)))
#write_csv(scaffolds_mapped_summary, "output/scaffolds_mapped.csv")
