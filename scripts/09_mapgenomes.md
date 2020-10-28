---
title: "09_mapgenomes.Rmd"
---

Map metagenome reads to reference reference genomes accessed from IMG.

## Mask RNA genes in reference genomes

1. Download reference genomes and annotations from IMG genome cart

2. Extract useful files (.fna and .gff) from downloaded folders:
- Copy .fna files to new folder with command line: cp */*.fna polynucleobacter_genomes/
- Copy .gff files to new folder with command line: cp */*.gff polynucleobacter_gff/

3. Identify RNA genes and their coordinates in reference genomes

Create .bed files of gene coordinates on scaffolds

```r
setwd("paleocapture/")

# Load libraries
library(tidyverse)
library(janitor)

# Import gff files
file_path <- "polynucleobacter/polynucleobacter_gff/"

(gff_files <- dir(path = file_path, pattern = "*.gff"))

# Create a list of gff file contents
gff <- tibble(file_name = gff_files) %>%
  mutate(genome = str_remove_all(file_name, ".gff")) %>%
  mutate(file_contents = map(file_name, ~ read_tsv(file.path(file_path, .),
                                                   col_names = c("scaffold_id", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "attributes"),
                                                   skip = 1))) %>%
  unnest(c(file_contents))

# Inspect unique gff types
unique(gff$type)

# Filter rRNA and tRNA genes
rna <- gff %>%
  select(-file_name) %>%
  filter(type == "rRNA" | type == "tRNA") %>%
  group_split(genome)

# Name lists by genome
rna %>%
  purrr::map(~pull(.,genome)) %>%  # Pull out genome variable
  purrr::map(~as.character(.)) %>%  # Convert factor to character
  purrr::map(~unique(.)) -> names(rna)

# Write .bed file containing gene coordinates on scaffolds
# Select desired columns
rna <- rna %>%
  purrr::map(~select(.,scaffold_id, start_coord, end_coord))

# Define function for writing files with the desired file names and into the right path
write_bed <- function(data, names){
  dir_path <- "polynucleobacter/polynucleobacter_rna/"
  write_tsv(data, paste0(dir_path, names, "_rna.bed"), col_names = FALSE)
}

Generate individual .bed files for each genome (RNA)
list(data = rna, names = names(rna)) %>%
  purrr::pmap(write_bed)
```

4. Mask DNA sequences based on bed file with bedtools maskfasta

Convert RNA gene sequences to N characters.

```shell
cd polynucleobacter_genomes/
for filename in *.fna
do
  prefix=`basename $filename .fna`
  bedtools maskfasta -fi $filename -bed polynucleobacter_rna/${prefix}_rna.bed -fo polynucleobacter_maskrna/${prefix}_maskrna.fasta
done
```

5. Concatenate rRNA and tRNA gene-masked sequence files

```shell
# 1. Polynucleobacter
cd polynucleobacter_maskrna/
cat *.fasta > burkholderiales_maskrna/polynucleobacter_cat_maskrna.fasta

# 2 Limnohabitans
cd limnohabitans_maskrna/
cat *.fasta > burkholderiales_maskrna/limnohabitans_cat_maskrna.fasta
```

## Map unassembled metagenomes to reference Burkholderiales genomes with BBMap

```shell
# There are separate concatenated reference genomes for Polynucleobacter and Limnohabitans
# Build references for individual genomes
cd burkholderiales_maskrna/
for filename in *_maskrna.fasta
do
	# Define prefix by removing _maskrna.fasta suffix from genome fasta files
	prefix=`basename $filename _maskrna.fasta`

	# Build genome reference with BBMap
	bbmap.sh ref=$filename path=burkholderiales_ref/${prefix}_ref
done


# Map unassembled metagenome to individual reference genomes (at 70% identity threshold)
cd burkholderiales_ref/
for filename in *_ref
do
	# cd into _ref directory
	cd burkholderiales_ref/${filename}
	
	# Define prefix by removing _ref from directory names
	prefix=`basename $filename _ref`
	
	# Recruit reads with BBMap
	bbmap.sh maxlen=500 ambiguous=best mappedonly=t idtag=t pigz=t minid=0.70 covstats=burkholderiales_bbmap/${prefix}_recruit_lacpaula-watercolumn_reads_70id_covstats.txt idhist=burkholderiales_bbmap/${prefix}_recruit_lacpaula-watercolumn_reads_70id_idhist.txt in=lacpaula-watercolumn_R1_p_trimmed.fastq.gz in2=lacpaula-watercolumn_R2_p_trimmed.fastq.gz outm=burkholderiales_bbmap/${prefix}_recruit_lacpaula-watercolumn_reads_70id.sam
	bbmap.sh maxlen=500 ambiguous=best mappedonly=t idtag=t pigz=t minid=0.70 covstats=burkholderiales_bbmap/${prefix}_recruit_lacpaula-top_reads_70id_covstats.txt idhist=burkholderiales_bbmap/${prefix}_recruit_lacpaula-top_reads_70id_idhist.txt in=lacpaula-top_R1_p_trimmed.fastq.gz in2=lacpaula-top_R2_p_trimmed.fastq.gz outm=burkholderiales_bbmap/${prefix}_recruit_lacpaula-top_reads_70id.sam
	bbmap.sh maxlen=500 ambiguous=best mappedonly=t idtag=t pigz=t minid=0.70 covstats=burkholderiales_bbmap/${prefix}_recruit_lacpaula-bottom_reads_70id_covstats.txt idhist=burkholderiales_bbmap/${prefix}_recruit_lacpaula-bottom_reads_70id_idhist.txt in=lacpaula-bottom_R1_p_trimmed.fastq.gz in2=lacpaula-bottom_R2_p_trimmed.fastq.gz outm=burkholderiales_bbmap/${prefix}_recruit_lacpaula-bottom_reads_70id.sam
	bbmap.sh maxlen=500 ambiguous=best mappedonly=t idtag=t pigz=t minid=0.70 covstats=burkholderiales_bbmap/${prefix}_recruit_eightmilelake-watercolumn_reads_70id_covstats.txt idhist=burkholderiales_bbmap/${prefix}_recruit_eightmilelake-watercolumn_reads_70id_idhist.txt in=eightmilelake-watercolumn_R1_p_trimmed.fastq.gz in2=eightmilelake-watercolumn_R2_p_trimmed.fastq.gz outm=burkholderiales_bbmap/${prefix}_recruit_eightmilelake-watercolumn_reads_70id.sam
	bbmap.sh maxlen=500 ambiguous=best mappedonly=t idtag=t pigz=t minid=0.70 covstats=burkholderiales_bbmap/${prefix}_recruit_eightmilelake-top_reads_70id_covstats.txt idhist=burkholderiales_bbmap/${prefix}_recruit_eightmilelake-top_reads_70id_idhist.txt in=eightmilelake-top_R1_p_trimmed.fastq.gz in2=eightmilelake-top_R2_p_trimmed.fastq.gz outm=burkholderiales_bbmap/${prefix}_recruit_eightmilelake-top_reads_70id.sam
	bbmap.sh maxlen=500 ambiguous=best mappedonly=t idtag=t pigz=t minid=0.70 covstats=burkholderiales_bbmap/${prefix}_recruit_eightmilelake-bottom_reads_70id_covstats.txt idhist=burkholderiales_bbmap/${prefix}_recruit_eightmilelake-bottom_reads_70id_idhist.txt in=eightmilelake-bottom_R1_p_trimmed.fastq.gz in2=eightmilelake-bottom_R2_p_trimmed.fastq.gz outm=burkholderiales_bbmap/${prefix}_recruit_eightmilelake-bottom_reads_70id.sam
	bbmap.sh maxlen=500 ambiguous=best mappedonly=t idtag=t pigz=t minid=0.70 covstats=burkholderiales_bbmap/${prefix}_recruit_grandlactouradi-watercolumn_reads_70id_covstats.txt idhist=burkholderiales_bbmap/${prefix}_recruit_grandlactouradi-watercolumn_reads_70id_idhist.txt in=grandlactouradi-watercolumn_R1_p_trimmed.fastq.gz in2=grandlactouradi-watercolumn_R2_p_trimmed.fastq.gz outm=burkholderiales_bbmap/${prefix}_recruit_grandlactouradi-watercolumn_reads_70id.sam
	bbmap.sh maxlen=500 ambiguous=best mappedonly=t idtag=t pigz=t minid=0.70 covstats=burkholderiales_bbmap/${prefix}_recruit_grandlactouradi-top_reads_70id_covstats.txt idhist=burkholderiales_bbmap/${prefix}_recruit_grandlactouradi-top_reads_70id_idhist.txt in=grandlactouradi-top_R1_p_trimmed.fastq.gz in2=grandlactouradi-top_R2_p_trimmed.fastq.gz outm=burkholderiales_bbmap/${prefix}_recruit_grandlactouradi-top_reads_70id.sam
	bbmap.sh maxlen=500 ambiguous=best mappedonly=t idtag=t pigz=t minid=0.70 covstats=burkholderiales_bbmap/${prefix}_recruit_grandlactouradi-bottom_reads_70id_covstats.txt idhist=burkholderiales_bbmap/${prefix}_recruit_grandlactouradi-bottom_reads_70id_idhist.txt in=grandlactouradi-bottom_R1_p_trimmed.fastq.gz in2=grandlactouradi-bottom_R2_p_trimmed.fastq.gz outm=burkholderiales_bbmap/${prefix}_recruit_grandlactouradi-bottom_reads_70id.sam
	cd burkholderiales_bbmap/
	pigz *
done
```

## Visualize metagenome read mapping to Burkholderiales reference genomes

```r
setwd("paleocapture/")

# Load libraries
library(tidyverse)
library(readxl)
library(scales)


#### Import idhist.txt files ####
file_path <- "burkholderiales/burkholderiales_bbmap/"

(idhist_files <- dir(path = file_path, pattern = "*idhist.txt"))


# Create a list of gff file contents
idhist <- tibble(file_name = idhist_files) %>%
  mutate(genome = str_remove(file_name, "_.*"),
         lake_id = str_extract(file_name, "\\d\\d\\-\\d\\d\\d"),
         metagenome = case_when(grepl("watercolumn", file_name) ~ "surface water",
                                grepl("top", file_name) ~ "top sediment",
                                grepl("bottom", file_name) ~ "bottom sediment")) %>%
  mutate(lake = case_when(lake_id == "06-126" ~ "Lac Paula",
                          lake_id == "17-054" ~ "Eightmile Lake",
                          lake_id == "17-067" ~ "Grand lac Touradi")) %>%
  mutate(file_contents = map(file_name, ~ read_tsv(file.path(file_path, .),
                                                   col_names = c("identity", "reads", "bases"),
                                                   comment = "#"))) %>%
  unnest(c(file_contents)) %>%
  filter(identity >= 70)


#### Normalize n reads mapped based on unassembled metagenome size ####
# Import read tracking data
read_tracking <- read_xlsx("data/lp2017_metagenomes_qaqc_read_tracking.xlsx", col_names = TRUE)

# Normalize reads mapped
idhist$reads_normalized <- NA
for (i in 1:nrow(idhist)) {
  idhist$reads_normalized[i] <- idhist$reads[i]/read_tracking$n_reads_R1_p_trimmed[which(read_tracking$lake == idhist$lake[i] & read_tracking$env_fraction == idhist$metagenome[i])] * 100
}

# Define function for generating scatter plots
palette_lakes <- c("#ff9f1a", "#18dcff", "#7d5fff")
names(palette_lakes) <- c("Lac Paula", "Eightmile Lake", "Grand lac Touradi")


#### Plot proportion of reads recruited vs. percent identity ####
(idhist_plot_06126_swa_refs <- idhist %>%
   filter(lake == "Lac Paula",
          metagenome == "surface water") %>%
   ggplot(aes(x = identity, y = reads_normalized, colour = lake, linetype = genome)) +
   geom_smooth(method = "loess", colour = "#ff9f1a", fill = NA, size = 2) +
   scale_x_reverse() +
   scale_y_continuous(labels = number_format(accuracy = 0.01)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),
         axis.title = element_blank(), legend.position = "none") +
   coord_cartesian(ylim = c(0,0.25), expand = FALSE))

(idhist_plot_06126_tsa_refs <- idhist %>%
    filter(lake == "Lac Paula",
           metagenome == "top sediment") %>%
    ggplot(aes(x = identity, y = reads_normalized, colour = lake, linetype = genome)) +
    geom_smooth(method = "loess", colour = "#ff9f1a", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank(), legend.position = "none") +
    coord_cartesian(ylim = c(0,0.105), expand = FALSE))

(idhist_plot_06126_bsa_refs <- idhist %>%
    filter(lake == "Lac Paula",
           metagenome == "bottom sediment") %>%
    ggplot(aes(x = identity, y = reads_normalized, colour = lake, linetype = genome)) +
    geom_smooth(method = "loess", colour = "#ff9f1a", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank(), legend.position = "none") +
    coord_cartesian(ylim = c(0,0.037), expand = FALSE))

(idhist_plot_17054_swa_refs <- idhist %>%
    filter(lake == "Eightmile Lake",
           metagenome == "surface water") %>%
    ggplot(aes(x = identity, y = reads_normalized, colour = lake, linetype = genome)) +
    geom_smooth(method = "loess", colour = "#18dcff", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank(), legend.position = "none") +
    coord_cartesian(ylim = c(0,0.33), expand = FALSE))

(idhist_plot_17054_tsa_refs <- idhist %>%
    filter(lake == "Eightmile Lake",
           metagenome == "top sediment") %>%
    ggplot(aes(x = identity, y = reads_normalized, colour = lake, linetype = genome)) +
    geom_smooth(method = "loess", colour = "#18dcff", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank(), legend.position = "none") +
    coord_cartesian(ylim = c(0,0.03), expand = FALSE))

(idhist_plot_17054_bsa_refs <- idhist %>%
    filter(lake == "Eightmile Lake",
           metagenome == "bottom sediment") %>%
    ggplot(aes(x = identity, y = reads_normalized, colour = lake, linetype = genome)) +
    geom_smooth(method = "loess", colour = "#18dcff", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank(), legend.position = "none") +
    coord_cartesian(ylim = c(0,0.009), expand = FALSE))

(idhist_plot_17067_swa_refs <- idhist %>%
    filter(lake == "Grand lac Touradi",
           metagenome == "surface water") %>%
    ggplot(aes(x = identity, y = reads_normalized, colour = lake, linetype = genome)) +
    geom_smooth(method = "loess", colour = "#7d5fff", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank(), legend.position = "none") +
    coord_cartesian(ylim = c(0,0.25), expand = FALSE))

(idhist_plot_17067_tsa_refs <- idhist %>%
    filter(lake == "Grand lac Touradi",
           metagenome == "top sediment") %>%
    ggplot(aes(x = identity, y = reads_normalized, colour = lake, linetype = genome)) +
    geom_smooth(method = "loess", colour = "#7d5fff", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank(), legend.position = "none") +
    coord_cartesian(ylim = c(0,0.085), expand = FALSE))

(idhist_plot_17067_bsa_refs <- idhist %>%
    filter(lake == "Grand lac Touradi",
           metagenome == "bottom sediment") %>%
    ggplot(aes(x = identity, y = reads_normalized, colour = lake, linetype = genome)) +
    geom_smooth(method = "loess", colour = "#7d5fff", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank(), legend.position = "none") +
    coord_cartesian(ylim = c(0,0.026), expand = FALSE))


#### Create composite plot of all lakes and metagenomes ####
(idhist_plot_all <- grid.arrange(idhist_plot_06126_swa_refs,
                                 idhist_plot_06126_tsa_refs,
                                 idhist_plot_06126_bsa_refs,
                                 
                                 idhist_plot_17054_swa_refs,
                                 idhist_plot_17054_tsa_refs,
                                 idhist_plot_17054_bsa_refs,
                                 
                                 idhist_plot_17067_swa_refs,
                                 idhist_plot_17067_tsa_refs,
                                 idhist_plot_17067_bsa_refs,
                                 
                                 nrow = 3, ncol = 3,
                                 layout_matrix = rbind(c(1, 4, 7),
                                                       c(2, 5, 8),
                                                       c(3, 6, 9)),
                                 left = "Mapped (%)",
                                 bottom = "Sequence identity (%)"))
#ggsave(filename = "figures/figs05_paleocapture_idhist_refgenomes.pdf", plot = idhist_plot_all, device = "pdf", units = "in", width = 8, height = 8, scale = 1)
```
