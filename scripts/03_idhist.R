# 03_idhist.R
# Read mapping vs. sequence identity in captured metagenomes

# Set working directory
setwd("paleocapture/")

# Load libraries
library(tidyverse)
library(readxl)
library(scales)
library(gridExtra)

#### Import read tracking data ####
# Import file containing the number of reads in each metagenome
read_tracking <- read_xlsx("data/lp2017_metagenomes_qaqc_read_tracking.xlsx", col_names = TRUE)

nreads_06126_bottom <- read_tracking %>% filter(sample_id == "06-126_bottom") %>% pull(n_reads_R1_p_trimmed)
nreads_06126_top <- read_tracking %>% filter(sample_id == "06-126_top") %>% pull(n_reads_R1_p_trimmed)

nreads_17054_bottom <- read_tracking %>% filter(sample_id == "17-054_bottom") %>% pull(n_reads_R1_p_trimmed)
nreads_17054_top <- read_tracking %>% filter(sample_id == "17-054_top") %>% pull(n_reads_R1_p_trimmed)

nreads_17067_bottom <- read_tracking %>% filter(sample_id == "17-067_bottom") %>% pull(n_reads_R1_p_trimmed)
nreads_17067_top <- read_tracking %>% filter(sample_id == "17-067_top") %>% pull(n_reads_R1_p_trimmed)


#### Import idhist.txt files ####
# Filter reads recruited at percent identity >= 70 %
# Normalize read recruitment by n unassembled reads
idhist_06126_swa_tsr <- read_tsv(file = "data/bbmap/idhist/lp2017_06-126-swa_tsr_70id_idhist.txt", col_names = c("identity", "reads", "bases"), comment = "#") %>% filter(identity >= 70) %>% mutate(reads_normalized = reads/nreads_06126_top * 100, lake = "Lac Paula", metagenome = "swa_tsr")
idhist_06126_swa_bsr <- read_tsv(file = "data/bbmap/idhist/lp2017_06-126-swa_bsr_70id_idhist.txt", col_names = c("identity", "reads", "bases"), comment = "#") %>% filter(identity >= 70) %>% mutate(reads_normalized = reads/nreads_06126_bottom * 100, lake = "Lac Paula", metagenome = "swa_bsr")
idhist_06126_tsa_bsr <- read_tsv(file = "data/bbmap/idhist/lp2017_06-126-tsa_bsr_70id_idhist.txt", col_names = c("identity", "reads", "bases"), comment = "#") %>% filter(identity >= 70) %>% mutate(reads_normalized = reads/nreads_06126_bottom * 100, lake = "Lac Paula", metagenome = "tsa_bsr")

idhist_17054_swa_tsr <- read_tsv(file = "data/bbmap/idhist/lp2017_17-054-swa_tsr_70id_idhist.txt", col_names = c("identity", "reads", "bases"), comment = "#") %>% filter(identity >= 70) %>% mutate(reads_normalized = reads/nreads_17054_top * 100, lake = "Eightmile Lake", metagenome = "swa_tsr")
idhist_17054_swa_bsr <- read_tsv(file = "data/bbmap/idhist/lp2017_17-054-swa_bsr_70id_idhist.txt", col_names = c("identity", "reads", "bases"), comment = "#") %>% filter(identity >= 70) %>% mutate(reads_normalized = reads/nreads_17054_bottom * 100, lake = "Eightmile Lake", metagenome = "swa_bsr")
idhist_17054_tsa_bsr <- read_tsv(file = "data/bbmap/idhist/lp2017_17-054-tsa_bsr_70id_idhist.txt", col_names = c("identity", "reads", "bases"), comment = "#") %>% filter(identity >= 70) %>% mutate(reads_normalized = reads/nreads_17054_bottom * 100, lake = "Eightmile Lake", metagenome = "tsa_bsr")

idhist_17067_swa_tsr <- read_tsv(file = "data/bbmap/idhist/lp2017_17-067-swa_tsr_70id_idhist.txt", col_names = c("identity", "reads", "bases"), comment = "#") %>% filter(identity >= 70) %>% mutate(reads_normalized = reads/nreads_17067_top * 100, lake = "Grand lac Touradi", metagenome = "swa_tsr")
idhist_17067_swa_bsr <- read_tsv(file = "data/bbmap/idhist/lp2017_17-067-swa_bsr_70id_idhist.txt", col_names = c("identity", "reads", "bases"), comment = "#") %>% filter(identity >= 70) %>% mutate(reads_normalized = reads/nreads_17067_bottom * 100, lake = "Grand lac Touradi", metagenome = "swa_bsr")
idhist_17067_tsa_bsr <- read_tsv(file = "data/bbmap/idhist/lp2017_17-067-tsa_bsr_70id_idhist.txt", col_names = c("identity", "reads", "bases"), comment = "#") %>% filter(identity >= 70) %>% mutate(reads_normalized = reads/nreads_17067_bottom * 100, lake = "Grand lac Touradi", metagenome = "tsa_bsr")

idhist_all <- bind_rows(idhist_06126_swa_tsr,
                        idhist_06126_swa_bsr,
                        idhist_06126_tsa_bsr,
                        
                        idhist_17054_swa_tsr,
                        idhist_17054_swa_bsr,
                        idhist_17054_tsa_bsr,
                        
                        idhist_17067_swa_tsr,
                        idhist_17067_swa_bsr,
                        idhist_17067_tsa_bsr)


#### Plot proportion of reads recruited vs. percent identity ####
(idhist_plot_06126_swa_tsr <- idhist_all %>%
   filter(lake == "Lac Paula",
          metagenome == "swa_tsr") %>%
   ggplot(aes(x = identity, y = reads_normalized)) +
   geom_smooth(method = "loess", colour = "#ff9f1a", fill = NA, size = 2) +
   scale_x_reverse() +
   scale_y_continuous(labels = number_format(accuracy = 0.01)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),
         axis.title = element_blank()) +
   coord_cartesian(ylim = c(0,0.27), expand = FALSE))

(idhist_plot_06126_swa_bsr <- idhist_all %>%
    filter(lake == "Lac Paula",
           metagenome == "swa_bsr") %>%
    ggplot(aes(x = identity, y = reads_normalized)) +
    geom_smooth(method = "loess", colour = "#ff9f1a", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()) +
    coord_cartesian(ylim = c(0,0.14), expand = FALSE))

(idhist_plot_06126_tsa_bsr <- idhist_all %>%
    filter(lake == "Lac Paula",
           metagenome == "tsa_bsr") %>%
    ggplot(aes(x = identity, y = reads_normalized)) +
    geom_smooth(method = "loess", colour = "#ff9f1a", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()) +
    coord_cartesian(ylim = c(0,4.3), expand = FALSE))

(idhist_plot_17054_swa_tsr <- idhist_all %>%
    filter(lake == "Eightmile Lake",
           metagenome == "swa_tsr") %>%
    ggplot(aes(x = identity, y = reads_normalized)) +
    geom_smooth(method = "loess", colour = "#18dcff", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()) +
    coord_cartesian(ylim = c(0,0.12), expand = FALSE))

(idhist_plot_17054_swa_bsr <- idhist_all %>%
    filter(lake == "Eightmile Lake",
           metagenome == "swa_bsr") %>%
    ggplot(aes(x = identity, y = reads_normalized)) +
    geom_smooth(method = "loess", colour = "#18dcff", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()) +
    coord_cartesian(ylim = c(0,0.082), expand = FALSE))

(idhist_plot_17054_tsa_bsr <- idhist_all %>%
    filter(lake == "Eightmile Lake",
           metagenome == "tsa_bsr") %>%
    ggplot(aes(x = identity, y = reads_normalized)) +
    geom_smooth(method = "loess", colour = "#18dcff", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()) +
    coord_cartesian(ylim = c(0,0.75), expand = FALSE))

(idhist_plot_17067_swa_tsr <- idhist_all %>%
    filter(lake == "Grand lac Touradi",
           metagenome == "swa_tsr") %>%
    ggplot(aes(x = identity, y = reads_normalized)) +
    geom_smooth(method = "loess", colour = "#7d5fff", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()) +
    coord_cartesian(ylim = c(0,0.189), expand = FALSE))

(idhist_plot_17067_swa_bsr <- idhist_all %>%
    filter(lake == "Grand lac Touradi",
           metagenome == "swa_bsr") %>%
    ggplot(aes(x = identity, y = reads_normalized)) +
    geom_smooth(method = "loess", colour = "#7d5fff", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()) +
    coord_cartesian(ylim = c(0,0.081), expand = FALSE))

(idhist_plot_17067_tsa_bsr <- idhist_all %>%
    filter(lake == "Grand lac Touradi",
           metagenome == "tsa_bsr") %>%
    ggplot(aes(x = identity, y = reads_normalized)) +
    geom_smooth(method = "loess", colour = "#7d5fff", fill = NA, size = 2) +
    scale_x_reverse() +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_blank()) +
    coord_cartesian(ylim = c(0,1), expand = FALSE))


#### Create composite plot of all lakes and metagenomes ####
(idhist_plot_all <- grid.arrange(idhist_plot_06126_swa_tsr,
                                 idhist_plot_06126_swa_bsr,
                                 idhist_plot_06126_tsa_bsr,
                                 
                                 idhist_plot_17054_swa_tsr,
                                 idhist_plot_17054_swa_bsr,
                                 idhist_plot_17054_tsa_bsr,
                                 
                                 idhist_plot_17067_swa_tsr,
                                 idhist_plot_17067_swa_bsr,
                                 idhist_plot_17067_tsa_bsr,
                                 
                                 nrow = 3, ncol = 3,
                                 layout_matrix = rbind(c(1, 4, 7),
                                                       c(2, 5, 8),
                                                       c(3, 6, 9)),
                                 left = "Mapped (%)",
                                 bottom = "Sequence identity (%)"))
#ggsave(filename = "figures/02_paleocapture_idhist.pdf", plot = idhist_plot_all, device = "pdf", units = "in", width = 8, height = 8, scale = 1)
