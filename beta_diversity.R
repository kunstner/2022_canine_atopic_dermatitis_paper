# General Information -----------------------------------------------------

# Authors: Mirja Thomsen 
# Date:    2022-06-13

# Beta diversity ----------------------------------------------------------

# Load libraries ----------------------------------------------------------

library(phyloseq)
library(tidyverse)
library(microbiome)
library(vegan)
library(RVAideMemoire)
library(ggplot2)

# Load data ---------------------------------------------------------------

ps <- readRDS(file = "phyloseq.OTU.RDS")

# Exemplarily for stool samples

ps_stool <- subset_samples(ps, location == "stool")

# Aitchison's distance

ps_clr <- microbiome::transform(ps_stool, "clr")
otu.table_clr <- otu_table(ps_clr) %>% t()
ps_clr_dist <- dist(otu.table_clr, method="euclidean")

# Hypothesis testing

metadata <- data.frame(sample_data(ps_clr))
vegan::adonis2(formula = ps_clr_dist ~ group + sex, data = metadata, permutations = 99999)
RVAideMemoire::pairwise.perm.manova(ps_clr_dist, metadata$group, nperm=99999)


# Constrained analysis of principal coordinates

ps_clr_ord <- phyloseq::ordinate(ps_clr, "CAP", distance = "euclidean", formula = ~group + Condition(sex))

phyloseq::plot_ordination( physeq = ps_clr, ordination = ps_clr_ord, color="group")+ 
    scale_color_manual( values=c( "orange1", "#d7191c", "steelblue1", "steelblue4" ) ) +
    geom_point(size=3) +
    theme_bw()+
    theme(legend.title = element_text(size=14),
          legend.text = element_text(size=14),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14)) +
    labs(color='Group') +
    ggtitle( paste0("CAP of Aitchison distance") )

ggsave(filename = "plots/beta.aitchison.cap.stool.C(Sex).pdf")

