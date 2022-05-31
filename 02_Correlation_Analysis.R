
# General Information -----------------------------------------------------

# Authors: Axel Künstner
# Date:    2022-05-31

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(phyloseq)

seedID <- 1253
set.seed(seedID) 

# Get data ----------------------------------------------------------------

ps <- readRDS(file = "data/phyloseq.OTU.RDS")

# selbal selected sites
ps <- subset_samples(physeq = ps, location %in% c("abdomen", "frontPaw", "lips", "stool", "tail"))

sample_data(ps) %>% colnames
sample_data(ps) %>% .$location %>% unique %>% sort
sample_data(ps) %>% .$group %>% unique

ps_ra  <- transform_sample_counts(ps, function(x){x / sum(x)})
ps_clr <- microbiome::transform(ps, "clr")

metadata <- sample_data(ps) %>% data.frame()

# Correlations PVAS -------------------------------------------------------

# Phylum
corr_df <- phylosmith::variable_correlation( phyloseq_obj = ps_clr %>% subset_samples(., !is.na(PVAS) & location == 'stool'),
                                             treatment = c('location', 'group'), 
                                             variables = c('PVAS'), # 'CADESI', 'CADESIall'),
                                             classification = "Phylum", method = 'spearman')

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

plot_corr_phyla <- corr_df %>%
    tidyr::separate(col = Treatment, into = c('Location', 'Group'), remove = FALSE) %>% 
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::mutate(X = gsub(pattern = "Candidatus", replacement = "Cand.", X)) %>%
    dplyr::mutate(Location = str_to_title(Location)) %>% 
    dplyr::mutate(Location = gsub(pattern = "Frontpaw", replacement = "Front paw", x = Location)) %>% 
    dplyr::mutate(Location = factor(x = Location, 
                                    levels = c("Abdomen", "Front paw", "Lips", "Tail", "Stool"))) %>%
    dplyr::mutate(Group = gsub(pattern = "AD", replacement = "AD ", x = Group)) %>% 
    dplyr::mutate(Group = factor(x = Group, 
                                 levels = c("AD pre", "AD post", "AD treatment"))) %>%
    ggplot(data = ., aes(x=Group, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Group, X,
                  label = p_clean),
              color = "black", size = 4) +
    facet_wrap(Location~., scales = "free_x", ncol = 5) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("") + ggtitle("Correlation with PVAS")
plot_corr_phyla

# Genera
( ps.sel <- phyloseq::tax_glom(physeq = ps_clr %>% subset_samples(., !is.na(PVAS) & location == "stool"),
                     taxrank = "Genus") %>%
        phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/15)
        )

corr_df <- phylosmith::variable_correlation( phyloseq_obj = ps.sel,
                                             treatment = c('location', 'group'), 
                                             variables = c('PVAS'), # 'CADESI', 'CADESIall'),
                                             classification = "Genus", method = 'spearman')

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""
corr_df <- corr_df[corr_df$X != "Unknown", ]
# keep only genera with p < 0.05
corr_genera <- corr_df$X[corr_df$p < 0.05]

plot_corr_genera <- corr_df %>%
    dplyr::filter(X %in% corr_genera) %>%
    tidyr::separate(col = Treatment, into = c('Location', 'Group'), remove = FALSE) %>% 
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::mutate(Location = str_to_title(Location)) %>% 
    dplyr::mutate(Location = gsub(pattern = "Frontpaw", replacement = "Front paw", x = Location)) %>% 
    dplyr::mutate(Location = factor(x = Location, 
                                    levels = c("Abdomen", "Front paw", "Lips", "Tail", "Stool"))) %>%
    dplyr::mutate(Group = gsub(pattern = "AD", replacement = "AD ", x = Group)) %>% 
    dplyr::mutate(Group = factor(x = Group, 
                                 levels = c("AD pre", "AD post", "AD treatment"))) %>%
    ggplot(data = ., aes(x=Group, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Group, X,
                  label = p_clean),
              color = "black", size = 4) +
    facet_wrap(Location~., scales = "free_x", ncol = 5) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("") + ggtitle("Correlation with PVAS")
plot_corr_genera

plot_corr_phyla
ggsave(filename = paste0( "plots/Correlations_PVAS_stool_Phylum.pdf"), width = 5, height = 10, plot = plot_corr_phyla)
plot_corr_genera
ggsave(filename = paste0( "plots/Correlations_PVAS_stool_Genus.pdf"), width = 5, height = 10, plot = plot_corr_genera)

# Correlations stool CADESI total -----------------------------------------

# Phylum
corr_df <- phylosmith::variable_correlation( phyloseq_obj = ps_clr %>% subset_samples(., !is.na(PVAS) & location == 'stool'),
                                             treatment = c('location', 'group'), 
                                             variables = c('CADESIall'), # 'CADESI', 'CADESIall'),
                                             classification = "Phylum", method = 'spearman')

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

plot_corr_phyla <- corr_df %>%
    tidyr::separate(col = Treatment, into = c('Location', 'Group'), remove = FALSE) %>% 
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::mutate(X = gsub(pattern = "Candidatus", replacement = "Cand.", X)) %>%
    dplyr::mutate(Location = str_to_title(Location)) %>% 
    dplyr::mutate(Location = gsub(pattern = "Frontpaw", replacement = "Front paw", x = Location)) %>% 
    dplyr::mutate(Location = factor(x = Location, 
                                    levels = c("Abdomen", "Front paw", "Lips", "Tail", "Stool"))) %>%
    dplyr::mutate(Group = gsub(pattern = "AD", replacement = "AD ", x = Group)) %>% 
    dplyr::mutate(Group = factor(x = Group, 
                                 levels = c("AD pre", "AD post", "AD treatment"))) %>%
    ggplot(data = ., aes(x=Group, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Group, X,
                  label = p_clean),
              color = "black", size = 4) +
    facet_wrap(Location~., scales = "free_x", ncol = 5) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("") + ggtitle("Correlation with total CADESI")
plot_corr_phyla

# Genera
( ps.sel <- phyloseq::tax_glom(physeq = ps_clr %>% subset_samples(., !is.na(PVAS) & location == "stool"),
                               taxrank = "Genus") %>%
        phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/15)
)

corr_df <- phylosmith::variable_correlation( phyloseq_obj = ps.sel,
                                             treatment = c('location', 'group'), 
                                             variables = c('CADESIall'), # 'CADESI', 'CADESIall'),
                                             classification = "Genus", method = 'spearman')

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""
corr_df <- corr_df[corr_df$X != "Unknown", ]
# keep only genera with p < 0.05
corr_genera <- corr_df$X[corr_df$p < 0.05]

plot_corr_genera <- corr_df %>%
    dplyr::filter(X %in% corr_genera) %>%
    tidyr::separate(col = Treatment, into = c('Location', 'Group'), remove = FALSE) %>% 
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::mutate(Location = str_to_title(Location)) %>% 
    dplyr::mutate(Location = gsub(pattern = "Frontpaw", replacement = "Front paw", x = Location)) %>% 
    dplyr::mutate(Location = factor(x = Location, 
                                    levels = c("Abdomen", "Front paw", "Lips", "Tail", "Stool"))) %>%
    dplyr::mutate(Group = gsub(pattern = "AD", replacement = "AD ", x = Group)) %>% 
    dplyr::mutate(Group = factor(x = Group, 
                                 levels = c("AD pre", "AD post", "AD treatment"))) %>%
    ggplot(data = ., aes(x=Group, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Group, X,
                  label = p_clean),
              color = "black", size = 4) +
    facet_wrap(Location~., scales = "free_x", ncol = 5) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("") + ggtitle("Correlation with total CADESI")
plot_corr_genera

plot_corr_phyla
ggsave(filename = paste0( "plots/Correlations_CADESI_stool_Phylum.pdf"), width = 5, height = 10, plot = plot_corr_phyla)
plot_corr_genera
ggsave(filename = paste0( "plots/Correlations_CADESI_stool_Genus.pdf"), width = 5, height = 10, plot = plot_corr_genera)

# Correlations CADESI -------------------------------------------------------

# Phylum
corr_df <- phylosmith::variable_correlation( phyloseq_obj = ps_clr %>% subset_samples(., !is.na(CADESI) & location != "stool"),
                                             treatment = c('location', 'group'), 
                                             variables = c('CADESI'), # 'CADESI', 'CADESIall'),
                                             classification = "Phylum", method = 'spearman')

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

plot_corr_phyla <- corr_df %>%
    tidyr::separate(col = Treatment, into = c('Location', 'Group'), remove = FALSE) %>% 
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::mutate(X = gsub(pattern = "Candidatus", replacement = "Cand.", X)) %>%
    dplyr::mutate(Location = str_to_title(Location)) %>% 
    dplyr::mutate(Location = gsub(pattern = "Frontpaw", replacement = "Front paw", x = Location)) %>% 
    dplyr::mutate(Location = gsub(pattern = "Lips", replacement = "Penlabial area", x = Location)) %>% 
    dplyr::mutate(Location = gsub(pattern = "Tail", replacement = "Ventral tail", x = Location)) %>% 
    dplyr::mutate(Location = factor(x = Location, 
                                    levels = c("Abdomen", "Front paw", "Penlabial area", "Ventral tail", "Stool"))) %>%
    dplyr::mutate(Group = gsub(pattern = "AD", replacement = "AD ", x = Group)) %>% 
    dplyr::mutate(Group = factor(x = Group, 
                                 levels = c("AD pre", "AD post", "AD treatment"))) %>%
    ggplot(data = ., aes(x=Group, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Group, X,
                  label = p_clean),
              color = "black", size = 4) +
    facet_wrap(Location~., scales = "free_x", ncol = 5) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("") + ggtitle("Correlation with CADESI")
plot_corr_phyla

# Genera
( ps.sel <- phyloseq::tax_glom(physeq = ps_clr %>% subset_samples(., !is.na(CADESI) & location != "stool"),
                               taxrank = "Genus") %>%
        phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/15)
)

corr_df <- phylosmith::variable_correlation( phyloseq_obj = ps.sel,
                                             treatment = c('location', 'group'), 
                                             variables = c('CADESI'), # 'CADESI', 'CADESIall'),
                                             classification = "Genus", method = 'spearman')

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""
corr_df <- corr_df[corr_df$X != "Unknown", ]
# keep only genera with p < 0.05
corr_genera <- corr_df$X[corr_df$p < 0.05]

plot_corr_genera <- corr_df %>%
    dplyr::filter(X %in% corr_genera) %>%
    tidyr::separate(col = Treatment, into = c('Location', 'Group'), remove = FALSE) %>% 
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::mutate(Location = str_to_title(Location)) %>% 
    dplyr::mutate(Location = gsub(pattern = "Frontpaw", replacement = "Front paw", x = Location)) %>% 
    dplyr::mutate(Location = gsub(pattern = "Lips", replacement = "Penlabial area", x = Location)) %>% 
    dplyr::mutate(Location = gsub(pattern = "Tail", replacement = "Ventral tail", x = Location)) %>% 
    dplyr::mutate(Location = factor(x = Location, 
                                    levels = c("Abdomen", "Front paw", "Penlabial area", "Ventral tail", "Stool"))) %>%
        dplyr::mutate(Group = gsub(pattern = "AD", replacement = "AD ", x = Group)) %>% 
    dplyr::mutate(Group = factor(x = Group, 
                                 levels = c("AD pre", "AD post", "AD treatment"))) %>%
    ggplot(data = ., aes(x=Group, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Group, X,
                  label = p_clean),
              color = "black", size = 4) +
    facet_wrap(Location~., scales = "free_x", ncol = 5) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("") + ggtitle("Correlation with CADESI")
plot_corr_genera

plot_corr_phyla
ggsave(filename = paste0( "plots/Correlations_CADESI_Phylum.pdf"), width = 12, height = 10, plot = plot_corr_phyla)
plot_corr_genera
ggsave(filename = paste0( "plots/Correlations_CADESI_Genus.pdf"), width = 12, height = 10, plot = plot_corr_genera)

# Session Info ------------------------------------------------------------

sessioninfo::session_info()

