
# General Information -----------------------------------------------------

# Authors: Mirja Thomsen & Axel Künstner
# Date:    2022-05-31

# Libraries ---------------------------------------------------------------

library(phyloseq)
library(tidyverse)

library(corncob)
library(ANCOMBC)

library(patchwork)

# Get data ----------------------------------------------------------------

ps <- readRDS(file = "data/phyloseq.OTU.RDS")

sample_data(ps)$group <- gsub(pattern = "healthy", 
                              replacement = "Healthy", 
                              x = sample_data(ps)$group)

# Setup -------------------------------------------------------------------

siglevel <- 0.1
TAXLEVEL <- 'Genus'

my_da_list <- list()

# stool 
# abdomen
# frontPaw
# lips -> perilabial area
# tail -> ventral tail

# Process data ------------------------------------------------------------

# genus level
ps_g <- ps %>% 
    microbiome::aggregate_taxa(x = ., level = TAXLEVEL) %>% 
    tax_glom(physeq = ., taxrank = TAXLEVEL, NArm = TRUE)

# subsets 
ps.ADpre_healthy   <- subset_samples(ps_g, group == "Healthy" | group == "ADpre")
ps.ADpost_ADpre    <- subset_samples(ps_g, group == "ADpost"  | group == "ADpre")
ps.ADtreat_healthy <- subset_samples(ps_g, group == "Healthy" | group == "ADtreatment")

sample_data(ps.ADpre_healthy)$group <- sample_data(ps.ADpre_healthy)$group %>% 
    fct_relevel(c("Healthy", "ADpre"))
sample_data(ps.ADpost_ADpre)$group <- sample_data(ps.ADpost_ADpre)$group %>% 
    fct_relevel(c("ADpre", "ADpost"))
sample_data(ps.ADtreat_healthy)$group <- sample_data(ps.ADtreat_healthy)$group %>% 
    fct_relevel(c("Healthy", "ADtreatment"))

# DA: stool ---------------------------------------------------------------

test_location <- "stool"

# ADpre vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpre_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.2) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = 0.1)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADpre vs Healthy",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da, 
    formula = "group+sex", group = "group", 
    p_adj_method = "BH", alpha = siglevel, 
    zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>% 
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADpre vs Healthy")]] <- da_results

# ADpost vs ADpre
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpost_ADpre %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.20) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = siglevel)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADpost vs ADpre",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da, 
    formula = "group+sex", group = "group", 
    p_adj_method = "BH", alpha = siglevel, 
    zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>% 
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADpost vs ADpre")]] <- da_results

# ADtreatment vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADtreat_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.20) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = siglevel)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADtreat vs Healthy",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da, 
    formula = "group+sex", group = "group", 
    p_adj_method = "BH", alpha = siglevel, 
    zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>% 
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADtreat vs Healthy")]] <- da_results

# DA: abdomen ---------------------------------------------------------------

test_location <- "abdomen"

# ADpre vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpre_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.2) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = 0.1)
da_corncob$significant_taxa
# da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
#     geom_point(size = 3)
# da_plot
# 
# da_results <- data.frame(
#     Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
#     # Taxa_veri = da_plot$data$taxa,
#     Location = test_location,
#     Comparison = "ADpre vs Healthy",
#     p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
#     Effect = da_plot$data$x,
#     Error_min = da_plot$data$xmin,
#     Error_max = da_plot$data$xmax, row.names = NULL,
#     ANCOMBC = FALSE
# )
# 
# # ANCOM-BC
# da_ancombc <- ancombc(
#     phyloseq = ps_da, 
#     formula = "group+sex", group = "group", 
#     p_adj_method = "BH", alpha = siglevel, 
#     zero_cut = 1, # no prev filtering necessary anymore 
#     lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
#     conserve = FALSE, # TRUE if small sample sizes
#     global = FALSE
# )
# da_ancombc <- da_ancombc$res$diff_abn %>% 
#     rownames_to_column("Feature")
# colnames(da_ancombc) <- c("Feature", "DA")
# da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE) 
# da_ancombc
# da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE
# 
# da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADpre vs Healthy")]] <- da_results

# ADpost vs ADpre
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpost_ADpre %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.20) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = siglevel)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADpost vs ADpre",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da, 
    formula = "group+sex", group = "group", 
    p_adj_method = "BH", alpha = siglevel, 
    zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>% 
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADpost vs ADpre")]] <- da_results

# ADtreatment vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADtreat_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.20) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = siglevel)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADtreat vs Healthy",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da, 
    formula = "group+sex", group = "group", 
    p_adj_method = "BH", alpha = siglevel, 
    zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>% 
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADtreat vs Healthy")]] <- da_results

# DA: frontPaw ---------------------------------------------------------------

test_location <- "frontPaw"

# ADpre vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpre_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.2) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = 0.1)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADpre vs Healthy",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da,
    formula = "group+sex", group = "group",
    p_adj_method = "BH", alpha = siglevel,
    zero_cut = 1, # no prev filtering necessary anymore
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000,
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>%
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE)
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADpre vs Healthy")]] <- da_results

# ADpost vs ADpre
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpost_ADpre %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.20) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = siglevel)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADpost vs ADpre",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da, 
    formula = "group+sex", group = "group", 
    p_adj_method = "BH", alpha = siglevel, 
    zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>% 
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADpost vs ADpre")]] <- da_results

# ADtreatment vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADtreat_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.20) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = siglevel)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADtreat vs Healthy",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da, 
    formula = "group+sex", group = "group", 
    p_adj_method = "BH", alpha = siglevel, 
    zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>% 
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADtreat vs Healthy")]] <- da_results

# DA: lips ---------------------------------------------------------------

test_location <- "lips"

# ADpre vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpre_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.2) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = 0.1)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADpre vs Healthy",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da,
    formula = "group+sex", group = "group",
    p_adj_method = "BH", alpha = siglevel,
    zero_cut = 1, # no prev filtering necessary anymore
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000,
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>%
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE)
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADpre vs Healthy")]] <- da_results

# ADpost vs ADpre
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpost_ADpre %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.20) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = siglevel)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADpost vs ADpre",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da, 
    formula = "group+sex", group = "group", 
    p_adj_method = "BH", alpha = siglevel, 
    zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>% 
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADpost vs ADpre")]] <- da_results

# ADtreatment vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADtreat_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.20) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = siglevel)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADtreat vs Healthy",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da, 
    formula = "group+sex", group = "group", 
    p_adj_method = "BH", alpha = siglevel, 
    zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>% 
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADtreat vs Healthy")]] <- da_results

# DA: tail ---------------------------------------------------------------

test_location <- "tail"

# ADpre vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpre_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.2) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = 0.1)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADpre vs Healthy",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da,
    formula = "group+sex", group = "group",
    p_adj_method = "BH", alpha = siglevel,
    zero_cut = 1, # no prev filtering necessary anymore
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000,
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>%
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE)
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADpre vs Healthy")]] <- da_results

# ADpost vs ADpre
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpost_ADpre %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.20) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = siglevel)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADpost vs ADpre",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da, 
    formula = "group+sex", group = "group", 
    p_adj_method = "BH", alpha = siglevel, 
    zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>% 
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADpost vs ADpre")]] <- da_results

# ADtreatment vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADtreat_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.20) # Prevalence filtering (20%)

da_corncob <- 
    differentialTest(formula = ~ group + sex,
                     phi.formula = ~ group + sex,
                     formula_null = ~ sex,
                     phi.formula_null = ~ group + sex,
                     test = "Wald", boot = FALSE, B = 0,
                     data = ps_da, verbose = F,
                     fdr = "BH", fdr_cutoff = siglevel)
da_corncob$significant_taxa
da_plot <- plot(da_corncob, level = c(TAXLEVEL)) +
    geom_point(size = 3)
da_plot

da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_da, level = TAXLEVEL),
    # Taxa_veri = da_plot$data$taxa,
    Location = test_location,
    Comparison = "ADtreat vs Healthy",
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = ps_da, 
    formula = "group+sex", group = "group", 
    p_adj_method = "BH", alpha = siglevel, 
    zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)
da_ancombc <- da_ancombc$res$diff_abn %>% 
    rownames_to_column("Feature")
colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature, DA) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

my_da_list[[paste0(test_location, " ADtreat vs Healthy")]] <- da_results

# Plotting ----------------------------------------------------------------

df_da <- bind_rows(my_da_list)

summary(df_da$Effect)
xmin <- min(df_da$Error_min) %>% floor
xmax <- max(df_da$Error_max) %>% ceiling

head(df_da)

table(df_da$Location, df_da$Comparison)

colnames(df_da)

p_stool <- df_da %>%
    dplyr::mutate(Comparison = factor(x = Comparison, 
                                      levels = c("ADpre vs Healthy", 
                                                 "ADtreat vs Healthy", 
                                                 "ADpost vs ADpre"))) %>% 
    dplyr::mutate(Feature = gsub(pattern = "_", replacement = " ", Feature)) %>%
    dplyr::filter(Location == "stool") %>% 
    ggplot(data = ., aes(x = Effect, 
                         y = factor(Feature, levels = rev(levels(factor(Feature)))))) +
    geom_point(size = 3, aes(colour = factor(Score))) +
    scale_color_manual(values = c("cornflowerblue", "orange")) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    facet_wrap(facets = Comparison~ ., scales = "free_x", ncol = 3) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=14),
          strip.text.x = element_text(angle = 0, size = 14),
          axis.text.y = element_text(size=12, face = "italic"),
          text = element_text(size = 14),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")
    ) +
    xlim(xmin, xmax) +
    xlab("Effect size") + ylab("")

p_abdomen <- df_da %>%
    dplyr::mutate(Comparison = factor(x = Comparison, 
                                      levels = c("ADpre vs Healthy", 
                                                 "ADtreat vs Healthy", 
                                                 "ADpost vs ADpre"))) %>% 
    dplyr::mutate(Feature = gsub(pattern = "_", replacement = " ", Feature)) %>%
    dplyr::filter(Location == "abdomen") %>% 
    ggplot(data = ., aes(x = Effect, 
                         y = factor(Feature, levels = rev(levels(factor(Feature)))))) +
    geom_point(size = 3, aes(colour = factor(Score))) +
    scale_color_manual(values = c("cornflowerblue", "orange")) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    facet_wrap(facets = Comparison~ ., scales = "free_x", ncol = 3) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=14),
          strip.text.x = element_text(angle = 0, size = 14),
          axis.text.y = element_text(size=12, face = "italic"),
          text = element_text(size = 14),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")
    ) +
    xlim(xmin, xmax) +
    xlab("Effect size") + ylab("")

p_frontPaw <- df_da %>%
    dplyr::mutate(Comparison = factor(x = Comparison, 
                                      levels = c("ADpre vs Healthy", 
                                                 "ADtreat vs Healthy", 
                                                 "ADpost vs ADpre"))) %>% 
    dplyr::mutate(Feature = gsub(pattern = "_", replacement = " ", Feature)) %>%
    dplyr::filter(Location == "frontPaw") %>% 
    ggplot(data = ., aes(x = Effect, 
                         y = factor(Feature, levels = rev(levels(factor(Feature)))))) +
    geom_point(size = 3, aes(colour = factor(Score))) +
    scale_color_manual(values = c("cornflowerblue", "orange")) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    facet_wrap(facets = Comparison~ ., scales = "free_x", ncol = 3) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=14),
          strip.text.x = element_text(angle = 0, size = 14),
          axis.text.y = element_text(size=12, face = "italic"),
          text = element_text(size = 14),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")
    ) +
    xlim(xmin, xmax) +
    xlab("Effect size") + ylab("")

p_lips <- df_da %>%
    dplyr::mutate(Comparison = factor(x = Comparison, 
                                      levels = c("ADpre vs Healthy", 
                                                 "ADtreat vs Healthy", 
                                                 "ADpost vs ADpre"))) %>% 
    dplyr::mutate(Feature = gsub(pattern = "_", replacement = " ", Feature)) %>%
    dplyr::filter(Location == "lips") %>% 
    ggplot(data = ., aes(x = Effect, 
                         y = factor(Feature, levels = rev(levels(factor(Feature)))))) +
    geom_point(size = 3, aes(colour = factor(Score))) +
    scale_color_manual(values = c("cornflowerblue", "orange")) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    facet_wrap(facets = Comparison~ ., scales = "free_x", ncol = 3) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=14),
          strip.text.x = element_text(angle = 0, size = 14),
          axis.text.y = element_text(size=12, face = "italic"),
          text = element_text(size = 14),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")
    ) +
    xlim(xmin, xmax) +
    xlab("Effect size") + ylab("")

p_tail <- df_da %>%
    dplyr::mutate(Comparison = factor(x = Comparison, 
                                      levels = c("ADpre vs Healthy", 
                                                 "ADtreat vs Healthy", 
                                                 "ADpost vs ADpre"))) %>% 
    dplyr::mutate(Feature = gsub(pattern = "_", replacement = " ", Feature)) %>%
    dplyr::filter(Location == "tail") %>% 
    ggplot(data = ., aes(x = Effect, 
                         y = factor(Feature, levels = rev(levels(factor(Feature)))))) +
    geom_point(size = 3, aes(colour = factor(Score))) +
    scale_color_manual(values = c("cornflowerblue", "orange")) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    facet_wrap(facets = Comparison~ ., scales = "free_x", ncol = 3) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=14),
          strip.text.x = element_text(angle = 0, size = 14),
          axis.text.y = element_text(size=12, face = "italic"),
          text = element_text(size = 14),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")
    ) +
    xlim(xmin, xmax) +
    xlab("Effect size") + ylab("")

layout <- "
AAX
AAX
BBB
BBB
CCC
CCC
CCC
DDD
DDD
DDD
"
p_abdomen +
    p_frontPaw +
    p_lips +
    p_tail +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'a')
ggsave(filename = paste0("plots/DA_testing_", TAXLEVEL, ".pdf"), height = 24, width = 12)

p_stool
ggsave(filename = paste0("plots/DA_testing_", TAXLEVEL, "_stool.pdf"), height = 5, width = 12)

# Write Excel -------------------------------------------------------------

WriteXLS::WriteXLS(x = my_da_list, 
                   ExcelFileName = paste0("tables/DA_testing_", TAXLEVEL, ".xlsx"), 
                   AdjWidth = T, BoldHeaderRow = T, FreezeRow = 1)

WriteXLS::WriteXLS(x = df_da, SheetNames = "DA",
                   ExcelFileName = paste0("tables/DA_testing_long_", TAXLEVEL, ".xlsx"), 
                   AdjWidth = T, BoldHeaderRow = T, FreezeRow = 1)

# Session Info ------------------------------------------------------------

sessioninfo::session_info()
