
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

ps <- readRDS(file = "phyloseq.OTU.RDS")

sample_data(ps)$group <- gsub(pattern = "healthy", 
                              replacement = "Healthy", 
                              x = sample_data(ps)$group)

# Staph
stap <- gdata::read.xls(xls = "../data_sintax_skin_stool/staphOTUs.xlsx")
rownames(stap) <- stap$OTU

# Add Staph to tax table --------------------------------------------------

ps_staph <- subset_taxa(physeq = ps, Genus == "Staphylococcus" )

rownames(tax_table(ps_staph)) == stap$OTU
stap <- stap[ rownames(tax_table(ps_staph)), ]
rownames(tax_table(ps_staph)) == stap$OTU

tt <- tax_table(ps_staph) %>% data.frame()
tt$Species <- stap$Species
tax_table(ps_staph) <- as.matrix(tt)

# Setup -------------------------------------------------------------------

siglevel <- 0.1
TAXLEVEL <- 'Species'

my_da_list <- list()

# abdomen
# frontPaw
# lips -> perilabial area
# tail -> ventral tail

# Process data ------------------------------------------------------------

# subsets 
ps.ADpre_healthy   <- subset_samples(ps_staph, group == "Healthy" | group == "ADpre")
ps.ADpost_ADpre    <- subset_samples(ps_staph, group == "ADpost"  | group == "ADpre")
ps.ADtreat_healthy <- subset_samples(ps_staph, group == "Healthy" | group == "ADtreatment")

sample_data(ps.ADpre_healthy)$group <- sample_data(ps.ADpre_healthy)$group %>% 
    fct_relevel(c("Healthy", "ADpre"))
sample_data(ps.ADpost_ADpre)$group <- sample_data(ps.ADpost_ADpre)$group %>% 
    fct_relevel(c("ADpre", "ADpost"))
sample_data(ps.ADtreat_healthy)$group <- sample_data(ps.ADtreat_healthy)$group %>% 
    fct_relevel(c("Healthy", "ADtreatment"))

# DA: abdomen ---------------------------------------------------------------

test_location <- "abdomen"

# ADpre vs healthy
da_plot <- NULL; da_results <- NULL; da_corncob <- NULL
ps_da <- ps.ADpre_healthy %>% 
    tax_glom(physeq = ., taxrank = "Species") %>% 
    tax_glom(physeq = ., taxrank = "Species")

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

da_results <- data.frame(
    Feature = rownames(da_ancombc$res$q_val[1]),
    Species = corncob::otu_to_taxonomy(OTU = rownames(da_ancombc$res$q_val[1]), data = ps_da, level = TAXLEVEL),
    Location = test_location,
    Comparison = "ADpre vs Healthy",
    p_fdr = da_ancombc$res$q_val[[1]],
    Effect = da_ancombc$res$beta[[1]],
    Error_min = da_ancombc$res$beta[[1]] - da_ancombc$res$se[[1]],
    Error_max = da_ancombc$res$beta[[1]] + da_ancombc$res$se[[1]], 
    W = da_ancombc$res$W[[1]],
    diff_abn = da_ancombc$res$diff_abn[[1]],
    row.names = NULL
)

my_da_list[[paste0(test_location, " ADpre vs Healthy")]] <- da_results

# ADpost vs ADpre
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpost_ADpre %>% 
    subset_samples(physeq = ., location == test_location) %>%      tax_glom(physeq = ., taxrank = "Species")

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

da_results <- data.frame(
    Feature = rownames(da_ancombc$res$q_val[1]),
    Species = corncob::otu_to_taxonomy(OTU = rownames(da_ancombc$res$q_val[1]), data = ps_da, level = TAXLEVEL),
    Location = test_location,
    Comparison = "ADpost vs ADpre",
    p_fdr = da_ancombc$res$q_val[[1]],
    Effect = da_ancombc$res$beta[[1]],
    Error_min = da_ancombc$res$beta[[1]] - da_ancombc$res$se[[1]],
    Error_max = da_ancombc$res$beta[[1]] + da_ancombc$res$se[[1]], 
    W = da_ancombc$res$W[[1]],
    diff_abn = da_ancombc$res$diff_abn[[1]],
    row.names = NULL
)

my_da_list[[paste0(test_location, " ADpost vs ADpre")]] <- da_results

# ADtreatment vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADtreat_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>%      tax_glom(physeq = ., taxrank = "Species")

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

da_results <- data.frame(
    Feature = rownames(da_ancombc$res$q_val[1]),
    Species = corncob::otu_to_taxonomy(OTU = rownames(da_ancombc$res$q_val[1]), data = ps_da, level = TAXLEVEL),
    Location = test_location,
    Comparison = "ADtreat vs Healthy",
    p_fdr = da_ancombc$res$q_val[[1]],
    Effect = da_ancombc$res$beta[[1]],
    Error_min = da_ancombc$res$beta[[1]] - da_ancombc$res$se[[1]],
    Error_max = da_ancombc$res$beta[[1]] + da_ancombc$res$se[[1]], 
    W = da_ancombc$res$W[[1]],
    diff_abn = da_ancombc$res$diff_abn[[1]],
    row.names = NULL
)

my_da_list[[paste0(test_location, " ADtreat vs Healthy")]] <- da_results

# DA: frontPaw ---------------------------------------------------------------

test_location <- "frontPaw"

# ADpre vs healthy
da_plot <- NULL; da_results <- NULL; da_corncob <- NULL
ps_da <- ps.ADpre_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>%      tax_glom(physeq = ., taxrank = "Species")

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

da_results <- data.frame(
    Feature = rownames(da_ancombc$res$q_val[1]),
    Species = corncob::otu_to_taxonomy(OTU = rownames(da_ancombc$res$q_val[1]), data = ps_da, level = TAXLEVEL),
    Location = test_location,
    Comparison = "ADpre vs Healthy",
    p_fdr = da_ancombc$res$q_val[[1]],
    Effect = da_ancombc$res$beta[[1]],
    Error_min = da_ancombc$res$beta[[1]] - da_ancombc$res$se[[1]],
    Error_max = da_ancombc$res$beta[[1]] + da_ancombc$res$se[[1]], 
    W = da_ancombc$res$W[[1]],
    diff_abn = da_ancombc$res$diff_abn[[1]],
    row.names = NULL
)

my_da_list[[paste0(test_location, " ADpre vs Healthy")]] <- da_results

# ADpost vs ADpre
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpost_ADpre %>% 
    subset_samples(physeq = ., location == test_location) %>%      tax_glom(physeq = ., taxrank = "Species")

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

da_results <- data.frame(
    Feature = rownames(da_ancombc$res$q_val[1]),
    Species = corncob::otu_to_taxonomy(OTU = rownames(da_ancombc$res$q_val[1]), data = ps_da, level = TAXLEVEL),
    Location = test_location,
    Comparison = "ADpost vs ADpre",
    p_fdr = da_ancombc$res$q_val[[1]],
    Effect = da_ancombc$res$beta[[1]],
    Error_min = da_ancombc$res$beta[[1]] - da_ancombc$res$se[[1]],
    Error_max = da_ancombc$res$beta[[1]] + da_ancombc$res$se[[1]], 
    W = da_ancombc$res$W[[1]],
    diff_abn = da_ancombc$res$diff_abn[[1]],
    row.names = NULL
)

my_da_list[[paste0(test_location, " ADpost vs ADpre")]] <- da_results

# ADtreatment vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADtreat_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>%      tax_glom(physeq = ., taxrank = "Species")

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

da_results <- data.frame(
    Feature = rownames(da_ancombc$res$q_val[1]),
    Species = corncob::otu_to_taxonomy(OTU = rownames(da_ancombc$res$q_val[1]), data = ps_da, level = TAXLEVEL),
    Location = test_location,
    Comparison = "ADtreat vs Healthy",
    p_fdr = da_ancombc$res$q_val[[1]],
    Effect = da_ancombc$res$beta[[1]],
    Error_min = da_ancombc$res$beta[[1]] - da_ancombc$res$se[[1]],
    Error_max = da_ancombc$res$beta[[1]] + da_ancombc$res$se[[1]], 
    W = da_ancombc$res$W[[1]],
    diff_abn = da_ancombc$res$diff_abn[[1]],
    row.names = NULL
)

my_da_list[[paste0(test_location, " ADtreat vs Healthy")]] <- da_results

# DA: lips ---------------------------------------------------------------

test_location <- "lips"

# ADpre vs healthy
da_plot <- NULL; da_results <- NULL; da_corncob <- NULL
ps_da <- ps.ADpre_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>%      tax_glom(physeq = ., taxrank = "Species")

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

da_results <- data.frame(
    Feature = rownames(da_ancombc$res$q_val[1]),
    Species = corncob::otu_to_taxonomy(OTU = rownames(da_ancombc$res$q_val[1]), data = ps_da, level = TAXLEVEL),
    Location = test_location,
    Comparison = "ADpre vs Healthy",
    p_fdr = da_ancombc$res$q_val[[1]],
    Effect = da_ancombc$res$beta[[1]],
    Error_min = da_ancombc$res$beta[[1]] - da_ancombc$res$se[[1]],
    Error_max = da_ancombc$res$beta[[1]] + da_ancombc$res$se[[1]], 
    W = da_ancombc$res$W[[1]],
    diff_abn = da_ancombc$res$diff_abn[[1]],
    row.names = NULL
)

my_da_list[[paste0(test_location, " ADpre vs Healthy")]] <- da_results

# ADpost vs ADpre
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpost_ADpre %>% 
    subset_samples(physeq = ., location == test_location) %>%      tax_glom(physeq = ., taxrank = "Species")

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

da_results <- data.frame(
    Feature = rownames(da_ancombc$res$q_val[1]),
    Species = corncob::otu_to_taxonomy(OTU = rownames(da_ancombc$res$q_val[1]), data = ps_da, level = TAXLEVEL),
    Location = test_location,
    Comparison = "ADpost vs ADpre",
    p_fdr = da_ancombc$res$q_val[[1]],
    Effect = da_ancombc$res$beta[[1]],
    Error_min = da_ancombc$res$beta[[1]] - da_ancombc$res$se[[1]],
    Error_max = da_ancombc$res$beta[[1]] + da_ancombc$res$se[[1]], 
    W = da_ancombc$res$W[[1]],
    diff_abn = da_ancombc$res$diff_abn[[1]],
    row.names = NULL
)

my_da_list[[paste0(test_location, " ADpost vs ADpre")]] <- da_results

# ADtreatment vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADtreat_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>%      tax_glom(physeq = ., taxrank = "Species")

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

da_results <- data.frame(
    Feature = rownames(da_ancombc$res$q_val[1]),
    Species = corncob::otu_to_taxonomy(OTU = rownames(da_ancombc$res$q_val[1]), data = ps_da, level = TAXLEVEL),
    Location = test_location,
    Comparison = "ADtreat vs Healthy",
    p_fdr = da_ancombc$res$q_val[[1]],
    Effect = da_ancombc$res$beta[[1]],
    Error_min = da_ancombc$res$beta[[1]] - da_ancombc$res$se[[1]],
    Error_max = da_ancombc$res$beta[[1]] + da_ancombc$res$se[[1]], 
    W = da_ancombc$res$W[[1]],
    diff_abn = da_ancombc$res$diff_abn[[1]],
    row.names = NULL
)

my_da_list[[paste0(test_location, " ADtreat vs Healthy")]] <- da_results

# DA: tail ---------------------------------------------------------------

test_location <- "tail"

# ADpre vs healthy
da_plot <- NULL; da_results <- NULL; da_corncob <- NULL
ps_da <- ps.ADpre_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>%      tax_glom(physeq = ., taxrank = "Species")

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

da_results <- data.frame(
    Feature = rownames(da_ancombc$res$q_val[1]),
    Species = corncob::otu_to_taxonomy(OTU = rownames(da_ancombc$res$q_val[1]), data = ps_da, level = TAXLEVEL),
    Location = test_location,
    Comparison = "ADpre vs Healthy",
    p_fdr = da_ancombc$res$q_val[[1]],
    Effect = da_ancombc$res$beta[[1]],
    Error_min = da_ancombc$res$beta[[1]] - da_ancombc$res$se[[1]],
    Error_max = da_ancombc$res$beta[[1]] + da_ancombc$res$se[[1]], 
    W = da_ancombc$res$W[[1]],
    diff_abn = da_ancombc$res$diff_abn[[1]],
    row.names = NULL
)

my_da_list[[paste0(test_location, " ADpre vs Healthy")]] <- da_results

# ADpost vs ADpre
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADpost_ADpre %>% 
    subset_samples(physeq = ., location == test_location) %>%      tax_glom(physeq = ., taxrank = "Species")

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

da_results <- data.frame(
    Feature = rownames(da_ancombc$res$q_val[1]),
    Species = corncob::otu_to_taxonomy(OTU = rownames(da_ancombc$res$q_val[1]), data = ps_da, level = TAXLEVEL),
    Location = test_location,
    Comparison = "ADpost vs ADpre",
    p_fdr = da_ancombc$res$q_val[[1]],
    Effect = da_ancombc$res$beta[[1]],
    Error_min = da_ancombc$res$beta[[1]] - da_ancombc$res$se[[1]],
    Error_max = da_ancombc$res$beta[[1]] + da_ancombc$res$se[[1]], 
    W = da_ancombc$res$W[[1]],
    diff_abn = da_ancombc$res$diff_abn[[1]],
    row.names = NULL
)

my_da_list[[paste0(test_location, " ADpost vs ADpre")]] <- da_results

# ADtreatment vs healthy
da_plot <- NULL; da_results <- NULL
ps_da <- ps.ADtreat_healthy %>% 
    subset_samples(physeq = ., location == test_location) %>%      tax_glom(physeq = ., taxrank = "Species")

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

da_results <- data.frame(
    Feature = rownames(da_ancombc$res$q_val[1]),
    Species = corncob::otu_to_taxonomy(OTU = rownames(da_ancombc$res$q_val[1]), data = ps_da, level = TAXLEVEL),
    Location = test_location,
    Comparison = "ADtreat vs Healthy",
    p_fdr = da_ancombc$res$q_val[[1]],
    Effect = da_ancombc$res$beta[[1]],
    Error_min = da_ancombc$res$beta[[1]] - da_ancombc$res$se[[1]],
    Error_max = da_ancombc$res$beta[[1]] + da_ancombc$res$se[[1]], 
    W = da_ancombc$res$W[[1]],
    diff_abn = da_ancombc$res$diff_abn[[1]],
    row.names = NULL
)

my_da_list[[paste0(test_location, " ADtreat vs Healthy")]] <- da_results

# Plotting ----------------------------------------------------------------

df_da <- bind_rows(my_da_list) %>% 
    dplyr::rename(Otu = Feature,
                  Feature = Species) %>% 
    dplyr::filter(diff_abn == TRUE) %>% 
    dplyr::filter(Feature != "") %>% 
    dplyr::mutate(Score = 1) %>% 
    dplyr::filter( abs(Effect) > 0.3 )

summary(df_da$Effect)
xmin <- min(df_da$Error_min) %>% floor
xmax <- max(df_da$Error_max) %>% ceiling

head(df_da)

table(df_da$Location, df_da$Comparison)

colnames(df_da)

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
    scale_color_manual(values = c("grey45", "grey75")) +
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
    xlab("beta coefficient") + ylab("")

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
    scale_color_manual(values = c("grey45", "grey75")) +
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
    xlab("beta coefficient") + ylab("")

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
    scale_color_manual(values = c("grey45", "grey75")) +
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
    xlab("beta coefficient") + ylab("")

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
    scale_color_manual(values = c("grey45", "grey75")) +
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
    xlab("beta coefficient") + ylab("")

layout <- "
AAA
AAA
AAA
AAA
AAA
AAA
AAA
AAA
AAA
AAA
BBB
BBB
BBB
BBB
BBB
BBB
BBB
CCC
CCC
CCC
CCC
DDD
DDD
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
ggsave(filename = paste0("plots/Fig3.pdf"), height = 10, width = 12)

# Session Info ------------------------------------------------------------

sessioninfo::session_info()
