
# General Information -----------------------------------------------------

# Author: Axel Künstner
# Date:   2023-08-08
# Analysis for revision

# Libraries ---------------------------------------------------------------

library(phyloseq)
library(tidyverse)
library(patchwork)
library(ggpubr)

library(lme4)

library(corncob)
library(ANCOMBC)

# devtools::install_github("davidgohel/flextable")
library(flextable)

theme_user <- 
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 0, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", linewidth=1.5, linetype="solid")) 

my_res <- list()

# Get data ----------------------------------------------------------

ps_all <- readRDS(file = "data/phyloseq.OTU.RDS")
haplo_data <- read.delim("data/haplogroup_data.txt")

# alpha diversity, the fast/easy way
sample_data(ps_all)$Shannon <- estimate_richness(physeq = ps_all, measures = "Shannon")$Shannon

location <- 'ear'

# Ear/medial pinnea---------------------------------------------------

ps <- subset_samples(physeq = ps_all, location == 'ear')

# alpha diversity DivNet
otu_table(ps) %>% rowSums() %>% sort(decreasing = T) %>% head
dv_ind <- DivNet::divnet(ps, base = "Otu2", ncores = 10,
                         network = "diagonal",
                         tuning = list(EMiter = 6, EMburn = 3, MCiter = 250, MCburn = 100))
sample_data(ps)$DivNet <- dv_ind$shannon %>% summary() %>% data.frame %>% .$estimate
sample_data(ps)$DivNet_lower <- dv_ind$shannon %>% summary() %>% data.frame %>% .$lower
sample_data(ps)$DivNet_upper <- dv_ind$shannon %>% summary() %>% data.frame %>% .$upper
sample_data(ps)$DivNet_error <- dv_ind$shannon %>% summary() %>% data.frame %>% .$error

table(haplo_data$Haplogroup)
unique( sample_data(ps)$dogID ) %in% haplo_data$dogID

cova <- sample_data(ps) %>% data.frame()
cova <- merge(x = cova, y = haplo_data, by = "dogID")
rownames(cova) <- cova$SequencingID

# rough chi squared calculation
table(cova$MajorGroup, cova$group)
chisq.test( table(cova$MajorGroup, cova$group) )
chisq.test( table(cova$MajorGroup, cova$group=="healthy") )
fisher.test( table(cova$MajorGroup, cova$group) )

table(cova$sex, cova$group)
chisq.test( table(cova$sex, cova$group) )

sample_data(ps) <- cova
rm(cova)

ps <- subset_samples(physeq = ps, group != 'ADpre')
ps <- subset_samples(physeq = ps, Haplogroup != 'B1')
ps <- subset_samples(physeq = ps, Haplogroup != 'A64')

ps_clr <- microbiome::transform(x = ps, transform = "clr")
ps_ra  <- transform_sample_counts(ps, function(x){x / sum(x)})

cova <- sample_data(ps) %>% data.frame()

colv2 <- c("steelblue", "orange")

table(cova$MajorGroup)

# Alpha diversity ---------------------------------------------------------

# https://yury-zablotski.netlify.app/post/mixed-models/

# m1 <- lmer(Shannon ~ 1 + MajorGroup + (1|group) , data = cova)
# summary(m1)

# Random effects model (both intercept and slope are random)
m1 <- lmerTest::lmer(DivNet ~ MajorGroup + (MajorGroup|sex) + (MajorGroup|group) , 
                     data = cova, REML = T)
summary(m1)
coef(m1)$group # show random intercepts and random slopes
lattice::dotplot(ranef(m1, condVar=T))

# no random slope
m2 <-lmerTest::lmer(DivNet ~ MajorGroup + (1|sex) + (1|group) , 
                    data = cova, REML = T)
# no random slope, no sex
m3 <-lmerTest::lmer(DivNet ~ MajorGroup + (1|group) , 
                    data = cova, REML = T)

anova(m1, m2, m3)

summary(m3)

lattice::dotplot(ranef(m1, condVar=T))
lattice::dotplot(ranef(m2, condVar=T))
lattice::dotplot(ranef(m3, condVar=T))

p1 <- cova %>% 
    dplyr::mutate(group = gsub(pattern = "AD", replacement = "AD ", x = group)) %>% 
    ggplot(data = ., aes(y = DivNet, x = MajorGroup)) +
    geom_point(mapping = aes(color = MajorGroup, alpha = 0.75)) +
    scale_color_manual(values = colv2) +
    xlab("") +
    ylab("DivNet estimate of Shannon") +
    theme_classic(base_size = 15) +
    stat_smooth(method = "lm", formula = 'y ~ x', se=F, fullrange = T) +
    facet_wrap(~group) +
    ylim(0,6) +
    theme_user +
    theme(legend.position = "none")
p1
p2 <- cova %>% 
    ggpubr::ggviolin(data = ., 
                     x = "MajorGroup", y = "DivNet", 
                     add = c("median_iqr", "jitter"), 
                     trim = TRUE, width = 0.5,
                     fill = "MajorGroup", palette = colv2, alpha = 0.75,
                     add.params = list(color = "grey35",
                                       size = 1)) +
    xlab('') + ylab("DivNet estimate of Shannon") +
    ylim(0,6) +
    theme_user +
    theme(legend.position = "none") # + stat_compare_means()
p2

breakaway::betta(chats = cova$DivNet,
                 ses = cova$DivNet_error,
                 X = model.matrix(~MajorGroup, data = cova) )$table

layout <- "
AAAB
CCCX"

m3.coeff <- coef(summary(m3, ddf="Satterthwaite")) %>% data.frame() %>% 
    rownames_to_column()
colnames(m3.coeff) <- c("Variable", "Estimate", "Std. error", "df", "t-value", "p-value")
p3 <- flextable( data = m3.coeff )

p1 + p2 +
    grid::rasterGrob(flextable::as_raster(x = p3)) +
    plot_layout(design = layout) +
    plot_annotation(title = 'Alpha diversity')
ggsave(filename = paste0("plots/haplogroups_", location, "_Alpha.pdf"), height = 6, width = 9)

# Beta diversity ----------------------------------------------------------

otu.table_clr <- otu_table(ps_clr) %>% t()
ps_clr_dist <- dist(otu.table_clr, method="euclidean")
ps_clr_ord <- phyloseq::ordinate(physeq = ps_clr, 
                                 method = "CAP", formula = ~ MajorGroup + Condition(group),
                                 distance = "euclidean")
p_beta <- plot_ordination( 
    physeq = ps_clr, 
    ordination = ps_clr_ord, color='MajorGroup')  + 
    scale_color_manual(values=colv2) +
    geom_point( mapping = aes(color = MajorGroup, alpha = 0.75), size = 2 ) + 
    # geom_text(mapping = aes(label = dogID), size = 3, nudge_x = 0.25) +
    ggtitle("CCA of Aitchison distance") +
    theme_user +
    theme(legend.position = "none") 

df1 <- vegan::anova.cca(ps_clr_ord, permutations = 99999) %>% 
    data.frame() %>% 
    rownames_to_column()
df2 <- vegan::anova.cca(ps_clr_ord, by="term", permutations = 99999, parallel=8) %>% 
    data.frame() %>% 
    rownames_to_column()
df3 <- vegan::anova.cca(ps_clr_ord, by="axis", permutations = 99999, parallel=8) %>% 
    data.frame() %>% 
    rownames_to_column()

colnames(df1) <- c("Variable", "df", "Variance", "F-statistic", "P-Value")
colnames(df2) <- c("Variable", "df", "Variance", "F-statistic", "P-Value")
colnames(df3) <- c("Variable", "df", "Variance", "F-statistic", "P-Value")

p_beta2 <- flextable( data = df1 )
p_beta3 <- flextable( data = df2 )
p_beta4 <- flextable( data = df3 )

layout <- "
AAB
AAC
AAD"

p_beta +
    grid::rasterGrob(flextable::as_raster(x = p_beta2)) +
    grid::rasterGrob(flextable::as_raster(x = p_beta3)) +
    grid::rasterGrob(flextable::as_raster(x = p_beta4)) +
    plot_layout(design = layout) +
    plot_annotation(title = 'Beta diversity')
ggsave(filename = paste0("plots/haplogroups_", location, "_Beta.pdf"), height = 6, width = 9)

# Differential abundance analysis -----------------------------------------

sigLevel <- 0.1
TAXRANK <- "Phylum"

corncob_dat <- ps %>%
    subset_samples(physeq = ., location == location) %>% 
    microbiome::aggregate_taxa(x = ., level = TAXRANK) %>% 
    tax_glom(physeq = ., taxrank = TAXRANK, NArm = TRUE) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.2) # Prevalence filtering (20%)

da_analysis <-
    differentialTest(formula = ~ MajorGroup + group,
                     phi.formula = ~ MajorGroup + group, # model to be fitted to the dispersion
                     formula_null = ~ group, # Formula for mean under null, without response
                     phi.formula_null = ~ MajorGroup + group, # Formula for overdispersion under null, without response
                     test = "Wald", boot = FALSE, B = 0,
                     data = corncob_dat,
                     fdr_cutoff = sigLevel)
( da_plot <- plot(da_analysis, level = c("Phylum"))  )
da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = corncob_dat, level = TAXRANK),
    # Taxa_veri = da_plot$data$taxa,
    Location = location,
    Comparison = "ADpre vs Healthy",
    p_fdr = da_analysis$p_fdr[da_analysis$p_fdr < sigLevel & !is.na(da_analysis$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = corncob_dat, 
    formula = "MajorGroup+group", group = "group", 
    p_adj_method = "BH", alpha = sigLevel, 
    # zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)

da_ancombc <- da_ancombc$res$diff_abn # %>% rownames_to_column("Feature")
# colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature = taxon, DA = MajorGroupCn) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

p_da_phylum <- da_results %>%
    dplyr::mutate(Feature = gsub(pattern = "_", replacement = " ", Feature)) %>%
    ggplot(data = ., aes(x = Effect, 
                         y = factor(Feature, levels = rev(levels(factor(Feature)))))) +
    geom_point(size = 3, aes(colour = factor(Score))) +
    scale_color_manual(values = c("cornflowerblue", "orange")) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    # facet_wrap(facets = Comparison~ ., scales = "free_x", ncol = 3) +
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
          strip.background = element_rect(colour="white", fill="white", linewidth =1.5, linetype="solid")
    ) +
    xlim(-10, 10) +
    xlab("Effect size") + ylab("")
p_da_phylum
my_res[[paste(TAXRANK, location, sep = " ")]] <- da_results

TAXRANK <- "Genus"

corncob_dat <- ps %>%
    subset_samples(physeq = ., location == location) %>% 
    microbiome::aggregate_taxa(x = ., level = TAXRANK) %>% 
    tax_glom(physeq = ., taxrank = TAXRANK, NArm = TRUE) %>% 
    metagMisc::phyloseq_filter_prevalence(physeq = ., prev.trh = 0.2) # Prevalence filtering (20%)

da_analysis <-
    differentialTest(formula = ~ MajorGroup + group ,
                     phi.formula = ~ MajorGroup + group, # model to be fitted to the dispersion
                     formula_null = ~ group, # Formula for mean under null, without response
                     phi.formula_null = ~ MajorGroup + group, # Formula for overdispersion under null, without response
                     test = "Wald", boot = FALSE, B = 0,
                     data = corncob_dat,
                     fdr_cutoff = sigLevel)
( da_plot <- plot(da_analysis, level = c("Phylum", "Genus"))  )
da_results <- data.frame(
    Feature = corncob::otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = corncob_dat, level = TAXRANK),
    # Taxa_veri = da_plot$data$taxa,
    Location = location,
    Comparison = "ADpre vs Healthy",
    p_fdr = da_analysis$p_fdr[da_analysis$p_fdr < sigLevel & !is.na(da_analysis$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax, row.names = NULL,
    ANCOMBC = FALSE
)

# ANCOM-BC
da_ancombc <- ancombc(
    phyloseq = corncob_dat, 
    formula = "MajorGroup+group", group = "group", 
    p_adj_method = "BH", alpha = sigLevel, 
    # zero_cut = 1, # no prev filtering necessary anymore 
    lib_cut = 0, struc_zero = TRUE, neg_lb = FALSE, tol = 1e-5, max_iter = 10000, 
    conserve = FALSE, # TRUE if small sample sizes
    global = FALSE
)

da_ancombc <- da_ancombc$res$diff_abn # %>% rownames_to_column("Feature")
# colnames(da_ancombc) <- c("Feature", "DA")
da_ancombc <- da_ancombc %>% dplyr::select(Feature = taxon, DA = MajorGroupCn) %>% dplyr::filter(DA == TRUE) 
da_ancombc
da_results$ANCOMBC[da_results$Feature %in% da_ancombc$Feature] <- TRUE

da_results$Score <- 1 + da_results$ANCOMBC

p_da_genus <- da_results %>%
    dplyr::mutate(Feature = gsub(pattern = "_", replacement = " ", Feature)) %>%
    ggplot(data = ., aes(x = Effect, 
                         y = factor(Feature, levels = rev(levels(factor(Feature)))))) +
    geom_point(size = 3, aes(colour = factor(Score))) +
    scale_color_manual(values = c("cornflowerblue", "orange")) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    # facet_wrap(facets = Comparison~ ., scales = "free_x", ncol = 3) +
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
    xlim(-10, 10) +
    xlab("Effect size") + ylab("")
p_da_genus
my_res[[paste(TAXRANK, location, sep = " ")]] <- da_results

layout <- "
AB
XB
XB"

p_da_phylum + p_da_genus + 
    plot_layout(design = layout) +
    plot_annotation(title = 'Differential Abundance Analysis')
ggsave(filename = paste0("plots/haplogroups_", location, "_da.pdf"), height = 6, width = 12)

WriteXLS::WriteXLS(x = my_res, ExcelFileName = paste0("tables/Results_Haplogroups_", location, ".xlsx"), row.names = F, AdjWidth = T, BoldHeaderRow = T, FreezeRow = 1)

# Pretty plotting ---------------------------------------------------------

layout <- "
AAB
AAB
CCD
CCE
CCE
"

p1 + 
    p2 +
    p_beta + ggtitle("") +
    p_da_phylum +
    p_da_genus + 
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'a')
ggsave(filename = paste0("plots/Fig5_haplogroups_", location, ".pdf"), height = 9, width = 12)

# Session Info ------------------------------------------------------------

sessioninfo::session_info()
