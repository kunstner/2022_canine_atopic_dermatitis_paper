
# General Information -----------------------------------------------------

# Authors: Axel Künstner
# Date:    2022-05-31

# Libraries ---------------------------------------------------------------

library(phyloseq)
library(tidyverse)
library(patchwork)
library(ggpubr)

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
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) 

layout <- "
AABB
DDCC
"

layout2 <- "
ABC
"

# Staphylococcus ----------------------------------------------------------

ps <- readRDS(file = "data/phyloseq.OTU.RDS")

sample_data(ps)$group <- gsub(pattern = "healthy", 
                              replacement = "Healthy", 
                              x = sample_data(ps)$group)

ps_clr <- microbiome::transform(ps, "clr")
ps_ra  <- transform_sample_counts(ps, function(x){x / sum(x)})

ps.sta <- phyloseq::tax_glom(physeq = ps_ra %>% subset_samples(., group != "ADpost"),
                             taxrank = "Genus") 

table(sample_data(ps.sta)$group)
table(sample_data(ps.sta)$location)

genus_table <- otu_table(ps.sta) %>% data.frame() %>% t()
colnames(genus_table) <- data.frame(tax_table(ps.sta))$Genus

cova <- sample_data(ps.sta) %>% data.frame()
sum( rownames(cova) == rownames(genus_table) ) == nrow(cova)

sort( colnames(genus_table) )

cova$Staphylococcus <- genus_table[, 'Staphylococcus']
cova$Fusobacterium <- genus_table[, 'Fusobacterium']

cova$Test <- genus_table[, 'Megamonas']

cova$group <- factor(cova$group, levels = c("ADpre", "ADtreatment", "Healthy" ))

cova$CADESI[is.na(cova$CADESI)] <- 0
cova$PVAS[is.na(cova$PVAS)] <- 0

m1 <- lm(formula = Staphylococcus ~ CADESI + group, data = cova %>% dplyr::filter(location == "abdomen")) 
m2 <- lm(formula = Staphylococcus ~ CADESI + group, data = cova %>% dplyr::filter(location == "axilla")) 
m3 <- lm(formula = Staphylococcus ~ CADESI + group, data = cova %>% dplyr::filter(location == "ear")) 
m4 <- lm(formula = Staphylococcus ~ CADESI + group, data = cova %>% dplyr::filter(location == "elbow")) 
m5 <- lm(formula = Staphylococcus ~ CADESI + group, data = cova %>% dplyr::filter(location == "flank")) 
m6 <- lm(formula = Staphylococcus ~ CADESI + group, data = cova %>% dplyr::filter(location == "frontPaw")) 
m7 <- lm(formula = Staphylococcus ~ CADESI + group, data = cova %>% dplyr::filter(location == "groin")) 
m8 <- lm(formula = Staphylococcus ~ CADESI + group, data = cova %>% dplyr::filter(location == "hindPaw")) 
m9 <- lm(formula = Staphylococcus ~ CADESI + group, data = cova %>% dplyr::filter(location == "lips")) 
m10<- lm(formula = Staphylococcus ~ CADESI + group, data = cova %>% dplyr::filter(location == "palm")) 
m11<- lm(formula = Staphylococcus ~ CADESI + group, data = cova %>% dplyr::filter(location == "perineum")) 
m12<- lm(formula = Staphylococcus ~ CADESI + group, data = cova %>% dplyr::filter(location == "tail")) 
m13<- lm(formula = Test ~ PVAS + group, data = cova %>% dplyr::filter(location == "stool")) 

m14<- lm(formula = Fusobacterium ~ PVAS + group, data = cova %>% dplyr::filter(location == "stool")) 

summary(m1) # 
summary(m2) # 
summary(m3)
summary(m4) # 
summary(m5) #
summary(m6) #
summary(m7)
summary(m8) # 
summary(m9) # 
summary(m10) #
summary(m11)
summary(m12)
summary(m13)

# double check sign. results
summary(m1)
summary(m2)
summary(m4)
summary(m5)
summary(m6)
summary(m8)
summary(m9)
summary(m10)

sjPlot::plot_model(m10, type = "pred", terms = "CADESI") + theme_user + ggtitle('') +  ylab('Staphylococcus') 

# abdomen
m1.group <- update(m1, Staphylococcus ~ group)
p1 <- sjPlot::plot_model(m1, show.values = T, sort.est = FALSE, vline.color = "grey45") + theme_user + ggtitle('')
p2 <- ggplot(cova %>% dplyr::filter(location == "abdomen"), aes(CADESI, Staphylococcus)) +
    geom_point() + geom_smooth(method='lm') + theme_user + ggtitle('') +  ylab('Staphylococcus') + ylim(-0.1,1.4)
p3 <- sjPlot::plot_model(m1.group, type = "pred", terms = "group", show.values = T, ) + theme_user + xlab('') + ylab('Staphylococcus') + ggtitle('')
emmeans::emmeans(m1.group, pairwise ~ group)$contrasts
p4 <- flextable::flextable( emmeans::emmeans(m1.group, pairwise ~ group)$contrasts %>% as.data.frame() )

p1 + p2 + p3 + 
    grid::rasterGrob(flextable::as_raster(x = p4)) +
    plot_layout(design = layout) +
    plot_annotation(title = 'Abdomen')
ggsave(filename = "plots/model_Staph_abdomen.pdf", height = 9, width = 9)
effectsize::eta_squared(car::Anova(m1, type = 2), partial = FALSE, ci = 0.95)

plot_abdomen <- p1 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16)) +
    p2 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    p3 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    plot_layout(design = layout2) 

# Axilla
m2.group <- update(m2, Staphylococcus ~ group)
p1 <- sjPlot::plot_model(m2, show.values = T, show.p = T, sort.est = FALSE, vline.color = "grey45") + theme_user + ggtitle('')
p2 <- ggplot(cova %>% dplyr::filter(location == "axilla"), aes(CADESI, Staphylococcus)) +
    geom_point() + geom_smooth(method='lm') + theme_user + ggtitle('') +  ylab('Staphylococcus') + ylim(-0.1,1.4)
p3 <- sjPlot::plot_model(m2.group, type = "pred", terms = "group", show.values = T, ) + theme_user + xlab('') + ylab('Staphylococcus') + ggtitle('')
emmeans::emmeans(m2.group, pairwise ~ group)$contrasts
p4 <- flextable::flextable( emmeans::emmeans(m2.group, pairwise ~ group)$contrasts %>% as.data.frame() )

p1 + p2 + p3 + 
    grid::rasterGrob(flextable::as_raster(x = p4)) +
    plot_layout(design = layout) +
    plot_annotation(title = 'Axilla')
ggsave(filename = "plots/model_Staph_axilla.pdf", height = 9, width = 9)
effectsize::eta_squared(car::Anova(m2, type = 2), partial = FALSE, ci = 0.95)

plot_axilla <- p1 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16)) +
    p2 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    p3 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    plot_layout(design = layout2) 

# Elbow
m4.group <- update(m4, Staphylococcus ~ group)
p1 <- sjPlot::plot_model(m4, show.values = T, sort.est = FALSE, transform = NULL, vline.color = "grey45") + theme_user + ggtitle('')
p2 <- ggplot(cova %>% dplyr::filter(location == "elbow"), aes(CADESI, Staphylococcus)) +
    geom_point() + geom_smooth(method='lm') + theme_user + ggtitle('') +  ylab('Staphylococcus') + ylim(-0.1,1.4)
p3 <- sjPlot::plot_model(m4.group, type = "pred", terms = "group", show.values = T, ) + theme_user + xlab('') + ylab('Staphylococcus') + ggtitle('')
emmeans::emmeans(m4.group, pairwise ~ group)$contrasts
p4 <- flextable::flextable( emmeans::emmeans(m4.group, pairwise ~ group)$contrasts %>% as.data.frame() )

p1 + p2 + p3 + 
    grid::rasterGrob(flextable::as_raster(x = p4)) +
    plot_layout(design = layout) +
    plot_annotation(title = 'Elbow')
ggsave(filename = "plots/model_Staph_elbow.pdf", height = 9, width = 9)
effectsize::eta_squared(car::Anova(m4, type = 2), partial = FALSE, ci = 0.95)

plot_elbow <- p1 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16)) +
    p2 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    p3 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    plot_layout(design = layout2) 

# Flank
m5.group <- update(m5, Staphylococcus ~ group)
p1 <- sjPlot::plot_model(m5, show.values = T, show.p = T, sort.est = FALSE, transform = NULL, vline.color = "grey45") + theme_user + ggtitle('')
p2 <- ggplot(cova %>% dplyr::filter(location == "flank"), aes(CADESI, Staphylococcus)) +
    geom_point() + geom_smooth(method='lm') + theme_user + ggtitle('') +  ylab('Staphylococcus') + ylim(-0.1,1.4) 
p3 <- sjPlot::plot_model(m5.group, type = "pred", terms = "group", show.values = T, ) + theme_user + xlab('') + ylab('Staphylococcus') + ggtitle('')
emmeans::emmeans(m5.group, pairwise ~ group)$contrasts
p4 <- flextable::flextable( emmeans::emmeans(m5.group, pairwise ~ group)$contrasts %>% as.data.frame() )

p1 + p2 + p3 + 
    grid::rasterGrob(flextable::as_raster(x = p4)) +
    plot_layout(design = layout) +
    plot_annotation(title = 'Flank')
ggsave(filename = "plots/model_Staph_flank.pdf", height = 9, width = 9)
effectsize::eta_squared(car::Anova(m5, type = 2), partial = FALSE, ci = 0.95)

plot_flank <- p1 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16)) +
    p2 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    p3 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    plot_layout(design = layout2) 

# Front Paw
m6.group <- update(m6, Staphylococcus ~ group)
p1 <- sjPlot::plot_model(m6, show.values = T, show.p = T, sort.est = FALSE, transform = NULL, vline.color = "grey45") + theme_user + ggtitle('')
p2 <- ggplot(cova %>% dplyr::filter(location == "frontPaw"), aes(CADESI, Staphylococcus)) +
    geom_point() + geom_smooth(method='lm') + theme_user + ggtitle('') +  ylab('Staphylococcus') + ylim(-0.1,1.4)
p3 <- sjPlot::plot_model(m6.group, type = "pred", terms = "group", show.values = T, ) + theme_user + xlab('') + ylab('Staphylococcus') + ggtitle('')
emmeans::emmeans(m6.group, pairwise ~ group)$contrasts
p4 <- flextable::flextable( emmeans::emmeans(m6.group, pairwise ~ group)$contrasts %>% as.data.frame() )

p1 + p2 + p3 + 
    grid::rasterGrob(flextable::as_raster(x = p4)) +
    plot_layout(design = layout) +
    plot_annotation(title = 'Front Paw')
ggsave(filename = "plots/model_Staph_frontPaw.pdf", height = 9, width = 9)
effectsize::eta_squared(car::Anova(m6, type = 2), partial = FALSE, ci = 0.95)

plot_frontpaw <- p1 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16)) +
    p2 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    p3 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    plot_layout(design = layout2) 

# hindPaw
m8.group <- update(m8, Staphylococcus ~ group)
p1 <- sjPlot::plot_model(m8, show.values = T, sort.est = FALSE, transform = NULL, vline.color = "grey45") + theme_user + ggtitle('')
p2 <- ggplot(cova %>% dplyr::filter(location == "hindPaw"), aes(CADESI, Staphylococcus)) +
    geom_point() + geom_smooth(method='lm') + theme_user + ggtitle('') +  ylab('Staphylococcus') + ylim(-0.1,1.4) 
p3 <- sjPlot::plot_model(m8.group, type = "pred", terms = "group", show.values = T, ) + theme_user + xlab('') + ylab('Staphylococcus') + ggtitle('')
emmeans::emmeans(m8.group, pairwise ~ group)$contrasts
p4 <- flextable::flextable( emmeans::emmeans(m8.group, pairwise ~ group)$contrasts %>% as.data.frame() )

p1 + p2 + p3 + 
    grid::rasterGrob(flextable::as_raster(x = p4)) +
    plot_layout(design = layout) +
    plot_annotation(title = 'Hind Paw')
ggsave(filename = "plots/model_Staph_hindPaw.pdf", height = 9, width = 9)
effectsize::eta_squared(car::Anova(m8, type = 2), partial = FALSE, ci = 0.95)

plot_hindpaw <- p1 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16)) +
    p2 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    p3 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    plot_layout(design = layout2) 


# Lips
m9.group <- update(m9, Staphylococcus ~ group)
p1 <- sjPlot::plot_model(m9, show.values = T, sort.est = FALSE, transform = NULL, vline.color = "grey45") + theme_user + ggtitle('')
p2 <- ggplot(cova %>% dplyr::filter(location == "lips"), aes(CADESI, Staphylococcus)) +
    geom_point() + geom_smooth(method='lm') + theme_user + ggtitle('') +  ylab('Staphylococcus') + ylim(-0.1,1.4)
p3 <- sjPlot::plot_model(m9.group, type = "pred", terms = "group", show.values = T, ) + theme_user + xlab('') + ylab('Staphylococcus') + ggtitle('')
emmeans::emmeans(m9.group, pairwise ~ group)$contrasts
p4 <- flextable::flextable( emmeans::emmeans(m9.group, pairwise ~ group)$contrasts %>% as.data.frame() )

p1 + p2 + p3 + 
    grid::rasterGrob(flextable::as_raster(x = p4)) +
    plot_layout(design = layout) +
    plot_annotation(title = 'Lips')
ggsave(filename = "plots/model_Staph_lips.pdf", height = 9, width = 9)
effectsize::eta_squared(car::Anova(m9, type = 2), partial = FALSE, ci = 0.95)

plot_lips <- p1 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16)) +
    p2 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    p3 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    plot_layout(design = layout2) 


# palm
m10.group <- update(m10, Staphylococcus ~ group)
p1 <- sjPlot::plot_model(m10, show.values = T, show.p = T, sort.est = FALSE, transform = NULL, vline.color = "grey45") + theme_user + ggtitle('')
p2 <- ggplot(cova %>% dplyr::filter(location == "palm"), aes(CADESI, Staphylococcus)) +
    geom_point() + geom_smooth(method='lm') + theme_user + ggtitle('') +  ylab('Staphylococcus') + ylim(-0.1,1.4)
p3 <- sjPlot::plot_model(m10.group, type = "pred", terms = "group", show.values = T, ) + theme_user + xlab('') + ylab('Staphylococcus') + ggtitle('')
emmeans::emmeans(m10.group, pairwise ~ group)$contrasts
p4 <- flextable::flextable( emmeans::emmeans(m10.group, pairwise ~ group)$contrasts %>% as.data.frame() )

p1 + p2 + p3 + 
    grid::rasterGrob(flextable::as_raster(x = p4)) +
    plot_layout(design = layout) +
    plot_annotation(title = 'Palm')
ggsave(filename = "plots/model_Staph_palm.pdf", height = 9, width = 9)
effectsize::eta_squared(car::Anova(m10, type = 2), partial = FALSE, ci = 0.95)

plot_palm <- p1 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16)) +
    p2 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    p3 + geom_point(size = 4) + theme(axis.text.y = element_text(size = 16), axis.title.y = element_text(face = "italic", size = 16)) +
    plot_layout(design = layout2) 


# Plot interaction effect
# sjPlot::plot_model(m5, type = "int")

# Fusobacterium
summary(m14)
m14.group <- update(m14, Fusobacterium ~ group)
summary(m14.group)
emmeans::emmeans(m14.group, pairwise ~ group)$contrasts


# Pretty plotting ---------------------------------------------------------

layout3 <- "
AAA
BBB
CCC
DDD
EEE
FFF
GGG
HHH
"

p_model <- plot_abdomen /
    plot_axilla /
    plot_elbow /
    plot_flank /
    plot_frontpaw /
    plot_hindpaw /
    plot_lips /
    plot_palm +
    plot_layout(design = layout3) +
    plot_annotation(tag_levels = 'a')
ggsave(filename = "plots/FigS6.pdf", height = 27, width = 18, plot = p_model)

# Staphylococcus vs Neisseria ---------------------------------------------

ps <- readRDS(file = "data/phyloseq.OTU.RDS")
ps_clr <- microbiome::transform(ps, "clr")
ps_ra  <- transform_sample_counts(ps, function(x){x / sum(x)})

# ps.genus <- phyloseq::tax_glom(physeq = ps_ra,  taxrank = "Genus") 
ps.genus <- microbiome::aggregate_taxa(x = ps_ra, level = "Genus")

table(sample_data(ps.genus)$group)
table(sample_data(ps.genus)$location)

genus_table <- otu_table(ps.genus) %>% data.frame() %>% t()
colnames(genus_table) <- data.frame(tax_table(ps.genus))$unique

cova <- sample_data(ps.genus) %>% data.frame()
sum( rownames(cova) == rownames(genus_table) ) == nrow(cova)

sort(colnames(genus_table))

cova$Staphylococcus <- genus_table[, 'Staphylococcus']
cova$Fusobacterium <- genus_table[, 'Fusobacterium']
cova$Neisseria <- genus_table[, 'Neisseria']
cova$Neisseria_unk <- genus_table[, 'Bacteria_Proteobacteria_Betaproteobacteria_Neisseriales_Neisseriaceae_']

summary(cova$Staphylococcus)
summary(cova$Neisseria)
summary(cova$Neisseria_unk)

cova$log_Staph <- log1p(cova$Staphylococcus)
cova$log_Neiss <- log1p(cova$Neisseria)
cova$log_Neiss_un <- log1p(cova$Neisseria_unk)
cova$log_Neiss_total <- log1p(cova$Neisseria + cova$Neisseria_unk)

cova$Staph_Neis <- cova$log_Staph - cova$log_Neiss

cova$group <- factor(cova$group, levels = c("ADpre", "ADpost", "ADtreatment", "Healthy" ))


comp <- list(
    c('ADpre', 'ADpost'),
    c('ADpre', 'ADtreatment'),
    c('ADpre', 'Healthy'),
    c('ADpost', 'ADtreatment'),
    c('ADpost', 'Healthy'),
    c('ADtreatment', 'Healthy')
)

cova %>% 
    dplyr::filter(location != "stool") %>% 
    ggviolin(data = ., x = 'group', y = 'Staph_Neis', 
             trim = T, facet.by = 'location', 
             draw_quantiles = c(0.1, 0.5, 0.9), add = 'jitter') +
#    facet_wrap(~location, ncol = 3) +
    stat_compare_means() +
    stat_compare_means(comparisons = comp, method = "wilcox.test") +
    ylim(-1.,2) + ylab('log ratio Staphylococcus vs Neisseria') +
    xlab('') +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
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
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) 
ggsave(filename = "plots/Staph_Neisseria.pdf", height = 15, width = 15)

# Session Info ------------------------------------------------------------

sessioninfo::session_info()
