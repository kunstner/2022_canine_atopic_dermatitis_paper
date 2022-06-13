# General Information -----------------------------------------------------

# Authors: Mirja Thomsen 
# Date:    2022-06-13

# Alpha diversity ---------------------------------------------------------

# Load libraries ----------------------------------------------------------

library(phyloseq)
library(tidyverse)
library(DivNet)
library(breakaway)
library(stats)


# Load data ---------------------------------------------------------------

ps <- readRDS(file = "phyloseq.OTU.RDS")


# DivNet estimate of Shannon for stool samples ----------------------------

ps_stool <- subset_samples(ps, location == "stool")

dv_stool <- DivNet::divnet(ps_stool, base = "Otu606") #using the most abundant OTU as base

combined_shannon <- ps_stool %>% sample_data %>% data.frame %>% 
    mutate(sample_names = rownames(.)) %>% 
    left_join(dv_stool$shannon %>% summary,
              by = "sample_names") %>% 
    dplyr::filter( !is.na(estimate) )

combined_shannon %>% group_by(group) %>% summarize(mean = mean(estimate), median = median(estimate))


# Hypothesis testing

breakaway::betta(chats = combined_shannon$estimate,
                 ses = combined_shannon$error,
                 X = model.matrix(~group, data = combined_shannon))$table


# Plot alpha diversity

combined_shannon %>% 
    ggpubr::ggviolin(data = ., x = "group", y = "estimate", width = 1,
                     trim = T, add = c("boxplot", "dotplot"), add.params = list(color="grey2",fill = "grey75", size=0.2, width = 0.4),
                     color = "group", palette = c("orange1","#d7191c", "steelblue1", "steelblue4")) +
    ylim(-0.15, 4) +
    ylab("DivNet estimate of Shannon") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=14),
          axis.text.y = element_text(size=14),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=14, vjust=2),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid"),
          legend.position = "none") +
    scale_color_manual(values= c("orange1","#d7191c", "steelblue1", "steelblue4"))

ggsave(filename="plots/Shannon_groups_stool.pdf", height = 5, width=4.5)


# Community-wise estimate

dv_stool_group <- ps_stool %>% DivNet::divnet( X = "group", base = "Otu606") #using the most abundant OTU as base


combined_shannon_group <- ps_stool %>% sample_data %>% data.frame %>% 
    mutate(sample_names = rownames(.)) %>% 
    left_join(dv_stool_group$shannon %>% summary,
              by = "sample_names") %>% 
    dplyr::filter( !is.na(estimate) )


# Hypothesis testing

DivNet::testDiversity(dv_stool_group, h0= "shannon")


# Plot community-wise alpha diversity

plot <- combined_shannon_group %>% group_by(group) %>% summarize(estimate = mean(estimate), error = mean(error))

ggplot(plot, aes(x=group, y=estimate, color=group)) + 
    geom_pointrange(aes(ymin=estimate-error, ymax=estimate+error))+
    ylim(-0.15,4) + 
    ylab("DivNet estimate of Shannon") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=14),
          axis.text.y = element_text(size=14),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=14, vjust=2),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid"),
          legend.position = "none") +
    scale_color_manual(values=c("orange1","#d7191c", "steelblue1", "steelblue4"))

ggsave(filename="plots/group-wise_Shannon_stool.pdf", height = 5, width=4)


# Shannon diversity for skin samples --------------------------------------

# Exemplarily for skin site abdomen

sample_data(ps)$Shannon <- phyloseq::estimate_richness(physeq = ps, measures = "Shannon")$Shannon

ps_abdomen <- subset_samples(ps, location == "abdomen")

sample_data(ps_abdomen) %>% as_tibble() %>% group_by(group) %>% summarise(Shannondiv = mean(Shannon), SD = sd(Shannon))

# Hypothesis testing

stats::kruskal.test(Shannon ~ group, data = as_tibble(sample_data(ps_abdomen)))

wilcox <- sample_data(ps_abdomen) %>% as_tibble %>% filter (group == "ADpre" | group == "healthy") #likewise for other group comparisons
stats::wilcox.test(Shannon ~ group, data = wilcox)


# Plot alpha diversity

sample_data(ps_abdomen) %>% as_tibble %>%
    ggpubr::ggviolin(data = ., x = "group", y = "Shannon", width = 1,
                     trim = T, add = c("boxplot", "dotplot"), add.params = list(color="grey2",fill = "grey75", size=0.3, width = 0.4),
                     color = "group", palette = c("orange1","#d7191c", "steelblue1", "steelblue4")) +
    ylim(-0.15, 6) +
    ggpubr::stat_compare_means() +
    ggtitle("Abdomen") +
    ylab("Estimated Shannon diversity") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=14),
          axis.text.y = element_text(size=14),
          axis.ticks.x=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=14, vjust=2),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid"),
          legend.position = "none") +
    scale_color_manual(values= c("orange1","#d7191c", "steelblue1", "steelblue4"))

ggsave(filename="plots/Shannon_groups_abdomen.pdf")
