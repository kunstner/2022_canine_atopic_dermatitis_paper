
# General Information -----------------------------------------------------

# Authors: Axel Künstner
# Date:    2022-05-31

# Libraries ---------------------------------------------------------------

# for data handling
library(tidyverse)
library(phyloseq)

# for estimates and analysis
library(selbal)

# Get data ----------------------------------------------------------------

ps <- readRDS(file = "data/phyloseq.OTU.RDS")

sample_data(ps) %>% .$location %>% unique %>% sort
sample_data(ps) %>% .$group %>% unique

my_selbal <- list()

# Selbal on stool ---------------------------------------------------------

run_selbal <- function(Location, TaxRank) {
    
    print( paste0("Run selbal on ", Location, " ", TaxRank) )
    
    ps_sel <- ps %>% phyloseq::subset_samples(physeq = ., location == Location) %>%
        phyloseq::subset_samples(physeq = ., group != "ADpost") %>% 
        phyloseq::tax_glom(taxrank = TaxRank) %>% 
        phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/3) 
    
    selbal_dat <- ps_sel %>%
        otu_table() %>% t() %>% data.frame()
    colnames(selbal_dat) <- ps_sel %>% 
        tax_table() %>% data.frame() %>% .[, TaxRank]
    
    y <- ps_sel %>%
        phyloseq::sample_data() %>% .$PVAS
    y[is.na(y)] <- 0
    
    selbal <- selbal.cv(x = as.matrix(selbal_dat), y = y, maxV = 12,
                        n.fold = 5, n.iter = 10,
                        zero.rep = "one",
                        logit.acc = 'AUC')
    pdf(file = paste0("plots/check_selbal_", Location, "_", TaxRank, ".pdf"), width = 8, height = 6)
    grid.draw( selbal$global.plot )
    dev.off()
    # plot.tab(selbal$cv.tab)
    # selbal$glm
    # summary(selbal$glm)$coefficients %>% data.frame()
    # selbal$cv.tab %>% data.frame()
    
    return(list(
        summary(selbal$glm)$coefficients %>% data.frame(), 
        selbal$cv.tab %>% data.frame())
    )
}

for( i in sort( unique( sample_data(ps)$location) ) ) {
    Location <- i
    my_selbal_tmp <- run_selbal(Location = i, TaxRank = "Phylum")
    my_selbal[[paste(Location, "Phylum", 'GLM',  sep = " ")]] <- my_selbal_tmp[[1]]
    my_selbal[[paste(Location, "Phylum", 'TAB',  sep = " ")]] <- my_selbal_tmp[[2]]
    
    my_selbal_tmp <- run_selbal(Location = i, TaxRank = "Genus")
    my_selbal[[paste(Location, "Genus", 'GLM',  sep = " ")]] <- my_selbal_tmp[[1]]
    my_selbal[[paste(Location, "Genus", 'TAB',  sep = " ")]] <- my_selbal_tmp[[2]]
    # stop()
}

WriteXLS::WriteXLS(x = my_selbal, ExcelFileName = "Selbal_check.xlsx", row.names = TRUE)

# Session Info ------------------------------------------------------------

sessioninfo::session_info()
