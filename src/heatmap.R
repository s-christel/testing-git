pdf("sweASS_heatmaps.pdf", paper="a4")

library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(RColorBrewer)

#in terminal: less sweASS.taxonomy | cut -f 1,3 | tr ';' '\t' | tail -n +2 > sweASS.taxonomy_mod

#set environment
setwd("~/Documents/xiaofen_ass/phyloseq")
otufile="sweASS.otutable"
taxfile="sweASS.taxonomy_mod"
sampleinfo="sweASS_sample_info.clean.csv"




###DATA IMPORT

#make phyloseq OTU table
uppOTU <- read.table(otufile, header=T, row.names=1, sep="\t")
phyOTU <- otu_table(as.matrix(uppOTU), taxa_are_rows=T)

#make phyloseq TAX table
uppTAX <- read.table(taxfile, row.names=1, fill=T, sep="\t")
colnames(uppTAX) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
phyTAX <- tax_table(as.matrix(uppTAX))

#make phyloseq sample data table
uppSAMPLEDATA <- read.table(sampleinfo,row.names=1, header=T, sep=";", dec=".")
phySAMPLEDATA <- sample_data(uppSAMPLEDATA)

#make phyloseq object
phySWEASS <- phyloseq(phyOTU,phyTAX,phySAMPLEDATA)
phySWEASS




###GLOBAL SUBSETS

ss_ind <- phySWEASS %>%
            subset_samples(
            extraction=="indirect") %>%
            prune_taxa(taxa_sums(.) > 0, .)
            
ss_ind_bac <- ss_ind %>%
                subset_taxa(
                Kingdom == "Bacteria" &
                Family  != "mitochondria" &
                Class   != "Chloroplast"
                )
                
ss_ind_loc <- phySWEASS %>%
            subset_samples(
            extraction=="indirect" &
            (location=="Aneset" | location=="Flarkback")) %>%
            prune_taxa(taxa_sums(.) > 0, .)
            
ss_ind_loc_bac <- ss_ind_loc %>%
                subset_taxa(
                Kingdom == "Bacteria" &
                Family  != "mitochondria" &
                Class   != "Chloroplast"
                )

             
                
                

##group by phyla, filter phyla below 2 percent
#ss_ind_bac1 <- ss_ind_bac %>%
#                tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
#                transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
#                psmelt() %>%                                         # Melt to long format
#                filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
#                arrange(Phylum)                                      # Sort data frame alphabetically by phylum

#group by phyla, filter phyla below 1 percent
ss_ind_bac_phylum <- ss_ind_bac %>%
                    tax_glom(taxrank = "Phylum") %>%                     
                    transform_sample_counts(function(x) {x/sum(x)} ) %>% 
                    psmelt() %>%                                         
                    filter(Abundance > 0.01) %>%                        
                    arrange(Phylum)
                                    
#group by family, filter below 1 percent
ss_ind_bac_family <- ss_ind_bac %>%
                    tax_glom(taxrank = "Family") %>%                     
                    transform_sample_counts(function(x) {x/sum(x)} ) %>% 
                    psmelt() %>%                                         
                    filter(Abundance > 0.01) %>%                        
                    arrange(Phylum) 




###PLOTS

p.heatmap_phylum <- ggplot(data=ss_ind_bac_phylum, aes(x=sample_name, y=Phylum)) +
                        facet_grid(.~location, drop=T, scales="free_x") +
                        scale_x_discrete(expand=c(0,0))+
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Abundance") +
                        geom_tile(aes(fill=Abundance)) +
                        scale_fill_gradient(low="#FFFFCC", high="#FF0000", na.value="white") +
                        theme_bw(base_size=6) +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p.heatmap_phylum



p.heatmap_family <- ggplot(data=ss_ind_bac_family, aes(x=sample_name, y=Family)) +
                        facet_grid(.~location, drop=T, scales="free_x") +
                        scale_x_discrete(expand=c(0,0))+
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Family") +
                        geom_tile(aes(fill=Abundance), col="white") +
                        scale_fill_gradient(low="white", high="steelblue", na.value="white") +
                        theme_bw(base_size=6) +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p.heatmap_family










dev.off()