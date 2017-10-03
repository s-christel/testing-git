pdf("sweASS_stackedbars.pdf", paper="a4")

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
            extraction=="indirect" &
            (location=="Aneset" | location=="Flarkback")) %>%
            prune_taxa(taxa_sums(.) > 0, .)

ss_ind_bac <- ss_ind %>%
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

#group by genus, filter below 1 percent
ss_ind_bac_genus <- ss_ind_bac %>%
                    tax_glom(taxrank = "Genus") %>%                     
                    transform_sample_counts(function(x) {x/sum(x)} ) %>% 
                    psmelt() %>%                                         
                    filter(Abundance > 0.01) %>%                        
                    arrange(Genus) 




###PLOTS

#colours <- c("#9DCC00", "#0075DC","#F0A3FF",
#"#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00", 
#"#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00")  # 21



colourCount = length(unique(ss_ind_bac_phylum$Phylum))      #calculate number of colors needed for this plot
getPalette = colorRampPalette(brewer.pal(9, "Set1"))        #create color gradient with x colors
p.stackedbar_phyla <- ggplot(data=ss_ind_bac_phylum, aes(x=sample_name, y=Abundance, fill=Phylum)) +
                        coord_cartesian(ylim=c(0,1)) +
                        scale_y_continuous(breaks=seq(0,1,0.2)) +
                        scale_fill_manual(values=getPalette(colourCount)) +
                        facet_grid(.~location, drop=T, scales="free_x") +
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Relative Abundance of Phyla > 1%") +
                        geom_bar(stat="identity") +
                        theme_bw(base_size=6) +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                        theme(legend.position="bottom")                                     
p.stackedbar_phyla



colourCount = length(unique(ss_ind_bac_family$Family))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
p.stackedbar_family <- ggplot(data=ss_ind_bac_family, aes(x=sample_name, y=Abundance, fill=Family)) +
                        coord_cartesian(ylim=c(0,1)) +
                        scale_y_continuous(breaks=seq(0,1,0.2)) +
                        scale_fill_manual(values=getPalette(colourCount)) +
                        facet_grid(.~location, drop=T, scales="free_x") +
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Relative Abundance of Families > 1%") +
                        geom_bar(stat="identity") +
                        theme_bw(base_size=6) +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                        theme(legend.position="bottom")                                         
p.stackedbar_family



colourCount = length(unique(ss_ind_bac_genus$Genus))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
p.stackedbar_genus <- ggplot(data=ss_ind_bac_genus, aes(x=sample_name, y=Abundance, fill=Genus)) +
                        coord_cartesian(ylim=c(0,1)) +
                        scale_y_continuous(breaks=seq(0,1,0.2)) +
                        scale_fill_manual(values=getPalette(colourCount)) +
                        facet_grid(.~location, drop=T, scales="free_x") +
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Relative Abundance of Genus > 1%") +
                        geom_bar(stat="identity") +
                        theme_bw(base_size=6) +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                        theme(legend.position="bottom")                                         
p.stackedbar_genus



dev.off()