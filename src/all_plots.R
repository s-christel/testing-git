pdf("sweASS_plots.pdf", paper="a4")

suppressMessages(library(ggplot2))
suppressMessages(library(vegan))
suppressMessages(library(dplyr))
suppressMessages(library(scales))
suppressMessages(library(grid))
suppressMessages(library(reshape2))
suppressMessages(library(phyloseq))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ape))
suppressMessages(library(genefilter))


#in terminal: less sweASS.taxonomy | cut -f 1,3 | tr ';' '\t' | tail -n +2 > sweASS.taxonomy_mod

#set environment
setwd("~/Documents/xiaofen_ass/phyloseq")
otufile="sweASS.otutable"
taxfile="sweASS.taxonomy_mod"
treefile="newick_tree.tree"
sampleinfo="sweASS_sample_info.clean.csv"




###DATA IMPORT

#make phyloseq OTU table
uppOTU <- read.table(otufile, header=T, row.names=1, sep="\t")
phyOTU <- otu_table(as.matrix(uppOTU), taxa_are_rows=T)

#make phyloseq TAX table
uppTAX <- read.table(taxfile, row.names=1, fill=T, sep="\t")
colnames(uppTAX) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
phyTAX <- tax_table(as.matrix(uppTAX))

#make sequence/tree object
phyTREE <- read.tree(treefile)

#make phyloseq sample data table
uppSAMPLEDATA <- read.table(sampleinfo,row.names=1, header=T, sep=";", dec=".")
phySAMPLEDATA <- sample_data(uppSAMPLEDATA)

#make phyloseq object
phySWEASS <- phyloseq(phyOTU,phyTAX,phySAMPLEDATA,phyTREE)
phySWEASS




###GLOBAL SUBSETS

phy_bac <- phySWEASS %>%
                subset_taxa(
                Kingdom == "Bacteria" &
                Family  != "mitochondria" &
                Class   != "Chloroplast") %>%
                prune_taxa(taxa_sums(.) > 0, .)
                phy_bac

phy_bac_ind <- phy_bac %>%
                subset_samples(
                extraction=="indirect") %>%
                prune_taxa(taxa_sums(.) > 0, .)
                phy_bac_ind               
               
phy_bac_ind_loc <- phy_bac_ind %>%
                subset_samples(
                (location=="Aneset" | location=="Flarkback")) %>%
                prune_taxa(taxa_sums(.) > 0, .)
                phy_bac_ind_loc
            


###LEVELS
percentage=0.01

ss_bac_ind_loc_phylum_per <- phy_bac_ind_loc %>%
                tax_glom(taxrank = "Phylum") %>%
                transform_sample_counts(function(x) {x/sum(x)} ) %>%
                filter_taxa(filterfun(kOverA(1, percentage)), TRUE) %>%
                psmelt %>% arrange(Phylum)

ss_bac_ind_loc_family_per <- phy_bac_ind_loc %>%
                tax_glom(taxrank = "Family") %>%
                transform_sample_counts(function(x) {x/sum(x)} ) %>%
                filter_taxa(filterfun(kOverA(1, percentage)), TRUE) %>%
                psmelt %>% arrange(Family)

ss_bac_ind_phylum_per <- phy_bac_ind %>%
                tax_glom(taxrank = "Phylum") %>%
                transform_sample_counts(function(x) {x/sum(x)} ) %>%
                filter_taxa(filterfun(kOverA(1, percentage)), TRUE) %>%
                psmelt %>% arrange(Phylum)

ss_bac_ind_family_per <- phy_bac_ind %>%
                tax_glom(taxrank = "Family") %>%
                transform_sample_counts(function(x) {x/sum(x)} ) %>%
                filter_taxa(filterfun(kOverA(1, percentage)), TRUE) %>%
                psmelt %>% arrange(Family)

#to remove small percentages from stacked bars:
#ss_bac_ind_loc_phylum_per_filt <- ss_bac_ind_loc_phylum_per %>% filter(Abundance > percentage)

list_ind_top100 <- phy_bac_ind %>% 
                transform_sample_counts(function(x) {x/sum(x)} ) %>%
                prune_taxa(names(sort(taxa_sums(.),TRUE)[1:100]),.) %>%
                psmelt() %>% .[1] %>% unique() %>% arrange(OTU)
lapply(list_ind_top100, write, file="ind_top100.otulist", append=F, ncolumns=1)


###PLOTS

#colours <- c("#9DCC00", "#0075DC","#F0A3FF",
#"#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00", 
#"#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00")  # 21

plot_tree(phy_bac_ind, size="abundance", color="location", label.tips="taxa_names")


colourCount = length(unique(ss_bac_ind_loc_phylum_per$Phylum))      #calculate number of colors needed for this plot
getPalette = colorRampPalette(brewer.pal(9, "Set1"))        #create color gradient with x colors
p <- ggplot(data=ss_bac_ind_loc_phylum_per, aes(x=sample_name, y=Abundance, fill=Phylum))                               +
                        coord_cartesian(ylim=c(0,1))                                                                    +
                        scale_y_continuous(breaks=seq(0,1,0.2))                                                         +
                        scale_fill_manual(values=getPalette(colourCount))                                               +
                        facet_grid(.~location, drop=T, scales="free_x")                                                 +
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Relative Abundance of Taxon > 1%")       +
                        geom_bar(stat="identity")                                                                       +
                        theme_bw(base_size=6)                                                                           +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))                                        +
                        theme(legend.position="bottom")                                     
                        p

colourCount = length(unique(ss_bac_ind_loc_family_per$Family))      #calculate number of colors needed for this plot
getPalette = colorRampPalette(brewer.pal(9, "Set1"))        #create color gradient with x colors
p <- ggplot(data=ss_bac_ind_loc_family_per, aes(x=sample_name, y=Abundance, fill=Family))                               +
                        coord_cartesian(ylim=c(0,1))                                                                    +
                        scale_y_continuous(breaks=seq(0,1,0.2))                                                         +
                        scale_fill_manual(values=getPalette(colourCount))                                               +
                        facet_grid(.~location, drop=T, scales="free_x")                                                 +
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Relative Abundance of Taxon > 1%")       +
                        geom_bar(stat="identity")                                                                       +
                        theme_bw(base_size=6)                                                                           +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))                                        +
                        theme(legend.position="bottom")                                     
                        p

colourCount = length(unique(ss_bac_ind_phylum_per$Phylum))      #calculate number of colors needed for this plot
getPalette = colorRampPalette(brewer.pal(9, "Set1"))        #create color gradient with x colors
p <- ggplot(data=ss_bac_ind_phylum_per, aes(x=sample_name, y=Abundance, fill=Phylum))                               +
                        coord_cartesian(ylim=c(0,1))                                                                    +
                        scale_y_continuous(breaks=seq(0,1,0.2))                                                         +
                        scale_fill_manual(values=getPalette(colourCount))                                               +
                        facet_grid(.~location, drop=T, scales="free_x")                                                 +
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Relative Abundance of Taxon > 1%")       +
                        geom_bar(stat="identity")                                                                       +
                        theme_bw(base_size=6)                                                                           +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))                                        +
                        theme(legend.position="bottom")                                     
                        p
                        
colourCount = length(unique(ss_bac_ind_family_per$Family))      #calculate number of colors needed for this plot
getPalette = colorRampPalette(brewer.pal(9, "Set1"))        #create color gradient with x colors
p <- ggplot(data=ss_bac_ind_family_per, aes(x=sample_name, y=Abundance, fill=Family))                               +
                        coord_cartesian(ylim=c(0,1))                                                                    +
                        scale_y_continuous(breaks=seq(0,1,0.2))                                                         +
                        scale_fill_manual(values=getPalette(colourCount))                                               +
                        facet_grid(.~location, drop=T, scales="free_x")                                                 +
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Relative Abundance of Taxpn > 1%")       +
                        geom_bar(stat="identity")                                                                       +
                        theme_bw(base_size=6)                                                                           +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))                                        +
                        theme(legend.position="bottom")                                     
                        p
                        
                        
                        
                        
#heatmap
p <- ggplot(data=ss_bac_ind_loc_phylum_per, aes(x=sample_name, y=Phylum))                           +
                        facet_grid(.~location, drop=T, scales="free_x")                             +
                        scale_x_discrete(expand=c(0,0))                                             +
                        ylim(rev(levels(ss_bac_ind_loc_phylum_per$Phylum)))                         +
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Taxon")              +
                        geom_tile(aes(fill=Abundance))                                              +
                        scale_fill_gradient(low="#FFFFCC", high="#FF0000", na.value="blue")         +
                        theme_bw(base_size=6)                                                       +
                        theme(panel.background=element_rect(fill="white"))                          +
                        theme(panel.grid=element_blank())                                           +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))                    
                        p

p <- ggplot(data=ss_bac_ind_loc_family_per, aes(x=sample_name, y=Family))                           +
                        facet_grid(.~location, drop=T, scales="free_x")                             +
                        scale_x_discrete(expand=c(0,0))                                             +
                        ylim(rev(levels(ss_bac_ind_loc_family_per$Family)))                         +
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Taxon")              +
                        geom_tile(aes(fill=Abundance))                                              +
                        scale_fill_gradient(low="#FFFFCC", high="#FF0000", na.value="blue")         +
                        theme_bw(base_size=6)                                                       +
                        theme(panel.background=element_rect(fill="white"))                          +
                        theme(panel.grid=element_blank())                                           +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))                    
                        p

p <- ggplot(data=ss_bac_ind_phylum_per, aes(x=sample_name, y=Phylum))                           +
                        facet_grid(.~location, drop=T, scales="free_x")                             +
                        scale_x_discrete(expand=c(0,0))                                             +
                        ylim(rev(levels(ss_bac_ind_phylum_per$Phylum)))                         +
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Taxon")              +
                        geom_tile(aes(fill=Abundance))                                              +
                        scale_fill_gradient(low="#FFFFCC", high="#FF0000", na.value="blue")         +
                        theme_bw(base_size=6)                                                       +
                        theme(panel.background=element_rect(fill="white"))                          +
                        theme(panel.grid=element_blank())                                           +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))                    
                        p
                        
p <- ggplot(data=ss_bac_ind_family_per, aes(x=sample_name, y=Family))                           +
                        facet_grid(.~location, drop=T, scales="free_x")                             +
                        scale_x_discrete(expand=c(0,0))                                             +
                        ylim(rev(levels(ss_bac_ind_family_per$Family)))                         +
                        labs(title="Swedish Acid Sulfate Soil", x="Sample", y="Taxon")              +
                        geom_tile(aes(fill=Abundance))                                              +
                        scale_fill_gradient(low="#FFFFCC", high="#FF0000", na.value="blue")         +
                        theme_bw(base_size=6)                                                       +
                        theme(panel.background=element_rect(fill="white"))                          +
                        theme(panel.grid=element_blank())                                           +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1))                    
                        p


























dev.off()