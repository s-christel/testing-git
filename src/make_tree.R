setwd("~/Documents/xiaofen_ass/phyloseq")
pdf("tree.pdf")

suppressMessages(library(msa))
suppressMessages(library(seqinr))
suppressMessages(library(phangorn))

                 

seqfile="test.seqs"


#read in sequences
uppSEQ <- readDNAStringSet(seqfile, format="fasta", use.names=T)

#align sequences (MSA)
uppSEQ_aligned <- msa(uppSEQ, "Muscle")

#convert alignment for other packages
seqinr_alignment <-msaConvert(uppSEQ_aligned, type="seqinr::alignment")
phangorn_alignment <- msaConvert(uppSEQ_aligned, type="phangorn::phyDat")
ape_alignment <- msaConvert(uppSEQ_aligned, type="ape::DNAbin")

#create Neighbor-joining tree as a template for max likelihood (APE)
distance_matrix <- dist.alignment(seqinr_alignment)
NJtree <- nj(distance_matrix)
plot(NJtree)

#fit the NJ tree with MLH model (PHANGORN)
MLH_fit <- pml(NJtree, phangorn_alignment)
plot(MLH_fit)

#optimise fit according to model
MLH_fit_JC <- optim.pml(MLH_fit, model="JC", rearrangement="stochastic")
plot(MLH_fit_JC)

#create bootstraps
bs <- bootstrap.pml(MLH_fit_JC, bs=100, optNni=T, multicore=T, control=pml.control(trace=0))

#plot tree with bootstraps
MLH_fit_JC_BS <- plotBS(midpoint(MLH_fit_JC$tree), bs, p=50, type="p")

#write out as newick
write.tree(MLH_fit_JC_BS, file="newick_tree.tree")

dev.off()
