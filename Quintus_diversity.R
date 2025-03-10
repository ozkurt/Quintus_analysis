library(phyloseq)
library(microbiomeutilities)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(biomeUtils)
library(microbiome)  # Required for evenness calculation

# Load phyloseq object
load("phyloseq.Rdata")

# Filter samples with more than 1000 reads and remove unidentified taxa
physeq_filtered <- prune_samples(sample_sums(physeq) > 1000, physeq)
physeq_filtered <- subset_taxa(physeq_filtered, Domain != "?")

# Normalize sample counts (relative abundance transformation)
physeq_norm <- transform_sample_counts(physeq_filtered, function(x) x / sum(x))

# Ordination using PCoA with Bray-Curtis distance
ord_norm <- ordinate(physeq_norm, "PCoA", "bray")

ggplot_ord <- plot_ordination(physeq_norm, ord_norm, shape='Time', color='HostType') +
  geom_point(size=5, alpha=0.6) +
  scale_color_manual(values = c("#00FA9A", "black", "#FFD700")) +  
  scale_shape_manual(values = c(3, 16, 15, 17)) +
  theme_classic()

ggplot_ord

# Rarefaction to minimum sample depth
rf <- min(sample_sums(physeq_filtered))
physeq_rrf <- rarefy_even_depth(physeq_filtered, rf, replace=TRUE, rngseed=699)

# Calculate Faith's Phylogenetic Diversity
meta_tib <- calculatePD(physeq_norm, justDF=TRUE, include_root=FALSE)

# Define pairwise comparisons
my_comparisons <- list(
  c("seed", "T1_dom"), c("T1_dom", "T2_dom"), c("T2_dom", "T3_dom"), 
  c("T1_dom", "T1_wild"), c("T3_dom", "T3_wild")
)

# Plot Phylogenetic Diversity
pdf("Phylogenetic_Diversity.pdf", useDingbats=FALSE, width=8, height=4)
ggplot(meta_tib, aes(x=Plant_TimeX, y=PD, fill=Plant_TimeX)) +
  geom_boxplot(alpha=0.4) + scale_fill_brewer(palette="Dark2") +
  xlab("") + ylab("Phylogenetic Diversity") + theme_bw() +
  stat_compare_means(label = "p.signif")
dev.off()

# Alpha Diversity Analysis
p <- plot_richness(physeq_rrf, 'Plant_TimeX', measures=c("Observed", "Shannon", "InvSimpson", "Chao1"))

# Generate Alpha Diversity Plots for Different Indices
alpha_diversity_plots <- function(data, metric, filename) {
  pdf(filename, useDingbats=FALSE, width=8, height=4)
  ggplot(data[data$variable == metric,], aes(x=Plant_TimeX, y=value, fill=Plant_TimeX)) +
    geom_boxplot(alpha=0.4) + scale_fill_brewer(palette="Dark2") +
    xlab("") + ylab(metric) + theme_bw() +
    stat_compare_means(label = "p.signif")
  dev.off()
}

alpha_diversity_plots(p$data, "InvSimpson", "InvSimpson.pdf")
alpha_diversity_plots(p$data, "Chao1", "Chao1.pdf")
alpha_diversity_plots(p$data, "Shannon", "Shannon.pdf")


# Beta Diversity Analysis: Distance to Centroid
physeq_dom <- subset_samples(physeq_filtered, HostType != "Wild")
dist_BC_dom <- distance(physeq_norm, "bray")
mod <- betadisper(dist_BC_dom, data.frame(sample_data(physeq_dom))$Time)

pdf("DistancetoCentroid_Time.pdf", useDingbats=FALSE, width=8, height=4)
boxplot(mod)
dev.off()

# Evenness Calculation using Pielou's Index
ev <- evenness(physeq_rrf, 'pielou')
data_merged <- data.frame(variable="Pielou", value=ev$pielou, meta=p$data$Time)

ggplot(data_merged, aes(x=meta, y=value, fill=meta)) +
  geom_boxplot() + scale_fill_brewer(palette="Pastel1") +
  xlab("") + ylab("Pielou's Evenness") + theme_bw()

dev.off()

# PERMANOVA Analysis (Bray-Curtis Dissimilarity)
metadata <- as(sample_data(physeq_norm), "data.frame")
aov <- adonis(distance(physeq_norm, method="bray") ~ Tissue, data=metadata)
aov$aov.tab

# Core Taxa Analysis
asv_ps <- core(physeq_norm, detection=0.0001, prevalence=0.50)
asv_ps <- format_to_besthit(asv_ps)

set.seed(14285)
p_v <- plot_abund_prev(asv_ps, label.core=TRUE, color="Phylum", mean.abund.thres=0.01,
                        mean.prev.thres=0.99, dot.opacity=0.7, label.size=3, label.opacity=1.0,
                        nudge.label=-0.15, bs.iter=999, size=20, replace=TRUE, label.color="#5f0f40")

pdf("SuppFig_CoreTaxa.pdf", useDingbats=FALSE, width=8, height=6)
p_v +
  geom_vline(xintercept=0.95, lty="dashed", alpha=0.7) +
  geom_hline(yintercept=0.01, lty="dashed", alpha=0.7) +
  scale_color_brewer(palette="Dark2")
dev.off()