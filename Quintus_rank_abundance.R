library(phyloseq)
library(ggplot2)
library(cowplot)  # For plot_grid

# Normalize the counts
physeq_norm <- transform_sample_counts(physeq2, function(x) x / sum(x))
list_p <- list()

# Loop through different plant timepoints
for (plant_timex in c("seed", "T1_dom", "T2_dom", "T3_dom", "T1_wild", "T3_wild")) {
  
  # Subset phyloseq object by plant timepoint
  physeq_norm_sub <- subset_samples(physeq_norm, Plant_TimeX == plant_timex)
  
  # Convert phyloseq components to data frames
  otu_table_df <- as.data.frame(otu_table(physeq_norm_sub))
  tax_table_df <- as.data.frame(tax_table(physeq_norm_sub))
  sample_data_df <- as.data.frame(sample_data(physeq_norm_sub))
  
  # Function to rank taxa within each sample and keep the top 7
  rank_taxa <- function(sample_data) {
    ranked <- sort(sample_data, decreasing = TRUE)
    top_7 <- ranked[1:7]
    return(data.frame(OTU = names(top_7), Rank = 1:7, Abundance = top_7))
  }
  
  # Apply ranking function to OTU table
  ranked_abundance_df <- apply(otu_table_df, 2, rank_taxa)
  
  # Combine ranked taxa into a single data frame
  ranked_abundance_df <- do.call(rbind, lapply(names(ranked_abundance_df), function(sample) {
    cbind(Sample = sample, ranked_abundance_df[[sample]])
  }))
  
  # Merge with taxonomy information
  ranked_abundance_df <- merge(
    ranked_abundance_df, tax_table_df[, c("Family", "Genus", "Species")], 
    by.x = "OTU", by.y = "row.names", all.x = TRUE
  )
  
  # Add plant metadata from sample data
  ranked_abundance_df <- merge(
    ranked_abundance_df, sample_data_df[, "Plant", drop = FALSE], 
    by.x = "Sample", by.y = "row.names", all.x = TRUE
  )
  
  # Generate boxplot of ranked OTUs
  plot <- ggplot(ranked_abundance_df, aes(x = factor(Rank), y = Abundance)) +
    geom_boxplot() +
    labs(x = "Rank", y = "Relative Abundance") +
    ggtitle(plant_timex) +
    theme_bw()
  
  # Annotate p-values using Wilcoxon tests for consecutive ranks
  for (i in 1:6) {
    abundance_rank_i <- ranked_abundance_df$Abundance[ranked_abundance_df$Rank == i]
    abundance_rank_next <- ranked_abundance_df$Abundance[ranked_abundance_df$Rank == i + 1]
    
    # Skip if next rank is empty
    if (length(abundance_rank_next) == 0) next
    
    wilcox_test <- wilcox.test(abundance_rank_i, abundance_rank_next)
    p_value <- wilcox_test$p.value
    
    # Add p-value annotation to plot
    plot <- plot + 
      annotate("text", 
               x = i + 0.5,  # Position between two compared ranks
               y = max(ranked_abundance_df$Abundance, na.rm = TRUE) * 1.1 - (i * 0.05),  
               label = paste("p =", format(p_value, digits = 2)), 
               size = 3, color = "red")
  }
  
  list_p[[plant_timex]] <- plot
}

# Combine individual plots into a grid layout
p <- plot_grid(plotlist = list_p, ncol = 2)
print(p)

# Save plot as PDF
pdf("RankAbun.pdf", useDingbats = FALSE)
print(p)
dev.off()


# Identify dominant OTUs in T1_dom samples
t1 <- which(sample_data(physeq2)[, "Plant_TimeX"] == "T1_dom")
dom_otu_ind <- table(apply(otu_table(physeq2)[, t1], 2, which.max))
dom_otu_ind <- as.integer(names(dom_otu_ind))
dom_otus <- rownames(otu_table(physeq2)[dom_otu_ind, ])

# Retrieve taxonomy information for dominant OTUs
for (i in dom_otus) {
  print(tax_table(physeq2)[rownames(tax_table(physeq2)) == i, c("Genus", "Species")])
}
