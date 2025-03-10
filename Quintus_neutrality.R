library(phyloseq)
library(cowplot)
library(vegan)
library(ggplot2)
library(data.table)
library(reltools)
library(stats4)
library(minpack.lm)
library(Hmisc)
library(tidyr)

# Convert phyloseq objects to data tables
DT_tab = data.table(otu_table(physeq.rrf.dom), keep.rownames = TRUE, key = "rn")
DT.taxa = as.data.table(tax_table(physeq.rrf.dom), keep.rownames = TRUE, key = "rn")

# Merge OTU and taxonomy data
DT.m = merge(DT_tab, DT.taxa)
OTU_B = t(otu_table(physeq.rrf.dom))
TAXA_B = tax_table(physeq.rrf.dom)

# Extract unique time points from metadata
time = names(table(sample_data(physeq.rrf.dom)[, "Time"]))

# Create list of phyloseq objects split by time
phy_list = list()
for (i in seq_along(time)) {
  phy_list[[i]] = subset_samples(physeq.rrf.dom, Time == time[i])
}

# Define consistent color palettes
Palette.phylum <- c(
  "Actinobacteriota" = "#1f77b4",
  "Firmicutes" = "#ff7f0e",
  "Proteobacteria" = "#d62728",
  "Other" = "gray"
)
Palette.shape.1 <- c(21, 16)

# Define plot theme
custom_theme <- theme(
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_rect(),
  axis.title.y = element_blank(),
  axis.title.x = element_blank(),
  panel.background = element_rect(fill = "white", colour = "black")
)

# Initialize data frames to store results
fit_stat_quintus = data.frame()
sncm.out_quintus = data.frame()

# Loop through each time point and fit the Sloan Neutral Model
for (i in seq_along(time)) {
  
  print(paste0("Computing SNCM for the ", time[[i]], " tissue ..."))
  
  # Filter non-zero mean taxa
  DT_tab$mean = rowMeans(DT_tab[, -1], na.rm = TRUE)
  LIST = DT_tab[DT_tab$mean != 0]$rn    
  
  # Fit Sloan Neutral Community Model (SNCM)
  sncm.out = fit_sncm(
    spp = OTU_B[rownames(sample_data(phy_list[[i]])), LIST],
    pool = OTU_B[rownames(sample_data(phy_list[c[1]])), LIST],
    taxon = TAXA_B[LIST, ]
  )
  
  # Save results to CSV
  write.csv(sncm.out, file = paste0("sncm.out_", time[[i]], ".csv"), sep = ",")
  
  # Process SNCM output
  DT.sncm = data.table(sncm.out$predictions, keep.rownames = TRUE, key = "rn")
  
  # Assign phylum categories
  ColPhylum = c("Actinobacteriota", "Firmicutes", "Proteobacteria")
  DT.sncm$ColPhylum = ifelse(
    DT.sncm$Phylum %in% ColPhylum, DT.sncm$Phylum, "Other"
  )
  
  # Classify taxa based on prediction fit
  DT.sncm$ShpClass = ifelse(
    DT.sncm$fit_class %in% c("Below prediction", "Above prediction"), "Div", "Prd"
  )
  
  # Store fit statistics
  fit.stat = sncm.out$fitstats
  fit_stat_quintus = rbind(fit_stat_quintus, fit.stat)
  sncm.out_quintus = rbind(
    sncm.out_quintus, 
    table(sncm.out$predictions$fit_class) / sum(table(sncm.out$predictions$fit_class))
  )
  colnames(fit_stat_quintus) = colnames(sncm.out$fitstats)
  
  # Generate plot
  pp1 = ggplot(DT.sncm, aes(x = log(p), y = freq)) +
    geom_jitter(aes(color = ColPhylum, shape = ShpClass), size = 1) +
    theme_bw() + custom_theme +
    scale_color_manual(values = Palette.phylum) +
    scale_shape_manual(values = c("Prd" = 1, "Div" = 16)) +
    geom_line(aes(x = log(p), y = freq.pred), color = "gray", alpha = 1) +
    geom_line(aes(x = log(p), y = pred.lwr), color = "gray", linetype = "longdash", alpha = 1) +
    geom_line(aes(x = log(p), y = pred.upr), color = "gray", linetype = "longdash", alpha = 1)
  
  # Add text labels for deviating taxa
  pp1 = pp1 + geom_text(
    data = DT.sncm[DT.sncm$ShpClass == "Div", ],
    aes(label = paste(Family, Genus, Species, sep = "_")),
    vjust = -0.5, size = 1
  )
  
  # Save plot to PDF
  pdf(paste0("NeutPlot", time[[i]], ".pdf"), useDingbats = FALSE, height = 3, width = 4)
  print(pp1)
  dev.off()
  
  print("######################################################################")
}
