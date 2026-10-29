# workflow/scripts/plot_rarefaction_basta.R

# Setup Logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(vegan)

# --- Inputs & Outputs ---
table_file <- snakemake@input[["table"]]
pdf_out    <- snakemake@output[["pdf"]]

# --- Read Data ---
message("Reading Count Table: ", table_file)
# Check.names=FALSE prevents R from replacing spaces in taxa names with dots
otu_table_raw <- read.table(table_file, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

# --- Pre-processing ---
# Vegan expects: Rows = Samples, Columns = Species/Taxa
# Our input is: Rows = Taxa, Columns = Samples
# So we transpose it (t)
otu_table <- as.data.frame(t(otu_table_raw))

# Ensure numeric
otu_table[] <- lapply(otu_table, as.numeric)

# Filter out empty samples (row sums > 0)
otu_table <- otu_table[rowSums(otu_table) > 0, , drop = FALSE]

# Check dimensions
n_samples <- nrow(otu_table)
n_taxa <- ncol(otu_table)

message("Samples found (non-empty): ", n_samples)
message("Taxa found: ", n_taxa)

if (n_samples == 0) {
    message("WARNING: No samples with valid BASTA hits found. Creating empty PDF.")
    pdf(file = pdf_out)
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "No BASTA hits found")
    text(1, 1, "No taxa identified in any sample.")
    dev.off()
    quit(save = "no")
}

# --- Plotting ---
message("Generating Rarefaction Curve...")
pdf(file = pdf_out, width=10, height=8)

# Calculate Curves
# step=20 makes it faster for large datasets
rarefactionCurve <- rarecurve(otu_table,
                              step = 50,
                              col = "black",
                              lty = "solid",
                              label = FALSE,
                              xlab = "Total Reads (Sequencing Depth)",
                              ylab = "Richness (Distinct Lineages)",
                              main = "BASTA Rarefaction Curves (LCA)")

# --- Labeling Logic (Top 20 steepest curves) ---
# This identifies samples with the highest diversity to label them
labelCutoff <- min(20, nrow(otu_table))
slope <- vector()
SampleID <- rownames(otu_table)

for (i in seq_along(rarefactionCurve)) {
  curve_i <- rarefactionCurve[[i]]
  
  # Calculate simple slope (richness / max_depth)
  if (length(curve_i) > 0) {
      max_y <- max(curve_i)
      max_x <- attr(curve_i, "Subsample")[length(curve_i)]
      slope_val <- max_y / max_x
  } else {
      slope_val <- 0
  }
  slope <- c(slope, slope_val)
}

# Sort by slope (diversity rate)
curvedf <- data.frame(SampleID, slope)
ordered_indices <- order(curvedf$slope, decreasing = TRUE)

# Highlight and Label top samples
for (i in 1:labelCutoff) {
  idx <- ordered_indices[i]
  # Get the coordinates for this specific curve
  curve_values <- rarefactionCurve[[idx]]
  curve_steps <- attr(curve_values, "Subsample")
  
  # Draw red line
  lines(curve_steps, curve_values, col="red", lwd=1.5)
  
  # Add text label at the end of the curve
  text(max(curve_steps), max(curve_values), 
       SampleID[idx], 
       cex=0.7, pos=4, col="red")
}

dev.off()
message("Done.")