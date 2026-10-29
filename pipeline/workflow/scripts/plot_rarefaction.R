# workflow/scripts/plot_rarefaction.R

# Setup Logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(biomformat)
library(stringr)
library(GUniFrac)
library(vegan)

# --- Inputs & Outputs ---
biom_file <- snakemake@input[["biom"]]
pdf_out   <- snakemake@output[["pdf"]]
table_out <- snakemake@output[["table"]]

# --- Read Data ---
message("Reading BIOM file: ", biom_file)
new <- read_biom(biom_file)
df <- as(biom_data(new), "matrix")
df <- as.data.frame(df)

# Adjust OTU ID (Last col -> First col logic)
df$OTU_id <- row.names(df)
# Reorder: Put OTU_id first, then samples
df_max <- df[, c(ncol(df), 1:(ncol(df)-1))]
row.names(df_max) <- NULL

# --- Taxonomy Parsing ---
META <- observation_metadata(new)

if(!is.null(META)) {
    tax_cols <- grep("taxonomy", names(META), value=TRUE, ignore.case=TRUE)
    if(length(tax_cols) > 0) {
        clean_tax <- function(x) {
            x <- str_remove(as.character(x), pattern = ".+ ?\\_") 
            x[is.na(x)] <- ""
            return(x)
        }
        META[tax_cols] <- lapply(META[tax_cols], clean_tax)
        if(length(tax_cols) > 1) {
            taxon <- do.call(paste, c(META[tax_cols], sep=";"))
        } else {
            taxon <- META[[tax_cols]]
        }
        df_max$taxonomy <- taxon
    } else {
        df_max$taxonomy <- "Unclassified"
    }
} else {
    df_max$taxonomy <- "Unclassified"
}

# Save Raw Table
write.table(df_max, file = table_out, sep = "\t", quote = F, row.names = TRUE, col.names = NA)

# --- Data Preparation ---
# Prepare numeric matrix for analysis (Remove ID and Taxonomy)
numeric_cols <- names(df_max)[!names(df_max) %in% c("OTU_id", "taxonomy")]

# [FIX 1] Protect single sample structure with drop=FALSE
otu_table <- df_max[, numeric_cols, drop = FALSE]
rownames(otu_table) <- df_max$OTU_id

# Convert to numeric
otu_table[] <- lapply(otu_table, as.numeric)
rownames(otu_table) <- df_max$OTU_id

# [FIX 2] Clean empty rows safely
otu_table <- otu_table[rowSums(otu_table) > 0, , drop=FALSE]

# --- Normalization Logic (Conditional) ---
method <- 1
level <- 0
normCutoff <- 1000 

if (level == 0) {
  min_sum <- min(colSums(otu_table))
} else {
  min_sum <- normCutoff
}

message("Target Normalization Depth: ", min_sum)

# [FIX 3] Smart Normalization Check
# Only normalize if we have multiple samples AND enough depth
if (ncol(otu_table) > 1 && min_sum > 0) {
    message("Multiple samples detected. Applying Normalization...")
    
    if (method == 0) {
        # Proportion-based (Result usually float, round for vegan)
        norm_otu_table <- t(min_sum * t(otu_table) / colSums(otu_table))
        norm_otu_table <- round(norm_otu_table) 
    } else {
        # GUniFrac Rarefy
        # Must transpose: GUniFrac expects Rows=Samples
        input_mat <- as.matrix(t(otu_table))
        
        tryCatch({
            r_out <- Rarefy(input_mat, depth = min_sum)
            # Transpose back: We want Rows=OTUs for internal consistency
            norm_otu_table <- t(as.data.frame(r_out$otu.tab.rff))
            message("Rarefaction successful.")
        }, error = function(e) {
            message("WARNING: Rarefaction failed (Error: ", e$message, "). Using raw data.")
            norm_otu_table <<- otu_table
        })
    }
} else {
    message("Single sample or zero depth detected. Skipping normalization to avoid crash.")
    norm_otu_table <- otu_table
}

# --- Plotting (Multi-Page PDF) ---
message("Generating Plots...")
pdf(file = pdf_out, width=10, height=8)

# Helper function for labeling curves
add_labels <- function(curves, otu_tab) {
    labelCutoff <- min(21, ncol(otu_tab))
    slope <- vector()
    SampleID <- vector()
    
    for (i in seq_along(curves)) {
        curve_i <- curves[[i]]
        richness <- 0
        if (length(curve_i) > 6) {
            y_end <- curve_i[length(curve_i)]
            y_prev <- curve_i[length(curve_i) - 5]
            x_end <- attr(curve_i, "Subsample")[length(curve_i)]
            x_prev <- attr(curve_i, "Subsample")[length(curve_i) - 5]
            richness <- (y_end - y_prev) / (x_end - x_prev)
        }
        slope <- c(slope, richness)
        SampleID <- c(SampleID, colnames(otu_tab)[i])
    }
    
    curvedf <- data.frame(SampleID, slope)
    ordered_vector <- order(as.numeric(curvedf$slope), decreasing = TRUE)
    
    for (i in 1:labelCutoff) {
        idx <- ordered_vector[i]
        N <- attr(curves[[idx]], "Subsample")
        lines(N, curves[[idx]], col="red")
        text(max(N), max(curves[[idx]]), SampleID[idx], cex=0.6, pos=4)
    }
}

# --- PAGE 1: RAW DATA ---
message("Plotting Page 1: Raw Data")
# Check if valid for plotting
if (ncol(otu_table) > 0 && nrow(otu_table) > 0) {
    # Vegan expects Rows=Samples
    curves_raw <- rarecurve(data.frame(t(otu_table)), step = 20, col = "black", lty = "solid", 
                            label = FALSE, xlab = "Total Reads", ylab = "Richness (OTUs)", 
                            main = "Rarefaction Curves (Raw Data)")
    add_labels(curves_raw, otu_table)
    abline(v = min_sum, col="blue", lty=2) # Show where normalization WOULD happen
    legend("bottomright", legend=c("Normalization Cutoff"), col="blue", lty=2)
}

# --- PAGE 2: NORMALIZED DATA ---
message("Plotting Page 2: Normalized Data")
if (ncol(norm_otu_table) > 0 && nrow(norm_otu_table) > 0) {
    # Vegan expects Rows=Samples
    curves_norm <- rarecurve(data.frame(t(norm_otu_table)), step = 20, col = "darkgrey", lty = "solid", 
                             label = FALSE, xlab = "Total Reads (Normalized)", ylab = "Richness (OTUs)", 
                             main = paste0("Rarefaction Curves (Normalized to ", min_sum, ")"))
    add_labels(curves_norm, norm_otu_table)
}

dev.off()
message("Done.")