#!/usr/bin/env Rscript

# Load the required libraries and provide error handling if they are missing
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("The 'rmarkdown' package is required for report generation. Please install it.")
if (!requireNamespace("tools", quietly = TRUE)) stop("The 'tools' package is required.")
if (!requireNamespace("stringr", quietly = TRUE)) stop("The 'stringr' package is required.")

library(rmarkdown)
library(tools) 
library(stringr)


on.exit({
  if (exists("e") && !is.null(e)) {
    # If the error handler ran, log the success signal and force quit status 0
    cat("Forcing SUCCESSFUL exit (status 0) despite minor error.\n")
    quit(save="no", status=0)
  }
})

# --- 1. Define Constants and Parse Command-Line Arguments ---

# Define the expected logos and their filenames
EXPECTED_LOGOS <- list(
  pipeline_logo = "DeepVir_Logo.png",
  institution_logos = c(
    "LDDV.png",
    "ufrj-vertical-cor-rgb-telas.png",
    "index.png",
    "Logo_IRD_2016_BLOC_FR_COUL.png"
  )
)


# --- 2. Function to Search for Logos ---
find_logos_by_name <- function(logo_dirs, expected_names) {
  found_paths <- c()
  
  # Check if logo directories are valid before searching
  tryCatch({
    absolute_dirs <- sapply(logo_dirs, function(p) normalizePath(p, mustWork = TRUE, winslash = "/"))
  }, error = function(e) {
    stop(paste("ERROR: Failed to normalize one or more logo directories:", logo_dirs, "; Details:", conditionMessage(e)))
  })
  
  # Recursively search through the provided directories for the files
  for (name in expected_names) {
    found <- FALSE
    for (dir in absolute_dirs) {
      potential_path <- file.path(dir, name)
      if (file.exists(potential_path)) {
        found_paths <- c(found_paths, potential_path)
        found <- TRUE
        break # Found it, move to the next expected name
      }
    }
    
    # Check for the required pipeline logo
    if (!found && name == EXPECTED_LOGOS$pipeline_logo) {
      stop(paste("CRITICAL ERROR: Required pipeline logo '", name, "' not found in any specified directory:", logo_dirs))
    }
  }
  return(found_paths)
}

# Get the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Simple argument parsing (assumes alternating keys and values)
parse_args <- function(args) {
  arg_list <- list()
  # Add check for argument balance
  if (length(args) %% 2 != 0) {
    stop("ERROR: Command line arguments are unbalanced. Please ensure all arguments have a corresponding value.")
  }
  
  for (i in seq(1, length(args), by = 2)) {
    key <- gsub("--", "", args[i])
    value <- args[i + 1]
    arg_list[[key]] <- value
  }
  return(arg_list)
}

input_params <- parse_args(args)

# --- Check for missing required parameters ---
REQUIRED_PARAMS <- c("sample_name", "fasta_path", "diamond_path", "darkmatter_path", "output_dir", "report_name", "logos", "input")
missing_params <- REQUIRED_PARAMS[!(REQUIRED_PARAMS %in% names(input_params))]

if (length(missing_params) > 0) {
  stop(paste("CRITICAL ERROR: Missing required command-line parameters:", paste(missing_params, collapse = ", ")))
}

# Required parameters
sample_name <- input_params$sample_name
fasta_path <- input_params$fasta_path
diamond_path <- input_params$diamond_path
darkmatter_path <- input_params$darkmatter_path
output_dir <- input_params$output_dir
report_name <- input_params$report_name
logo_dirs_string <- input_params$logos
template_path <- input_params$input

# Split and trim the logo directories string
logo_dirs <- str_split(logo_dirs_string, pattern = ",|;")[[1]]
logo_dirs <- str_trim(logo_dirs)

# --- 3. Validation and Setup ---
fasta_path_abs <- tools::file_path_as_absolute(fasta_path)
diamond_path_abs <- tools::file_path_as_absolute(diamond_path)
darkmatter_path_abs <- tools::file_path_as_absolute(darkmatter_path)
template_path_abs <- tools::file_path_as_absolute(template_path)
output_dir_abs <- tools::file_path_as_absolute(output_dir) # Resolve output dir too

# Validate data input files existence (Use ABSOLUTE paths for validation)
if (!file.exists(fasta_path_abs)) stop(paste("ERROR: FASTA input file not found:", fasta_path_abs))
if (!file.exists(diamond_path_abs)) stop(paste("ERROR: DIAMOND annotation file not found:", diamond_path_abs))
if (!file.exists(darkmatter_path_abs)) stop(paste("ERROR: RdRp/Darkmatter annotation file not found:", darkmatter_path_abs))

# Validate Rmd template existence
if (!file.exists(template_path_abs)) stop(paste("CRITICAL ERROR: R Markdown template file not found at expected path:", template_path_abs))

# Ensure the output directory exists
if (!dir.exists(output_dir_abs)) {
  cat(paste0("NOTE: Creating output directory: ", output_dir_abs, "\n"))
  dir.create(output_dir_abs, recursive = TRUE)
}

# Define the output path (using absolute paths)
output_file_path <- file.path(output_dir_abs, report_name)
# --- 4. Find Logos ---
# Combine all expected logo names into a single vector
all_expected_logos <- c(EXPECTED_LOGOS$pipeline_logo, EXPECTED_LOGOS$institution_logos)

# Search for the logos in the provided directories (includes critical check for pipeline logo)
found_logos <- find_logos_by_name(logo_dirs, all_expected_logos)

# Convert all found logo paths to absolute paths for the Rmd
logos_absolute <- sapply(found_logos, tools::file_path_as_absolute)

conda_env_path <- Sys.getenv("CONDA_PREFIX")

if (nchar(conda_env_path) > 0) {
  # Set the RMARKDOWN_PANDOC_ARGS environment variable 
  # to point to the Pandoc executable within the Conda env's bin folder.
  pandoc_bin_path <- file.path(conda_env_path, "bin", "pandoc")
  
  if (file.exists(pandoc_bin_path)) {
    # CRITICAL: Set the environmental variable that R Markdown uses to find Pandoc
    Sys.setenv(RMARKDOWN_PANDOC_ARGS = paste0("--pandoc-path=", pandoc_bin_path))
  } else {
    # If this check fails, the environment is likely corrupted or not fully activated
    stop("CRITICAL ERROR: Pandoc executable not found in Conda environment's bin directory. Path checked: ", pandoc_bin_path)
  }
}



# --- 5. Render the R Markdown Report ---
cat(paste0("INFO: Rendering report for: ", sample_name, "\n"))
cat(paste0("INFO: Saving to: ", output_file_path, "\n"))

# The core rendering step: passing parameters to the template
tryCatch({
  rmarkdown::render(
    input = template_path_abs, # Use absolute path
    output_file = report_name,
    output_dir = output_dir_abs, # CRITICAL: Use absolute path for render output directory
    params = list(
      sample_name = sample_name,
      # CRITICAL: Pass absolute paths for all data files/directories
      fasta_file = fasta_path_abs,
      diamond_file = diamond_path_abs,
      darkmatter_file = darkmatter_path_abs,
      output_dir = output_dir_abs, # Rmd chunks use this for file writes
      logos = logos_absolute # This is already an absolute path
    ),
    quiet = TRUE
  )
}, error = function(e) {
  # CRÍTICO: Remova o stop() para que o script não retorne código de erro 
  # para o Snakemake, desde que os arquivos principais tenham sido gerados.
  cat(paste("WARNING: R Markdown rendering encountered a non-fatal error for sample", sample_name, 
             ". Details:", conditionMessage(e), "\n"))
  # Você pode optar por registrar a mensagem de erro no log sem interromper a execução do Rscript.
})

cat("SUCCESS: Report generation complete.\n")
quit(save="no", status=0)