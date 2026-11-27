# 000_config.R
# Central configuration for the GSE333704 project scripts
# Edit these values to customize behavior for your environment.

# Series / dataset identifier
series <- 'GSE333704'

# Verbose messaging
verbose <- TRUE

# Project root (default: project working directory)
project_root <- normalizePath('.', winslash = '/', mustWork = TRUE)

# Paths (relative to project_root)
raw_dir <- file.path(project_root, 'data', 'raw', series)
processed_dir <- file.path(project_root, 'data', 'processed', series)
results_dir <- file.path(project_root, 'results')
figures_dir <- file.path(project_root, 'figures')

# Download settings
download_max_attempts <- 2
download_methods <- c('auto','libcurl')

# QC thresholds (tweak as needed)
qc_min_features <- 200
qc_max_features <- 7500
qc_max_mt <- 20

# Marker finding parameters
marker_min_pct <- 0.25
marker_logfc_threshold <- 0.25

# Manifest file
manifest_path <- file.path(project_root, 'data', paste0(series, '_download_manifest.tsv'))

# Helper message function
msg <- function(...) if (isTRUE(verbose)) message(...)

# End of config
