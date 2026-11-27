
# 002_download_GSE333704.R
# Download supplementary data files for the GSE333704 dataset
# Run from the project root in RStudio:
# source("scripts/002_download_GSE333704.R")


# Load project configuration
config_path <- file.path("scripts", "000_config.R")

if (file.exists(config_path)) {
  source(config_path)
} else {
  verbose <- TRUE
  msg <- function(...) {
    if (isTRUE(verbose)) message(...)
  }
}


# Set default project and dataset information if not already defined
if (!exists("series")) {
  series <- "GSE333704"
}

if (!exists("project_root")) {
  project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
}


# Keep both naming styles available for compatibility with older scripts
if (exists("projectRoot") && !exists("project_root")) {
  project_root <- projectRoot
}

if (exists("project_root") && !exists("projectRoot")) {
  projectRoot <- project_root
}


# Set up directories for raw and processed data
if (!exists("raw_dir")) {
  if (exists("rawDir")) {
    raw_dir <- rawDir
  } else {
    raw_dir <- file.path(project_root, "data", "raw", series)
  }
}

if (!exists("processed_dir")) {
  if (exists("processedDir")) {
    processed_dir <- processedDir
  } else {
    processed_dir <- file.path(project_root, "data", "processed", series)
  }
}

# Maintain the camelCase names for scripts that still use them
rawDir <- raw_dir
processedDir <- processed_dir


# Create the required directories if they do not exist
if (!dir.exists(raw_dir)) {
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(processed_dir)) {
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
}


# Build the GEO supplementary files URL
base <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  paste0(substr(series, 1, 6), "nnn"),
  "/", series, "/suppl/"
)

msg("Checking supplementary files at: ", base)


# Retrieve the file listing, with a few retries in case of network issues
resp <- NULL
maxAttempts <- if (exists("download_max_attempts")) {
  download_max_attempts
} else {
  2
}

for (attempt in seq_len(maxAttempts)) {
  
  try({
    html <- readLines(base, warn = FALSE)
    resp <- paste(html, collapse = "\n")
    break
  }, silent = TRUE)
  
  msg("Could not retrieve the file listing (attempt ", attempt, "). Retrying...")
  Sys.sleep(1 * attempt)
}


# Create an empty manifest when no files are available to download
create_empty_manifest <- function(path) {
  
  manifest <- data.frame(
    filename = character(0),
    url = character(0),
    dest = character(0),
    status = character(0),
    message = character(0),
    timestamp = character(0),
    size_bytes = numeric(0),
    stringsAsFactors = FALSE
  )
  
  if (!dir.exists(dirname(path))) {
    dir.create(dirname(path), recursive = TRUE)
  }
  
  write.table(
    manifest,
    file = path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  return(manifest)
}


# Define where the download report will be saved
manifest_path <- if (exists("manifest_path")) {
  manifest_path
} else {
  file.path(projectRoot, "data", paste0(series, "_download_manifest.tsv"))
}


# Stop here if the GEO file listing could not be retrieved
if (is.null(resp)) {
  
  msg("Unable to access the supplementary file listing at: ", base)
  msg("The GEO server may be temporarily unavailable or access may be blocked.")
  msg("No files were downloaded.")
  
  manifest <- create_empty_manifest(manifest_path)
  
  msg("Empty download manifest saved to: ", manifest_path)
  
  return(invisible(manifest))
}


# Extract file links from the GEO directory listing
matches <- regmatches(
  resp,
  gregexpr('href="[^"]+"', resp, perl = TRUE)
)[[1]]

links <- character(0)

if (length(matches) > 0) {
  links <- gsub('href="', "", gsub('"$', "", matches))
  
  # Ignore links pointing to the parent directory
  links <- links[links != "../"]
}


# Handle datasets with no supplementary files
if (length(links) == 0) {
  
  msg("No supplementary files were found for ", series, ".")
  msg("The data may not have been uploaded to GEO or may be restricted.")
  
  manifest <- create_empty_manifest(manifest_path)
  
  msg("Empty download manifest saved to: ", manifest_path)
  
  return(invisible(manifest))
}


msg("Found ", length(links), " file link(s). Starting download process...")


# Identify files that are likely to contain processed data
is_processed_like <- function(filename) {
  
  grepl(
    "\\.(rds|rda)(\\.gz)?$|matrix|counts|processed|seurat|\\.h5$|\\.h5ad$|\\.loom$|\\.mtx$",
    filename,
    ignore.case = TRUE
  )
}


# Initialize a manifest to keep track of downloaded files
manifest <- data.frame(
  filename = character(0),
  url = character(0),
  dest = character(0),
  status = character(0),
  message = character(0),
  timestamp = character(0),
  size_bytes = numeric(0),
  stringsAsFactors = FALSE
)


# Download each supplementary file
for (link in links) {
  
  fileUrl <- if (grepl("^https?://", link)) {
    link
  } else {
    paste0(base, link)
  }
  
  fileName <- basename(fileUrl)
  
  
  # Skip links that point to directories rather than actual files
  if (grepl("/$", link) || fileName == series || !grepl("\\.", fileName)) {
    
    msg("Skipping directory or invalid file link: ", link)
    
    manifest <- rbind(manifest, data.frame(
      filename = fileName,
      url = fileUrl,
      dest = NA_character_,
      status = "skipped_dir",
      message = "Directory link or missing file extension",
      timestamp = as.character(Sys.time()),
      size_bytes = NA_real_,
      stringsAsFactors = FALSE
    ))
    
    next
  }
  
  
  # Save processed files separately from raw data
  dest <- if (is_processed_like(fileName)) {
    file.path(processedDir, fileName)
  } else {
    file.path(rawDir, fileName)
  }
  
  
  # Avoid downloading files that are already present
  if (file.exists(dest)) {
    
    msg("File already exists, skipping: ", fileName)
    
    manifest <- rbind(manifest, data.frame(
      filename = fileName,
      url = fileUrl,
      dest = dest,
      status = "exists",
      message = "File already exists",
      timestamp = as.character(Sys.time()),
      size_bytes = file.info(dest)$size,
      stringsAsFactors = FALSE
    ))
    
    next
  }
  
  
  msg("Downloading ", fileName, "...")
  
  downloaded <- FALSE
  last_err <- NULL
  
  
  # Try different download methods in case the default method fails
  methods_to_try <- if (exists("download_methods")) {
    download_methods
  } else {
    c("auto", "libcurl")
  }
  
  if (nzchar(Sys.which("wget"))) {
    methods_to_try <- c(methods_to_try, "wget")
  }
  
  if (nzchar(Sys.which("curl"))) {
    methods_to_try <- c(methods_to_try, "curl")
  }
  
  
  # Try each available method until the file downloads successfully
  for (m in methods_to_try) {
    
    try({
      
      utils::download.file(
        fileUrl,
        destfile = dest,
        mode = "wb",
        quiet = TRUE,
        method = m
      )
      
      downloaded <- TRUE
      
      msg("Successfully downloaded ", fileName, " using ", m, ".")
      
      break
      
    }, silent = TRUE)
    
    last_err <- paste0("Download failed using method: ", m)
  }
  
  
  # Record the result of the download attempt
  if (!downloaded) {
    
    msg("Download failed for ", fileName, ". Tried: ",
        paste(methods_to_try, collapse = ", "))
    
    manifest <- rbind(manifest, data.frame(
      filename = fileName,
      url = fileUrl,
      dest = dest,
      status = "failed",
      message = last_err,
      timestamp = as.character(Sys.time()),
      size_bytes = NA_real_,
      stringsAsFactors = FALSE
    ))
    
  } else {
    
    manifest <- rbind(manifest, data.frame(
      filename = fileName,
      url = fileUrl,
      dest = dest,
      status = "downloaded",
      message = "",
      timestamp = as.character(Sys.time()),
      size_bytes = file.info(dest)$size,
      stringsAsFactors = FALSE
    ))
  }
}


# Save the download report for future reference
if (!dir.exists(dirname(manifest_path))) {
  dir.create(dirname(manifest_path), recursive = TRUE)
}

write.table(
  manifest,
  file = manifest_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# Print the final download summary
msg("\nDownload process completed.")
msg("Download report saved to: ", manifest_path)
msg("Raw data directory: ", rawDir)
msg("Processed data directory: ", processedDir)