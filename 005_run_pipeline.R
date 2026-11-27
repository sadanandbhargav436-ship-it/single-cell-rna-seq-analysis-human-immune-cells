# 005_run_pipeline.R
# Convenience script to run the full pipeline: setup, download, process, render
# Usage: source('scripts/005_run_pipeline.R')

config_path <- file.path('scripts', '000_config.R')
if (file.exists(config_path)) {
  source(config_path)
} else {
  verbose <- TRUE
  msg <- function(...) if (isTRUE(verbose)) message(...)
}

if (!exists('series')) series <- 'GSE333704'

msg('Running full pipeline for series: ', series)

# 1) Setup
tryCatch({
  msg('1) Running setup (install packages)')
  source(file.path('scripts', '001_setup.R'))
}, error = function(e) msg('Setup failed: ', e$message))

# 2) Download
tryCatch({
  msg('2) Running downloader')
  source(file.path('scripts', '002_download_GSE333704.R'))
}, error = function(e) msg('Downloader failed: ', e$message))

# 3) Process (Fixed file name: 003_import_and_process_GSE333704.R)
tryCatch({
  msg('3) Running import and processing')
  source(file.path('scripts', '003_import_and_process_GSE333704.R'))
}, error = function(e) msg('Processing failed: ', e$message))

# 4) Render report
tryCatch({
  msg('4) Rendering report')
  source(file.path('scripts', '004_render_report_GSE333704.R'))
}, error = function(e) msg('Report render failed: ', e$message))

msg('Pipeline complete. Check results/, figures/, and data/ for outputs.')