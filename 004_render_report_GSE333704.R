# 004_render_report_GSE333704.R
# Render an HTML report from notebooks/004_report_GSE333704.Rmd
# Usage: source('scripts/004_render_report_GSE333704.R')

# Source config if present
config_path <- file.path('scripts', '000_config.R')
if (file.exists(config_path)) {
  source(config_path)
} else {
  verbose <- TRUE
  msg <- function(...) if (isTRUE(verbose)) message(...)
}

# Ensure rmarkdown package is available
if (!requireNamespace('rmarkdown', quietly = TRUE)) {
  msg('rmarkdown is required to render the report. Run scripts/001_setup.R to install it. Skipping report render.')
} else {
  rmd_in <- file.path('notebooks', '004_report_GSE333704.Rmd')
  
  if (!file.exists(rmd_in)) {
    msg('R Markdown template not found: ', rmd_in, '. Skipping report render.')
  } else {
    project_root <- if (exists('project_root')) project_root else normalizePath('.', winslash = '/', mustWork = TRUE)
    results_dir <- file.path(project_root, 'results')
    
    # Ensure project results directory exists
    if (!dir.exists(results_dir)) {
      dir.create(results_dir, recursive = TRUE)
      msg('Created directory: ', results_dir)
    }
    
    report_filename <- paste0('GSE333704_report_', Sys.Date(), '.html')
    
    # Render report, forcing knit_root_dir to project_root
    tryCatch({
      rmarkdown::render(
        input = rmd_in,
        output_file = report_filename,
        output_dir = results_dir,
        knit_root_dir = project_root,
        quiet = FALSE
      )
      msg('Rendered report to: ', file.path(results_dir, report_filename))
    }, error = function(e) {
      msg('Report render failed: ', e$message)
    })
  }
}