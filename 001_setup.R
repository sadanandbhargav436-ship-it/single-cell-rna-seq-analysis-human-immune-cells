# 001_setup.R
# setup script for scRNA-seq stuff
# Reminder: Run this from the project root in RStudio -> source('scripts/001_setup.R')
# TODO: clean this up later or split into parts if it gets too long

# let's grab the config file if it's there
config_file <- "scripts/000_config.R"
if (file.exists(config_file)) {
  source(config_file)
} else {
  # just print a warning, not a big deal
  print("No config file found, using defaults...")
}

# quick logging thing so I don't have to type message() every time
if (!exists('msg')) {
  msg <- function(...) {
    message(...)
  }
}

# packages we absolutely need
my_cran_pkgs <- c('tidyverse', 'Seurat', 'Matrix', 'rmarkdown')
maybe_cran <- c('SeuratDisk') # sometimes acts up on Windows
bioc_pkgs <- c('SingleCellExperiment', 'scater', 'scran', 'rhdf5')

# setting CRAN mirror just in case
options(repos = c(CRAN = 'https://cran.rstudio.com'))

# Helper function to install stuff
install_packages_manually <- function(pkg_list, custom_installer = install.packages) {
  for (pkg in pkg_list) {
    # check if installed
    already_here <- requireNamespace(pkg, quietly = TRUE)
    
    if (already_here == FALSE) {
      msg('Hey, installing ', pkg, ' for you...')
      
      # wrap in tryCatch because packages fail all the time
      tryCatch({
        custom_installer(pkg)
      }, error = function(err) {
        msg('Oops, failed to install ', pkg, ': ', err$message)
      })
      
    } else {
      msg(pkg, ' is already here, skipping.')
    }
  }
}

# 1. get CRAN packages first
msg('--- Installing CRAN packages ---')
install_packages_manually(my_cran_pkgs, install.packages)

# 2. Check Optional CRAN packages (manually looped because why not)
for (p in maybe_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    msg(paste0("Warning: Optional package '", p, "' is not installed."))
    msg("Note: I remember that this is used for converting between HDF5/AnnData formats (h5seurat/h5ad).")
    msg('If you are on Windows, remember you might need Rtools to build it from source!')
    msg('Check out: https://cran.r-project.org/bin/windows/Rtools/ for the right version.')
    msg(paste0("To install later, try: install.packages('", p, "') or use remotes::install_github('mojaveazure/seurat-disk')\n"))
  } else {
    msg(paste0("Awesome, optional package '", p, "' is already installed."))
  }
}

# 3. Bioconductor stuff
# Need BiocManager first
if (!requireNamespace('BiocManager', quietly = TRUE)) {
  msg('BiocManager missing. Trying to install it...')
  tryCatch({
    install.packages('BiocManager')
  }, error = function(e) {
    msg('Oh no, BiocManager install failed: ', e$message)
  })
}

# Now install the Bioc packages using BiocManager
if (requireNamespace('BiocManager', quietly = TRUE)) {
  msg('--- Installing Bioconductor packages ---')
  # using an anonymous function here just to match the style
  install_packages_manually(bioc_pkgs, function(x) {
    BiocManager::install(x, ask = FALSE, update = FALSE)
  })
} else {
  msg('Still no BiocManager... skipping Bioconductor packages. Install it manually!')
}

# 4. renv check
has_renv <- requireNamespace('renv', quietly = TRUE)
if (has_renv == FALSE) {
  msg('\nTip: You really should use renv for managing your package versions!')
  msg('Run: install.packages("renv") and then renv::init()')
} else {
  msg('\nrenv is detected. Good job keeping things reproducible!')
}

msg('\nAll done! If you used renv, don\'t forget to run renv::snapshot() now.')