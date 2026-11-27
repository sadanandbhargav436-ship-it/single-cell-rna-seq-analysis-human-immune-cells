# 006_import_and_process_GSE333704.R
# Improved, safer import & processing pipeline for GSE333704
# Usage: source('scripts/006_import_and_process_GSE333704.R')

# Load central config if present
config_path <- file.path('scripts', '000_config.R')
if (file.exists(config_path)) {
  source(config_path)
} else {
  verbose <- TRUE
  msg <- function(...) if (isTRUE(verbose)) message(...)
}

# Set series and paths
if (!exists('series')) series <- 'GSE333704'
project_root <- if (exists('project_root')) project_root else normalizePath('.', winslash = '/', mustWork = TRUE)
processed_dir <- if (exists('processed_dir')) processed_dir else file.path(project_root, 'data', 'processed', series)
raw_dir <- if (exists('raw_dir')) raw_dir else file.path(project_root, 'data', 'raw', series)
results_dir <- if (exists('results_dir')) results_dir else file.path(project_root, 'results')
figures_dir <- if (exists('figures_dir')) figures_dir else file.path(project_root, 'figures')

# Create directories
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Check required packages
required <- c('Seurat', 'Matrix', 'tidyverse')
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop('Missing required packages: ', paste(missing_pkgs, collapse = ', '), '\nRun scripts/001_setup.R to install them.')
}

library(Seurat)
library(tidyverse)

proc_files <- list.files(processed_dir, full.names = TRUE)
raw_files <- list.files(raw_dir, full.names = TRUE)
seu <- NULL

# Helper: convert SingleCellExperiment if needed
try_as_seurat <- function(obj) {
  if (inherits(obj, 'Seurat')) return(obj)
  if (inherits(obj, 'SingleCellExperiment')) {
    msg('Converting SingleCellExperiment to Seurat')
    return(tryCatch(as.Seurat(obj, counts = 'counts'), error = function(e) { msg('Conversion failed: ', e$message); NULL }))
  }
  return(NULL)
}

# 1) R serialized files
r_files <- proc_files[grepl('\\.(rds|rda|RData)(\\.gz)?$', proc_files, ignore.case = TRUE)]
if (length(r_files) > 0) {
  f <- r_files[1]
  msg('Trying serialized file: ', f)
  if (grepl('\\.rds(\\.gz)?$', f, ignore.case = TRUE)) {
    obj <- tryCatch(readRDS(f), error = function(e) { msg('readRDS failed: ', e$message); NULL })
  } else {
    e <- new.env()
    ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(e) { msg('load() failed: ', e$message); FALSE })
    obj <- NULL
    if (ok) {
      objs <- ls(envir = e)
      seurat_obj <- objs[sapply(objs, function(n) inherits(get(n, envir = e), 'Seurat'))]
      if (length(seurat_obj) > 0) obj <- get(seurat_obj[1], envir = e)
      else if (length(objs) == 1) obj <- get(objs[1], envir = e)
    }
  }
  if (!is.null(obj)) {
    if (inherits(obj, 'Seurat')) seu <- obj
    else if (inherits(obj, 'SingleCellExperiment')) seu <- try_as_seurat(obj)
  }
}

# 2) 10X HDF5
if (is.null(seu)) {
  h5s <- proc_files[grepl('\\.(h5|hdf5)$', proc_files, ignore.case = TRUE)]
  if (length(h5s) > 0) {
    f <- h5s[1]
    msg('Trying HDF5: ', f)
    counts <- tryCatch(Read10X_h5(f), error = function(e) { msg('Read10X_h5 failed: ', e$message); NULL })
    if (!is.null(counts)) seu <- tryCatch(CreateSeuratObject(counts = counts, project = series), error = function(e) NULL)
  }
}

# 3) 10X mtx or directory (with automatic file standardization)
if (is.null(seu)) {
  all_files <- c(proc_files, raw_files)
  mtxs <- all_files[grepl('matrix\\.mtx', all_files, ignore.case = TRUE)]
  bcs <- all_files[grepl('barcodes(\\.tsv|\\.txt)', all_files, ignore.case = TRUE)]
  
  target_dir <- NULL
  if (length(mtxs) > 0) target_dir <- dirname(mtxs[1])
  else if (length(bcs) > 0) target_dir <- dirname(bcs[1])
  
  if (!is.null(target_dir)) {
    msg('Found 10X directory: ', target_dir)
    
    # Standardize filenames if GEO prefixed
    dir_f <- list.files(target_dir, full.names = TRUE)
    for (f in dir_f) {
      fn <- basename(f)
      if (grepl('barcodes', fn, ignore.case = TRUE) && !fn %in% c('barcodes.tsv.gz', 'barcodes.tsv')) {
        file.rename(f, file.path(target_dir, if(grepl('\\.gz$', fn)) 'barcodes.tsv.gz' else 'barcodes.tsv'))
      }
      if (grepl('features|genes', fn, ignore.case = TRUE) && !fn %in% c('features.tsv.gz', 'features.tsv', 'genes.tsv.gz', 'genes.tsv')) {
        file.rename(f, file.path(target_dir, if(grepl('\\.gz$', fn)) 'features.tsv.gz' else 'features.tsv'))
      }
      if (grepl('matrix', fn, ignore.case = TRUE) && !fn %in% c('matrix.mtx.gz', 'matrix.mtx')) {
        file.rename(f, file.path(target_dir, if(grepl('\\.gz$', fn)) 'matrix.mtx.gz' else 'matrix.mtx'))
      }
    }
    
    counts <- tryCatch(Read10X(data.dir = target_dir), error = function(e) { msg('Read10X failed: ', e$message); NULL })
    if (!is.null(counts)) {
      seu <- tryCatch(CreateSeuratObject(counts = counts, project = series), error = function(e) { msg('CreateSeuratObject failed: ', e$message); NULL })
    }
  }
}

# Strict Guard Check
if (is.null(seu) || !inherits(seu, 'Seurat')) {
  stop('No valid Seurat object constructed. Check data/processed/', series, ' for matrix.mtx.gz, barcodes.tsv.gz, and features.tsv.gz.')
}

msg('Seurat object ready. Proceeding with QC and analysis')

# Detect mitochondrial genes
feat_names <- rownames(seu)
if (is.null(feat_names) || length(feat_names) == 0) {
  stop('No feature names found in Seurat object.')
}

mt_pattern <- if (any(grepl('^MT-', feat_names))) '^MT-' else if (any(grepl('^mt-', feat_names))) '^mt-' else '^MT-'
msg('Using mitochondrial pattern: ', mt_pattern)
seu[['percent.mt']] <- PercentageFeatureSet(seu, pattern = mt_pattern)

# QC thresholds
min_features <- if (exists('qc_min_features')) qc_min_features else 200
max_features <- if (exists('qc_max_features')) qc_max_features else 7500
max_mt <- if (exists('qc_max_mt')) qc_max_mt else 20

n_before <- ncol(seu)
seu <- subset(seu, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < max_mt)
msg('Cells kept after QC: ', ncol(seu), ' (before: ', n_before, ')')

# Processing pipeline
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu, verbose = FALSE)

pca_ok <- tryCatch({ seu <- RunPCA(seu, verbose = FALSE); TRUE }, error = function(e) { msg('RunPCA failed: ', e$message); FALSE })

n_pcs <- 10
if (pca_ok && 'pca' %in% names(seu@reductions)) {
  avail <- ncol(Embeddings(seu, 'pca'))
  n_pcs <- min(n_pcs, avail)
}

umap_ok <- tryCatch({ seu <- RunUMAP(seu, dims = 1:n_pcs, verbose = FALSE); TRUE }, error = function(e) FALSE)
tryCatch({ seu <- FindNeighbors(seu, dims = 1:n_pcs); seu <- FindClusters(seu, resolution = 0.5) }, error = function(e) NULL)

# Save Seurat object
seu_path <- file.path(processed_dir, paste0(series, '_seurat_object.rds'))
saveRDS(seu, seu_path)
msg('Saved Seurat object: ', seu_path)

cluster_tab <- tibble(cell = colnames(seu), cluster = as.character(Idents(seu)))
write_tsv(cluster_tab, file.path(results_dir, paste0(series, '_clusters.tsv')))

# Markers
top_markers <- NULL
markers <- tryCatch({ FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) }, error = function(e) NULL)
if (!is.null(markers) && nrow(markers) > 0) {
  write_tsv(markers, file.path(results_dir, paste0(series, '_all_markers.tsv')))
  lfc_col <- c('avg_log2FC','avg_logFC','avg_log')[which(c('avg_log2FC','avg_logFC','avg_log') %in% names(markers))[1]]
  if (!is.na(lfc_col)) {
    top_markers <- markers %>% group_by(cluster) %>% slice_max(order_by = .data[[lfc_col]], n = 10)
    write_tsv(top_markers, file.path(results_dir, paste0(series, '_top10_markers_per_cluster.tsv')))
  }
}

# Plots
tryCatch({
  vln <- VlnPlot(seu, features = c('nFeature_RNA','nCount_RNA','percent.mt'), ncol = 3)
  ggsave(file.path(figures_dir, paste0(series, '_QC_violin.png')), plot = vln, width = 14, height = 6, dpi = 150)
}, error = function(e) NULL)

if (!is.null(top_markers) && nrow(top_markers) > 0) {
  top_genes <- unique(top_markers$gene)
  if (length(top_genes) > 50) top_genes <- top_genes[1:50]
  tryCatch({
    hm <- DoHeatmap(seu, features = top_genes) + ggtitle('Top markers heatmap')
    ggsave(file.path(figures_dir, paste0(series, '_markers_heatmap.png')), plot = hm, width = 14, height = 10, dpi = 150)
  }, error = function(e) NULL)
}

if (umap_ok) {
  tryCatch({
    up <- DimPlot(seu, reduction = 'umap', label = TRUE) + ggtitle(paste(series, 'UMAP'))
    ggsave(file.path(figures_dir, paste0(series, '_UMAP.png')), plot = up, width = 12, height = 9, dpi = 150)
  }, error = function(e) NULL)
}

writeLines(capture.output(sessionInfo()), con = file.path(results_dir, paste0(series, '_sessionInfo.txt')))
msg('Processing complete. Check results/ and figures/ for outputs.')  