# Single-cell RNA-seq Analysis of Human Immune Cells

This RStudio project provides a reproducible pipeline and example analysis for single-cell RNA-seq focused on human immune cells (Series: GSE333704). The repository contains numbered scripts, a downloader, an importer/analysis script that builds a Seurat object, and an RMarkdown report template.

Important: this README was updated after a project-wide review to accurately reflect the current scripts and their expected behavior. No other files were changed.

Project layout (important folders)
- scripts/: numbered pipeline scripts and helpers. Key scripts:
  - 000_config.R        : central configuration (series, paths, thresholds, manifest path)
  - 001_setup.R         : installs required R packages (SeuratDisk treated as optional)
  - 002_download_GSE333704.R : downloader (R) that writes a manifest and saves files to data/raw/ or data/processed/
  - 003_import_and_process_GSE333704.R : legacy importer (kept for reference)
  - 004_render_report_GSE333704.R : renders notebooks/004_report_GSE333704.Rmd to HTML
  - 005_run_pipeline.R  : convenience runner (setup -> download -> process -> render). Updated to call 006.
  - 006_import_and_process_GSE333704.R : improved, hardened importer/analysis (preferred by the runner)
- notebooks/: RMarkdown report templates (notebooks/004_report_GSE333704.Rmd)
- data/raw/: downloaded supplementary files per-series (e.g. data/raw/GSE333704)
- data/processed/: processed matrices and saved Seurat objects per-series (e.g. data/processed/GSE333704)
- results/: TSVs, marker tables, sessionInfo, and rendered HTML report(s)
- figures/: generated plots (UMAP, QC, heatmaps, dotplots)
- docs/: notes, usage guides

Quickstart (RStudio, project root)
1. Open the project (.Rproj) in RStudio.
2. Install required R packages (one-time):
   source('scripts/001_setup.R')
   - Note: SeuratDisk is configured as optional because on Windows it often requires Rtools to build from source. If you need SeuratDisk (h5seurat/h5ad conversion), install Rtools (https://cran.r-project.org/bin/windows/Rtools/) and then run install.packages('SeuratDisk') or use remotes::install_github().
3. Download the GEO supplementary files for GSE333704:
   source('scripts/002_download_GSE333704.R')
   - The downloader inspects the GEO supplement directory and writes a manifest at: data/GSE333704_download_manifest.tsv (or the manifest path in scripts/000_config.R). The downloader classifies files as "processed-like" (saved to data/processed/<series>/) or raw (saved to data/raw/<series>/).
4. Process / import and run the analysis (recommended):
   source('scripts/005_run_pipeline.R')
   - 005_run_pipeline.R runs the sequence: 001_setup.R -> 002_download_GSE333704.R -> 006_import_and_process_GSE333704.R -> 004_render_report_GSE333704.R
   - If you prefer to run steps manually, use the exact scripts above. After successful download, you can run only the importer:
     source('scripts/006_import_and_process_GSE333704.R')

What the importer expects (so it can build a Seurat object)
- One of the following in data/processed/<series>/ (or specified processed_dir in scripts/000_config.R):
  - An R-serialized Seurat or SingleCellExperiment object (.rds, .rda, .RData)
  - A 10X HDF5 file (.h5 or .hdf5) readable by Read10X_h5
  - 10X-style matrix files (matrix.mtx.gz or matrix.mtx plus features.tsv and barcodes.tsv)
- If only FASTQ are available, external quantification (CellRanger, STARsolo, kallisto|bustools, etc.) is required beforehand to create counts/matrices.

Typical outputs produced by the pipeline (paths relative to project root)
- Seurat object: data/processed/GSE333704/GSE333704_seurat_object.rds
- Cluster assignments: results/GSE333704_clusters.tsv
- All markers and top markers per cluster: results/GSE333704_all_markers.tsv and results/GSE333704_top10_markers_per_cluster.tsv
- Figures: figures/GSE333704_QC_violin.png, GSE333704_QC_feature_scatter.png, GSE333704_markers_heatmap.png, GSE333704_markers_dotplot.png, GSE333704_UMAP.png
- Rendered report (HTML): results/GSE333704_report_<DATE>.html
- Download manifest: data/GSE333704_download_manifest.tsv

Recent project observations (from a local run used to validate these scripts)
- The GEO supplemental archive GSE333704_RAW.tar contained processed 10X-style matrices (filtered_matrix.mtx.gz, barcodes, features) for multiple samples. The downloader marks these as processed-like and places them in data/processed/GSE333704/ when present.
- A truncated download can occur (partial tar). If you see tar/truncation errors, re-download the archive (scripts/002_download_GSE333704.R will skip existing files by default). Manual re-download or deleting the partial file then re-running the downloader is an effective fix.

Troubleshooting & common fixes
- "No suitable processed files found": check data/processed/GSE333704/ and data/raw/GSE333704/. Inspect the manifest: data/GSE333704_download_manifest.tsv to see what was downloaded or skipped.
- SeuratDisk install warnings on Windows: install Rtools (link above) before installing SeuratDisk. SeuratDisk is optional for most steps; only required for h5seurat/h5ad conversions.
- Read10X_h5 errors: ensure rhdf5 is installed and the system has HDF5 libs. On Windows, installing rhdf5 via BiocManager typically handles R-side pieces but system HDF5 may be needed.
- Truncated tar archive: delete the partial file in data/raw/GSE333704/ and re-run the downloader. The downloader writes a manifest so you can confirm file sizes.

Checklist to run the pipeline in RStudio
1. Open project in RStudio.
2. Optionally edit scripts/000_config.R to change series, paths, or thresholds.
3. Run: source('scripts/001_setup.R')  # installs packages (optional package warnings are informational)
4. Run: source('scripts/002_download_GSE333704.R')  # fetches supplement files and writes manifest
5. Inspect manifest: data/GSE333704_download_manifest.tsv
6. Run: source('scripts/005_run_pipeline.R')  # full pipeline (or source 006 directly after download)
7. Inspect outputs in results/ and figures/

Reproducibility
- Consider using renv::init() and renv::snapshot() to capture package versions. Commit renv.lock so collaborators can restore exact package versions.

Notes & limitations
- The scripts attempt to be "soft": they log errors and return rather than aborting the R session. This is intentional for interactive RStudio use.
- System-level dependencies (Rtools on Windows, libhdf5) are not installed by these scripts and must be provided by the user.
- The downloader uses a simple HTML parsing approach (readLines + regex). For very large or complex GEO pages, a more robust approach using httr/xml2 could be added (this would add package dependencies).

If you want further README polish
- I updated README.md only. If you want, I can add a short one-page usage file in docs/ (USAGE_GSE333704.md) summarizing the exact RStudio commands and the most common fixes — say the word and I will add it.

Maintainer / Contact
- Project scripts are intended to be runnable in RStudio. For help running the pipeline locally, provide the R console output (copy-paste) and I will assist.