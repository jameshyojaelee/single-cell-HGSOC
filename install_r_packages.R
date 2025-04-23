# Set library path
.libPaths("/gpfs/commons/home/jameslee/R/x86_64-conda-linux-gnu-library/4.4")

# Set repositories
options(repos = c(
    CRAN = "https://cloud.r-project.org",
    BioCsoft = "https://bioconductor.org/packages/3.20/bioc",
    BioCann = "https://bioconductor.org/packages/3.20/data/annotation",
    BioCexp = "https://bioconductor.org/packages/3.20/data/experiment"
))

# Install BiocManager if not installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = "/gpfs/commons/home/jameslee/R/x86_64-conda-linux-gnu-library/4.4")

# Install required packages with force=TRUE to ensure installation
BiocManager::install("escape", lib = "/gpfs/commons/home/jameslee/R/x86_64-conda-linux-gnu-library/4.4", force = TRUE)
BiocManager::install("dittoSeq", lib = "/gpfs/commons/home/jameslee/R/x86_64-conda-linux-gnu-library/4.4", force = TRUE)
BiocManager::install("S4Vectors", lib = "/gpfs/commons/home/jameslee/R/x86_64-conda-linux-gnu-library/4.4", force = TRUE)
