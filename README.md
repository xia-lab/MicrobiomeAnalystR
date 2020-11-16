<p align="center">
  <img src="https://github.com/xia-lab/MicrobiomeAnalystR/blob/master/docs/microbiomeanalystr_logo.png" width="600">
</p>

### Note - MicrobiomeAnalystR is still under development - we cannot guarantee full functionality ### 

## Description 

**MicrobiomeAnalystR** contains the R functions and libraries underlying the popular MicrobiomeAnalyst web server, including > 200 functions for statistical, functional, and visual analysis of microbiome data. The package is synchronized with the MicrobiomeAnalyst web server. After installing and loading the package, users will be able to reproduce the same results from their local computers using the corresponding R command history downloaded from MicrobiomeAnalyst, thereby achieving maximum flexibility and reproducibility. With this R package we also aim to address an important gap left in its web version. Raw sequence data processing. 

To demonstrate this new functionality, we provide the "MicrobiomeAnalystR Workflow" vignette. In this vignette, we perform end-to-end microbiome data analysis on a subset of clinical IBD samples.   

## Getting Started 

### Step 1. Install package dependencies 

If you are using RStudio, ensure that it has been updated to the latest version for smoother installation of MicrobiomeAnalystR and overall better compatibility with all other R packages.

To use MicrobiomeAnalystR , first install all package dependencies. Ensure that you are able to download packages from Bioconductor - the Bioconductor package ("BiocManager") and RTools should be pre-installed. To install package dependencies, one can use the pacman R package (for those with >R 3.5.1). Note that some of these packages may require additional library dependencies that need to be installed prior to their own successful installation. For users who wish to perform raw sequence data processing, dada2 will also need to be installed. 

**MicrobiomeAnalystR will not be successfully installed until all package dependencies and their associated dependencies are also installed.** An example message signaling the R package installation failure is "non-zero exit status". The most common reason is that not all R package dependencies were successfully installed. If you are unable to run the pacman function, you will have to install each R package dependency one by one using **install.packages("x", dependencies = TRUE)** if the package is from CRAN or **BiocManager::install("x")** if the package is from Bioconductor. Note to know where the package is deposited, simply google the R package - i.e. "phyloseq R" will return the Bioconductor page where you can follow the installation instructions for that R package.

Note about Tax4Fun: We are in the progress of migrating to Tax4Fun2 - the older version is no longer supported on CRAN but can still be installed from nick-youngblut/Tax4Fun by ensuring the following dependencies are also installed - rhdf5; qiimer; joey711/biom. 

```R
install.packages("pacman")

library(pacman)

pacman::p_load(phyloseq, metacoder, pryr, biomformat, RColorBrewer, ggplot2, gplots, Cairo, igraph, BiocParallel, randomForest, metagenomeSeq, MASS, DESeq2, vegan, RJSONIO, ggfortify, pheatmap, xtable, genefilter, data.table, reshape, stringr, ape, grid, gridExtra, splitstackshape, edgeR, globaltest, R.utils, viridis, ggrepel, ppcor, qs)
```
### Step 2. Install the package 

MicrobiomeAnalystR is freely available from GitHub. Note that the Rpackage is currently under construction. The package documentation, including the vignettes for each module and user manual is available within the downloaded R package file. If all package dependencies were installed, you will be able to install the MicrobiomeAnalystR. There are three options, A) using the R package devtools, B) cloning the github, C) manually downloading the .tar.gz file. Note that the MicrobiomeAnalystR github will have the most up-to-date version of the package. 

#### Option A) Install the package directly from github using the *devtools* package. Open R and enter:

Due to issues with Latex, some users may find that they are only able to install MicrobiomeAnalystR without any documentation (i.e. vignettes). 

```R
# Step 1: Install devtools
install.packages("devtools")
library(devtools)

# Step 2: Install MicrobiomeAnalystR WITHOUT documentation
devtools::install_github("xia-lab/MicrobiomeAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))

# Step 2: Install MicrobiomeAnalystR WITH documentation
devtools::install_github("xia-lab/MicrobiomeAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))

```

#### Option B) Clone Github and install locally

The * must be replaced by what is actually downloaded and built. For instance, check your Downloads folder to see what tar.gz file was downloaded. So if you download MicrobiomeAnalystR_1.0.1.tar.gz, replace the * with the downloaded version number.  

```R
git clone https://github.com/xia-lab/MicrobiomeAnalystR.git
R CMD build MicrobiomeAnalystR
R CMD INSTALL MicrobiomeAnalystR_*.tar.gz

```

#### Option C) Manual download of MicrobiomeAnalystR.tar.gz and install locally - not yet available, stable release to come soon!!

Manually download the .tar.gz file from [here](https://www.dropbox.com/s/wk43rs9hswzypgt/MicrobiomeAnalystR_0.0.0.9000.tar.gz?dl=0). The * must be replaced by what is actually downloaded and built.  

```R
cd ~/Downloads
R CMD INSTALL MicrobiomeAnalystR_*.tar.gz

```

## Case Studies

### MicrobiomeAnalyst Workflow

To showcase how to utilize MicrobiomeAnalystR , we provide a detailed tutorial to perform a comprehensive end-to-end workflow from raw sequence data preprocessing to knowledge-based analysis. The dataset showcased in the tutorial consists of a subset of pediatric IBD stool samples obtained from the Integrative Human Microbiome Project Consortium (https://ibdmdb.org/). The tutorial is available inside the R package as a vignette.

## Tutorials

For detailed tutorials on how to use MicrobiomeAnalystR, please refer to the R package vignettes. These vignettes include a comprehensive tutorial introducing MicrobiomeAnalystR, four detailed step-by-step tutorials with example data for each of the main MetaboAnalytR  modules, and a case-study showcasing the end-to-end functionality of MicrobiomeAnalystR. Note, the functions below work only if the R package vignettes were built. 

Within R:
```R
vignette(package="MicrobiomeAnalystR")
```

Within a web-browser:
```R
browseVignettes("MicrobiomeAnalystR")
```

## Citation

MicrobiomeAnalystR has been developed by the [XiaLab](http://xialabresearch.com/) at McGill University. The original manuscript (web-based version) can be found [here](https://www.ncbi.nlm.nih.gov/pubmed/28449106). 

We encourage users to further develop the package to suit their needs. If you use the R package, please cite us: 

Dhariwal A, Chong J, Habib S, King IL, Agellon LB, Xia J. MicrobiomeAnalyst: a web-based tool for comprehensive statistical, visual and meta-analysis of microbiome data. Nucleic acids research. 2017 Jul 3;45(W1):W180-8.

*From within R:*

```R
citation("MicrobiomeAnalystR")
```

## Bugs or feature requests

To inform us of any bugs or requests, please open a new issue or send an email to #jasmine.chong@mail.mcgill.ca.

## MicrobiomeAnalystR History & Updates

11-16-2020 - Code update w. web + change files from .rds to .qs - users need to install qs R package now
02-24-2020 - Code update w. web + added note about usage
09-05-2019 - Bug fixing w. web
08-07-2019 - Added function to import SILVA annotated biom files (handling Domain in taxonomy)
07-11-2019 - Added volcano + dot plots for RNAseq analysis
07-08-2019 - Testing R code for local use + creating vignettes
07-03-2019 - Updating R code + documentation
06-22-2019 - Prepping R package for stable release

## MicrobiomeAnalystR TO DO

Add 3D visualizations using plotly

Add function to make multiple feature box plots

