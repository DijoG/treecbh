# treecbh 🌳

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue.svg)](https://www.r-project.org/)
[![platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)]()

The **treecbh** package provides functions to detect individual tree level Crown Base Height (CBH) using high-resolution LiDAR data. 

Individual tree segmentation must be conducted prior. The package is meant to be used within the framework of the **lidR** package. 

***Update (29.10.2025):*** Cumulative percentile-based density method implemented instead of the original kernel density approach.

 - **Robust**: works with various tree shapes
 - **Conservative**: 5th percentile avoids false positives from low outliers
 - **Simple**: much easier to understand and modify (compared to the original 2D kernel method)

## Original paper
*Testing treecbh in Central European forests: an R package for crown base height detection using high-resolution aerial laser-scanned data*

https://doi.org/10.1093/forestry/cpae044

## System Requirements

- **R** version 4.0.0 or higher
- **CloudCompare** (for tree isolation functionality)
- **Recommended**: 8GB+ RAM for processing large point clouds

## Required packages

tidyverse, lidR, RCSF, sf, data.table, crayon, fpc, geometry, gtools

## Installation

```r
devtools::install_github("DijoG/treecbh")
```
## Example
This demonstration uses low point-density data and shows how to use **treecbh** combined with **lidR**. 

### Data preparation

```r
require(lidR);require(tidyverse)

# Forest point cloud (low resolution example data)
LAS <- system.file("extdata", "MixedConifer.laz", package = "lidR")
Alas <- lidR::readLAS(LAS, filter = "-drop_z_below 0") 

# Black, white and green color palette for visualizing CHM
bgcol <- function(x)
{
  col <- grDevices::colorRampPalette(c("grey1", "white", "forestgreen"))
  return(col(x))
}
```
Computing canopy height model using the pitfree algorithm.
```r
CHM <- lidR::rasterize_canopy(Alas, 0.5, pitfree(subcircle = 0.25))
plot(CHM, main = "CHM 0.5 pitfree", col = bgcol(50))
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/01_chm_pitfree.png">

Computing treetops using a constant windows size of 5 m.
```r
ws <- 5
treetops <- lidR::locate_trees(CHM, lmf(ws))
plot(CHM, main = "CHM 0.5 pitfree", col = bgcol(50))
plot(sf::st_geometry(treetops), add = T, pch = "+", col = "firebrick3")
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/02_chm_pitfree_treetops.png">

Point-based crown segmentation using the Dalponte algorithm.
```r
algo_crowns <- lidR::dalponte2016(CHM, treetops)
las_crowns <- lidR::segment_trees(Alas, algo_crowns, attribute = "ID")
```
Obtaining Individual Tree Segments (ITS) as polygons.
```r
ITS <- lidR::crown_metrics(las_crowns, attribute = "ID", geom = "concave", func = NULL)
plot(sf::st_geometry(ITS), reset = FALSE, col = "forestgreen", border = "grey80")
plot(sf::st_geometry(treetops), add = T, pch = "+", col = "firebrick3")
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/03_its_treetops.png">

Some trees were not segmented.

Extracting tree point clouds to ITS (classifying ground points and saving tree point cloud files by default).
```r
oudir <- "<path/to/directory>"
treecbh::get_3DTREE(Alas, ITS, output_dir = oudir, FEATURE= "ID")
```

### CBH detection 

The following functions and functionality apply to the usage of high point-density LiDAR data:

```r
# List las files
las_l <- list.files(oudir, pattern = ".las", full.names = T) %>%
  gtools::mixedsort()
  
# Get automated parameters
tictoc::tic()
params <- treecbh::get_PARAMS(oudir)
tictoc::toc()
# 25.84 sec

# Run
# 1) Optimization activated (default cbh_ONLY = 3): executing cbh detection only
O_CBH <- treecbh::get_CBH(list_LAS_char = las_l[1:5],
                          # Not necessary (default):
                          cbh_ONLY = params$parameters$cbh_ONLY,
                          # Not necessary (default):
                          kM = params$parameters$kM)
```
User R Console interaction:

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_O_table.png">

User R Plots with R Console interaction:

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_001.png">
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_002.png">
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_003.png">
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_004.png">
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_005.png">

Let's check O_CBH.
```r
O_CBH
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_O_tableO.png">

Two other modes of **get_CBH()**.
```r
# 2) Optimization deactivated (cbh_ONLY = 2): executing treeiso only (PRE-PROCESSING step!) 
# SENSIBLE above 20 points/m², it skips input points clouds with smaller than 20 points/point cloud (4-7 points/m²)

outdi_treeiso <- "<path/to/dir_treeiso>"    # treeiso isolated tree cloud
outdi_filtered <- "<path/to/dir_filtered>"  # filtered tree cloud (stem plus first leaved branch)
cc_dir <- "<path/to/CloudCompare.exe>"      

treecbh::get_CBH(list_LAS_char = las_l,
                 outdir1 = outdi_treeiso,
                 outdir2 = outdi_filtered,
                 # Run only tree isolation:
                 cbh_ONLY = 2,
                 cc_dir = cc_dir)

# 3) Optimization deactivated (cbh_ONLY = 1): performing treeiso plus CBH detection
# EXECUTE ONLY if params$density_stats$mean_density > 20
D_CBH <- treecbh::get_CBH(list_LAS_char = las_l,
                          outdir1 = outdi_treeiso,
                          outdir2 = outdi_filtered,
                          # Run both tree isolation and CBH detection:
                          cbh_ONLY = 1,
                          cc_dir = cc_dir)
```
## Common issues

  - "No ground points found" error

    Solution: Ensure your point cloud includes ground points or set normalize = FALSE in get_3DTREE()

  - CloudCompare not found

    Solution: Verify the path to CloudCompare.exe in the cc_dir parameter

  - Empty output directories

    Solution: Check that input LAS files have sufficient points (>20 points per tree)

  - Interactive mode not working

    Solution: Ensure you're running R in interactive mode (RStudio recommended)

  - "Cannot write a file with 0 point" error

    Solution: The function automatically skips files with insufficient points. Check input file  quality.
