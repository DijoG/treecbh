# treecbh 
==========

The **treecbh** package provides functions to detect individual tree level Crown Base Height (CBH) using high-resolution LiDAR data. 

Individual tree segmentation must be conducted prior. The package is meant to be used within the framework of the **lidR** package. 

***Update (29.10.2025):*** Cumulative percentile-based density method implemented instead of the original Kernel density approach.

==========

## Original paper
*Testing treecbh in Central European forests: an R package for crown base height detection using high-resolution aerial laser-scanned data*

https://doi.org/10.1093/forestry/cpae044

## Required packages

tidyverse, lidR, RCSF, sf, data.table, crayon, fpc, geometry, gtools

## Installation

```r
devtools::install_github("DijoG/treecbh")
```

# Example
Demonstration of how to use **treecbh** combined with **lidR**. 

### Data preparation

```r
require(lidR);require(tidyverse)

# Forest point cloud (low resolution example data)
LAS <- system.file("extdata", "MixedConifer.laz", package = "lidR")
Alas <- readLAS(LAS, filter = "-drop_z_below 0") 

# Black, white and green color palette for visualizing CHM
bgcol <- function(x)
{
  col <- grDevices::colorRampPalette(c("grey1", "white", "forestgreen"))
  return(col(x))
}
```
Computing canopy height model using the pitfree algorithm.

```r
CHM <- rasterize_canopy(Alas, 0.5, pitfree(subcircle = 0.25))
plot(CHM, main = "CHM 0.5 pitfree", col = bgcol(50))
```

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/01_chm_pitfree.png">

### Computing treetops
Using a constant windows size of 5 m.

```r
ws <- 5
treetops <- locate_trees(CHM, lmf(ws))
plot(CHM, main = "CHM 0.5 pitfree", col = bgcol(50))
plot(sf::st_geometry(treetops), add = T, pch = "+", col = "firebrick3")
```

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/02_chm_pitfree_treetops.png">

Point-based crown segmentation using the Dalponte algorithm.

```r
algo_crowns <- dalponte2016(CHM, treetops)
las_crowns <- segment_trees(Alas, algo_crowns, attribute = "ID")
```

`las_crowns` is a las object storing the ids of individual trees (ID attribute). This object can also be fed into the `treecbh::get_CBH()` function after storing each individual point cloud (tree) in a directory. Since there are other point- and CHM-based its algorithms outside the liDR framework, `trecbh::get_CBH()` accepts las files that can previously be isolated based both on points and CHMs. 

Obtaining individual tree segments (its) as polygons.

```r
Apoly <- crown_metrics(las_crowns, attribute = "ID", geom = "concave", func = NULL)
plot(sf::st_geometry(Apoly), reset = FALSE, col = "forestgreen", border = "grey80")
plot(sf::st_geometry(treetops), add = T, pch = "+", col = "firebrick3")
```

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/03_its_treetops.png">

Some trees were not segmented.

Extracting point clouds to its (classifying ground points and saving las files by default)
This step is only necessary if you have different its shape file coming from an outer source.

```r
oudir <- "<path to directory>"
treecbh::get_3DTREE(Alas, Apoly, output_dir = oudir, FEATURE= "ID", RETURN = F)
```

### CBH detection 

```r
# List las files
las_l <- list.files(oudir, pattern = ".las", full.names = T) %>%
  gtools::mixedsort()
  
# Get automated parameters
tictoc::tic()
params <- get_PARAMS(oudir)
tictoc::toc()
# 25.84 sec

# Run
# 1) Optimization activated, executing cbh detection only (default)
O_CBH <- get_CBH(list_LAS_char = las_l[1:5],
                 # Not necessary (default):
                 cbh_ONLY = params$parameters$cbh_ONLY,
                 # Not necessary (default):
                 kM = params$parameters$kM)
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_O_seq.png">

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_O_02.png">
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_O_03.png">
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_O_04.png">
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_O_05.png">

Let's check 'O_CBH':

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/treecbh_O_table.png">

```r
# 2) Optimization deactivated, executing treeiso only (PRE-PROCESSING step!) 
# Sensible above 20 points/m², it skips las trees with smaller than 20 points/point cloud (4-7 points/m²)

outdi_treeiso <- "D:/test/treecbh/treeiso"    # treeiso isolated tree cloud
outdi_filtered <- "D:/test/treecbh/filtered"  # filtered tree cloud (stem plus first leaved branch)
cc_dir <- "D:/CloudCompare/CloudCompare.exe"

get_CBH(list_LAS_char = las_l,
        outdir1 = outdi_treeiso,
        outdir2 = outdi_filtered,
        # Run only tree isolation:
        cbh_ONLY = 2,
        cc_dir = cc_dir)

# 3) Optimization deactivated, performing treeiso plus CBH detection
# EXECUTE ONLY if params$density_stats$mean_density > 20
D_CBH <- get_CBH(list_LAS_char = las_l,
                 outdir1 = outdi_treeiso,
                 outdir2 = outdi_filtered,
                 # Run both tree isolation and CBH detection
                 cbh_ONLY = 1,
                 cc_dir = cc_dir)
```



