# treecbh 
==========

The **treecbh** package provides functions to detect individual tree level Canopy Base Height (CBH) using high-resolution LiDAR data. Individual tree segmentation must be conducted prior. The package is meant to be used within the framework of the **lidR** package. 

==========

### Required packages

Make sure that the following R packages are installed before using **treecbh**:

tidyverse,
lidR,
sf,
data.table,
crayon,
factoextra,
fpc,
geometry,
gtools.

### Installation

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

### Computing canopy height model
Using the pitfree algorithm.

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

### Point-based crown segmentation
Using the Dalponte algorithm.

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

### Extracting point clouds to its (classifying ground points and saving las files by default)

This step is only necessary if you have different its shape file coming from an outer source.

```r
oudir <- "<path to directory>"
treecbh::get_3DTREE(Alas, Apoly, output_dir = oudir, FEATURE= "ID", RETURN = F)
```

### CBH detection

Optimization deactivated, performing treeiso plus cbh detection (`cbh_ONLY = 1`). Let's proceed with five trees that have sufficient points.

<img align="right" src="https://raw.githubusercontent.com/DijoG/storage/main/README/3D_tree_example.png">

```r
?treecbh::get_CBH()

# list of individual las tree clouds and the selection of five
its_l <- list.files(oudit, pattern = ".las", full.names = T) %>%
  gtools::mixedsort()
las_l <- its_l[c(7,9,72,78,131)]

# visualizing one
plot(readLAS(las_l[5]), bg = "white", size = 5, axis = T)

# output directories
outdi1 <- "<path to directory>" # treeiso isolated tree cloud
outdi2 <- "<path to directory>" # treeiso segmented clouds of tree parts

# CC executable 
cc_dir <- "<path to /CloudCompare.exe>"

# running 
A_CBH <- treecbh::get_CBH(las_l,
                          # run tree isolation and cbh detection
                          cbh_ONLY = 1,
                          # its point cloud directory
                          outdir1 = outdi1,
                          # point cloud segments (stem plus first leaved branch) directory
                          outdir2 = outdi2,
                          # CC executeable directory
                          cc_dir = cc_dir)
A_CBH
```

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/04_A_CBH.png">

Optimization deactivated, executing **treeiso** only (`cbh_ONLY = 2`).

```r
treecbh::get_CBH(its_l,
                 outdir1 = outdi1,
                 outdir2 = outdi2,
                 # run only tree isolation (disabling cbh detection)
                 cbh_ONLY = 2,
                 cc_dir = cc_dir)
``` 

Optimization activated with `kM = TRUE` (k-Means Clustering set to TRUE), executing cbh only (`cbh_ONLY = 3`). The interactive optimization process starts. 'Do you accept k?'. User answers, if the answer is no, the user is asked to type a number for `k` after 'Enter k:'.

<img align="right" src="https://raw.githubusercontent.com/DijoG/storage/main/README/sugg_k_no.png">

```r
O_CBH <- treecbh::get_CBH(its_l,
                          outdir1 = outdi1,
                          outdir2 = outdi2,
                          # run only cbh detection (disabling treeiso)
                          cbh_ONLY = 3,
                          # Optimization activator
                          # kM = FALSE is default
                          kM = TRUE,
                          # inactive CC executeable directory
                          cc_dir = cc_dir)
O_CBH
```

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/05_O_CBH.png">

