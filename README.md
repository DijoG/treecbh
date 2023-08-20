
The treecbh package provides functions to detect individual tree level Canopy Base Height (CBH) using high-resolution LiDAR data. Individual tree segmentation must be conducted prior. The package is meant to be used within the framework of the lidR package. 

### Required packages

tidyverse, lidR, sf, data.table, crayon, factoextra, fpc, geometry

### Installation

```r
devtools::install_github("DijoG/treecbh")
```
### Example

```r
#>>>>>> Data preparation

# Forest point cloud
Alas <- lidR::readLAS("<file.las>") 
  
# Individual tree segments (its)
Apoly <- sf::st_read("<file.shp>") 

# Extracting point clouds to its, saving las files
oudir <- "<path to directory>"
treecbh::get_3DTREE(Alas, Apoly, output_dir = oudir, FEATURE= "Species")

#>>>>>> CBH detection
#>>> 1> Optimization deactivated (treeiso plus cbh detection)

?treecbh::get_CBH()

its_l <- list.files(oudir, pattern = ".las", full.names = T)
outdi1 <- "<path to directory>"
outdi2 <- "<path to directory>"
cc_dir <- "<path to /CloudCompare.exe>"

A_CBH <- treecbh::get_CBH(its_l,
                          # run tree isolation and cbh detection
                          cbh_ONLY = 1,
                          # its point cloud directory
                          outdir1 = outdi1,
                          # point cloud segments (stem plus first leaved branch) directory 
                          outdir2 = outdi2,
                          cc_dir = cc_dir)
                 
#>>> 2> Optimization deactivated (treeiso only)

treecbh::get_CBH(its_l,
                 outdir1 = outdi1,
                 outdir2 = outdi2,
                 # run only tree isolation (disabling cbh detection)
                 cbh_ONLY = 2,
                 cc_dir = cc_dir) 
                           
#>>> 3> Optimization activated (cbh only)
# Run after executing #>>> 1> or >>> 2>, ideally after >>> 2>

O_CBH <- treecbh::get_CBH(its_l,
                          outdir1 = outdi1,
                          outdir2 = outdi2,
                          # run only cbh detection (disabling treeiso)
                          cbh_ONLY = 3,
                          # Activator:
                          kM = T,  
                          cc_dir = cc_dir)
```


