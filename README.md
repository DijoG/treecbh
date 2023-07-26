
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
Alas <- lidR::readLAS("<file.las>") %>%
  lidR::normalize_height(., tin())
  
# Individual tree segments
Apoly <- sf::st_read("<file.shp>") 

# Extracting point clouds to individual tree segments
Aits <- treecbh::get_3DTREE(Alas, Apoly, "Species")

# If needed, point cloud files can be saved
for (i in 1:length(Aits)) {
  setwd("<path to directory>")
  lidR::writeLAS(Aits[[i]], str_c("A_0", i, ".las"))
}

#>>>>>> CBH detection
#>>> 1> Optimization deactivated (default)

?treecbh::get_CBH()

outdi1 <- "<path to directory>"
outdi2 <- "<path to directory>"
cc_dir <- "<path to /CloudCompare.exe>"

A_CBH <- treecbh::get_CBH(Aits,
                          outdir1 = outdi1,
                          outdir2 = outdi2,
                          cc_dir = cc_dir)
                 
#>>> 2> Optimization activated
# Run after executing #>>> 1>

A_OCBH <- treecbh::get_CBH(Aits,
                           outdir1 = outdi1,
                           outdir2 = outdi2,
                           # Activator:
                           kM = T,  
                           # Disabling treeiso:
                           ONLY = T,
                           cc_dir = cc_dir)
```


