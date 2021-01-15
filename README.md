# Regional Smoothing of Insecticide Use Data


## Depends

To run this analysis two files are required:  
(i) `rat_ins.tif` raster file produced in a similar way to the method reported in the repository **[R_SpatialAnalysis_Mapping](https://github.com/eginamalaj/R_SpatialAnalysis_Mapping)**. This is a raster file of insecticide use across the Prairie Pothole Region of Canada;
(ii) Folder `catchment` with a vector file containing deliniations of 5000+ river basins. 

R version 4.0.3

Packages: mgcv, tidyverse, rgdal, raster, proj4, spdep, mgcv, viridis, gridExtra

Code and analysis was largely based on **[this blog](https://fromthebottomoftheheap.net/2017/10/19/first-steps-with-mrf-smooths/)** created by Gavin Simpson.

## Description

Here by smoothing the insecticide distribution, it would be possible to identify the “hotspot” areas with regard to insecticide use across the Prairies.

Initially, I produced the raster `rat_ins.tif`, which represents insecticide use (kg/ha) across the Prairie Pothole Region of Canada at 1-km resolution. This raster file was overlaid to the river basin vector file `catchment.shp`. I calculated mean insecticide distribution for each river basin. I fitted generalized additive models (GAM) to the data found in the packages `mgcv`. The level of smoothness is defined by the Markov Random Field (MRF), which allows modeling of spatial data with an intrinsic  Gaussian Markov random field (GMRF). As explained by **[G.Simpson](https://fromthebottomoftheheap.net/2017/10/19/first-steps-with-mrf-smooths/)**: 'MRFs are quite flexible as you can think about them as representing an undirected graph whose nodes are your samples and the connections between the nodes are specified via a neighbourhood structure'.


Below I'll show two sets of figures produced by this analysis. Figure 1 shows the distribution when mean insecticide use is used (basically a choropleth map) and and with full rank MRF, which is pretty close to the mean map.


Figure 2 shows three sets of GAM models with 30-1000 MRF ranks. You can visually see the smoothing happening more in 30 ranks than in 1000 MRF. In comparison to Figure 1, in Figure 2 it is easier to visually identify areas which are different from neighboring locations, and the regional hotspots are more evident as a result of smoothing. 


## References: 
https://www.fromthebottomoftheheap.net/2017/10/19/first-steps-with-mrf-smooths/

https://pudding.cool/process/regional_smoothing/