########################################################################################## 
#  
# Smoothing of insecticide areas across the Prairies
#
# contact: eginamalaj@gmail.com
#
# GAM models and plotting adaptated from: 
# https://fromthebottomoftheheap.net/2017/10/19/first-steps-with-mrf-smooths/
########################################################################################## 
#
#
require(raster) 
require(rgdal)
require(proj4)
require(spdep)
require(mgcv)
require(tidyverse)
require(viridis)
require(gridExtra)
#
#
########################################################################################## 
# Data
########################################################################################## 
#
# 1. Raster file with insecticide use
ras<-raster("rat_ins.tif")
#
# 2. Polygon file with river catchments
pol<-readOGR('.', 'catchment_red')
#
# 3. Load results of GAM models from the .RData below or run them yourself (part 1)
# Depending on your computer, they take can take a long time to run
#
load(file="m1.RData")
load(file="m2.RData")
load(file="m3.RData")
load(file="m4.RData")
#
# 4. Themes for plotting in ggplot
#
theme_map <- function(...) {
  theme_minimal() +
    theme(...,
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank())
}
#
myTheme <- theme_map(legend.position = 'bottom')
myScale <- scale_fill_viridis(name = '[kg/ha]', option = 'plasma',
                              limits = c(0.001, 2),
                              guide = guide_colorbar(direction = "horizontal",
                                                     barheight = unit(2, units = "mm"),
                                                     barwidth = unit(75, units = "mm"),
                                                     title.position = 'left',
                                                     title.hjust = 0.5,
                                                     label.hjust = 0.5))
#
#
########################################################################################## 
#
# Polygon Prep
#
# Polygons re-projections
#
aea.proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96"
pol <- spTransform(pol, CRS(aea.proj))
#
x11(); plot(pol)
dfs <- droplevels(as(pol, 'data.frame'))
#
# Fortify
pol_fort <- fortify(pol, region = 'OBJECTID_1')
#
# The same amount of regions as polygons. This creates a neighbourhood structure: 
# find which regions (polygons) are neighbours by checking if they share a border.
#
nb<- poly2nb(pol, row.names = dfs$OBJECTID_1)
names(nb) <- attr(nb, "region.id")
#
#
########################################################################################## 
# GAM modeling
########################################################################################## 
#
# Note the running times below
# Use provided files for first visualizations
#
ctrl <- gam.control(nthreads = 6)
#
# relate PestApl(from MEAN) with each region - Full rank takes 6h to run
system.time(m1 <- gam(MEAN~ s(OBJECTID_1, bs = 'mrf', xt = list(nb = nb)), # define MRF smooth
          data = dfs,
          method = 'REML', # fast version of REML smoothness selection
          control = ctrl,
          family = gaussian)) # linear model
#summary(m1)
#
#save(m1, file="m1.RData") #- file is 700MB
#load(file="m1.RData")
#
#
# Smooth it some more - 6 minutes to run
#
system.time(m2 <- gam(MEAN ~ s(OBJECTID_1, bs = 'mrf', k = 30, xt = list(nb = nb)),
          data = dfs,
          method = 'REML',
          control = ctrl,
          family = gaussian()))
summary(m2)
#
#save(m2, file="m2.RData")
#
m2_bas<-gam(MEAN+MEAN_1 ~ s(OBJECTID_1, bs = 'mrf', k = 30, xt = list(nb = nb)),
            data = dfs,
            method = 'REML',
            control = ctrl,
            family = gaussian())
#
#
# Smooth it again 
system.time(m3 <- gam(MEAN ~ s(OBJECTID_1, bs = 'mrf', k = 300, xt = list(nb = nb)),
                      data = dfs,
                      method = 'REML', 
                      control = ctrl,
                      family = gaussian()))
summary(m3)
#save(m3, file="m3.RData")
#
# Keep smoothing 
system.time(m4 <- gam(MEAN ~ s(OBJECTID_1, bs = 'mrf', k = 1000, xt = list(nb = nb)),
                      data = dfs,
                      method = 'REML', 
                      control = ctrl,
                      family = gaussian()))
summary(m4)
#save(m4, file="m4.RData")
#
#
########################################################################################## 
##### Plotting
########################################################################################## 
#
# 
dfs <- transform(dfs,
                mrfFull= predict(m1, type = 'response'),
                mrfRrank30= predict(m2, type = 'response'),
                mrfRrank300= predict(m3, type = 'response'),
                mrfRrank1000= predict(m4, type = 'response'))
#
as.factor(pol_fort$id)->pol_fort$id
poldt <- left_join(pol_fort, dfs, by = c('id' = 'OBJECTID_1'))
#
#
# Raster data - Mean Pesticide application
# For plotting in ggplot RasterLayer needs to be transformed into a dataframe
#
# Raw - Mean Pesticide application
#
praw<- ggplot(poldt, aes(x = long, y = lat, group = id)) +
  geom_polygon(aes(fill = MEAN)) +
  geom_path(col = NA, alpha = 0.5, size = 0.1) +
  coord_equal() +
  labs(x = NULL, y = NULL, title = 'Mean Insecticide Use per catchment')+
  myTheme + myScale
#
# Full MRF ranks
pfull<-ggplot(poldt, aes(x = long, y = lat, group = id)) +
  geom_polygon(aes(fill = mrfFull)) +
  geom_path(col = NA, alpha = 0.5, size = 0.1) +
  coord_equal() +
  labs(x = NULL, y = NULL, title = 'GAM model with full rank MRF')+
  myTheme + myScale
#
# 30 MRF ranks
pmrf30<-ggplot(poldt, aes(x = long, y = lat, group = id)) +
  geom_polygon(aes(fill = mrfRrank30)) +
  geom_path(col = NA, alpha = 0.5, size = 0.1) +
  coord_equal() +
  labs(x = NULL, y = NULL, title = 'GAM model with 30 MRF ranks')+
  myTheme + myScale
#
# 300 MRF ranks
pmrf300<-ggplot(poldt, aes(x = long, y = lat, group = id)) +
  geom_polygon(aes(fill = mrfRrank300)) +
  geom_path(col = NA, alpha = 0.5, size = 0.1) +
  coord_equal() +
  labs(x = NULL, y = NULL, title = 'GAM model with 300 MRF ranks')+
  myTheme + myScale
#
# 1000 MRF ranks
pmrf1000<-ggplot(poldt, aes(x = long, y = lat, group = id)) +
  geom_polygon(aes(fill = mrfRrank1000)) +
  geom_path(col = NA, alpha = 0.5, size = 0.1) +
  coord_equal() +
  labs(x = NULL, y = NULL, title = 'GAM model with 1000 MRF ranks')+
  myTheme + myScale
#
#
x11(width =10, height=4)
raw = grid.arrange(praw, pfull, nrow = 1)
ggsave(filename = "Mean_GAM.png", dpi=500, width =10, height=4, plot = raw)
#
#
x11(width =13, height=4)
pmrf <- grid.arrange(pmrf30, pmrf300,pmrf1000, nrow = 1)
ggsave(filename = "PMRF_range.png", dpi=500, width =13, height=4, plot = pmrf)
#
#
# References: 
# https://www.fromthebottomoftheheap.net/2017/10/19/first-steps-with-mrf-smooths/
# https://pudding.cool/process/regional_smoothing/
#
#
  
  