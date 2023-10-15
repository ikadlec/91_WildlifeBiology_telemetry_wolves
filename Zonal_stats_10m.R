require(ctmm)
require(raster)
require(rgdal)
library(stars)
library(tidyverse)

# Loading the necessary data
load("UDS_10m.RData")
load("UDS_month_wolves_10m.RData")

rasfn <- "ESA_clip_proj.tif"
ESA <- raster(rasfn)
name <- read.table("rast_table_ESA.csv", header = T, sep = ";") # table containing legend to ESA landcover

habitat <- data.frame(Value = seq(10, 110, by = 10),
                   Habitat = c("Tree cover", "Shrubland", "Grassland", "Cropland",
                          "Built-up", "Bare/sparse vegetation", 
                          "Snow and Ice", "Pernament water bodies",
                          "Herbaceous wetland", "Mangroves", "Moss and lichen"))

# UDS - res 10m
UDS_36927_10m <- akde(tel_36927, FITS_36927, level = 0.95, grid = 10) # UD estimate
UDS_36928_10m <- akde(tel_36928, FITS_36928, level = 0.95, grid = 10) # UD estimate
UDS_36930_10m <- akde(tel_36930, FITS_36930, level = 0.95, grid = 10) # UD estimate
UDS_36931_10m <- akde(tel_36931, FITS_36931, level = 0.95, grid = 10) # UD estimate
UDS_36932_10m <- akde(tel_36932, FITS_36932, level = 0.95, grid = 10) # UD estimate
UDS_29456_10m <- akde(tel_29456, FITS_29456, level = 0.95, grid = 10) # UD estimate
UDS_47016_10m <- akde(tel_47016, FITS_47016, level = 0.95, grid = 10) # UD estimate
UDS_30776_10m <- akde(tel_30776, FITS_30776, level = 0.95, grid = 10) # UD estimate
UDS_30777_10m <- akde(tel_30777, FITS_30777, level = 0.95, grid = 10) # UD estimate
UDS_30779_10m <- akde(tel_30779, FITS_30779, level = 0.95, grid = 10) # UD estimate
UDS_46090_10m <- akde(tel_46090, FITS_46090, level = 0.95, grid = 10) # UD estimate
UDS_46088_10m <- akde(tel_46088, FITS_46088, level = 0.95, grid = 10) # UD estimate
UDS_29456_AUT_10m <- akde(tel_29456_AUT, FITS_29456_AUT, level = 0.95, grid = 10) # UD estimate


UDS_10m <- list(UDS_36927_10m, UDS_36928_10m ,UDS_36930_10m, UDS_36931_10m, UDS_36932_10m,
                UDS_29456_10m, UDS_47016_10m, 
                UDS_30776_10m, UDS_30777_10m, UDS_30779_10m, UDS_46088_10m,
                UDS_46090_10m, UDS_29456_AUT_10m)

names(UDS_10m) <- c("36927", "36928", "36930", "36931", "36932", "29456", "47016",
                "30776", "30777", "30779", "46088", "46090", "29456_AUT")

save(UDS_10, file = "UDS_10m_wolves.RData")


# The analysis of wolf use of different landscape --------

# export UDs to raster
rs_UDS_10m <- lapply(UDS_10m, function(b) ctmm::raster(b, level.UD = 0.95, level = 0.95, DF = "PDF") %>% 
                       st_as_stars %>% 
                       st_warp(crs = st_crs(st_as_stars(ESA)), cellsize = res(ESA)) %>% 
                       as("Raster") %>% 
                       "/"(sum(cellStats(.,"sum"))))

# ESA landocver cropping
ESA_crop <- lapply(1:13, function(i) crop(ESA, rs_UDS_10m[[i]]) %>% 
                       st_as_stars %>% 
                       st_warp(st_as_stars(rs_UDS_10m[[i]])) %>% 
                       as("Raster")) 

# zonal statistics
zonal_stats_ESA <- lapply(1:13, function(i) 
  zonal(rs_UDS_10m[[i]], ESA_crop[[i]], fun = "sum") %>%
    as.data.frame %>% 
    left_join(zonal(rs_UDS_10m[[i]], ESA_crop[[i]], fun = "mean") %>% 
                as.data.frame) %>% 
    mutate(count = sum/mean,
           percent = sum/sum(sum),
           uniform = count*sum(sum)/sum(count, na.rm = T),
           relative = sum/uniform,
           ID = names(UDS_10m)[i])) %>% 
  bind_rows() %>% # joining into one data frame
  mutate(habitat = recode(zone,
                          "10" = "Tree cover",
                          "20" = "Shrubland",
                          "30" = "Grassland",
                          "40" = "Cropland",
                          "50" = "Built-up", 
                          "60" = "Bare/sparse vegetation",
                          "70" = "Snow and Ice",
                          "80" = "Pernament water bodies",
                          "90" = "Herbaceous wetland",
                          "100" = "Mangroves",
                          "110" = "Moss and lichen",
  ))


save(zonal_stats_ESA, file = "zonal_10m.RData")

ggplot(zonal_stats_ESA, 
       aes(y=factor(ID, 
                    levels = c("47015", "36927", "36930", "36931", "47016", 
                               "29456",
                               "36928", "36932", "30776", "30777", "30779", "46088", "46090", "29456_AUT")), 
           x = sum, fill = habitat)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(legend.position = "bottom", legend.key.width = unit(0.5, 'cm')) 




# Analysis of the use of different landscape areas by wolves during the year - month by month --------


# Zonal stats_ESA_10m by month --------------------------------------------
UDS_ESA_month <- list(UDS_1, UDS_2, UDS_3, UDS_4, UDS_5, UDS_6, UDS_8,
                  UDS_9, UDS_10, UDS_12, UDS_13, UDS_14, UDS_15)
names(UDS_ESA_month) <- c("36927", "36928", "36930", "36931", "36932", "29456", "47016",
                      "30776", "30777", "30779", "46090", "46088", "29456_AUT")

rs_UDS_ESA_month <- vector("list",13)

UDS_ESA_month <- UDS_month



for (i in 1:13) {
  rs_UDS_ESA_month[[i]] <- lapply(UDS_ESA_month[[i]], function(b) ctmm::raster(b, level.UD = 0.95, level = 0.95, DF = "PDF") %>% 
                                st_as_stars %>% 
                                st_warp(crs = st_crs(st_as_stars(ESA))) %>% 
                                as("Raster") %>% 
                                "/"(sum(cellStats(.,"sum"))))
  
}



ESA_36927 <- lapply(rs_UDS_ESA_month[[1]], function(i) crop(ESA, i) %>%
                        st_as_stars %>%
                        st_warp(st_as_stars(i)) %>%
                        as("Raster")) 
ESA_36928 <- lapply(rs_UDS_ESA_month[[2]], function(i) crop(ESA, i) %>%
                        st_as_stars %>%
                        st_warp(st_as_stars(i)) %>%
                        as("Raster"))
ESA_36930 <- lapply(rs_UDS_ESA_month[[3]], function(i) crop(ESA, i) %>%
                        st_as_stars %>%
                        st_warp(st_as_stars(i)) %>%
                        as("Raster"))
ESA_36931 <- lapply(rs_UDS_ESA_month[[4]], function(i) crop(ESA, i) %>%
                        st_as_stars %>%
                        st_warp(st_as_stars(i)) %>%
                        as("Raster"))
ESA_36932 <- lapply(rs_UDS_ESA_month[[5]], function(i) crop(ESA, i) %>%
                        st_as_stars %>%
                        st_warp(st_as_stars(i)) %>%
                        as("Raster"))
ESA_29456 <- lapply(rs_UDS_ESA_month[[6]], function(i) crop(ESA, i) %>%
                        st_as_stars %>%
                        st_warp(st_as_stars(i)) %>%
                        as("Raster"))
ESA_47016 <- lapply(rs_UDS_ESA_month[[7]], function(i) crop(ESA, i) %>%
                        st_as_stars %>%
                        st_warp(st_as_stars(i)) %>%
                        as("Raster"))

ESA_30776 <- lapply(rs_UDS_ESA_month[[8]], function(i) crop(ESA, i) %>%
                        st_as_stars %>%
                        st_warp(st_as_stars(i)) %>%
                        as("Raster"))
ESA_30777 <- lapply(rs_UDS_ESA_month[[9]], function(i) crop(ESA, i) %>%
                        st_as_stars %>%
                        st_warp(st_as_stars(i)) %>%
                        as("Raster"))
ESA_30779 <- lapply(rs_UDS_ESA_month[[10]], function(i) crop(ESA, i) %>%
                        st_as_stars %>%
                        st_warp(st_as_stars(i)) %>%
                        as("Raster"))
ESA_46090 <- lapply(rs_UDS_ESA_month[[11]], function(i) crop(ESA, i) %>%
                        st_as_stars %>%
                        st_warp(st_as_stars(i)) %>%
                        as("Raster"))
ESA_46088 <- lapply(rs_UDS_ESA_month[[12]], function(i) crop(ESA, i) %>%
                        st_as_stars %>%
                        st_warp(st_as_stars(i)) %>%
                        as("Raster"))
ESA_29456_AUT <- lapply(rs_UDS_ESA_month[[13]], function(i) crop(ESA, i) %>%
                            st_as_stars %>%
                            st_warp(st_as_stars(i)) %>%
                            as("Raster"))


ESA_ESA_month <- list(ESA_36927, ESA_36928, ESA_36930, ESA_36931, ESA_36932,
                    ESA_29456, ESA_47016,
                    ESA_30776, ESA_30777, ESA_30779, ESA_46090,
                    ESA_46088, ESA_29456_AUT)

names(ESA_ESA_month) <- c("36927", "36928", "36930", "36931", "36932", "29456", "47016",
                        "30776", "30777", "30779", "46090", "46088", "29456_AUT")


# zonal statistics
## 36927
zonal_month_36927 <- lapply(1:11, function (b) zonal(rs_UDS_ESA_month[[1]][[b]], ESA_36927[[b]], fun = "sum") %>%
                        as.data.frame %>% 
                        left_join(zonal(rs_UDS_ESA_month[[1]][[b]], ESA_36927[[b]], fun = "mean") %>% 
                                    as.data.frame) %>% 
                        mutate(count = sum/mean,
                               percent = sum/sum(sum),
                               uniform = count*sum(sum)/sum(count, na.rm = T),
                               relative = sum/uniform,
                               ID = "36927",
                               month = names(ESA_36927)[b])) %>% 
  bind_rows() #  joining into one data frame

# adding missing months in the year (filling NA with values)
Oct <- c(10, NA, NA, NA, NA, NA, NA, "36927", "Oct") %>%
  rbind(c(30, NA, NA, NA, NA, NA, NA, "36927", "Oct"),
        c(20, NA, NA, NA, NA, NA, NA, "36927", "Oct"),
        c(40, NA, NA, NA, NA, NA, NA, "36927", "Oct"),
        c(50, NA, NA, NA, NA, NA, NA, "36927", "Oct"),
        c(60, NA, NA, NA, NA, NA, NA, "36927", "Oct"),
        c(80, NA, NA, NA, NA, NA, NA, "36927", "Oct"),
        c(100, NA, NA, NA, NA, NA, NA, "36927", "Oct"),
        c(90, NA, NA, NA, NA, NA, NA, "36927", "Oct")) %>%
  as.data.frame()
names(Oct) <- names(zonal_month_36927)
Oct$zone <- as.numeric(Oct$zone)
Oct$sum <- as.numeric(Oct$sum)
Oct$mean <- as.numeric(Oct$mean)
Oct$count <- as.numeric(Oct$count)
Oct$percent <- as.numeric(Oct$percent)
Oct$uniform <- as.numeric(Oct$uniform)
Oct$relative <- as.numeric(Oct$relative)

NA_36927 <- Oct 
zonal_month_36927 <- rbind(zonal_month_36927, NA_36927)

levels(factor(NA_36927$zone))
levels(factor(zonal_month_36927$zone))

## 36928
zonal_month_36928 <- lapply(1:8, function (b) zonal(rs_UDS_ESA_month[[2]][[b]], ESA_36928[[b]], fun = "sum") %>%
                        as.data.frame %>% 
                        left_join(zonal(rs_UDS_ESA_month[[2]][[b]], ESA_36928[[b]], fun = "mean") %>% 
                                    as.data.frame) %>% 
                        mutate(count = sum/mean,
                               percent = sum/sum(sum),
                               uniform = count*sum(sum)/sum(count, na.rm = T),
                               relative = sum/uniform,
                               ID = "36928",
                               month = names(ESA_36928)[b])) %>% 
  bind_rows() 

## 36930
zonal_month_36930 <- lapply(1:8, function (b) zonal(rs_UDS_ESA_month[[3]][[b]], ESA_36930[[b]], fun = "sum") %>%
                        as.data.frame %>% 
                        left_join(zonal(rs_UDS_ESA_month[[3]][[b]], ESA_36930[[b]], fun = "mean") %>% 
                                    as.data.frame) %>% 
                        mutate(count = sum/mean,
                               percent = sum/sum(sum),
                               uniform = count*sum(sum)/sum(count, na.rm = T),
                               relative = sum/uniform,
                               ID = "36930",
                               month = names(ESA_36930)[b])) %>% 
  bind_rows() 

# adding missing months in the year (filling NA with values)
Aug <- c(10, NA, NA, NA, NA, NA, NA, "36930", "Aug") %>%
  rbind(c(20, NA, NA, NA, NA, NA, NA, "36930", "Aug"),
        c(30, NA, NA, NA, NA, NA, NA, "36930", "Aug"),
        c(40, NA, NA, NA, NA, NA, NA, "36930", "Aug"),
        c(50, NA, NA, NA, NA, NA, NA, "36930", "Aug"),
        c(60, NA, NA, NA, NA, NA, NA, "36930", "Aug"),
        c(80, NA, NA, NA, NA, NA, NA, "36930", "Aug"),
        c(90, NA, NA, NA, NA, NA, NA, "36930", "Aug")) %>%
  as.data.frame()
names(Aug) <- names(zonal_month_36930)

Sep <- c(10, NA, NA, NA, NA, NA, NA, "36930", "Sep") %>% 
  rbind(c(20, NA, NA, NA, NA, NA, NA, "36930", "Sep"),
        c(30, NA, NA, NA, NA, NA, NA, "36930", "Sep"),
        c(40, NA, NA, NA, NA, NA, NA, "36930", "Sep"),
        c(50, NA, NA, NA, NA, NA, NA, "36930", "Sep"),
        c(60, NA, NA, NA, NA, NA, NA, "36930", "Sep"),
        c(80, NA, NA, NA, NA, NA, NA, "36930", "Sep"),
        c(90, NA, NA, NA, NA, NA, NA, "36930", "Sep")) %>%
  as.data.frame()
names(Sep) <- names(zonal_month_36930)

Oct <- c(10, NA, NA, NA, NA, NA, NA, "36930", "Oct") %>%
  rbind(c(20, NA, NA, NA, NA, NA, NA, "36930", "Oct"),
        c(30, NA, NA, NA, NA, NA, NA, "36930", "Oct"),
        c(40, NA, NA, NA, NA, NA, NA, "36930", "Oct"),
        c(40, NA, NA, NA, NA, NA, NA, "36930", "Oct"),
        c(50, NA, NA, NA, NA, NA, NA, "36930", "Oct"),
        c(60, NA, NA, NA, NA, NA, NA, "36930", "Oct"),
        c(80, NA, NA, NA, NA, NA, NA, "36930", "Oct"),
        c(90, NA, NA, NA, NA, NA, NA, "36930", "Oct")) %>%
  as.data.frame()
names(Oct) <- names(zonal_month_36930)

Nov <- c(10, NA, NA, NA, NA, NA, NA, "36930", "Nov") %>%
  rbind(c(20, NA, NA, NA, NA, NA, NA, "36930", "Nov"),
        c(30, NA, NA, NA, NA, NA, NA, "36930", "Nov"),
        c(40, NA, NA, NA, NA, NA, NA, "36930", "Nov"),
        c(40, NA, NA, NA, NA, NA, NA, "36930", "Nov"),
        c(50, NA, NA, NA, NA, NA, NA, "36930", "Nov"),
        c(60, NA, NA, NA, NA, NA, NA, "36930", "Nov"),
        c(70, NA, NA, NA, NA, NA, NA, "36930", "Nov"),
        c(80, NA, NA, NA, NA, NA, NA, "36930", "Nov"),
        c(90, NA, NA, NA, NA, NA, NA, "36930", "Nov"),
        c(100, NA, NA, NA, NA, NA, NA, "36930", "Nov")) %>%
  as.data.frame()
names(Nov) <- names(zonal_month_36930)


NA_36930 <- rbind(Aug, Sep, Oct, Nov) 
NA_36930$zone <- as.numeric(NA_36930$zone)
NA_36930$sum <- as.numeric(NA_36930$sum)
NA_36930$mean <- as.numeric(NA_36930$mean)
NA_36930$count <- as.numeric(NA_36930$count)
NA_36930$percent <- as.numeric(NA_36930$percent)
NA_36930$uniform <- as.numeric(NA_36930$uniform)
NA_36930$relative <- as.numeric(NA_36930$relative)

zonal_month_36930 <- rbind(zonal_month_36930, NA_36930)

levels(factor(NA_36930$zone))
levels(factor(zonal_month_36930$zone))

## 36931
zonal_month_36931 <- lapply(1:4, function (b) zonal(rs_UDS_ESA_month[[4]][[b]], ESA_36931[[b]], fun = "sum") %>%
                        as.data.frame %>% 
                        left_join(zonal(rs_UDS_ESA_month[[4]][[b]], ESA_36931[[b]], fun = "mean") %>% 
                                    as.data.frame) %>% 
                        mutate(count = sum/mean,
                               percent = sum/sum(sum),
                               uniform = count*sum(sum)/sum(count, na.rm = T),
                               relative = sum/uniform,
                               ID = "36931",
                               month = names(ESA_36931)[b])) %>% 
  bind_rows() 

# adding missing months in the year (filling NA with values)
Feb <- c(10, NA, NA, NA, NA, NA, NA, "36931", "Feb") %>%
  rbind(c(20, NA, NA, NA, NA, NA, NA, "36931", "Feb"),
        c(30, NA, NA, NA, NA, NA, NA, "36931", "Feb"),
        c(40, NA, NA, NA, NA, NA, NA, "36931", "Feb"),
        c(50, NA, NA, NA, NA, NA, NA, "36931", "Feb"),
        c(60, NA, NA, NA, NA, NA, NA, "36931", "Feb"),
        c(80, NA, NA, NA, NA, NA, NA, "36931", "Feb"),
        c(100, NA, NA, NA, NA, NA, NA, "36931", "Feb"),
        c(90, NA, NA, NA, NA, NA, NA, "36931", "Feb")) %>%
  as.data.frame()
names(Feb) <- names(zonal_month_36931)

Mar <- c(10, NA, NA, NA, NA, NA, NA, "36931", "Mar") %>%
  rbind(c(20, NA, NA, NA, NA, NA, NA, "36931", "Mar"),
        c(30, NA, NA, NA, NA, NA, NA, "36931", "Mar"),
        c(40, NA, NA, NA, NA, NA, NA, "36931", "Mar"),
        c(50, NA, NA, NA, NA, NA, NA, "36931", "Mar"),
        c(60, NA, NA, NA, NA, NA, NA, "36931", "Mar"),
        c(80, NA, NA, NA, NA, NA, NA, "36931", "Mar"),
        c(100, NA, NA, NA, NA, NA, NA, "36931", "Mar"),
        c(90, NA, NA, NA, NA, NA, NA, "36931", "Mar")) %>%
  as.data.frame()
names(Mar) <- names(zonal_month_36931)

Apr <- c(10, NA, NA, NA, NA, NA, NA, "36931", "Apr") %>%
  rbind(c(20, NA, NA, NA, NA, NA, NA, "36931", "Apr"),
        c(30, NA, NA, NA, NA, NA, NA, "36931", "Apr"),
        c(40, NA, NA, NA, NA, NA, NA, "36931", "Apr"),
        c(50, NA, NA, NA, NA, NA, NA, "36931", "Apr"),
        c(60, NA, NA, NA, NA, NA, NA, "36931", "Apr"),
        c(80, NA, NA, NA, NA, NA, NA, "36931", "Apr"),
        c(100, NA, NA, NA, NA, NA, NA, "36931", "Apr"),
        c(90, NA, NA, NA, NA, NA, NA, "36931", "Apr")) %>%
  as.data.frame()
names(Apr) <- names(zonal_month_36931)

May <- c(10, NA, NA, NA, NA, NA, NA, "36931", "May") %>%
  rbind(c(20, NA, NA, NA, NA, NA, NA, "36931", "May"),
        c(30, NA, NA, NA, NA, NA, NA, "36931", "May"),
        c(40, NA, NA, NA, NA, NA, NA, "36931", "May"),
        c(50, NA, NA, NA, NA, NA, NA, "36931", "May"),
        c(60, NA, NA, NA, NA, NA, NA, "36931", "May"),
        c(80, NA, NA, NA, NA, NA, NA, "36931", "May"),
        c(100, NA, NA, NA, NA, NA, NA, "36931", "May"),
        c(90, NA, NA, NA, NA, NA, NA, "36931", "May")) %>%
  as.data.frame()
names(May) <- names(zonal_month_36931)

Jun <- c(10, NA, NA, NA, NA, NA, NA, "36931", "Jun") %>%
  rbind(c(20, NA, NA, NA, NA, NA, NA, "36931", "Jun"),
        c(30, NA, NA, NA, NA, NA, NA, "36931", "Jun"),
        c(40, NA, NA, NA, NA, NA, NA, "36931", "Jun"),
        c(50, NA, NA, NA, NA, NA, NA, "36931", "Jun"),
        c(60, NA, NA, NA, NA, NA, NA, "36931", "Jun"),
        c(80, NA, NA, NA, NA, NA, NA, "36931", "Jun"),
        c(100, NA, NA, NA, NA, NA, NA, "36931", "Jun"),
        c(90, NA, NA, NA, NA, NA, NA, "36931", "Jun")) %>%
  as.data.frame()
names(Jun) <- names(zonal_month_36931)

Jul <- c(10, NA, NA, NA, NA, NA, NA, "36931", "Jul") %>%
  rbind(c(20, NA, NA, NA, NA, NA, NA, "36931", "Jul"),
        c(30, NA, NA, NA, NA, NA, NA, "36931", "Jul"),
        c(40, NA, NA, NA, NA, NA, NA, "36931", "Jul"),
        c(50, NA, NA, NA, NA, NA, NA, "36931", "Jul"),
        c(60, NA, NA, NA, NA, NA, NA, "36931", "Jul"),
        c(80, NA, NA, NA, NA, NA, NA, "36931", "Jul"),
        c(100, NA, NA, NA, NA, NA, NA, "36931", "Jul"),
        c(90, NA, NA, NA, NA, NA, NA, "36931", "Jul")) %>%
  as.data.frame()
names(Jul) <- names(zonal_month_36931)

Aug <- c(10, NA, NA, NA, NA, NA, NA, "36931", "Aug") %>%
  rbind(c(20, NA, NA, NA, NA, NA, NA, "36931", "Aug"),
        c(30, NA, NA, NA, NA, NA, NA, "36931", "Aug"),
        c(40, NA, NA, NA, NA, NA, NA, "36931", "Aug"),
        c(50, NA, NA, NA, NA, NA, NA, "36931", "Aug"),
        c(60, NA, NA, NA, NA, NA, NA, "36931", "Aug"),
        c(80, NA, NA, NA, NA, NA, NA, "36931", "Aug"),
        c(100, NA, NA, NA, NA, NA, NA, "36931", "Aug"),
        c(90, NA, NA, NA, NA, NA, NA, "36931", "Aug")) %>%
  as.data.frame()
names(Aug) <- names(zonal_month_36931)

Sep <- c(10, NA, NA, NA, NA, NA, NA, "36931", "Sep") %>% 
  rbind(c(20, NA, NA, NA, NA, NA, NA, "36931", "Sep"),
        c(30, NA, NA, NA, NA, NA, NA, "36931", "Sep"),
        c(40, NA, NA, NA, NA, NA, NA, "36931", "Sep"),
        c(50, NA, NA, NA, NA, NA, NA, "36931", "Sep"),
        c(60, NA, NA, NA, NA, NA, NA, "36931", "Sep"),
        c(80, NA, NA, NA, NA, NA, NA, "36931", "Sep"),
        c(100, NA, NA, NA, NA, NA, NA, "36931", "Sep"),
        c(90, NA, NA, NA, NA, NA, NA, "36931", "Sep")) %>%
  as.data.frame()
names(Sep) <- names(zonal_month_36931)

NA_36931 <- rbind(Feb, Mar, Apr, May, Jun, Jul, Aug, Sep) 
NA_36931$zone <- as.numeric(NA_36931$zone)
NA_36931$sum <- as.numeric(NA_36931$sum)
NA_36931$mean <- as.numeric(NA_36931$mean)
NA_36931$count <- as.numeric(NA_36931$count)
NA_36931$percent <- as.numeric(NA_36931$percent)
NA_36931$uniform <- as.numeric(NA_36931$uniform)
NA_36931$relative <- as.numeric(NA_36931$relative)

zonal_month_36931 <- rbind(zonal_month_36931, NA_36931)

levels(factor(NA_36931$zone))
levels(factor(zonal_month_36931$zone))

##9 36932
zonal_month_36932 <- lapply(1:14, function (b) zonal(rs_UDS_ESA_month[[5]][[b]], ESA_36932[[b]], fun = "sum") %>%
                        as.data.frame %>% 
                        left_join(zonal(rs_UDS_ESA_month[[5]][[b]], ESA_36932[[b]], fun = "mean") %>% 
                                    as.data.frame) %>% 
                        mutate(count = sum/mean,
                               percent = sum/sum(sum),
                               uniform = count*sum(sum)/sum(count, na.rm = T),
                               relative = sum/uniform,
                               ID = "36932",
                               month = names(ESA_36932)[b])) %>% 
  bind_rows()

# splitting the months into their appropriate years
zonal_month_36932[c("month", "year")] <- str_split_fixed(zonal_month_36932$month, "_",2)
zonal_month_36932_22 <- zonal_month_36932 %>%
  filter(year == 22) %>%
  mutate(ID = recode(ID, "36932" = "36932 - 2022"))
zonal_month_36932_23 <- zonal_month_36932 %>%
  filter(year == 23) %>%
  mutate(ID = recode(ID, "36932" = "36932 - 2023"))

## 29456
zonal_month_29456 <- lapply(1:17, function (b) zonal(rs_UDS_ESA_month[[6]][[b]], ESA_29456[[b]], fun = "sum") %>%
                        as.data.frame %>% 
                        left_join(zonal(rs_UDS_ESA_month[[6]][[b]], ESA_29456[[b]], fun = "mean") %>% 
                                    as.data.frame) %>% 
                        mutate(count = sum/mean,
                               percent = sum/sum(sum),
                               uniform = count*sum(sum)/sum(count, na.rm = T),
                               relative = sum/uniform,
                               ID = "29456",
                               month = names(ESA_29456)[b])) %>% 
  bind_rows() 

# splitting the months into their appropriate years
zonal_month_29456[c("month", "year")] <- str_split_fixed(zonal_month_29456$month, "_",2)
zonal_month_29456_20 <- zonal_month_29456 %>%
  filter(year == 20) %>%
  mutate(ID = recode(ID, "29456" = "29456 - 2020"))
zonal_month_29456_21 <- zonal_month_29456 %>%
  filter(year == 21) %>%
  mutate(ID = recode(ID, "29456" = "29456 - 2021"))

## 47016
zonal_month_47016 <- lapply(1:10, function (b) zonal(rs_UDS_ESA_month[[7]][[b]], ESA_47016[[b]], fun = "sum") %>%
                        as.data.frame %>% 
                        left_join(zonal(rs_UDS_ESA_month[[7]][[b]], ESA_47016[[b]], fun = "mean") %>% 
                                    as.data.frame) %>% 
                        mutate(count = sum/mean,
                               percent = sum/sum(sum),
                               uniform = count*sum(sum)/sum(count, na.rm = T),
                               relative = sum/uniform,
                               ID = "47016",
                               month = names(ESA_47016)[b])) %>% 
  bind_rows() 

# adding missing months in the year (filling NA with values)
Aug <- c(10, NA, NA, NA, NA, NA, NA, "47016", "Aug") %>%
  rbind(c(30, NA, NA, NA, NA, NA, NA, "47016", "Aug"),
        c(20, NA, NA, NA, NA, NA, NA, "47016", "Aug"),
        c(40, NA, NA, NA, NA, NA, NA, "47016", "Aug"),
        c(50, NA, NA, NA, NA, NA, NA, "47016", "Aug"),
        c(60, NA, NA, NA, NA, NA, NA, "47016", "Aug"),
        c(80, NA, NA, NA, NA, NA, NA, "47016", "Aug"),
        c(90, NA, NA, NA, NA, NA, NA, "47016", "Aug")) %>%
  as.data.frame()
names(Aug) <- names(zonal_month_47016)

Sep <- c(10, NA, NA, NA, NA, NA, NA, "47016", "Sep") %>%
  rbind(c(30, NA, NA, NA, NA, NA, NA, "47016", "Sep"),
        c(20, NA, NA, NA, NA, NA, NA, "47016", "Sep"),
        c(40, NA, NA, NA, NA, NA, NA, "47016", "Sep"),
        c(50, NA, NA, NA, NA, NA, NA, "47016", "Sep"),
        c(60, NA, NA, NA, NA, NA, NA, "47016", "Sep"),
        c(80, NA, NA, NA, NA, NA, NA, "47016", "Sep"),
        c(90, NA, NA, NA, NA, NA, NA, "47016", "Sep")) %>%
  as.data.frame()
names(Sep) <- names(zonal_month_47016)

NA_47016 <- rbind(Aug, Sep) 
NA_47016$zone <- as.numeric(NA_47016$zone)
NA_47016$sum <- as.numeric(NA_47016$sum)
NA_47016$mean <- as.numeric(NA_47016$mean)
NA_47016$count <- as.numeric(NA_47016$count)
NA_47016$percent <- as.numeric(NA_47016$percent)
NA_47016$uniform <- as.numeric(NA_47016$uniform)
NA_47016$relative <- as.numeric(NA_47016$relative)

zonal_month_47016 <- rbind(zonal_month_47016, NA_47016)

levels(factor(NA_47016$zone))
levels(factor(zonal_month_47016$zone))

## 30776
zonal_month_30776 <- lapply(1:3, function (b) zonal(rs_UDS_ESA_month[[8]][[b]], ESA_30776[[b]], fun = "sum") %>%
                        as.data.frame %>% 
                        left_join(zonal(rs_UDS_ESA_month[[8]][[b]], ESA_30776[[b]], fun = "mean") %>% 
                                    as.data.frame) %>% 
                        mutate(count = sum/mean,
                               percent = sum/sum(sum),
                               uniform = count*sum(sum)/sum(count, na.rm = T),
                               relative = sum/uniform,
                               ID = "30776",
                               month = names(ESA_30776)[b])) %>% 
  bind_rows() 

## 30777
zonal_month_30777 <- lapply(1:6, function (b) zonal(rs_UDS_ESA_month[[9]][[b]], ESA_30777[[b]], fun = "sum") %>%
                        as.data.frame %>% 
                        left_join(zonal(rs_UDS_ESA_month[[9]][[b]], ESA_30777[[b]], fun = "mean") %>% 
                                    as.data.frame) %>% 
                        mutate(count = sum/mean,
                               percent = sum/sum(sum),
                               uniform = count*sum(sum)/sum(count, na.rm = T),
                               relative = sum/uniform,
                               ID = "30777",
                               month = names(ESA_30777)[b])) %>% 
  bind_rows() 

## 30779
zonal_month_30779 <- lapply(1:7, function (b) zonal(rs_UDS_ESA_month[[10]][[b]], ESA_30779[[b]], fun = "sum") %>%
                        as.data.frame %>% 
                        left_join(zonal(rs_UDS_ESA_month[[10]][[b]], ESA_30779[[b]], fun = "mean") %>% 
                                    as.data.frame) %>% 
                        mutate(count = sum/mean,
                               percent = sum/sum(sum),
                               uniform = count*sum(sum)/sum(count, na.rm = T),
                               relative = sum/uniform,
                               ID = "30779",
                               month = names(ESA_30779)[b])) %>% 
  bind_rows() 

# adding missing months in the year (filling NA with values)
Jun <- c(10, NA, NA, NA, NA, NA, NA, "30779", "Jun") %>%
  rbind(c(30, NA, NA, NA, NA, NA, NA, "30779", "Jun"),
        c(20, NA, NA, NA, NA, NA, NA, "30779", "Jun"),
        c(40, NA, NA, NA, NA, NA, NA, "30779", "Jun"),
        c(50, NA, NA, NA, NA, NA, NA, "30779", "Jun"),
        c(60, NA, NA, NA, NA, NA, NA, "30779", "Jun"),
        c(80, NA, NA, NA, NA, NA, NA, "30779", "Jun"),
        c(90, NA, NA, NA, NA, NA, NA, "30779", "Jun")) %>%
  as.data.frame()
names(Jun) <- names(zonal_month_30779)

Jul <- c(10, NA, NA, NA, NA, NA, NA, "30779", "Jul") %>%
  rbind(c(30, NA, NA, NA, NA, NA, NA, "30779", "Jul"),
        c(20, NA, NA, NA, NA, NA, NA, "30779", "Jul"),
        c(40, NA, NA, NA, NA, NA, NA, "30779", "Jul"),
        c(50, NA, NA, NA, NA, NA, NA, "30779", "Jul"),
        c(60, NA, NA, NA, NA, NA, NA, "30779", "Jul"),
        c(80, NA, NA, NA, NA, NA, NA, "30779", "Jul"),
        c(90, NA, NA, NA, NA, NA, NA, "30779", "Jul")) %>%
  as.data.frame()
names(Jul) <- names(zonal_month_30779)

Aug <- c(10, NA, NA, NA, NA, NA, NA, "30779", "Aug") %>%
  rbind(c(30, NA, NA, NA, NA, NA, NA, "30779", "Aug"),
        c(20, NA, NA, NA, NA, NA, NA, "30779", "Aug"),
        c(40, NA, NA, NA, NA, NA, NA, "30779", "Aug"),
        c(50, NA, NA, NA, NA, NA, NA, "30779", "Aug"),
        c(60, NA, NA, NA, NA, NA, NA, "30779", "Aug"),
        c(80, NA, NA, NA, NA, NA, NA, "30779", "Aug"),
        c(90, NA, NA, NA, NA, NA, NA, "30779", "Aug")) %>%
  as.data.frame()
names(Aug) <- names(zonal_month_30779)

Sep <- c(10, NA, NA, NA, NA, NA, NA, "30779", "Sep") %>% 
  rbind(c(30, NA, NA, NA, NA, NA, NA, "30779", "Sep"),
        c(20, NA, NA, NA, NA, NA, NA, "30779", "Sep"),
        c(40, NA, NA, NA, NA, NA, NA, "30779", "Sep"),
        c(50, NA, NA, NA, NA, NA, NA, "30779", "Sep"),
        c(60, NA, NA, NA, NA, NA, NA, "30779", "Sep"),
        c(80, NA, NA, NA, NA, NA, NA, "30779", "Sep"),
        c(90, NA, NA, NA, NA, NA, NA, "30779", "Sep")) %>%
as.data.frame()
names(Sep) <- names(zonal_month_30779)



Oct <- c(10, NA, NA, NA, NA, NA, NA, "30779", "Oct") %>%
  rbind(c(30, NA, NA, NA, NA, NA, NA, "30779", "Oct"),
        c(20, NA, NA, NA, NA, NA, NA, "30779", "Oct"),
        c(40, NA, NA, NA, NA, NA, NA, "30779", "Oct"),
        c(50, NA, NA, NA, NA, NA, NA, "30779", "Oct"),
        c(60, NA, NA, NA, NA, NA, NA, "30779", "Oct"),
        c(80, NA, NA, NA, NA, NA, NA, "30779", "Oct"),
        c(90, NA, NA, NA, NA, NA, NA, "30779", "Oct")) %>%
  as.data.frame()
names(Oct) <- names(zonal_month_30779)



NA_30779 <- rbind(Jun, Jul, Aug, Sep, Oct) 
NA_30779$zone <- as.numeric(NA_30779$zone)
NA_30779$sum <- as.numeric(NA_30779$sum)
NA_30779$mean <- as.numeric(NA_30779$mean)
NA_30779$count <- as.numeric(NA_30779$count)
NA_30779$percent <- as.numeric(NA_30779$percent)
NA_30779$uniform <- as.numeric(NA_30779$uniform)
NA_30779$relative <- as.numeric(NA_30779$relative)

levels(factor(NA_30779$zone))
levels(factor(zonal_month_30779$zone))

zonal_month_30779 <- rbind(zonal_month_30779, NA_30779)

## 46090
zonal_month_46090 <- lapply(1:6, function (b) zonal(rs_UDS_ESA_month[[11]][[b]], ESA_46090[[b]], fun = "sum") %>%
                        as.data.frame %>% 
                        left_join(zonal(rs_UDS_ESA_month[[11]][[b]], ESA_46090[[b]], fun = "mean") %>% 
                                    as.data.frame) %>% 
                        mutate(count = sum/mean,
                               percent = sum/sum(sum),
                               uniform = count*sum(sum)/sum(count, na.rm = T),
                               relative = sum/uniform,
                               ID = "46090",
                               month = names(ESA_46090)[b])) %>% 
  bind_rows() 

# 46088
zonal_month_46088 <- lapply(1:4, function (b) zonal(rs_UDS_ESA_month[[12]][[b]], ESA_46088[[b]], fun = "sum") %>%
                        as.data.frame %>% 
                        left_join(zonal(rs_UDS_ESA_month[[12]][[b]], ESA_46088[[b]], fun = "mean") %>% 
                                    as.data.frame) %>% 
                        mutate(count = sum/mean,
                               percent = sum/sum(sum),
                               uniform = count*sum(sum)/sum(count, na.rm = T),
                               relative = sum/uniform,
                               ID = "46088",
                               month = names(ESA_46088)[b])) %>% 
  bind_rows() 

## 29456_AUT
zonal_month_29456_AUT <- lapply(1:8, function (b) zonal(rs_UDS_ESA_month[[13]][[b]], ESA_29456_AUT[[b]], fun = "sum") %>%
                            as.data.frame %>% 
                            left_join(zonal(rs_UDS_ESA_month[[13]][[b]], ESA_29456_AUT[[b]], fun = "mean") %>% 
                                        as.data.frame) %>% 
                            mutate(count = sum/mean,
                                   percent = sum/sum(sum),
                                   uniform = count*sum(sum)/sum(count, na.rm = T),
                                   relative = sum/uniform,
                                   ID = "29456_AUT",
                                   month = names(ESA_29456_AUT)[b])) %>% 
  bind_rows() 


# merging results from zonal statis together into one data frame
zonal_month_ESA <- rbind(zonal_month_36927, zonal_month_36928, zonal_month_36930, zonal_month_36931,
                         zonal_month_36932_22[,1:9], zonal_month_36932_23[,1:9],
                     zonal_month_29456_20[,1:9], zonal_month_29456_21[,1:9],
                     zonal_month_47016,
                     zonal_month_30776, zonal_month_30777, zonal_month_30779, zonal_month_46090,
                     zonal_month_46088, zonal_month_29456_AUT)

zonal_ESA_month <- zonal_month_ESA %>%
  mutate(month = factor(month, 
                        levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep","Oct", "Nov", "Dec") 
  ))

# smazat
ESA_mont <- zonal_ESA_month %>%
  mutate(Value = zone) %>%
  left_join(name)

zonal_ESA_month_sum <- zonal_ESA_month %>%
  group_by(month, zone) %>% 
  summarise(total = sum(sum),
            mean = mean(sum),
            relative = mean(relative %>%  na.omit())) %>%
  mutate(Value = zone) %>%
  left_join(name)


save(ESA_month, file = "ESA_zonal_month.RData")
