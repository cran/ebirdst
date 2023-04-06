## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      comment = "#>",
                      out.width = "100%",
                      fig.height = 4, 
                      fig.width = 7, 
                      fig.align = "center")
# only build vignettes locally and not for R CMD check
knitr::opts_chunk$set(eval = nzchar(Sys.getenv("BUILD_VIGNETTES")))

## ----access-------------------------------------------------------------------
#  library(ebirdst)
#  library(terra)
#  library(sf)
#  library(dplyr)
#  
#  # download a simplified example dataset for Yellow-bellied Sapsucker in Michigan
#  path <- ebirdst_download(species = "example_data")

## ----species------------------------------------------------------------------
#  glimpse(ebirdst_runs)

## ----types_dir----------------------------------------------------------------
#  # for non-example data use the species code or name instead of "example_data"
#  path <- get_species_path("example_data")

## ----types_weekly-------------------------------------------------------------
#  # weekly, low res, median occurrence
#  occ_lr <- load_raster(path, product = "occurrence", resolution = "lr")
#  occ_lr
#  # use parse_raster_dates() to get the date associated which each raster layer
#  parse_raster_dates(occ_lr)
#  
#  # weekly, low res, abundance confidence intervals
#  abd_lower <- load_raster(path, product = "abundance", metric = "lower",
#                           resolution = "lr")
#  abd_upper <- load_raster(path, product = "abundance", metric = "upper",
#                           resolution = "lr")

## ----types_seasonal-----------------------------------------------------------
#  # seasonal, low res, mean relative abundance
#  abd_seasonal_mean <- load_raster(path, product = "abundance",
#                                   period = "seasonal", metric = "mean",
#                                   resolution = "lr")
#  # season that each layer corresponds to
#  names(abd_seasonal_mean)
#  # just the breeding season layer
#  abd_seasonal_mean[["breeding"]]
#  
#  # seasonal, low res, max occurrence
#  occ_seasonal_max <- load_raster(path, product = "occurrence",
#                                  period = "seasonal", metric = "max",
#                                  resolution = "lr")

## ----types_ranges-------------------------------------------------------------
#  # seasonal, low res, smoothed ranges
#  ranges <- load_ranges(path, resolution = "lr")
#  ranges
#  
#  # subset to just the breeding season range using dplyr
#  range_breeding <- filter(ranges, season == "breeding")

