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
# library(dplyr)
# library(sf)
# library(terra)
# library(ebirdst)
# 
# # download a simplified example dataset for Yellow-bellied Sapsucker in Michigan
# ebirdst_download_status(species = "yebsap-example", download_all = TRUE)

## ----species------------------------------------------------------------------
# glimpse(ebirdst_runs)

## ----review-------------------------------------------------------------------
# ebirdst_runs |>
#   filter(species_code == "yebsap-example") |>
#   glimpse()

## ----types_weekly-------------------------------------------------------------
# # weekly, 27km res, median relative abundance
# abd_lr <- load_raster("yebsap-example", product = "abundance",
#                       resolution = "27km")
# 
# # weekly, 27km res, median proportion of population
# prop_pop_lr <- load_raster("yebsap-example", product = "proportion-population",
#                       resolution = "27km")
# 
# # weekly, 27km res, abundance confidence intervals
# abd_lower <- load_raster("yebsap-example", product = "abundance", metric = "lower",
#                          resolution = "27km")
# abd_upper <- load_raster("yebsap-example", product = "abundance", metric = "upper",
#                          resolution = "27km")

## ----types_weekly_dates-------------------------------------------------------
# as.Date(names(abd_lr))

## ----types_seasonal-----------------------------------------------------------
# # seasonal, 27km res, mean relative abundance
# abd_seasonal_mean <- load_raster("yebsap-example", product = "abundance",
#                                  period = "seasonal", metric = "mean",
#                                  resolution = "27km")
# # season that each layer corresponds to
# names(abd_seasonal_mean)
# # just the breeding season layer
# abd_seasonal_mean[["breeding"]]
# 
# # seasonal, 27km res, max occurrence
# occ_seasonal_max <- load_raster("yebsap-example", product = "occurrence",
#                                 period = "seasonal", metric = "max",
#                                 resolution = "27km")

## ----types_fullyear-----------------------------------------------------------
# # full year, 27km res, maximum relative abundance
# abd_fy_max <- load_raster("yebsap-example", product = "abundance",
#                           period = "full-year", metric = "max",
#                           resolution = "27km")

## ----types_ranges-------------------------------------------------------------
# # seasonal, 27km res, smoothed ranges
# ranges <- load_ranges("yebsap-example", resolution = "27km")
# ranges
# 
# # subset to just the breeding season range using dplyr
# range_breeding <- filter(ranges, season == "breeding")

## ----types_regional-----------------------------------------------------------
# regional <- load_regional_stats("yebsap-example")
# glimpse(regional)

## ----types_ppms---------------------------------------------------------------
# pr_auc <- load_ppm("yebsap-example", ppm = "occ_pr_auc_normalized")
# print(pr_auc)

## ----types_ppms_plot----------------------------------------------------------
# plot(trim(pr_auc[[26]]))

## ----coverage_dl, eval=FALSE--------------------------------------------------
# ebirdst_download_data_coverage()

## ----coverage_load, echo=-1---------------------------------------------------
# par(mar = c(0, 0, 0, 0))
# site_sel <- load_data_coverage("selection-probability", week = "05-10")
# plot(site_sel, axes = FALSE)

