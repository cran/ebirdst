## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, 
                      message = FALSE,
                      collapse = TRUE,
                      comment = "#>",
                      out.width = "100%", 
                      fig.height = 5, 
                      fig.width = 7, 
                      fig.align = "center")
# only build vignettes local and not for R CMD check
knitr::opts_chunk$set(eval = nzchar(Sys.getenv("BUILD_VIGNETTES")))

## ----load_pipd----------------------------------------------------------------
#  library(ebirdst)
#  library(dplyr)
#  
#  # Because the non-raster data is large, there is a parameter on the
#  # ebirdst_download function that defaults to downloading only the raster data.
#  # To access the non-raster data, set tifs_only = FALSE.
#  sp_path <- ebirdst_download(species = "example_data", tifs_only = FALSE)
#  
#  # predictor importance
#  pis <- load_pis(sp_path)
#  glimpse(pis)
#  
#  # partial dependence
#  pds <- load_pds(sp_path)
#  glimpse(pds)

## ----extent-------------------------------------------------------------------
#  # define a spatiotemporal extent
#  lp_extent <- ebirdst_extent(c(xmin = -86, xmax = -83, ymin = 42, ymax = 45),
#                              t = c(0.425, 0.475))
#  print(lp_extent)
#  
#  # subset to this extent
#  pis_ss <- ebirdst_subset(pis, ext = lp_extent)
#  nrow(pis)
#  nrow(pis_ss)

## ----plot_extents-------------------------------------------------------------
#  par(mfrow = c(1, 1), mar = c(0, 0, 0, 6))
#  footprint <- stixel_footprint(sp_path, ext = lp_extent)
#  plot(footprint)

## ----plot_pis-----------------------------------------------------------------
#  # with all classes
#  plot_pis(pis, ext = lp_extent, by_cover_class = FALSE, n_top_pred = 15)
#  
#  # aggregating fragstats for cover classes
#  plot_pis(pis, ext = lp_extent, by_cover_class = TRUE, n_top_pred = 15)

## ----plot_pds-----------------------------------------------------------------
#  # in the interest of speed, run with 5 bootstrap iterations
#  # in practice, best to run with the default number of iterations (100)
#  pd_smooth <- plot_pds(pds, "solar_noon_diff", ext = lp_extent, n_bs = 5)
#  dplyr::glimpse(pd_smooth)
#  
#  # deciduous broadleaf forest
#  plot_pds(pds, "mcd12q1_lccs1_fs_c14_1500_pland", ext = lp_extent, n_bs = 5)

## ----ppms---------------------------------------------------------------------
#  ppms <- ebirdst_ppms(sp_path, ext = lp_extent)
#  plot(ppms)

## ----binary_by_time-----------------------------------------------------------
#  ppms_monthly <- ebirdst_ppms_ts(sp_path, ext = lp_extent, summarize_by = "months")
#  # plot binary kappa
#  plot(ppms_monthly, type = "binary", metric = "kappa")
#  # plot occurrence probability auc
#  plot(ppms_monthly, type = "occurrence", metric = "auc")

