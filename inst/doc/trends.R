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

## ----runs---------------------------------------------------------------------
# library(dplyr)
# library(ggplot2)
# library(rnaturalearth)
# library(sf)
# library(terra)
# library(ebirdst)
# 
# trends_runs <- ebirdst_runs |>
#   filter(has_trends) |>
#   select(species_code, common_name,
#          trends_season, trends_region,
#          trends_start_year, trends_end_year,
#          trends_start_date, trends_end_date,
#          rsquared, beta0, trends_version_year)
# glimpse(trends_runs)

## ----crossing-----------------------------------------------------------------
# trends_runs |>
#   filter(common_name == "Canvasback") |>
#   select(trends_start_year, trends_end_year,
#          trends_start_date, trends_end_date)

## ----download-----------------------------------------------------------------
# ebirdst_download_trends("Sage Thrasher")

## ----load---------------------------------------------------------------------
# trends_sagthr <- load_trends("Sage Thrasher")

## ----dates--------------------------------------------------------------------
# trends_runs |>
#   filter(common_name == "Sage Thrasher") |>
#   select(trends_start_year, trends_end_year,
#          trends_start_date, trends_end_date)

## ----spatial-points-----------------------------------------------------------
# trends_sf <- st_as_sf(trends_sagthr,
#                       coords = c("longitude", "latitude"),
#                       crs = 4326)
# print(trends_sf)

## ----spatial-points-export, eval=FALSE----------------------------------------
# # be sure to modify the path to the file to save the file to directory of
# # your choice on your hard drive
# write_sf(trends_sf, "ebird-trends_sagthr_2022.gpkg",
#          layer = "sagthr_trends")

## ----spatial-raster-----------------------------------------------------------
# # rasterize the percent per year trend with confidence limits (default)
# ppy_raster <- rasterize_trends(trends_sagthr)
# print(ppy_raster)
# # rasterize the cumulative trend estimate
# trends_raster <- rasterize_trends(trends_sagthr, layers = "abd_trend")
# print(trends_raster)

## ----spatial-raster-export, eval=FALSE----------------------------------------
# writeRaster(trends_raster, filename = "ebird-trends_sagthr_2021.tif")

## ----spatial-raster-map-------------------------------------------------------
# # define breaks and palettes similar to those on status and trends website
# breaks <- seq(-4, 4)
# breaks[1] <- -Inf
# breaks[length(breaks)] <- Inf
# pal <- ebirdst_palettes(length(breaks) - 1, type = "trends")
# 
# # make a simple map
# plot(ppy_raster[["abd_ppy"]],
#      col = pal, breaks =  breaks,
#      main = "Sage Thrasher breeding trend 2012-2022 [% change per year]",
#      cex.main = 0.75,
#      axes = FALSE)

## ----uncertainty--------------------------------------------------------------
# trends_sagthr_folds <- load_trends("sagthr", fold_estimates = TRUE)
# print(trends_sagthr_folds)

## ----applications-regional----------------------------------------------------
# # boundaries of states in the united states
# states <- ne_states(iso_a2 = "US", returnclass = "sf") |>
#   filter(iso_a2 == "US", !postal %in% c("AK", "HI")) |>
#   transmute(state = iso_3166_2)
# 
# # convert fold-level trends estimates to sf format
# trends_sagthr_sf <-  st_as_sf(trends_sagthr_folds,
#                               coords = c("longitude", "latitude"),
#                               crs = 4326)
# 
# # attach state to the fold-level trends data
# trends_sagthr_sf <- st_join(trends_sagthr_sf, states, left = FALSE)
# 
# # abundance-weighted average trend by region and fold
# trends_states_folds <- trends_sagthr_sf |>
#   st_drop_geometry() |>
#   group_by(state, fold) |>
#   summarize(abd_ppy = sum(abd * abd_ppy) / sum(abd),
#             .groups = "drop")
# 
# # summarize across folds for each state
# trends_states <- trends_states_folds |>
#   group_by(state) |>
#   summarise(abd_ppy_median = median(abd_ppy, na.rm = TRUE),
#             abd_ppy_lower = quantile(abd_ppy, 0.10, na.rm = TRUE),
#             abd_ppy_upper = quantile(abd_ppy, 0.90, na.rm = TRUE),
#             .groups = "drop") |>
#   arrange(abd_ppy_median)

## ----applications-regional-plot-----------------------------------------------
# trends_states_sf <- left_join(states, trends_states, by = "state")
# ggplot(trends_states_sf) +
#   geom_sf(aes(fill = abd_ppy_median)) +
#   scale_fill_distiller(palette = "Reds",
#                        limits = c(NA, 0),
#                        na.value = "grey80") +
#   guides(fill = guide_colorbar(title.position = "top", barwidth = 15)) +
#   labs(title = "Sage Thrasher state-level breeding trends 2012-2022",
#        fill = "Relative abundance trend [% change / year]") +
#   theme_bw() +
#   theme(legend.position = "bottom")

## ----applications-multi-runs--------------------------------------------------
# sagebrush_species <- c("Brewer's Sparrow", "Sagebrush Sparrow", "Sage Thrasher")
# trends_runs |>
#   filter(common_name %in% sagebrush_species)

## ----applications-multi-dl----------------------------------------------------
# ebirdst_download_trends(sagebrush_species)

## ----applications-multi-------------------------------------------------------
# trends_sagebrush_species <- load_trends(sagebrush_species)
# 
# # calculate mean trend for each cell
# trends_sagebrush <- trends_sagebrush_species |>
#   group_by(srd_id, latitude, longitude) |>
#   summarize(n_species = n(),
#             abd_ppy = mean(abd_ppy, na.rm = TRUE),
#             .groups = "drop")
# print(trends_sagebrush)

## ----applications-multi-map---------------------------------------------------
# # convert the points to sf format
# all_species <- trends_sagebrush |>
#   filter(n_species == length(sagebrush_species)) |>
#   st_as_sf(coords = c("longitude", "latitude"),
#            crs = 4326)
# 
# # make a map
# ggplot(all_species) +
#   geom_sf(aes(color = abd_ppy), size = 2) +
#   scale_color_gradient2(low = "#CB181D", high = "#2171B5",
#                         limits = c(-4, 4),
#                         oob = scales::oob_squish) +
#   guides(color = guide_colorbar(title.position = "left", barheight = 15)) +
#   labs(title = "Sagebrush species breeding trends (2012-2022)",
#        color = "Relative abundance trend [% change / year]") +
#   theme_bw() +
#   theme(legend.title = element_text(angle = 90))

