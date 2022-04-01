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

## ----libraries----------------------------------------------------------------
#  library(ebirdst)
#  library(raster)
#  library(sf)
#  library(smoothr)
#  library(rnaturalearth)
#  library(dplyr)
#  library(tidyr)
#  library(stringr)
#  library(ggplot2)
#  # resolve namespace conflicts
#  select <- dplyr::select
#  sf_use_s2(FALSE)

## ----st-download--------------------------------------------------------------
#  sp_path <- ebirdst_download(species = "example_data")
#  # load the abundance data
#  # this automatically labels layers with their dates
#  abd <- load_raster(sp_path, "abundance")

## ----ne-data, results="hide"--------------------------------------------------
#  proj <- "+proj=eck4 +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
#  
#  # download natural earth data
#  temp_gpkg <- tempfile(fileext = ".gpkg")
#  file.path("https://github.com/CornellLabofOrnithology/ebirdst/raw/",
#            "main/data-raw/ebirdst_gis-data.gpkg") %>%
#    download.file(temp_gpkg)
#  
#  # land polygon
#  ne_land <- read_sf(temp_gpkg, layer = "ne_land") %>%
#    st_transform(crs = proj) %>%
#    st_geometry()
#  # country lines
#  ne_country_lines <- read_sf(temp_gpkg, layer = "ne_country_lines") %>%
#    st_transform(crs = proj) %>%
#    st_geometry()
#  # state lines
#  ne_state_lines <- read_sf(temp_gpkg, layer = "ne_state_lines") %>%
#    st_transform(crs = proj) %>%
#    st_geometry()
#  # rivers
#  ne_rivers <- read_sf(temp_gpkg, layer = "ne_rivers") %>%
#    st_transform(crs = proj) %>%
#    st_geometry()
#  # lakes
#  ne_lakes <- read_sf(temp_gpkg, layer = "ne_lakes") %>%
#    st_transform(crs = proj) %>%
#    st_geometry()

## ----season-defs--------------------------------------------------------------
#  # subset to the yellow-bellied sapsucker season definitions
#  run_review <- filter(ebirdst_runs, species_code == "yebsap")
#  yebsap_dates <- run_review %>%
#    # just keep the seasonal definition columns
#    select(setdiff(matches("(start)|(end)"), matches("year_round"))) %>%
#    # transpose
#    gather("label", "date") %>%
#    # spread data so start and end dates are in separate columns
#    separate(label, c("season", "start_end"), "_(?=s|e)") %>%
#    spread(start_end, date) %>%
#    select(season, start_dt = start, end_dt = end) %>%
#    filter(season != "resident")
#  # did the season pass review
#  quality_rating <- run_review[paste0(yebsap_dates$season, "_quality")]
#  yebsap_dates$pass <- as.integer(quality_rating) > 1
#  yebsap_dates

## ----season-assignment--------------------------------------------------------
#  # dates for each abundance layer
#  weeks <- parse_raster_dates(abd)
#  # assign to seasons
#  weeks_season <- rep(NA_character_, length(weeks))
#  for (i in seq_len(nrow(yebsap_dates))) {
#    s <- yebsap_dates[i, ]
#    # skip seasona assignment if season failed
#    if (!s$pass) {
#      next()
#    }
#    # handle seasons cross jan 1 separately
#    if (s$start_dt <= s$end_dt) {
#      in_season <- weeks >= s$start_dt & weeks <= s$end_dt
#    } else {
#      in_season <- weeks >= s$start_dt | weeks <= s$end_dt
#    }
#    weeks_season[in_season] <- s$season
#  }
#  table(weeks_season)

## ----seasonal-average---------------------------------------------------------
#  # drop weeks not assigned to season
#  week_pass <- !is.na(weeks_season)
#  abd <- abd[[which(week_pass)]]
#  weeks <- weeks[week_pass]
#  weeks_season <- weeks_season[week_pass]
#  # average over weeks in season
#  mean_season <- function(s) {
#    calc(abd[[which(weeks_season == s)]], mean, na.rm = TRUE)
#  }
#  seasons <- unique(weeks_season)
#  abd_season <- lapply(seasons, mean_season) %>%
#    stack() %>%
#    setNames(seasons)
#  abd_season

## ----split-seasons------------------------------------------------------------
#  migration_threshold <- 0.4
#  mig_seasons <- c("prebreeding_migration", "postbreeding_migration")
#  if (all(mig_seasons %in% names(abd_season))) {
#    # identify areas with abundance in only one season
#    abd_nz <- abd_season[[mig_seasons]] > 0
#    just_pre <- mask(abd_nz[["prebreeding_migration"]],
#                     abd_nz[["postbreeding_migration"]],
#                     maskvalue = 1)
#    just_post <- mask(abd_nz[["postbreeding_migration"]],
#                      abd_nz[["prebreeding_migration"]],
#                      maskvalue = 1)
#    # count the number of cells with abundance in only one season
#    n_just <- cellStats(stack(just_pre, just_post), sum)
#    n_all <- cellStats(abd_nz, sum)
#    # is the proportion of one season cells above the 40% threshold
#    split_migration <- max(n_just / n_all, na.rm = TRUE) >= migration_threshold
#  } else {
#    split_migration <- FALSE
#  }
#  n_just / n_all
#  split_migration

## ----show-yr------------------------------------------------------------------
#  threshold_yearround <- 0.01
#  # decide whether to show year-round layer
#  if (nlayers(abd_season) == 4) {
#    # annual abundance
#    abd_yr <- calc(abd, fun = mean, na.rm = TRUE)
#    # mask out cells that aren't occupied year-round
#    year_round <- calc(abd_season > 0, fun = sum, na.rm = TRUE) == 4
#    abd_yr_mask <- mask(abd_yr, year_round, maskvalue = 0)
#    # determine proportion of celss that are occupied year round
#    n_yr <- cellStats(abd_yr_mask > 0, sum)
#    n_an <- cellStats(abd_yr > 0, sum)
#    # only show year round abundance if it's above 1% of range threshold
#    show_yearround <- ((n_yr / n_an) >= threshold_yearround)
#  } else {
#    show_yearround <- FALSE
#  }
#  show_yearround

## ----breaks-------------------------------------------------------------------
#  vals <- abd_season %>%
#    getValues() %>%
#    na.omit() %>%
#    as.numeric()
#  bin_breaks <- quantile(vals[vals > 0], seq(0, 1, by = 0.05))
#  lbls <- signif(bin_breaks[c("5%", "50%", "95%")], 3)
#  rm(vals)

## ----abd-map------------------------------------------------------------------
#  # project the abundance data to mollweide
#  # use nearest neighbour resampling to preserve true zeros
#  abd_season_proj <- projectRaster(abd_season, crs = proj, method = "ngb")
#  # determine spatial extent for plotting
#  ext <- calc_full_extent(abd_season_proj)
#  # set the plotting order of the seasons
#  season_order <- c("postbreeding_migration", "prebreeding_migration",
#                    "nonbreeding", "breeding")
#  
#  # prediction region, cells with predicted value in at least one week
#  pred_region <- calc(abd_season_proj, mean, na.rm = TRUE)
#  # mask to land area
#  ne_land_buffer <- st_buffer(ne_land, dist = max(res(pred_region)) / 2)
#  pred_region <- mask(pred_region, as_Spatial(ne_land_buffer))
#  
#  # remove zeros from abundance layers
#  abd_no_zero <- subs(abd_season_proj, data.frame(from = 0, to = NA),
#                      subsWithNA = FALSE)
#  
#  # set up plot area
#  par(mar = c(0 , 0, 0, 0))
#  plot(ne_land, col = "#cfcfcf", border = NA,
#       xlim = c(ext@xmin, ext@xmax),
#       ylim = c(ext@ymin, ext@ymax))
#  # prediction region and explicit zeros
#  plot(pred_region, col = "#e6e6e6", maxpixels = raster::ncell(pred_region),
#       legend = FALSE, add = TRUE)
#  # lakes
#  plot(ne_lakes, col = "#ffffff", border =  "#444444", lwd = 0.5, add = TRUE)
#  # land border
#  plot(ne_land, col = NA, border = "#444444", lwd = 0.5, add = TRUE)
#  # seasonal layer
#  plot_seasons <- intersect(season_order, names(abd_no_zero))
#  for (s in plot_seasons) {
#    # handle splitting of migration seasons into different colors
#    if (!split_migration && s %in% c("prebreeding_migration",
#                                     "postbreeding_migration")) {
#      pal_season <- "migration"
#  
#    } else {
#      pal_season <- s
#    }
#    pal <- abundance_palette(length(bin_breaks) - 1, pal_season)
#    plot(abd_no_zero[[s]], col = pal, breaks = bin_breaks,
#         maxpixels = ncell(abd_no_zero[[s]]),
#         legend = FALSE, add = TRUE)
#  }
#  # year round
#  if (show_yearround) {
#    year_round_proj <- projectRaster(year_round, crs = mollweide, method = "ngb")
#    plot(year_round_proj,
#         col = abundance_palette(length(bin_breaks$bins) - 1, "year_round"),
#         breaks = bin_breaks$bins,
#         maxpixels = ncell(year_round_proj),
#         legend = FALSE, add = TRUE)
#  }
#  # linework
#  plot(ne_rivers, col = "#ffffff", lwd = 0.75, add = TRUE)
#  plot(ne_state_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
#  plot(ne_country_lines, col = "#ffffff", lwd = 2, add = TRUE)
#  
#  # legends
#  legend_seasons <- plot_seasons
#  if (!split_migration) {
#    legend_seasons[legend_seasons %in% c("prebreeding_migration",
#                                         "postbreeding_migration")] <- "migration"
#    legend_seasons <- unique(legend_seasons)
#  }
#  if (show_yearround) {
#    legend_seasons <- c(legend_seasons, "year_round")
#  }
#  # thin out labels
#  # plot legends
#  lbl_at <- seq(0, 1, length.out = length(lbls))
#  for (i in seq_along(legend_seasons)) {
#    pal <- abundance_palette(length(bin_breaks) - 1, legend_seasons[i])
#    if (i == 1) {
#      axis_args <- list(at = lbl_at, labels = lbls, line = -1,
#                        cex.axis = 0.75, lwd = 0)
#    } else {
#      axis_args <- list(at = lbl_at, labels = rep("", length(lbls)),
#                        cex.axis = 0.75, lwd = 0)
#    }
#    legend_title <- legend_seasons[i] %>%
#      str_replace_all("_", " ") %>%
#      str_to_title()
#    fields::image.plot(zlim = c(0, 1),
#                       legend.only = TRUE,
#                       breaks = seq(0, 1, length.out = length(bin_breaks)),
#                       col = pal,
#                       smallplot = c(0.05, 0.35, 0.01 + 0.06 * i, 0.03 + 0.06 * i),
#                       horizontal = TRUE,
#                       axis.args = axis_args,
#                       legend.args = list(text = legend_title, side = 3,
#                                          cex = 0.9, col = "black", line = 0.1))
#  }
#  title("Yellow-bellied Sapsucker Relative Abundance",
#        line = -1, cex.main = 1)

## ----raster-to-polygon--------------------------------------------------------
#  # aggregate
#  abd_season_agg <- aggregate(abd_season_proj, fact = 3)
#  # raster to polygon, one season at a time
#  range <- list()
#  pred_area <- list()
#  for (s in names(abd_season_agg)) {
#    # range
#    range[[s]] <- rasterToPolygons(abd_season_agg[[s]],
#                                   fun = function(y) {y > 0},
#                                   digits = 6) %>%
#      st_as_sfc() %>%
#      # combine polygon pieces into a single multipolygon
#      st_set_precision(1e6) %>%
#      st_union() %>%
#      st_sf() %>%
#      # tag layers with season
#      mutate(season = s, layer = "range")
#    # prediction area
#    pred_area[[s]] <- rasterToPolygons(abd_season_agg[[s]],
#                                       fun = function(y) {!is.na(y)},
#                                       digits = 6) %>%
#      st_as_sfc() %>%
#      # combine polygon pieces into a single multipolygon
#      st_set_precision(1e6) %>%
#      st_union() %>%
#      st_sf() %>%
#      # tag layers with season
#      mutate(season = s, layer = "prediction_area")
#  }
#  # combine the sf objects for all seasons
#  range <- rbind(do.call(rbind, range), do.call(rbind, pred_area))
#  row.names(range) <- NULL
#  print(range)

## ----smoothr------------------------------------------------------------------
#  # clean and smooth
#  cell_area <- (1.5 * prod(res(abd_season_agg)))
#  range_smooth <- range %>%
#    # drop fragment polygons smaller than 1.5 times the aggregated cell size
#    drop_crumbs(threshold = cell_area) %>%
#    # drop holes in polygons smaller than 1.5 times the aggregated cell size
#    fill_holes(threshold = cell_area) %>%
#    # smooth the polygon edges
#    smooth(method = "ksmooth", smoothness = 2)
#  # clip zeros to land border, range to buffered land to handle coastal species
#  range_split <- split(range_smooth, range_smooth$layer)
#  range_smooth <- rbind(
#    st_intersection(range_split$range, ne_land_buffer),
#    st_intersection(range_split$prediction_area, ne_land))

## ----range-map----------------------------------------------------------------
#  # range map color palette
#  range_palette <- c(nonbreeding = "#1d6996",
#                     prebreeding_migration = "#73af48",
#                     breeding = "#cc503e",
#                     postbreeding_migration = "#edad08",
#                     migration = "#edad08",
#                     year_round = "#6f4070")
#  
#  # set up plot area
#  par(mar = c(0 , 0, 0, 0))
#  plot(ne_land, col = "#cfcfcf", border = NA,
#       xlim = c(ext@xmin, ext@xmax),
#       ylim = c(ext@ymin, ext@ymax))
#  # prediction region and explicit zeros
#  annual_pred_area <- filter(range_smooth, layer == "prediction_area") %>%
#    st_union()
#  plot(annual_pred_area, col = "#e6e6e6", border = NA, add = TRUE)
#  # lakes
#  plot(ne_lakes, col = "#ffffff", border =  "#444444", lwd = 0.5, add = TRUE)
#  # land border
#  plot(ne_land, col = NA, border = "#444444", lwd = 0.5, add = TRUE)
#  # seasonal layer
#  for (s in intersect(season_order, unique(range_smooth$season))) {
#    # handle splitting of migration seasons into different colors
#    if (!split_migration && s %in% c("prebreeding_migration",
#                                     "postbreeding_migration")) {
#      col_season <- "migration"
#    } else {
#      col_season <- s
#    }
#    rng_season <- filter(range_smooth, season == s, layer == "range") %>%
#      st_geometry()
#    plot(rng_season, col = range_palette[col_season], border = NA, add = TRUE)
#  }
#  # year round
#  if (show_yearround) {
#    # find common area between all seasons
#    range_combined <- filter(range_smooth, layer == "range")
#    range_yearround <- range_combined[1, ]
#    range_combined <- sf::st_geometry(range_combined)
#    for (i in 2:length(range_combined)) {
#      range_yearround <- sf::st_intersection(range_yearround, range_combined[i])
#    }
#    plot(st_geometry(range_yearround),
#         col = range_palette["year_round"], border = NA,
#         add = TRUE)
#  }
#  # linework
#  plot(ne_rivers, col = "#ffffff", lwd = 0.75, add = TRUE)
#  plot(ne_state_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
#  plot(ne_country_lines, col = "#ffffff", lwd = 2, add = TRUE)
#  
#  # legend
#  rng_legend <- rev(range_palette[legend_seasons])
#  names(rng_legend) <- names(rng_legend) %>%
#    str_replace_all("_", " ") %>%
#    str_to_title()
#  legend("bottomleft", legend = names(rng_legend), fill = rng_legend)
#  title("Yellow-bellied Sapsucker Seasonal Range Map",
#        line = -1, cex.main = 1)

