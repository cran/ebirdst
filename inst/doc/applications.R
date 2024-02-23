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

## ----packages-----------------------------------------------------------------
#  library(dplyr)
#  library(ebirdst)
#  library(fields)
#  library(ggplot2)
#  library(lubridate)
#  library(rnaturalearth)
#  library(sf)
#  library(terra)
#  library(tidyr)
#  extract <- terra::extract

## ----map-load-----------------------------------------------------------------
#  # download the yellow-bellied sapsucker data
#  ebirdst_download_status("sagthr",
#                          pattern = "abundance_seasonal_mean")
#  
#  # load seasonal mean relative abundance at 3km resolution
#  abd_seasonal <- load_raster("sagthr",
#                              product = "abundance",
#                              period = "seasonal",
#                              metric = "mean",
#                              resolution = "3km")
#  
#  # extract just the breeding season relative abundance
#  abd_breeding <- abd_seasonal[["breeding"]]

## ----map-simple, echo=-1------------------------------------------------------
#  par(mar = c(0.25, 0.25, 0.25, 2))
#  plot(abd_breeding, axes = FALSE)

## ----map-extent, echo=-1------------------------------------------------------
#  par(mar = c(0.25, 0.25, 0.25, 0.25))
#  # wyoming boundary
#  region_boundary <- ne_states(iso_a2 = "US") |>
#    filter(name == "Wyoming")
#  
#  # project boundary to match raster data
#  region_boundary_proj <- st_transform(region_boundary, st_crs(abd_breeding))
#  
#  # crop and mask to boundary of wyoming
#  abd_breeding_mask <- crop(abd_breeding, region_boundary_proj) |>
#    mask(region_boundary_proj)
#  
#  # map the cropped data
#  plot(abd_breeding_mask, axes = FALSE)

## ----map-projection, echo=-1--------------------------------------------------
#  par(mar = c(0.25, 0.25, 0.25, 2))
#  # find the centroid of the region
#  region_centroid <- region_boundary |>
#    st_geometry() |>
#    st_transform(crs = 4326) |>
#    st_centroid() |>
#    st_coordinates() |>
#    round(1)
#  
#  # define projection
#  crs_laea <- paste0("+proj=laea +lat_0=", region_centroid[2],
#                     " +lon_0=", region_centroid[1])
#  
#  # transform to the custom projection using nearest neighbor resampling
#  abd_breeding_laea <- project(abd_breeding_mask, crs_laea, method = "near") |>
#    # remove areas of the raster containing no data
#    trim()
#  
#  # map the cropped and projected data
#  plot(abd_breeding_laea, axes = FALSE, breakby = "cases")

## ----map-bins, echo=-1--------------------------------------------------------
#  par(mar = c(0.25, 0.25, 0.25, 2))
#  # quantiles of non-zero values
#  v <- values(abd_breeding_laea, na.rm = TRUE, mat = FALSE)
#  v <- v[v > 0]
#  breaks <- quantile(v, seq(0, 1, by = 0.1))
#  # add a bin for 0
#  breaks <- c(0, breaks)
#  
#  # status and trends palette
#  pal <- ebirdst_palettes(length(breaks) - 2)
#  # add a color for zero
#  pal <- c("#e6e6e6", pal)
#  
#  # map using the quantile bins
#  plot(abd_breeding_laea, breaks = breaks, col = pal, axes = FALSE)

## ----map-basemap, echo=-1-----------------------------------------------------
#  par(mar = c(0.25, 0.25, 0.25, 0.25))
#  # natural earth boundaries
#  countries <- ne_countries(returnclass = "sf") |>
#    st_geometry() |>
#    st_transform(crs_laea)
#  states <- ne_states(iso_a2 = "US") |>
#    st_geometry() |>
#    st_transform(crs_laea)
#  
#  # define the map plotting extent with the region boundary polygon
#  region_boundary_laea <- region_boundary |>
#    st_geometry() |>
#    st_transform(crs_laea)
#  plot(region_boundary_laea)
#  # add basemap
#  plot(countries, col = "#cfcfcf", border = "#888888", add = TRUE)
#  # add relative abundance
#  plot(abd_breeding_laea,
#       breaks = breaks, col = pal,
#       maxcell = ncell(abd_breeding_laea),
#       legend = FALSE, add = TRUE)
#  # add boundaries
#  lines(vect(countries), col = "#ffffff", lwd = 3)
#  lines(vect(states), col =  "#ffffff", lwd = 1.5, xpd = TRUE)
#  lines(vect(region_boundary_laea), col = "#ffffff", lwd = 3, xpd = TRUE)
#  
#  # add legend using the fields package
#  # label the bottom, middle, and top
#  labels <- quantile(breaks, c(0, 0.5, 1))
#  label_breaks <- seq(0, 1, length.out = length(breaks))
#  image.plot(zlim = c(0, 1), breaks = label_breaks, col = pal,
#             smallplot = c(0.90, 0.93, 0.15, 0.85),
#             legend.only = TRUE,
#             axis.args = list(at = c(0, 0.5, 1),
#                              labels = round(labels, 2),
#                              col.axis = "black", fg = NA,
#                              cex.axis = 0.9, lwd.ticks = 0,
#                              line = -0.5))

## ----chron--------------------------------------------------------------------
#  region_boundary <- ne_states(iso_a2 = "US") |>
#    filter(name == "Kansas")

## ----chron-single-dl----------------------------------------------------------
#  # download data if they haven't already been downloaded
#  ebirdst_download_status("Snowy Plover",
#                          pattern = "abundance_(median|upper|lower)_3km")
#  
#  # load raster data
#  abd_median <- load_raster("snoplo5", product = "abundance", metric = "median")
#  abd_lower <- load_raster("snoplo5", product = "abundance", metric = "lower")
#  abd_upper <- load_raster("snoplo5", product = "abundance", metric = "upper")
#  
#  # project region boundary to match raster data
#  region_boundary_proj <- st_transform(region_boundary, st_crs(abd_median))

## ----chron-single-region------------------------------------------------------
#  # extract values within region and calculate the mean
#  abd_median_region <- extract(abd_median, region_boundary_proj,
#                               fun = "mean", na.rm = TRUE, ID = FALSE)
#  abd_lower_region <- extract(abd_lower, region_boundary_proj,
#                              fun = "mean", na.rm = TRUE, ID = FALSE)
#  abd_upper_region <- extract(abd_upper, region_boundary_proj,
#                              fun = "mean", na.rm = TRUE, ID = FALSE)
#  
#  # transform to data frame format
#  abd_median_region <- data.frame(week = as.Date(names(abd_median_region)),
#                                  median = as.numeric(abd_median_region[1, ]))
#  abd_lower_region <- data.frame(week = as.Date(names(abd_lower_region)),
#                                 lower = as.numeric(abd_lower_region[1, ]))
#  abd_upper_region <- data.frame(week = as.Date(names(abd_upper_region)),
#                                 upper = as.numeric(abd_upper_region[1, ]))
#  
#  # combine median and confidence intervals
#  chronology <- abd_median_region |>
#    inner_join(abd_lower_region, by = "week") |>
#    inner_join(abd_upper_region, by = "week")

## ----chron-single-chart-------------------------------------------------------
#  ggplot(chronology) +
#    aes(x = week, y = median) +
#    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
#    geom_line() +
#    scale_x_date(date_labels = "%b", date_breaks = "1 month") +
#    labs(x = "Week",
#         y = "Mean relative abundance in Kansas",
#         title = "Migration chronology for Snowy Plover in Kansas")

## ----chron-multi-chron--------------------------------------------------------
#  species <- c("American Avocet", "Snowy Plover", "Hudsonian Godwit",
#               "Semipalmated Sandpiper")
#  
#  chronologies <- NULL
#  for (species in species) {
#    # download propotion of population data
#    ebirdst_download_status(species,
#                            pattern = "proportion-population_median_3km")
#  
#    # load raster data
#    prop_pop <- load_raster(species, product = "proportion-population")
#  
#    # estimate proportion of population within the region
#    prop_pop_region <- extract(prop_pop, region_boundary_proj,
#                               fun = "sum", na.rm = TRUE, ID = FALSE)
#    prop_pop_region <- data.frame(species = species,
#                                  week = as.Date(names(prop_pop_region)),
#                                  median = as.numeric(prop_pop_region[1, ]))
#  
#    # combine with other species
#    chronologies <- bind_rows(chronologies, prop_pop_region)
#  }

## ----chron-multi-chart--------------------------------------------------------
#  ggplot(chronologies) +
#    aes(x = week, y = median, color = species) +
#    geom_line(linewidth = 1) +
#    scale_x_date(date_labels = "%b", date_breaks = "1 month") +
#    scale_y_continuous(labels = scales::label_percent()) +
#    scale_color_brewer(palette = "Set1") +
#    labs(x = NULL,
#         y = "Percent of population in Kansas",
#         title = "Migration chronologies for shorebirds in Kansas") +
#      theme(legend.position = "bottom")

## ----stats-dl-----------------------------------------------------------------
#  ebirdst_download_status("Golden Eagle")

## ----stats-seasonal-load------------------------------------------------------
#  # seasonal proportion of population
#  prop_pop_seasonal <- load_raster("goleag",
#                                   product = "proportion-population",
#                                   period = "seasonal")
#  
#  # state boundaries, excluding hawaii
#  states <- ne_states(iso_a2 = "US") |>
#    filter(name != "Hawaii") |>
#    select(state = name) |>
#    # transform to match projection of raster data
#    st_transform(crs = st_crs(prop_pop_seasonal))

## ----stats-seasonal-extract---------------------------------------------------
#  state_prop_pop <- extract(prop_pop_seasonal, states,
#                            fun = "sum", na.rm = TRUE, weights = TRUE,
#                            bind = TRUE) |>
#    as.data.frame() |>
#    # sort in descending order or breeding proportion of population
#    arrange(desc(breeding))
#  head(state_prop_pop)

## ----stats-relative-prep------------------------------------------------------
#  # seasonal relative abundance
#  abd_seasonal <- load_raster("goleag",
#                              product = "abundance",
#                              period = "seasonal")
#  
#  # load country polygon, union into a single polygon, and project
#  noram <- ne_countries(country = c("United States of America",
#                                    "Canada", "Mexico")) |>
#    st_union() |>
#    st_transform(crs = st_crs(abd_seasonal)) |>
#    # vect converts an sf object to terra format for mask()
#    vect()
#  
#  # mask seasonal abundance
#  abd_seasonal_noram <- mask(abd_seasonal, noram)
#  
#  # total north american relative abundance for each season
#  abd_noram_total <- global(abd_seasonal_noram, fun = "sum", na.rm = TRUE)
#  
#  # proportion of north american population
#  prop_pop_noram <- abd_seasonal_noram / abd_noram_total$sum

## ----stats-relative-calc------------------------------------------------------
#  state_prop_noram_pop <- extract(prop_pop_noram, states,
#                                  fun = "sum", na.rm = TRUE, weights = TRUE,
#                                  bind = TRUE) |>
#    as.data.frame() |>
#    # sort in descending order or breeding proportion of population
#    arrange(desc(breeding))
#  head(state_prop_noram_pop)

## ----stats-custom-calc--------------------------------------------------------
#  # weekly relative abundance, masked to north america
#  abd_weekly_noram <- load_raster("goleag",
#                                  product = "abundance",
#                                  resolution = "27km") |>
#    mask(noram)
#  
#  # total north american relative abundance for each week
#  abd_weekly_total <- global(abd_weekly_noram, fun = "sum", na.rm = TRUE)
#  
#  # proportion of north american population
#  prop_pop_weekly_noram <- abd_weekly_noram / abd_weekly_total$sum
#  
#  # proportion of weekly population in california
#  california <- filter(states, state == "California")
#  cali_prop_noram_pop <- extract(prop_pop_weekly_noram, california,
#                                 fun = "sum", na.rm = TRUE,
#                                 weights = TRUE, ID = FALSE)
#  prop_pop_weekly_noram <- data.frame(
#    week = as.Date(names(cali_prop_noram_pop)),
#    prop_pop = as.numeric(cali_prop_noram_pop[1, ]))
#  head(prop_pop_weekly_noram)

## ----stats-custom-month-------------------------------------------------------
#  prop_pop_weekly_noram |>
#    filter(month(week) == 1) |>
#    summarize(prop_pop = mean(prop_pop))

## ----stats-coastal-wrong------------------------------------------------------
#  # download only the season proportion of population layer
#  ebirdst_download_status("Surf Scoter",
#                          pattern = "proportion-population_seasonal_mean_3km")
#  
#  # breeding season proportion of population
#  abd_nonbreeding <- load_raster("Surf Scoter",
#                                 product = "proportion-population",
#                                 period = "seasonal") |>
#    subset("nonbreeding")
#  
#  # load a polygon for the boundary of Mexico
#  mexico <- ne_countries(country = "Mexico") |>
#    st_transform(crs = st_crs(abd_nonbreeding))
#  
#  # proportion in mexico
#  extract(abd_nonbreeding, mexico, fun = "sum", na.rm = TRUE)

## ----stats-coastal-buffer-----------------------------------------------------
#  # buffer by 5000m = 5km
#  mexico_buffer <- st_buffer(mexico, dist = 5000)
#  
#  # proportion in mexico
#  extract(abd_nonbreeding, mexico_buffer, fun = "sum", na.rm = TRUE,
#          touches = TRUE)
#  

