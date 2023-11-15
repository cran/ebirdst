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

## ----load---------------------------------------------------------------------
#  library(ebirdst)
#  library(terra)
#  library(sf)
#  library(dplyr)
#  library(tidyr)
#  library(rnaturalearth)
#  library(geodata)
#  library(ggplot2)
#  library(fields)
#  extract <- terra::extract
#  
#  # download the example yellow-bellied sapsucker data
#  # this simplified dataset doesn't require an access key
#  ebirdst_download_status("yebsap-example", download_ranges = TRUE)
#  
#  # load seasonal mean relative abundance at 27km resolution
#  abd_seasonal <- load_raster("yebsap-example",
#                              product = "abundance",
#                              period = "seasonal",
#                              metric = "mean",
#                              resolution = "27km")
#  
#  # get the seasons corresponding to each layer
#  names(abd_seasonal)
#  
#  # extract just the breeding season relative abundance
#  abd_breeding <- abd_seasonal[["breeding"]]

## ----seasons------------------------------------------------------------------
#  ebirdst_runs %>%
#    # note that the example data are for yellow-bellied sapsucker
#    filter(species_code == "yebsap-example") %>%
#    glimpse()

## ----map_simple, echo=-1------------------------------------------------------
#  par(mar = c(0.25, 0.25, 0.25, 2))
#  plot(abd_breeding, axes = FALSE)

## ----map_extent, echo=-1------------------------------------------------------
#  par(mar = c(0.25, 0.25, 0.25, 2))
#  # boundary polygon for michigan
#  mi <- ne_states(iso_a2 = "US", returnclass = "sf") %>%
#    filter(postal == "MI") %>%
#    # project to same coordinate reference system as the raster data
#    st_transform(st_crs(abd_seasonal))
#  
#  # crop data to michigan
#  abd_breeding_mi <- crop(abd_breeding, mi)
#  
#  # map the cropped data
#  plot(abd_breeding_mi, axes = FALSE)

## ----map_projection, echo=-1--------------------------------------------------
#  par(mar = c(0.25, 0.25, 0.25, 2))
#  # load the mapping parameters
#  fac_parameters <- load_fac_map_parameters("yebsap-example")
#  crs <- fac_parameters$custom_projection
#  
#  # transform to the custom projection using nearest neighbor resampling
#  abd_projected <- project(abd_breeding_mi, crs, method = "near")
#  
#  # map the cropped and projected data
#  plot(abd_projected, axes = FALSE)

## ----map_bins, echo=-1--------------------------------------------------------
#  par(mar = c(0.25, 0.25, 0.25, 2))
#  # quantiles of non-zero values
#  v <- values(abd_projected)
#  v <- v[!is.na(v) & v > 0]
#  bins <- quantile(v, seq(0, 1, by = 0.1))
#  # add a bin for 0
#  bins <- c(0, bins)
#  
#  # status and trends palette
#  pal <- ebirdst_palettes(length(bins) - 2)
#  # add a color for zero
#  pal <- c("#e6e6e6", pal)
#  
#  # map using the quantile bins
#  plot(abd_projected, breaks = bins, col = pal, axes = FALSE)

## ----map_basemap, echo=-1-----------------------------------------------------
#  par(mar = c(0.25, 0.25, 0.25, 0.25))
#  # natural earth boundaries
#  countries <- ne_countries(returnclass = "sf") %>%
#    st_geometry() %>%
#    st_transform(crs)
#  states <- ne_states(iso_a2 = "US", returnclass = "sf") %>%
#    st_geometry() %>%
#    st_transform(crs)
#  
#  # define the map extent with the michigan polygon
#  mi_ext <- mi %>%
#    st_geometry() %>%
#    st_transform(crs)
#  plot(mi_ext)
#  # add basemap
#  plot(countries, col = "#cfcfcf", border = "#888888", add = TRUE)
#  # add data
#  plot(abd_projected,
#       breaks = bins, col = pal,
#       axes = FALSE, legend = FALSE, add = TRUE)
#  # add boundaries
#  plot(countries, col = NA, border = "#888888", lwd = 3, add = TRUE)
#  plot(states, col = NA, border = "#888888", add = TRUE)
#  
#  # add legend using the fields package
#  # label the bottom, middle, and top
#  labels <- quantile(bins, c(0, 0.5, 1))
#  label_breaks <- seq(0, 1, length.out = length(bins))
#  image.plot(zlim = c(0, 1), breaks = label_breaks, col = pal,
#             smallplot = c(0.90, 0.93, 0.15, 0.85),
#             legend.only = TRUE,
#             axis.args = list(at = c(0, 0.5, 1),
#                              labels = round(labels, 2),
#                              col.axis = "black", fg = NA,
#                              cex.axis = 0.9, lwd.ticks = 0,
#                              line = -0.5))

## ----trajectories_load--------------------------------------------------------
#  abd_median <- load_raster("yebsap-example", product = "abundance",
#                            metric = "median", resolution = "27km")
#  abd_lower <- load_raster("yebsap-example", product = "abundance",
#                           metric = "lower", resolution = "27km")
#  abd_upper <- load_raster("yebsap-example", product = "abundance",
#                           metric = "upper", resolution = "27km")

## ----trajectories_extract-----------------------------------------------------
#  # set a point
#  pt <- st_point(c(-88.1, 46.7)) %>%
#    st_sfc(crs = 4326) %>%
#    st_transform(crs = st_crs(abd_median)) %>%
#    st_coordinates()
#  
#  # extract
#  traj_median <- as.matrix(extract(abd_median, pt))[1, ]
#  traj_upper <- as.matrix(extract(abd_upper, pt))[1, ]
#  traj_lower <- as.matrix(extract(abd_lower, pt))[1, ]
#  
#  # plot trajectories
#  plot_frame <- data.frame(x = seq_len(length(traj_median)),
#                           y = unname(traj_median),
#                           lower = unname(traj_lower),
#                           upper = unname(traj_upper))
#  ggplot(plot_frame, aes(x, y)) +
#    geom_line(data = plot_frame) +
#    geom_ribbon(data = plot_frame,
#                aes(ymin = lower, ymax = upper),
#                alpha = 0.3) +
#    ylab("Relative abundance") +
#    xlab("Week") +
#    theme_light()

## ----stats_counties-----------------------------------------------------------
#  mi_counties <- gadm(country = "USA", level = 2, path = tempdir()) %>%
#    st_as_sf() %>%
#    filter(NAME_1 == "Michigan") %>%
#    select(county = NAME_2, county_code = HASC_2) %>%
#    # remove lakes which aren't true counties
#    filter(county_code != "US.MI.WB")
#  # project to sinusoidal
#  mi_counties_proj <- st_transform(mi_counties, crs = st_crs(abd_median))

## ----stats_load---------------------------------------------------------------
#  pop_seasonal <- load_raster("yebsap-example", product = "proportion-population",
#                              period = "seasonal", resolution = "27km")
#  ranges <- load_ranges("yebsap-example", resolution = "27km", smoothed = FALSE)

## ----mean-rel-abd-------------------------------------------------------------
#  prop_pop <- extract(pop_seasonal, mi_counties_proj, fun = sum, na.rm = TRUE) %>%
#    # attach county attributes
#    mutate(county_code = mi_counties$county_code) %>%
#    # transpose to long format, one season per row
#    select(-ID) %>%
#    pivot_longer(cols = -county_code,
#                 names_to = "season",
#                 values_to = "proportion_population")
#  head(prop_pop)

## ----stats_pop----------------------------------------------------------------
#  # join back to county boundaries
#  prop_pop_proj <- prop_pop %>%
#    filter(season %in% c("breeding", "nonbreeding")) %>%
#    inner_join(mi_counties, ., by = "county_code") %>%
#    # transform to custom projection for plotting
#    st_transform(crs = crs)
#  
#  # plot
#  ggplot(prop_pop_proj) +
#    geom_sf(aes(fill = proportion_population)) +
#    scale_fill_viridis_c(trans = "sqrt") +
#    guides(fill = guide_colorbar(title.position = "top", barwidth = 15)) +
#    facet_wrap(~ season, ncol = 2) +
#    labs(title = "Seasonal proportion of population in MI counties",
#         fill = "Proportion of population") +
#    theme_bw() +
#    theme(legend.position = "bottom")

## ----stats_occ----------------------------------------------------------------
#  # add the area of each region
#  mi_counties$area <- st_area(mi_counties)
#  # for each season, intersect with the county boundaries and calculate area
#  range_pct_occupied <- NULL
#  for (s in ranges$season) {
#    range_pct_occupied <- ranges %>%
#      filter(season == s) %>%
#      st_intersection(mi_counties, .) %>%
#      mutate(proportion_occupied = as.numeric(st_area(.) / area)) %>%
#      select(season, county_code, proportion_occupied) %>%
#      st_drop_geometry() %>%
#      bind_rows(range_pct_occupied, .)
#  }
#  head(range_pct_occupied)

