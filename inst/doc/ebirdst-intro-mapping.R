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

## ----load_raster_stack--------------------------------------------------------
#  library(ebirdst)
#  library(raster)
#  library(sf)
#  library(rnaturalearth)
#  library(ggplot2)
#  library(viridisLite)
#  # handle namespace conflicts
#  extract <- raster::extract
#  
#  # Currently, example data is available on a public s3 bucket. The following
#  # ebirdst_download() function copies the species results to a selected path and
#  # returns the full path of the results. In this vignette, we'll use example data
#  # for Yellow-bellied Sapsucker
#  sp_path <- ebirdst_download(species = "example_data")
#  
#  # load median abundances and label dates
#  abd <- load_raster(sp_path, "abundance")
#  
#  # use parse_raster_dates() to get actual date objects for each layer
#  date_vector <- parse_raster_dates(abd)
#  print(date_vector)

## ----lowres-------------------------------------------------------------------
#  abd_lr <- load_raster(sp_path, "abundance", resolution = "lr")
#  # native resolution
#  res(abd)
#  # low resolution
#  res(abd_lr)

## ----mr-----------------------------------------------------------------------
#  abd <- load_raster(sp_path, "abundance", resolution = "mr")

## ----project_stack------------------------------------------------------------
#  # define an eckert iv projection centered on the western hemisphere
#  eck <- "+proj=eck4 +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84"
#  
#  # project single layer from stack to mollweide
#  week38_proj <- projectRaster(abd[[38]], crs = eck, method = "ngb")
#  week38_proj <- trim(week38_proj)
#  
#  # optionally, you can project an entire stack, but it takes much longer
#  # abd_proj <- projectRaster(abund, crs = eck, method = "ngb")
#  
#  # map single layer with full annual extent
#  par(mar = c(0, 0, 0, 2))
#  plot(week38_proj,
#       col = abundance_palette(10, season = "weekly"),
#       axes = FALSE, box = FALSE,
#       maxpixels = ncell(week38_proj))

## ----map_occurrence-----------------------------------------------------------
#  occ <- load_raster(sp_path, "occurrence", resolution = "mr")
#  
#  # select a week in the summer
#  occ_week26 <- occ[[26]]
#  
#  # create breaks every 0.1 from 0 to 1
#  occ_bins <- seq(0, 1, by = 0.1)
#  occ_proj <- projectRaster(occ_week26, crs = eck, method = "ngb")
#  occ_proj <- trim(occ_proj)
#  
#  par(mar = c(0, 0, 0, 2), cex = 0.9)
#  plot(occ_proj,
#       breaks = occ_bins,
#       col = abundance_palette(length(occ_bins) - 1, season = "weekly"),
#       axes = FALSE, box = FALSE,
#       maxpixels = ncell(occ_proj),
#       legend.width = 2, legend.shrink = 0.97)

## ----map_linear---------------------------------------------------------------
#  year_max <- max(maxValue(abd), na.rm = TRUE)
#  
#  week14_proj <- projectRaster(abd[[14]], crs = eck, method = "ngb")
#  week14_proj <- trim(week14_proj)
#  
#  # set graphical params
#  par(mfrow = c(1, 2), mar = c(0, 0, 0, 4))
#  
#  # use raster bounding box to set the spatial extent for the plot
#  bb <- st_as_sfc(st_bbox(week14_proj))
#  plot(bb, col = "white", border = "white")
#  # plot the abundance
#  plot(week38_proj, zlim = c(0, year_max),
#       col = abundance_palette(20, season = "weekly"),
#       maxpixels = ncell(week38_proj),
#       axes = FALSE, box = FALSE, legend = FALSE, add = TRUE)
#  
#  # do the same for week 14
#  par(mar = c(0, 0, 0, 4))
#  plot(bb, col = "white", border = "white")
#  plot(week14_proj, zlim = c(0, year_max),
#       col = abundance_palette(20, season = "weekly"),
#       maxpixels = ncell(week14_proj),
#       axes = FALSE, box = FALSE, legend.shrink = 0.97, add = TRUE)

## ----map_bins-----------------------------------------------------------------
#  # calculate ideal color bins for relative abundance
#  count <- load_raster(sp_path, "count", resolution = "mr")
#  year_bins <- calc_bins(abundance = abd, count = count)
#  
#  # plot
#  par(mfrow = c(1, 2), mar = c(0, 0, 0, 6))
#  plot(bb, col = "white", border = "white")
#  plot(week38_proj,
#       breaks = year_bins,
#       col = abundance_palette(length(year_bins) - 1, season = "weekly"),
#       maxpixels = ncell(week38_proj),
#       axes = FALSE, box = FALSE, legend = FALSE, add = TRUE)
#  par(mar = c(0, 0, 0, 6))
#  plot(bb, col = "white", border = "white")
#  plot(week14_proj,
#       breaks = year_bins,
#       col = abundance_palette(length(year_bins) - 1, season = "weekly"),
#       maxpixels = ncell(week14_proj),
#       axes = FALSE, box = FALSE, legend = FALSE, add = TRUE)
#  
#  # plot legend
#  lbls <- attr(year_bins, "labels")
#  plot(week14_proj, zlim = c(0, 1), legend.only = TRUE,
#       col = abundance_palette(length(year_bins) - 1, season = "weekly"),
#       breaks = seq(0, 1, length.out = length(year_bins)),
#       legend.shrink = 0.97, legend.width = 2,
#       axis.args = list(at = seq(0, 1, length.out = length(lbls)),
#                        labels = signif(lbls, 3),
#                        col.axis = "black", fg = NA,
#                        cex.axis = 0.9, lwd.ticks = 0,
#                        line = -0.5))

## ----map_w_es-----------------------------------------------------------------
#  # to add context, let's pull in some reference data to add
#  wh_states <- ne_states(country = c("United States of America", "Canada"),
#                      returnclass = "sf") %>%
#    st_transform(crs = eck) %>%
#    st_geometry()
#  
#  # we'll plot a week in the middle of summer
#  week26_proj <- projectRaster(abd[[26]], crs = eck, method = "ngb")
#  week26_proj <- trim(week26_proj)
#  
#  # set graphics params
#  par(mfrow = c(1, 1), mar = c(0, 0, 0, 6))
#  
#  # use the extent object to set the spatial extent for the plot
#  plot(st_as_sfc(st_bbox(week26_proj)), col = "white", border = "white")
#  
#  # add background spatial context
#  plot(wh_states, col = "#cfcfcf", border = NA, add = TRUE)
#  
#  # plot zeroes as gray
#  plot(week26_proj, col = "#e6e6e6",
#       maxpixels = ncell(week26_proj),
#       axes = FALSE, legend = FALSE, add = TRUE)
#  
#  # define color bins
#  qcol <- abundance_palette(length(year_bins) - 1, "weekly")
#  
#  # plot abundances
#  plot(week26_proj, breaks = year_bins, col = qcol,
#       maxpixels = ncell(week26_proj),
#       axes = FALSE, legend = FALSE, add = TRUE)
#  
#  # add zero to the bin labels
#  bin_labels <- signif(c(0, attr(year_bins, "labels")), 3)
#  
#  # create colors that include gray for 0
#  lcol <- c("#e6e6e6", qcol)
#  
#  # plot legend
#  plot(week14_proj, zlim = c(-0.05, 1), legend.only = TRUE,
#       col = lcol,
#       breaks = c(-0.05, seq(0, 1, length.out = length(lcol))),
#       legend.shrink = 0.97, legend.width = 2,
#       axis.args = list(at = c(-0.05, seq(0, 1, length.out = length(bin_labels) - 1)),
#                        labels = bin_labels,
#                        col.axis = "black", fg = NA,
#                        cex.axis = 0.9, lwd.ticks = 0,
#                        line = -0.5))
#  
#  # add state boundaries on top
#  plot(wh_states, border = "white", lwd = 1.5, add = TRUE)

## ----map_confidence_band------------------------------------------------------
#  # load lower and upper stacks
#  lower <- load_raster(sp_path, "abundance_lower", resolution = "mr")
#  upper <- load_raster(sp_path, "abundance_upper", resolution = "mr")
#  
#  # calculate band width
#  conf_band <- upper[[26]] - lower[[26]]
#  
#  conf_week26 <- projectRaster(conf_band, crs = eck, method = "ngb")
#  
#  par(mar = c(0, 0, 0, 2))
#  plot(trim(conf_week26), col = magma(20),
#       maxpixel = ncell(conf_week26),
#       axes = FALSE, box = FALSE)

## ----trajectories-------------------------------------------------------------
#  # set a point
#  pt <- st_point(c(-88.1, 46.7)) %>%
#    st_sfc(crs = 4326) %>%
#    st_transform(crs = projection(abd)) %>%
#    st_coordinates()
#  
#  # extract
#  abd_traj <- extract(abd, pt, fun = mean, na.rm = TRUE)[1, ]
#  upper_traj <- extract(upper, pt, fun = mean, na.rm = TRUE)[1, ]
#  lower_traj <- extract(lower, pt, fun = mean, na.rm = TRUE)[1, ]
#  
#  # plot trajectories
#  plot_frame <- data.frame(x = 1:length(abd_traj),
#                           y = unname(abd_traj),
#                           upper = unname(upper_traj),
#                           lower = unname(lower_traj))
#  
#  ggplot(plot_frame, aes(x, y)) +
#    geom_line(data = plot_frame) +
#    geom_ribbon(data = plot_frame,
#                aes(ymin = lower, ymax = upper),
#                alpha = 0.3) +
#    ylab("Expected Count") +
#    xlab("Week") +
#    theme_light()

## ----trajectories_region------------------------------------------------------
#  # set an extent based on polygon
#  mi <- ne_states(country = "United States of America", returnclass = "sf") %>%
#    st_transform(crs = projection(abd))
#  mi <- mi[mi$name == "Michigan"]
#  
#  # extract
#  # because we're using a region, we get lots of values that we need to average
#  abd_traj <- extract(abd, mi, fun = mean, na.rm = TRUE)
#  abd_traj <- apply(abd_traj, 2, mean, na.rm = TRUE)
#  
#  upper_traj <- extract(upper, mi, fun = mean, na.rm = TRUE)
#  upper_traj <- apply(upper_traj, 2, mean, na.rm = TRUE)
#  
#  lower_traj <- extract(lower, mi, fun = mean, na.rm = TRUE)
#  lower_traj <- apply(lower_traj, 2, mean, na.rm = TRUE)
#  
#  # plot trajectories
#  plot_frame <- data.frame(x = 1:length(abd_traj),
#                           y = unname(abd_traj),
#                           upper = unname(upper_traj),
#                           lower = unname(lower_traj))
#  
#  ggplot(plot_frame, aes(x, y)) +
#    geom_line(data = plot_frame) +
#    geom_ribbon(data = plot_frame,
#                aes(ymin = lower, ymax =upper),
#                alpha = 0.3) +
#    ylab("Expected Count") +
#    xlab("Week") +
#    theme_light()

