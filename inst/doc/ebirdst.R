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

## ----quick_start--------------------------------------------------------------
#  library(ebirdst)
#  
#  # download data
#  # download a simplified example dataset from aws s3
#  # example data are for Yellow-bellied Sapsucker in Michigan
#  # by default file will be stored in a persistent data directory:
#  # rappdirs::user_data_dir("ebirdst"))
#  sp_path <- ebirdst_download(species = "example_data")

## ----conversion---------------------------------------------------------------
#  library(raster)
#  
#  # load median abundances
#  abd <- load_raster("abundance", path = sp_path)
#  
#  # use parse_raster_dates() to get actual date objects for each layer
#  date_vector <- parse_raster_dates(abd)
#  
#  # to convert the data to a simpler geographic format and access tabularly	
#  # reproject into geographic (decimal degrees)
#  abd_ll <- projectRaster(abd[[26]], crs = "+init=epsg:4326", method = "ngb")
#  
#  # Convert raster object into a matrix
#  p <- rasterToPoints(abd_ll)
#  colnames(p) <- c("longitude", "latitude", "abundance_umean")
#  head(p)

## ----conversion_write, eval = FALSE-------------------------------------------
#  write.csv(p, file = "yebsap_week26.csv", row.names = FALSE)

