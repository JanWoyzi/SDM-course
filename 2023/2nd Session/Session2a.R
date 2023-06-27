
remotes::install_github("ropensci/scrubr")
library("scrubr")

install.packages(c("rgbif", "sdmpredictors", "fuzzySim"))
library(rgbif)
library(sdmpredictors)
library(fuzzySim)
library(terra)

dir.create("SDM-course-main/2023/2nd Session")
setwd("./SDM-course-main/2023/2nd Session/")
file.copy("~/nextcloud/5405.AG_Schnittler/Lectures/Lecture material/Seminar Plant Conversation/2023/2nd Session/data/", "./", recursive = T)


# DOWNLOAD SPECIES OCCURRENCE DATA ####

# here we'll download GBIF occurrence data for an example species; after running the script and 
# understanding how everything works, you can replace this with another species of your choice 
# (as long as it exists on GBIF) and run it again; but note things can be quite slow if there 
# are many occurrence points!

myspecies <- "Physarum albescens"

# mind that data download takes a long time when there are many occurrences!
gbif_data <- occ_data(scientificName = myspecies, hasCoordinate = TRUE, limit = 1000)

# if your species is widespread, you can download points only within a specified window of 
# longitude and latitude coordinates (otherwise, scripts may take too long to run during the course!):
countries <- vect("./data/countries/world_countries.shp")  # "../" means go up one level from the 
# current folder or working directory

plot(countries)  # see the coordinates in the plot axes

# set the xmin, xmax, ymin, ymax coordinates of your desired extent:
my_window <- ext(c(-14, 51, 30, 80))

plot(as.polygons(my_window), border = "red", lwd = 2, add = TRUE)  # check if it's where you want it on the map

# plot country borders within 'my_window' only:
plot(countries, xlim = my_window[1:2], ylim = my_window[3:4])

# if global data is too much, download GBIF data from this window only:
gbif_data <- occ_data(scientificName = myspecies, hasCoordinate = TRUE, limit = 1000, 
                      decimalLongitude = paste0(my_window[1:2], collapse = ", "), 
                      decimalLatitude = paste0(my_window[3:4], collapse = ", "))

gbif_data_all <- occ_data(scientificName = myspecies, hasCoordinate = TRUE, limit = 1000)

# fill in your gbif.org credentials 
user <- "woyzichovj" # your gbif.org username 
pwd <- "Stresemann1!GBIF" # your gbif.org password
email <- "jan.woyzichovski@uni-greifswald.de" # your email 

# New way:
query <- occ_download(pred("taxonKey", 3214840),
                          format = "SIMPLE_CSV",
                          user=user, pwd=pwd, email=email)

# Check status with
occ_download_wait(query[1])

# After it finishes, use:
gbif_down <- occ_download_get(query[1], overwrite = T) %>%
  occ_download_import()


gbif_data  # scroll up; if "Records found" is larger than "Records returned", you need to increase 
# the 'limit' argument above (or decrease the coordinate window size) -- see help(occ_data) for options
# and limitations

gbif_citation(gbif_data)  # NOTE: If you plan to use GBIF data in any report or publication, you need to 
# either cite all the references that you get here, or download the data directly from www.gbif.org 
# (then import the .csv to R) and note down the DOI and citation for that particular dataset. 
# It is very important to properly cite the data sources! GBIF is not a source, just a repository 
# for many people who put in very hard work to collect these data and make them available.

# check how the data are organized:
names(gbif_data)
names(gbif_data$meta)
names(gbif_data$data)

occurrences <- gbif_data$data #not neccessary if you used the function: occ_download_get()!

# map the occurrence records to see how they look:
plot(countries, xlim = range(my_window[1:2]), ylim = my_window[3:4])
points(occurrences[ , c("decimalLongitude", "decimalLatitude")], pch = 04, col = "red")  # compare e.g. with the range map of this species at https://www.iucnredlist.org to assess if the distribution is well represented

#for all occurencies:
plot(countries)  # see the coordinates in the plot axes
points(occurrences[ , c("decimalLongitude", "decimalLatitude")], pch = 04, col = "red")  # compare e.g. with the range map of this species at https://www.iucnredlist.org to assess if the distribution is well represented
# points(gbif_down[ , c("decimalLongitude", "decimalLatitude")], pch = 04, col = "blue")

# export the data as a .csv file:
# first, check if some columns are lists, which cannot be saved in a .csv
unique(sapply(occurrences, class))
which(sapply(occurrences, class) == "list")
# remove list columns:
occurrences[which(sapply(occurrences, class) == "list")] <- NULL

# create a folder for the output files of this course (if it doesn't already exist)
if (!file.exists("./outputs")) dir.create("./outputs")  # '../' goes up one level from the current working directory, so this creates the 'outputs' folder just outside the 'Rscripts' folder

# create a folder on disk and save the data there:
occurrences_dir <- "./outputs/species_occurrences"
if(!file.exists(occurrences_dir)) dir.create(occurrences_dir)

write.csv(occurrences, paste0(occurrences_dir, "/occurrences_", myspecies, "_raw.csv"), row.names = FALSE)

# from here on, you don't need to download these data again from GBIF - you can just import them from the .csv:
occurrences <- read.csv(paste0(occurrences_dir, "/occurrences_", myspecies, "_raw.csv"))

# CLEAN SPECIES OCCURRENCE DATA ####

# mind that data may contain many errors; careful mapping, inspection and cleaning are necessary!
# here we'll first remove records of absence or zero-abundance (if any):
names(occurrences)
sort(unique(occurrences$occurrenceStatus))  # check for different indications of "absent", which could be in different languages!
absence_rows <- which(occurrences$occurrenceStatus %in% c("absent", "Absent", "ABSENT", "ausente", "Ausente", "AUSENTE"))
nrow(occurrences)
if (length(absence_rows) > 0)  occurrences <- occurrences[-absence_rows, ]
nrow(occurrences)

# let's do some further data cleaning with functions of the 'scrubr' package (but note this cleaning is not exhaustive!)
occurrences <- coord_incomplete(coord_imprecise(coord_impossible(coord_unlikely(occurrences))))
nrow(occurrences)

# add the cleaned occurrence data to the map:
points(occurrences[ , c("decimalLongitude", "decimalLatitude")], pch = 04, col = "green")  # excluded points are not added in gree colour, so they remain visible in red

# also eliminate presences with reported coordinate uncertainty (location error, spatial resolution) larger than 10x10 km2:
max_uncertainty <- (10000 * sqrt(2)) / 2  # meters from centroid to corner (half the diagonal) of a square with 10-km side))
occurrences <- coord_uncertain(occurrences, coorduncertainityLimit = max_uncertainty)
nrow(occurrences)
# but note that this will only discard records where coordinate uncertainty is adequately reported in the dataset, which may not always be the case! Careful mapping and visual inspection are necessary

# add these less uncertain occurrence records with a different colour on top of the previous ones:
points(occurrences[ , c("decimalLongitude", "decimalLatitude")], pch = 04, col = "turquoise")  # excluded points are not added in turquoise colour, so they remain in red or blue

library(dplyr)
# remove specific issues stated by gbif
occurrences <- occurrences %>% filter(!grepl('COUNTRY_COORDINATE_MISMATCH', issues))
nrow(occurrences)

# remove specific collection IDs
names(occurrences)
occurrences <- filter(occurrences, !grepl('MGYA', occurrences$occurrenceID))
nrow(occurrences)

points(occurrences[ , c("decimalLongitude", "decimalLatitude")], pch = 04, col = "orange")

# save the cleaned data to disk as a .csv file:
write.csv(occurrences, paste0(occurrences_dir, "/occurrences_", myspecies, "_cleaned_all.csv"), row.names = FALSE)

# see the data you have on disk so far:
list.files(occurrences_dir)

# DOWNLOAD ENVIRONMENTAL VARIABLES ####

# we'll use functions of the 'sdmpredictors' package to access different online datasets
pred_datasets <- list_datasets(terrestrial = TRUE, marine = TRUE)
pred_datasets
names(pred_datasets)
pred_datasets[ , 1:4]  # these are the datasets currently available for download using the 'sdmpredictors' package; you can check their URLs for more info on their contents
pred_datasets[ , c("dataset_code", "citation")]  # remember to ALWAYS cite the actual data sources, not just the package you used for downloading!

pred_layers <- list_layers(datasets = pred_datasets)
unique(pred_layers$dataset_code)
unique(pred_layers[pred_layers$dataset_code == "WorldClim", ]$name)  # example of terrestrial variables dataset
unique(pred_layers[pred_layers$dataset_code == "MARSPEC", ]$name)  # example of marine variables dataset

# let's choose one dataset (e.g. WorldClim) and one particular set of variables (e.g. altitude and the bioclimatic ones, which are in rows 1 to 20):
layers_choice <- unique(pred_layers[pred_layers$dataset_code == "WorldClim", c("name", "layer_code")])
layers_choice
layers_choice <- layers_choice[1:20, ]
layers_choice

# define folder for downloading / fetching the variables' map layers:
options(sdmpredictors_datadir = "./outputs/sdmpredictors")
# load the layers to the current R session (downloading them if they aren't already in the folder defined above):
layers <- load_layers(layers_choice$layer_code, rasterstack = FALSE)  # rasterstack=TRUE gives error when there are layers with different extent
layers  # a list of raster maps
# see how many elements in 'layers':
length(layers)

# plot a couple of layers to see how they look:
names(layers)
# convert each layer to 'SpatRaster' class (from package 'terra'), which is much faster to process
layers <- lapply(layers, rast)
plot(layers[[1]], main = names(layers)[1])
plot(layers[[5]], main = names(layers)[5])

# find out if your layers have different extents or resolutions:
unique(pred_layers[pred_layers$dataset_code == "WorldClim", ]$cellsize_lonlat)  # in this case 0.08333333 - spatial resolution can then be coarsened as adequate for your species data and modelling region (see below)
sapply(layers, ext)  # if you get different extents (which doesn't happen with WorldClim, but may happen with other datasets), you'll have to crop all layers to the minimum common extent before proceeding

# for example, if the first layer has the smallest extent:
#layers <- lapply(layers, crop, extent(layers[[1]]))

# once all layers have the same extent and resolution, you can combine them in a single multi-layer raster map:
layers <- rast(layers)
layers
plot(layers)

# DELIMIT THE MODELLING REGION ####

# convert species occurrences table to a spatial object (like when you import a delimited text file into a GIS, you need to say which columns contain the spatial coordinates (point geometry) and what is the projection / coordinate reference system):
occurrence_points <- vect(occurrences, geom = c("decimalLongitude", "decimalLatitude"), crs = "+proj=longlat")
plot(occurrence_points, col = "darkblue")

# add countries map:
plot(countries, border = "tan", lwd = 2, add = TRUE)


# select the modelling region, e.g. using the countries where this species has occurrence points (which means the species was surveyed in those countries):

countries_with_points <- countries[occurrence_points, ]
plot(countries_with_points, border = "tan")
plot(occurrence_points, col = "darkblue", cex = 0.2, add = TRUE)
# judge if some countries are visibly insufficiently surveyed (in this dataset) for this species; compare with occurrence data from other sources, e.g. https://www.iucnredlist.org, national atlases, data papers

# for the example species, some european countries seem to be the only evenly surveyed countries in this GBIF dataset
# also, only the mainland should be included in the modelling region, as on islands species may be limited by factors other than climate variables (e.g. dispersal)
# so, assuming we can't complement the dataset with occurrences from other sources, let's select only the mainland for modelling:

# create a unique polygon identifier for 'countries_with_points' and add it as labels to the map, to see which polygon(s) we want to select:
countries_with_points$my_id <- 1:length(countries_with_points)
text(countries_with_points, labels = countries_with_points$my_id, col = "red", font = 2, halo = TRUE)
# select only the desired polygon(s) for the modelling region (polygons 1 and 2 FOR THE EXAMPLE DATA - CHANGE AS APPROPRIATE!!!):
selected_polygons <- c(2:9, 12:20)  # c(4, 12, 17)
mod_region <- subset(countries_with_points, countries_with_points$my_id %in% selected_polygons)
plot(mod_region, border = "green", lwd = 2, add = TRUE)  # check if desired polygon(s) outlined in green

# select the points in the modelling region:
occurrence_points_sel <- occurrence_points[mod_region, ]
plot(occurrence_points_sel, cex = 0.1, col = "green", add = TRUE)#

# if you can't select evenly surveyed countries (e.g. if you're working with marine species), 
# you can delimit the modelling region as a buffer of a given distance -- 
# e.g. 1 geographic degree, or 100 km, or the mean distance among points:


##############
# Save here! #
##############

mean_dist <- mean(terra::distance(occurrence_points_sel))  # takes time if there are many points; may give "x$.self$finalize()" errors, but check if object is successfully created nonetheless:
mean_dist
#pres_buff <- aggregate(buffer(occurrence_points, width = mean_dist))
pres_buff <- aggregate(buffer(occurrence_points_sel, width = 50000))
plot(pres_buff, lwd = 2) #, xlim = c(5, 8), ylim = c(44, 46) # -14, 51, 30, 80
plot(occurrence_points, col = "darkgreen", add = TRUE)
plot(countries, border = "tan", lwd = 2, add = TRUE)
plot(mod_region, border = "green", lwd = 2, add = TRUE)


# if the buffer is to be COMBINED WITH previously selected countries / polygon(s), the modelling region should be:
mod_region <- intersect(pres_buff, mod_region)
plot(mod_region, border = "lightblue", lwd = 2, add = TRUE)

# if you don't have any previously selected polygons and want to use only the buffer, the modelling region should be:
mod_region <- pres_buff

# IF YOU USED A LIMITED WINDOW OF COORDINATES to download the occurrence data, you need to intersect or crop with that too:
mod_region <- crop(mod_region, ext(my_window))
plot(mod_region, border = "green", lwd = 3, add = TRUE)

# aggregate modelling region into a single polygon:
mod_region <- aggregate(mod_region)
plot(mod_region)

# now import and cut (mask + trim) the variable maps to the extent of the modelling region:
layers <- rast(list.files("./outputs/sdmpredictors", pattern = "\\.tif$", full.names = TRUE))
#plot(layers)
layers_cut <- mask(crop(layers, mod_region), mod_region)
#plot(layers_cut)
names(layers_cut)
plot(layers_cut[[1]])
plot(countries, border = "darkgrey", add = TRUE)
plot(mod_region, add = TRUE)
plot(occurrence_points_sel, pch = 20, cex = 0.1, add = TRUE)
# check that everything overlaps correctly!

# SET THE APPROPRIATE SPATIAL RESOLUTION ####

# closely inspect your species data vs. the size of the variables' pixels:
plot(layers_cut[[1]], xlim = c(7, 13), ylim = c(59, 63)) #-14, 51, 35, 71
points(occurrence_points, cex = 0.5)

# plot within different x/y limits if necessary to see if your presence point resolution matches pixel resolution (i.e., if you don't have many evenly spaced points with pixels in between)
# notice that, in the example data, pixels have approximately the same spatial resolution as the presence points, but they don't match exactly (i.e., points are not at the centroids of these pixels); ideally, you should find the original grid over which (most of) these presences were sampled, and extract the raster values to that grid

# if necessary, you can aggregate the layers, to e.g. a 5-times coarser resolution (choose the 'fact' value that best matches your presence data resolution to your variables' resolution):
layers_aggr <- aggregate(layers_cut, fact = 5, fun = mean)
res(layers_aggr)
plot(layers_aggr[[1]], xlim = c(7, 13), ylim = c(59, 63))
points(occurrence_points, pch = ".")

# run the command below only if you did need to aggregate the layers:
#layers_cut <- layers_aggr


# save the cut layers to a folder on disk (we'll need them later!):
outdir_layers_cut <- paste0("./outputs/sdmpredictors/layers_cut_", myspecies)
if (!file.exists(outdir_layers_cut)) dir.create(outdir_layers_cut)
terra::writeRaster(layers_cut, filename = paste0(outdir_layers_cut, "/layers_cut.tif"), overwrite=T)

# now make a dataframe of the species occurrence data gridded to the resolution of the raster variables
# i.e., one row per pixel with the values of the variables and the presence/absence of species records:

##############
# Save here! #
##############


head(occurrences)
?gridRecords
gridded_data <- gridRecords(rst = layers_cut, pres.coords = occurrences[ , c("decimalLongitude", "decimalLatitude")])
head(gridded_data)

nrow(gridded_data)  # should be the same number as:
sum(!is.na(values(layers_cut[[5]])))
names(gridded_data)
myspecies

# plot the gridded records:
plot(layers_cut[[1]])
# plot the absences (pixels without presence records):
points(gridded_data[gridded_data[ , "presence"] == 0, c("x", "y")], col = "red", cex = 0.1)
# plot the presences (pixels with presence records):
points(gridded_data[gridded_data[ , "presence"] == 1, c("x", "y")], col = "blue", cex = 0.2)

# plot within a narrower coordinate range to see closer:
plot(occurrence_points, xlim = c(7, 13), ylim = c(59, 63), add = TRUE)  # see the coordinate range on the plot axes, and pick some limits within which to look closer
plot(layers_cut[[1]], xlim = c(7, 13), ylim = c(59, 63))
points(gridded_data[gridded_data[ , "presence"] == 0, c("x", "y")], col = "red", pch = 1, cex = 0.5)
points(gridded_data[gridded_data[ , "presence"] == 1, c("x", "y")], col = "blue", pch = 20, cex = 0.5)

# save the modelling dataframe to a .csv file on disk:
write.csv(gridded_data, paste0("./outputs/dat_gridded_", myspecies, ".csv"), row.names = FALSE)


