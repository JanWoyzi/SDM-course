'''
write me a function for R that accesses an excel sheets, within the sheet is a column with utm information 
that should be used to convert those utm units to lon lat units. 
The results should be added to the existing table as a new column on the right side.
'''

library(readxl)
library(sf)

convertUTMtoLatLon <- function(file, sheet, utm_column) {
  # Read the Excel file
  data <- read_excel(file, sheet = sheet)
  
  # Extract the UTM column
  utm <- data[[utm_column]]
  
  # Convert UTM coordinates to Lat/Lon
  utm_sf <- st_as_sf(data, coords = c(utm_column), crs = 32600) # Assuming UTM zone 32N (EPSG:32632)
  latlon_sf <- st_transform(utm_sf, 4326) # Convert to WGS84 (EPSG:4326)
  
  # Add the Lat/Lon column to the existing table
  data$lat <- st_coordinates(latlon_sf)[, 2]
  data$lon <- st_coordinates(latlon_sf)[, 1]
  
  return(data)
}

'''
can you improve the function script above by adding two parameters, where i write the character names of 
the easting and northing column from the excel sheets.
'''

convertUTMtoLatLon <- function(file, sheet, utm_column, easting_column, northing_column) {
  # Read the Excel file
  data <- read_excel(file, sheet = sheet)
  
  # Extract the UTM, easting, and northing columns
  utm <- data[[utm_column]]
  easting <- data[[easting_column]]
  northing <- data[[northing_column]]
  
  # Convert UTM coordinates to Lat/Lon
  utm_sf <- st_as_sf(data, coords = c(utm_column, easting_column, northing_column), crs = 32600) # Assuming UTM zone 32N (EPSG:32632)
  latlon_sf <- st_transform(utm_sf, 4326) # Convert to WGS84 (EPSG:4326)
  
  # Add the Lat/Lon column to the existing table
  data$lat <- st_coordinates(latlon_sf)[, 2]
  data$lon <- st_coordinates(latlon_sf)[, 1]
  
  return(data)
}


'''
get rid of the utm_column parameter

'''

convertUTMtoLatLon <- function(file, sheet, easting_column, northing_column) {
  # Read the Excel file
  data <- read_excel(file, sheet = sheet)
  
  # Extract the UTM column
  utm <- data[["UTM"]]
  
  # Extract the easting and northing columns
  easting <- data[[easting_column]]
  northing <- data[[northing_column]]
  
  # Convert UTM coordinates to Lat/Lon
  utm_sf <- st_as_sf(data, coords = c("UTM", easting_column, northing_column), crs = 32600) # Assuming UTM zone 32N (EPSG:32632)
  latlon_sf <- st_transform(utm_sf, 4326) # Convert to WGS84 (EPSG:4326)
  
  # Add the Lat/Lon column to the existing table
  data$lat <- st_coordinates(latlon_sf)[, 2]
  data$lon <- st_coordinates(latlon_sf)[, 1]
  
  return(data)
}

'''
This line:  "latlon_sf <- st_transform(utm_sf, 4326) # Convert to WGS84 (EPSG:4326)" does not work. 
afterwards i have empty fields


no, this did not change the issue

'''

# Check the current CRS
st_crs(utm_sf)

# Assign the correct CRS if needed
utm_sf <- st_set_crs(utm_sf, 32632)  # Replace 32632 with the appropriate EPSG code for your UTM zone

'''
Final and working product:
'''

library(readxl)
library(sf)

list.files("../nextcloud/5405.AG_Schnittler/Lectures/Lecture material/Seminar Plant Conversation/2023/1st Session/")

convertUTMtoLatLon <- function(file, sheet, easting_column, northing_column) {
  # Read the Excel file
  data <- readxl::read_excel(file, sheet = sheet)
  
  # Extract the easting and northing columns
  easting <- data[[easting_column]]
  northing <- data[[northing_column]]
  
  # Convert UTM coordinates to Lat/Lon
  utm_sf <- st_as_sf(data, coords = c(easting_column, northing_column), crs = 32632) # Assuming UTM zone 32N (EPSG:32632)
  latlon_sf <- st_transform(utm_sf, 4326) # Convert to WGS84 (EPSG:4326)
  
  # Add the Lat/Lon column to the existing table
  data$lat <- st_coordinates(latlon_sf)[, 2]
  data$lon <- st_coordinates(latlon_sf)[, 1]
  
  return(data)
}


file <- "../nextcloud/5405.AG_Schnittler/Lectures/Lecture material/Seminar Plant Conversation/2023/1st Session/UTM_lonlat.xlsx"

convertUTMtoLatLon(file, 1, "Easting", "Northing")