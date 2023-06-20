# My first steps with R

# data types:

# integers: 1, 2, 3, 4, ...
# basically all natural numbers

# two ways to asign a variable:
a = 1
# or
a <- 1

is.integer(a) # function to check if "a" is an integer (mode)
a <- as.integer(a) # formating the content of variable "a" to the type integer (class)
is.integer(a) # check again, to confirm our formation process

a = 1L # direct asignment of a variable with an integer value
is.integer(a)

# double: 0.01, 2.43, 43.002, ...
# a.k.a. float 
# as soon as a comma is involved in a number it is a float/double
# Does not mean the system recognize that! You need to declare and check it!

b = 2.3
is.double(b)
# is.numeric() works here too, but technically those two are differnt functions 
# working on different levels (S3 vs S4 methods)

# both types (integer and double) belong to the class and mode numeric 
# but if you test it for its class different results will come out!

class(a)
.class2(a)
mode(a)

class(b)
.class2(b)
mode(b)

# logic (boolean): TRUE (T), FALSE (F)

c = TRUE 

is.logical(c)

# characters (string): "abc defg"

d = "Hello World!"

is.character(d)


#### Level 2 - higher types or objects

# Lists: contain elements of different types

list_a <- list(12, 13, 14)

class(list_a)
mode(list_a)

# is.list()
# as.list()

# Call single elements in the list with [] and the respective index number
# R starts to count with 1! 
# Untypical usually most of the other program language start with 0

list_a[1]

list_a[1:2]

list_a[-1]

# Lists are quite powerful, often used and without any further packages 
# the fastest way to handle big data in loops or other complex processes

# You can put nearly everthing in to a list and even mix them together. 

list_mix <- list(1, "Hello", TRUE, c(1.2 ,2.34,3.2134123))

# Vectors: contains element of the same type

vec_a <- c(1, 2, 3, 4, 5)
vec_b <- c(1:5)

vec_c <- c(seq(from = 1, to = 5, by = 1)) # seq() is function to generate a sequence
vec_c <- c(seq(1, 5 , 1)) 

vec_d <- c("Hello", "World", "!")

is.vector(vec_d)

# useful function for vectors:
length(vec_d) # works for lists as well

# Factors: a data structure with predefined, finite number of values

fac_a <- factor(c("A", "B", "B", "C", "A"))
fac_b <- factor(c("A", "B", "B", "C", "A"), levels = c("A", "B", "C", "D"))

# useful function for vectors:
levels(fac_a)

# Matrices: a two dimensional data structure

mat_a <- matrix(1:9, nrow = 3, ncol = 3, byrow = T)
is.matrix(mat_a)

mat_zeros <- matrix(0, 3, 3) # zeros 3x3 matrix
mat_ones <- matrix(1, 3, 3) # ones 3x3 matrix

# call a certain cell:
mat_a[2,3] # mat_a[x,y] x=rows, y=columns

mat_a[,2] # call the entirety of column 2
mat_a[2,] # call the entirety of row 2

# useful function for matrices:
attributes(mat_a)
dim(mat_a)

# rename columns and rows:
dimnames(mat_a) = list(c("row1", "row2", "row2"), # row names 
                       c("col1", "col2", "col3")) # column names 

mat_at <- t(mat_a) # transpose the matrix
diag(mat_a) # show only the diagonal elements
diag(diag(mat_a)) # show only the diagonal elements in a matrix form


# Data frames: THE representative for tables (there are more but this is most often used)

data <- data.frame(x=1, y= 1:10)

data
View(data) # alternative way to see the data (for huge data frames not recommended)
head(data) # only the first 6 rows

# useful function for data frames:

attributes(data)
# too many to call them all! I will show them to you along the course.

# in most cases you don't create a data frame from scratch
# you load in a file like an excel or .csv file directly as a data frame
# the tricky part here is, how to load this data in and turn it correctly 
# in to your desired table?

# lets start by finding our file to import:

getwd() # helps you to find where you are in the system --> wd means "working directory"

list.dirs() # give you a list of all folders that are in the given path
# an empty field here means it looks in the working directory, same goes with the 
# next function

list.files("../") # presents you all files in the given path

list.dirs(, recursive = F) # "recursive" a boolean parameter restricts the output 
# to only the files/folders on the current given path level (F/FALSE) or to all 
# files/folders on and under! the given path level (T/TRUE)

list.files(, pattern = "next", full.names = T) # "pattern" is a string parameter that let you restrict the 
# output to only certain characters you look for

# If you need to search in an upper level of directory just use in the path parameter "../"
# one dot ("./") is still in the working directory
# two dots ("../") is on level above the working directory

# by setting the working directory to a known place helps to find all the files later on!
# But sometimes the place you want to put everything in doesn't exist yet

dir.create("./Seminar_test") # creates a folder, location and name of the folder is 
# given in the path parameter

setwd("./Seminar_test") # now let the new created folder be our new working directory

getwd() # test if it was successfull

# Alright everything is set up to receive our first data
# There are as always mutliple ways to receive data (raw, zip, xlsx, csv, txt, ...)
# I will show you two ways:
#  1) you have direct access to it via usb stick, on your computer, or external 
#     storage device
#  2) online, here we try via my github folder


# 1) direct access:

data01 <- read.csv2("../nextcloud/5405.AG_Schnittler/Lectures/Lecture material/Seminar Plant Conversation/2023/1st Session/Jonas20230605.csv")
# everthing is squashed up, what went wrong?

# We changed no parameters except the "file" one, so lets see what default values are on them!
# Select the function and press "F1"
# In most cases you need to play around with "sep=", "header=", "dec=" and "skip="

# sep --> describes the pattern your data is separated by
# header --> uses the first row to set the column names
# dec --> to determine whitch sign is used as the decimal sign
# skip --> helps you to jump/skip over the first rows 

data01 <- read.csv2("../nextcloud/5405.AG_Schnittler/Lectures/Lecture material/Seminar Plant Conversation/2023/1st Session/Jonas20230605.csv",
                    sep = ",", header = F)


# 2) online access:

download.file(url = "https://github.com/JanWoyzi/SDM-course/archive/refs/heads/main.zip",
              destfile = "Session01.zip")

unzip("Session01.zip") # be carefull, per default it will overwrite

data01 <- read.csv2("",
                    sep = ",", header = F)


# Here the data is conveniently all as csv-files but what if it is a xlsx-file?
# Without any additional help base R wont come here far.
# So lets start to "upgrade" R so that it is possible to handle those files as well.
# We do this by installing various packages from an online source. In most cases 
# we don't need to specificity the source.
# When we activate multiple packages the order of activating them does matter! Sometimes 
# different packages uses the same names for their functions (but they do something different!).
# R will inform you during package activation which functions get masked due some naming conflicts 
# with other packages. The mos recent activated package overwrites older ones.
# There is a way to use those masked function anyway by putting the name of the package in front 
# of the function separated only by "::" (see below for example)

install.packages("readxl") # download and install --> needs to be done only once, afterwards it can be commented
library(readxl) # activate --> needs to be done once every session

data02 <- read_xlsx("", skip = 0)
data02 <- readxl::read_xlsx("", skip = 0)

# Okay we have our data and now we can do some analyzing and ploting steps.
# for this lets load in more useful packages:

library(purrr) # reduce()
library(readr) # read_csv()
library(data.table) # plenty of functions for merging and quite effective with big data
library(tidyverse) # useful function for handling data frames

# And the challenge this time is that you have multiple files and need to merge them together.
dirs <- list.dirs(".", recursive = T, full.names = T)

files <- list.files("./Session01/Data01", recursive = T, full.names = T) # we take the result 
# of all files (with full path names) to a variable...
files <- grep("*.csv", files, value = T) # ...and take only (like a filter) files who ends with .csv

data <- files %>%
  map(read.csv2, sep=";") %>%    # read in all the files at once, using
  reduce(rbind)                  # the function read_csv() from the readr package
                                 # reduce with rbind() into one data frame

# more functions to combine data frames
rbind() # Combine R Objects by Rows
cbind() # Combine R Objects by Columns
merge() # Merge two data frames by common columns or row names, or do other versions of database join operations.

# Lets assume we need to change certain values column-wise, maybe as some kind of transformation or 
# scale factor. For this there are many ways but so far I will show you my favorite ones 
# (because they are easy to follow)


# 1) creating a new column based on the old ones and leave the old ones (with piping):
data <- data %>%
  mutate(
    Sepal.Length.trans = as.numeric(data$Sepal.Length) / 0.5*0.12, .keep = "all")

# Manipulate value in column A by certain condition on column B
data[data$Collector == "Ahmed Adoudi",]$Sepal.Length.trans <- data[data$Collector == "Ahmed Adoudi",]$Sepal.Length.trans+5

#  alternative with piping:
data <- data %>%
  dplyr::mutate(Sepal.Length.trans = ifelse(Collector == "Ahmed Adoudi", Sepal.Length.trans-5, Sepal.Length.trans))

# Okay maybe you need to create a total new column:
# with empty values:

data$Processed <- NA

# with specific values
data$Processed <- TRUE # "ABC", 1, 2.3...whatever

# Or delete an entire column
data$Processed <- NULL

# Filter rows based on a condition
filtered_data <- data %>%
  filter(Petal.Width > 0.6)

# Summarize data within each group
summary_data <- data %>%
  group_by(Species) %>%
  summarize(mean_value = mean(Sepal.Length),
            max_value = max(Sepal.Length))

# Sort the data by a variable
sorted_data <- data %>%
  arrange(Collector)

# Recode categorical variables
recoded_data <- data %>%
  mutate(Species = recode(Species, "setosa" = "Iris setosa", "virginica" = "Iris virginica", "versicolor" = "Iris versicolor"))

# If you are done you should save your data frame!
# save csv or excel files:

write.csv2(data, file ="./Session01/Data01/iris_all.csv")

library(xlsx)
xlsx::write.xlsx(data, file ="./Session01/Data01/iris_all.xlsx")


# Okay lets do graphics! with ggplot2

# Load the necessary libraries
library(ggplot2) # most often used for this job. Of course there are more out there but ggplot2 covers the most
# cases and only special cases needs special packages for generating special graphics (animation, networks, phylogeny trees) 

# View the structure of the data
str(data)

# Usually the task here is to find all parameters that will be accepted in journals. 
# For instance a complete white background without lines, black axis, larger fonts etc.

# Scatter plot
ggplot(data, aes(x = Sepal.Width, y = Sepal.Length)) +
  geom_point()

# Line plot
ggplot(data, aes(x = X, y = Sepal.Width)) +
  geom_line()

# Bar plot
ggplot(data, aes(x = Species, y = Petal.Length)) +
  geom_bar(stat = "identity")

# Histogram
ggplot(data, aes(x = Petal.Length)) +
  geom_histogram(binwidth = 0.5)

# Box plot
ggplot(data, aes(x = Species, y = Petal.Length)) +
  geom_boxplot()

# Customize plot aesthetics
ggplot(data, aes(x = Sepal.Width, y = Sepal.Length, color = Species)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "Title of the Plot", x = "X Axis", y = "Y Axis")

# Save the plot as an image file
ggsave("plot.png") # in this way it will saved in your working directory

# Show the plot if you work in R not on an IDE (integrated development environment, like RStudio)
print(ggplot(data, aes(x = Sepal.Width, y = Sepal.Length, color = Species)) +
        geom_point(size = 3) +
        theme_bw() +
        labs(title = "Title of the Plot", x = "X Axis", y = "Y Axis"))

data %>%  
  ggplot(aes(x= Species, y=Sepal.Width)) + 
  geom_boxplot(fill="gray") +
  stat_summary(fun=mean, geom="point", shape=23, size=4) +
  labs(title="Overview - Iris dataset",x="Species", y = "Sepal width (mm)") +
  theme_classic()

# Now we can as well plot some maps with ggplot2:
# here we need some more packages:
install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))

library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")

library(sf)

ggplot(data = world) +
  geom_sf() +
  coord_sf(expand = FALSE) + 
  # geom_sf_label(aes(label = name)) + # Chaos!
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("World map", subtitle = paste0("(", length(unique(world$name)), " countries)"))

# color country-wise the population:
ggplot(data = world) +
  geom_sf(aes(fill = pop_est)) +
  scale_fill_viridis_c(option = "plasma", trans = "sqrt")

# a play with different European Petroleum Survey Group (EPSG) projections:

ggplot(data = world) +
  geom_sf() +
  coord_sf(crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ")

ggplot(data = world) +
  geom_sf() +
  coord_sf(crs = st_crs(3035))
 
ggplot(data = world) +
  geom_sf() +
  coord_sf(crs = st_crs(4326))

ggplot(data = world) +
  geom_sf() +
  coord_sf(crs = st_crs(3857))

# show only a smaller window of the map

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-15.15, 45.12), ylim = c(30.65, 75), expand = FALSE)

#Okay enough of this, more in the next session! Lets focus on loops, if/else clauses and functions.

# Introduction to Loops and ifelse Clauses

# Looping with for loop
for (i in 1:5) {
  print(i)
}

# Looping with while loop
counter <- 1
while (counter <= 5) {
  print(counter)
  counter <- counter + 1
}

# Looping over a vector using a for loop
vector <- c("a", "b", "c", "d", "e")
for (element in vector) {
  print(element)
}

# Looping over a data frame using a for loop
data02 <- data.frame(x = 1:5, y = 6:10)
for (i in 1:nrow(data02)) {
  print(paste("x:", data02[i, "x"], "y:", data02[i, "y"]))
}

# If cases with multiple conditions
z <- 10
if (z > 0) {
  print("z is positive")
} else if (z == 0) {
  print("z is zero")
} else {
  print("z is negative")
}

# Nested if statements
a <- 7
b <- 5
if (a > b) {
  if (a %% 2 == 0) {
    print("a is greater than b and even")
  } else {
    print("a is greater than b but odd")
  }
} else {
  print("a is less than or equal to b")
}

# Using ifelse clause
x <- 10
result <- ifelse(x > 5, "Greater than 5", "Less than or equal to 5")
print(result)

# Looping with ifelse clause
vector <- c(2, 4, 6, 8, 10)
for (num in vector) {
  result <- ifelse(num > 5, "Greater than 5", "Less than or equal to 5")
  print(paste(num, ":", result))
}

# Looping with nested ifelse clauses
vector <- c(2, 4, 5, 6, 8, 10)
for (num in vector) {
  result <- ifelse(num > 5, "Greater than 5", ifelse(num == 5, "Equal to 5", "Less than 5"))
  print(paste(num, ":", result))
}

# And now the big finale: functions!

# Function to calculate average measurements for a given species in the iris dataset
calculate_species_avg <- function(species) {
  # Subset the iris dataset for the given species
  subset_data <- data[data$Species == species, ]
  
  # Calculate the average measurements
  avg_measurements <- colMeans(subset_data[, 1:4])
  
  # Print the average measurements
  cat("Average measurements for", species, ":\n")
  cat("Sepal Length:", avg_measurements[1], "\n")
  cat("Sepal Width:", avg_measurements[2], "\n")
  cat("Petal Length:", avg_measurements[3], "\n")
  cat("Petal Width:", avg_measurements[4], "\n")
}

# Finally:
calculate_species_avg("setosa")



#Function to convert utm to lon lat
library(rgdal)

convert_utm_to_latlon <- function(easting, northing, zone, hemisphere) {
  # Create a spatial reference object for the UTM zone
  utm_proj <- CRS(paste0("+proj=utm +zone=", zone, " +", hemisphere, " +ellps=WGS84"))
  
  # Create a data frame with UTM coordinates
  utm_data <- data.frame(Easting = easting, Northing = northing)
  
  # Convert UTM coordinates to latitude and longitude
  latlon_data <- spTransform(SpatialPoints(utm_data, proj4string = utm_proj), CRS("+proj=longlat +datum=WGS84"))
  
  # Extract the latitude and longitude values
  latitude <- slot(latlon_data, "coords")[, "Northing"]
  longitude <- slot(latlon_data, "coords")[, "Easting"]
  
  # Return the latitude and longitude values as a data frame
  data.frame(Longitude = longitude, Latitude = latitude)
}

# Usage example:
easting <- c(500000, 600000, 700000)
northing <- c(4000000, 4100000, 4200000)
zone <- 32
hemisphere <- "N"

result <- convert_utm_to_latlon(easting, northing, zone, hemisphere)
print(result)


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
