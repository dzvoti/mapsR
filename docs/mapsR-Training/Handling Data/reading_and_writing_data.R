##########################################################################
##########################################################################
#
# Basic R training
#
# 2.  Reading and writing data in R.
#
##########################################################################
##########################################################################
#
# The objective of this R script is to provide information on the topic under
# consideration, along with examples and exercises.  You should be able to
# work through it in R or R studio.  
#
# This script is concerned with how we get numbers into R, and then write 
# numbers out again.  By the end you should have a better understanding
# of the R scripts, and should be better-placed to start developing and editing scripts yourself. 
# The particular topics we shall cover are:
#
# 1. Reading simple data files into R
# 2. Writing data out of R
# 3. Reading data from Excel files
#
#
# In most applications of R we need to start by reading in data from a file
# which has been created elsewhere, perhaps using an editor, a spreadsheet
# or as output from a particular piece of equipment or software. In this 
# script we explore some options for doing this.

####################################################################################################

# SETTING A WORKING DIRECTORY

# Before you read data into R, you must ensure that R knows where to look
# for the data on your computer.  It is possible to give a full path to
# the file, but is usually most convenient to set the folder with the data
# as the working directory.  This can be done using the "Change directory.."
# option from the "Session" drop down menu when the R console is active.
# Or from the "File browser" three dots ... (click to navigate to folder), 
# and then "More" icon drop down to set working directory
# Alternatively you can use the "setwd" (set working directory) command as
# follows:

# Three ways of setting a working directory:

# 1.setwd: you need to know the working directory path, which I find complicated

setwd("C:/Users/sbzhp/OneDrive - The University of Nottingham/Documents/Basic R and Summary Statistics/Intro_to_R/Tanz_Training")

# 2. Use the Session tool on the Tool bar:
# Session etc.

# 3. Use the File Browser of RStudio:
# Navigate to the folder you want to set as working directory
# Then "click" on the "More" tab (has a cog/wheel next to it)
#Then "Click" on "Set as Working Directory"

#You can derive the full path to the file/working directory using "getwd" function:

getwd()

# Note that this is a path to the directory on my computer
# Edit it to give the correct path to the directory where you are
# working.  Also note that R uses forward slashes, /, in directory paths
# rather than the more common backslash.

##########################################################################
#
# 2.1  Reading simple data files into R
#
##########################################################################
#
# The read.table function is one of the most commonly used functions for 
# reading data. It can be used for reading data in simple ascii formats such
# as .txt files, or .dat files produced by many applications.  If you want to 
# read in data from a .csv file (comma-separated values), as is commonly output
# from an Excel spreadsheet, then you can use the variant read.csv command.
#
# read.table/read.csv assume that the data are organized in columns in the 
# file, one column for each variable, and one row corresponds to a single
# observation (e.g. an experimental plot or a soil core) for which we may have
# several factor labels, and several variables measured. 
#
#  read.table/read.csv has a few important arguments:

# file: the name of the file which you want to read.  At its simplest this 
# is just a file name, in double quotes, ("data.dat", for example), when this
# file is present in the working directory but there are more complex options.
# See ?read.table for details

# header: this is a logical argument so it takes value T (true) if the first
# row of the file is a header with names of the variables in columns, and F
# otherwise.

# These two arguments are often all that you will need, for small files.
# There are some other useful ones, for example:

# skip: the number of lines to skip from the beginning (e.g if you don't
# want to look at the first 100 lines of the file.

# stringsAsFactors: should character variables be coded as factors?  This
# will ensure that R treats any variables which have letters in them (e.g. "A1")
# as factors (i.e. labels for levels of a categorical variable such as
# variety or cultivation method).  Recent releases of R do not do this by
# default, so this argument can be very useful.  If you wanted all character
# variables to be treated as factors then include stringsAsFactors=T in the 
# command.


# We use read.table to create a data frame from the contents of a file 
# Variety_yileds.txt.

data.df <- read.delim("Variety_yields.txt",header=T,stringsAsFactors=T)


# We use read.table to create a data frame from the contents of a file 
# Cashmore_soil.dat.

data.df <- read.table("Cashmore_sol.dat",header=T,stringsAsFactors=T)

###########

##ERROR ALERT: matching names in editor and directory!!

##########

# We use read.csv to create a data frame from the contents of a file 
# yielddata.csv.

data.df <- read.csv("yielddata.csv",header=T,stringsAsFactors=T)

# The data frame is now ready for use.  It can be helpful to examine the data frame
# with various R tools before proceeding to any analyses.  For example, 

str(data.df) 

# str function describes the structure of the data frame.  It tells you the number of rows 
# (observations) in the data.frame and the number of columns (variables) that it contains.
# It also tells you the kind of data each variable comprises (recall that a data frame can
# contain different data structures).
#
# This example data frame contains data from a field experiment 
# It was a Split-split design carried out in 2016/2017 season, with 
# pd (planting date) as main plot; 
# variety (maize) as sub plot; 
# nrate (nitrogen fertilizer rate) as sub-sub plot; 
# x_100_grain_wt (100 seed grain weight)
#
# The following commands can also be useful:
# 
names(data.df)	# displays the names of the variables in the dataframe
head(data.df) 	# displays the first 6 rows of the dataframe
tail(data.df) 	# displays last 6 rows of the dataframe
head(data.df,10) 	# displays 10 rows of the dataframe
nrow(data.df) 	# displays the number of rows
ncol(data.df) 	# displays number of columns

#
# When you are happy with the content and structure of the data frame, then you can
# access variables from within it with R commands.  A very simple way to access a 
# variable in a dataframe is by the "dollar notation", given the name of the variable.  


# For example, to produce a set of summary plots of the 100-grain weights in 
# the data frame created above, using the "summaplot" command from CEPHaStat
# we do the following:
# load required library

source("CEPHaStat_2.R")

# Note: CEPHaStat is a script of functions prepared beforehand.
# The script is saved in the working directory from where it is "sourced"

summaplot(data.df$x_100_grain_wt)

# note that you refer to the vector of data using (i) the name of the dataframe, 
# (ii) a $ sign and (iii) the variable name, with no spaces.
#

# You can create a new variable within a data frame using the <- assign operator. 
# For example, to create a new variable which is the hundred grain weight in mg, 
# you just multiply the original values (in grams) by 1000

data.df$x_100_grain_wt_mg <- data.df$x_100_grain_wt*1000

#  you can now see that the data frame contains an additional variable

head(data.df)

#
###################################################################################
#
# EXERCISE 2.1
#
# Use read.csv to read the soil data in the file "Cashmore_soil.csv" into an R
# dataframe, making sure that character variables are read in as factors.
#

# The variables in the data file are:
#
#  GWC_T, GWC_S   gravimetric water content of the topsoil (_T) and subsoil.
#  pH_T,  pH_S   pH of the topsoil and subsoil
#  OM_T,  OM_S   organic matter content of the topsoil and subsoil
#  Soil_Series   the soil series (soil class)

# Examine the structure of the data set
str(data.df)
# Make summary plots of the continuous soil variables
# Make a new variable in the data frame which is the log of the gravimetric
# water content of the subsoil, and make a summary plot of it. Recall that
# the R function "log" will return the natural logarithm thus y<-log(x).


#####################################################################################

# Homework: The data are also provided in an ascii file "Cashmore_soil.dat".  Use the
# read.table command to convince yourself that these are the same data!
#
##########################################################################
#
#
#
#########################################################################
#
# 2.2  Writing data out of R
#
##########################################################################
#
#  After manipulating data in a data frame we might want to save the output, either
# to be used with other software, or to save for archiving or further use in R. 
# To produce a .txt file as output you can use the write.table command.

# First, we make a new data frame from random variables

Variety <- (rep(c("A","B","C"),each=10))
example.df <- data.frame(Variety,stringsAsFactors=T)# converting varieties to factors
example.df$Yield <- rnorm(30,10,5)# number of obsv, mean, sd
head(example.df)


#  First try the following

write.table(example.df,"Variety_yields.txt")

# Go to your working directory and look at the resulting file.  The first few
# lines will be something like this

# "Variety" "Yield"
# "1" "A" 14.9320801147097
# "2" "A" 13.2040635718118
# "3" "A" 17.6780697654481
# "4" "A" 10.3252632315969
# "5" "A" 7.79996535504178
# "6" "A" 18.0229009561374
# "7" "A" 10.2652234419097
# "8" "A" 4.44808130068076
# "9" "A" 6.62607203181109
# "10" "A" 9.59345532285773
# "11" "B" 13.4659787880458

# There are some features we might not like, first, the row numbers, which will often
# be a nuisance when reading the file back into R or other software.  Second, the 
# quotation marks around the variable names.  We can fix these as follows

write.table(example.df,file="Variety_yields.txt",row.names=F,quote=F)

# The first few lines of the file now look like this

# Variety Yield
# A 14.9320801147097
# A 13.2040635718118
# A 17.6780697654481
# A 10.3252632315969
# A 7.79996535504178
# A 18.0229009561374
# A 10.2652234419097
# A 4.44808130068076
# A 6.62607203181109
# A 9.59345532285773
# B 13.4659787880458


#  You can use the write.csv command to make your output a .csv file.

write.csv(example.df,file="Variety_yields.csv",row.names=F,quote=F)

# WRITING OUT A FILE WHILE ONE WITH THE SAME NAME IS OPEN IN THE WORKING DIRECTORY:

# If you write out a file with the same name as one that is open in the working directory
# you will get an error message. To try it:
# Go to your working directory and open the csv file "Variety_yields" you have just written
# Now, rewrite the the file and see what happens

write.csv(example.df,file="Variety_yields.csv",row.names=F,quote=F)


# OVERWRITING A FILE:
# HAPPENS WHEN YOU WRITE OUT A FILE USING THE SAME NAME AS ONE THAT IS CLOSED IN THE WORKING DIRECTORY:

# write out a file to the working directory. Do not open it. Or close it if open.
# now, write out another file with the same name to the same working directory 
# How many files of the same name appear in the working directory?
# Only one because the other has been overwritten

write.csv(example.df,file="Variety_yields.csv",row.names=F,quote=F)

# you can use help("write.csv") or ?write.csv to see the options provided by write.csv
#
###################################################################################
#
# EXERCISE 2.2
#
# Go back to the data frame you made with the soil data in the file 
# "Cashmore_soil.csv".
#
# The variables in the data file are:
#
#  GWC_T, GWC_S  gravimetric water content of the topsoil (_T) and subsoil.
#  pH_T,  pH_S   pH of the topsoil and subsoil
#  OM_T,  OM_S   organic matter content of the topsoil and subsoil
#  Soil_Series   the soil series (soil class)

# As before, add a variable to the dataframe which is the log of subsoil GWC.
# Then write this supplemented data file to a new .txt format file, and to 
# a new .csv file.  

##########################################################################
#
# 2.3 Reading data from Excel files. Homework
#
##########################################################################

# Often your data may be in an Excel file, and it is possible to read
# data directly from such a file into a data frame.  This requires a
# suitable package.  The readxl pacakge is particularly useful and can
# be used to read .xls or .xlsx .  You should install it on your computer
# using the command install.packages("readxl")
#
# When readxl is installed (package: install.readxl), then you only need to load it in future R
# sessions, using the command:

install.packages("readxl")
library(readxl) 

# use
help("readxl") # you can also read the R Documentation for readxl or Help Pages

# or 
?readxl # you can also read the R Documentation for readxl or Help Pages

# The readxl package makes it easy to get data out of Excel and into R. 
# Compared to many of the existing packages (e.g. gdata, xlsx, xlsReadWrite) 
# readxl has no external dependencies, so it's easy to install and use on all 
# operating systems. It is designed to work with tabular data.

#  You are provided with an Excel file called soil.xlsx.  This contains
# two sheets.  The first contains some data on soil clay content and on soil
# organic carbon content.  The second sheet contains some data on soil pH
# measured on the same sample in water or in calcium chloride.
#
# The excel_sheets() command in readxl will tell you the names of the 
# sheets in an Excel file available in your working directory:

excel_sheets("soil.xlsx") # this shows the sheet names within the "soil.xlsx" file

# Using this information, you can then read a particular sheet from
# the file into a data frame in R.  You can do this either by using the name
# of the sheet, or the number in the sequence of names:

data.df<-read_excel("soil.xlsx", sheet = "clay_SOC")
head(data.df)

# ... or, equivalently

data.df<-read_excel("soil.xlsx", sheet = 1)
head(data.df)

# The data are now in the dataframe, you can examine its contents, for example,
# with

str(data.df)

# and you can reference the variables using the $ notation, e.g.

plot(data.df$clay,data.df$SOC,xlab="Clay content /%",
ylab="Soil organic carbon /%")

# There are various ways to control which cells are read. 
# See the examples below

#... by specifying the maximum number of rows to read in:

data.df<-read_excel("soil.xlsx", n_max = 3) # this shows the specific number of rows using n_max
print(data.df)

# ... by specifying a range in the sheet, using the normal Excel
# notation (note that row 1 contains the variable names)

data.df<-read_excel("soil.xlsx", sheet=1, range = "A1:B8")
print(data.df)

#... by specifying a number of rows:

data.df<-read_excel("soil.xlsx", sheet=1, range = cell_rows(1:4))
print(data.df)

#... by specifying the cell columns
data.df<-read_excel("soil.xlsx", sheet=1,  range = cell_cols("A"))
print(data.df)

##########################################################################
#
#  EXERCISE 2.3: Homework
# 
#
# Use read.xl functions to read in the first 100 observations on soil pH in 
# soil.xlsx.  Then produce a plot of soil pH measured in water against soil 
# pH measured in calcium chloride. 
#
excel_sheets("soil.xlsx")
data.df<-read_excel("soil.xlsx", sheet = "pH", range = cell_rows(1:100))
head(data.df)
plot(data.df$pHw,data.df$pHCaCl, xlab ="Ph water", ylab="pHCalc")
#  EXERCISE 2.4
#
#  You have been provided with an Excel file called liempe_climatedata.xls
#  Use the appropriate command from read.xl to list the sheets that it contains
#  Next, read the first 360 rows of the monthly weather data into a data frame
#  and plot a graph of mean monthly temperature against rainfall.

#######################################################################################################
