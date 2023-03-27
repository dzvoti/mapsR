##########################################################################
##########################################################################
#
# CEPHaS project R training.  Going back to basics
#
#
##########################################################################
##########################################################################

# 1. Simple subsetting of data frames

##########################################################################
##########################################################################
#
# # Subsetting
# There are a number of operators that can be used to extract subsets of R objects.
# [ always returns an object of the same class as the original; can be used to select more than one element (there is one exception)
# [[ is used to extract elements of a list or a data frame; it can only be used to extract a single element
# and the class of the returned object will not necessarily be a list or data frame
# $ is used to extract elements of a list or data frame by name; semantics are similar to that of [[.


# ########################################################################
# Data Frames
# ########################################################################
# # Data frames are used to store tabular data
# They are represented as a special type of list where every element of the 
# list has to have the same length

# Each element of the list can be thought of as a column and the length of 
# each element of the list is the number of rows

# Unlike matrices, data frames can store different classes of objects 
# in each column (just like lists); matrices must have every element be 
# the same class

# Data frames also have a special attribute called row.names

# Data frames are usually created by calling read.table() or read.csv()

# Can be converted to a matrix by calling data.matrix()

# Data Frames
x <- data.frame(foo = 1:4, bar = c(T, T, F, F))
x
nrow(x)
ncol(x)

# # Data Frame Subsetting
# If provided with a single value, data frames assume you want to subset a column 
# or columns - multiple values then the data frame is treated as a matrix.

df = data.frame(a = 1:2, b = 3:4) # create a data frame

df[1] # returns an object of the same class as the original

df[[1]] # class of the returned object will not necessarily be a list or data frame

df[, "a"] # class of the returned object will not necessarily be a list or data frame

df["a"] # returns an object of the same class as the original

df[, "a", drop = FALSE] # returns an object of the same class as the original

df[1, ] 

df[c("a","b","a")]

##########################################################################
# # Taking a Subset of a Data Frame in R
#########################################################################
# Let's load in yielddata into a data frame
# This contains the 100 grain weight of maize grain at ZARI. 
# The data are on planting dates, nitrogen fertilizer rates (n1, n2, n3) & maize cultivars (v1, v2, v3)
# Using the str command, we find that there are 81 observations in this data frame

# load the yielddata dataset

data.df <- read.csv("yielddata.csv",header=T,stringsAsFactors=T)

head(data.df) # display the first 6 rows
str(data.df)  # check the structure of the data

# Let's try to get One Column after from the data frame named yielddata loaded into R, 
# we can take subsets of these 81 observations. 
# First, let's assume we just want to pull out the column of rep, pd, variety, nrate & x_100_grain_wt. 
# There are two ways we can do this: 
# specifying the column by name, or specifying 
# the column by its order of appearance. 
# 
# The general form for pulling information from data frames is 
# data.frame[rows,columns] 
# ... so you can get the first column in either of these two ways:

data.df[ ,3]         # get all rows, but only the third column
data.df[ ,c("pd")]   # get all rows, and only the column named "pd"

# # Getting Multiple Columns.
# If you want more than one column, you can specify the column numbers or the names of the variables that you want to extract. 
# If you want to get the pd, variety, & nrate columns, you would do this:
#   

data.df[,c(3,4)]                # get all rows, but only 3rd and 4th columns
data.df[,c("pd","variety")]    # get all rows, only "pd" & "variety" columns

# If you want more than one column and those columns are next to each other, you can do this:

data.df[,c(2:6)]

# # Getting One Row
# You can get the first row similarly to how you got the first column, 
# and any other row the same way:

data.df[1,]         # get first row, and all columns
data.df[72,]        # get 72nd row, and all columns

# # Getting Multiple Rows
# If you want more than one row, you can specify the row numbers you want like this:
#   
data.df[c(1:3,5,6),]

##########################################################################
# Exercise 4.1
# Try to get row (10 to 41) and columns (2 to 5) from the yielddata dataset




##########################################################################
# Subsetting Lists: Homework
##########################################################################
x <- list(foo = 1:4, bar = 0.6)
x[1]
x[[1]]
x$bar
x[["bar"]]
x["bar"]

# Subsetting Lists
x <- list(foo = 1:4, bar = 0.6, baz = "hello")
x
x[c(1, 3)]

# # Subsetting Lists
# The [[ operator can be used with computed indices; $ can only be used with literal names.
x <- list(foo = 1:4, bar = 0.6, baz = "hello")
name <- "foo"
x[[name]] ## computed index for 'foo'
x$name ## element 'name' doesn't exist!
x$foo

# # Subsetting Nested Elements of a List
# The [[ can take an integer sequence.
x <- list(a = list(10, 12, 14), b = c(3.14, 2.81))
x
x[[c(1, 3)]]
x[[1]][[3]]
x[[c(2, 1)]]

##########################################################################
# # Subsetting a Matrix: Homework
##########################################################################
# Matrices can be subsetted in the usual way with (i,j) type indices.
# Indices can also be missing.
x <- matrix(1:6, 2, 3)
x
x[1, 2] # R reads rows first and columns second
x[2, 1] # Reads both rows, columns

# Indices can also be missing.
x[1, ]  # Reads rows only
x[, 2]  # Reads columns only

# # Subsetting a Matrix
# By default, when a single element of a matrix is retrieved, it is returned as a vector of length 1 rather than a 1 ? 1 matrix. This behavior can be turned off by setting drop = FALSE.
x <- matrix(1:6, 2, 3)
x
x[1, 2]
x[1, 2, drop = FALSE]

# # Subsetting a Matrix
# Similarly, subsetting a single column or a single row will give you a vector, not a matrix (by default).
x <- matrix(1:6, 2, 3)
x
x[1, ] # displays all first row and all columns
x[1, , drop = FALSE]

# # Partial Matching
# Partial matching of names is allowed with [[ and $.
x <- list(aardvark = 1:5)
x
x$a
x[["a"]]
x[["a", exact = FALSE]]

##########################################################################
# # Subsetting a Vector: Homework
##########################################################################
# Vectors are basic objects in R and they can be subsetted using the [ operator.
x <- c("a", "b", "c", "c", "d", "a")

x[1] ## Extract the first element
x[2] ## Extract the second element

# The [ operator can be used to extract multiple elements of a vector by passing the operator an integer sequence. Here we extract the first four elements of the vector.

x[1:4]

# The sequence does not have to be in order; you can specify any arbitrary integer vector.
x[c(1, 3, 4)]

# We can also pass a logical sequence to the [ operator to extract elements of a vector that satisfy a given condition. For example, here we want the elements of x that come lexicographically after the letter a.

u <- x > "a"
u
x[u] # this displays all the elements fulfiling the condition 

# Another, more compact, way to do this would be to skip the creation of a logical vector and just subset the vector directly with the logical expression.

x[x > "a"]


##########################################################################
# It's frequently necessary to extract some of the elements of a larger vector.
# In R you can use square brackets to select an individual element or group of elements:
x <- c(5,9,2,14,-4)
x[3]

# note indexing starts from 1
x[c(2,3,5)]
x[1:3]
x[3:length(x)]

# There are two other methods for getting subvectors. The first is using a logical vector (i.e. containing TRUE and FALSE) of the same length:
x > 4
x[x > 4]

# or using negative indices to specify which elements should not be selected:
x[-1]
x[-c(1,4)]

# (Note that this is rather different to what other languages such as C or   Python would interpret negative indices to mean.)

##########################################################################
# 
# Exercise 4.2: Homework 
# 
# The built-in vector LETTERS contains the uppercase letters of the alphabet. Produce a vector of (i) the first 12 letters; (ii) the odd 'numbered' letters; (iii) the (English) consonants.


##########################################################################
# # Removing NA Values
##########################################################################
# Here we will use a complete.cases function in R
# # Find Complete Cases
# Return a logical vector indicating which cases are complete, i.e., have no missing values or NAs.

?complete.cases # read the help pages on how to use the complete.cases function

# A common task is to remove missing values (NAs).
x <- c(1, 2, NA, 4, NA, 5) 

bad <- is.na(x)
x[!bad]

# # Removing NA Values
# What if there are multiple things and you want to take the subset with no missing values?
x <- c(1, 2, NA, 4, NA, 5)
y <- c("a", "b", NA, "d", NA, "f")
good <- complete.cases(x, y)
x[good]
y[good]

# Removing NA Values
# Let's use airquality, one of R's built in datasets.

airquality[1:6, ] # this shows the first six row of the airquality data frame
good <- complete.cases(airquality)
airquality[good, ][1:6, ] # all rows containing an NA are removed completely

#Exercise: Write a script to remove NAs from airquality using columns


###########################################################################
# # Subsetting a factor: Homework
###########################################################################
# # R factors
# Factors are the data objects which are used to categorize the data and store it as levels.
# They can store both strings and integers. 
# They are useful in the columns which have a limited number of unique values. 
# Like "Male, "Female" and True, False etc. 
# They are useful in data analysis for statistical modeling.

# Factors are created using the factor() function by taking a vector as input.

# As an example, let's create a vector as input
data <- c("East","West","East","North","North","East","West","West","West","East","North")

print(data)
print(is.factor(data))

# Apply the factor function.
factor_data <- factor(data) # we are applying the factor() function to data
# save the data into the factor_data

print(factor_data)
print(is.factor(factor_data))

# When we execute the above code, it produces the following result ???
# 
# [1] "East"  "West"  "East"  "North" "North" "East"  "West"  "West"  "West"  "East" "North"
# [1] FALSE
# [1] East  West  East  North North East  West  West  West  East  North
# Levels: East North West
# [1] TRUE

# # Factors in Data Frame
# On creating any data frame with a column of text data, R treats the text column as categorical data and creates factors on it.

# Create the vectors for data frame.
height <- c(132,151,162,139,166,147,122)
weight <- c(48,49,66,53,67,52,40)
gender <- c("male","male","female","female","male","female","male")

# Create the data frame.
input_data <- data.frame(height, weight, gender)
print(input_data)

# Test if the gender column is a factor.
print(is.factor(input_data$gender))

# Print the gender column so see the levels.
print(input_data$gender)

# When we execute the above code, it produces the following result ???

#   height weight gender
# 1    132     48   male
# 2    151     49   male
# 3    162     66 female
# 4    139     53 female
# 5    166     67   male
# 6    147     52 female
# 7    122     40   male
# [1] TRUE
# [1] male   male   female female male   female male  
# Levels: female male

# # Changing the Order of Levels
# The order of the levels in a factor can be changed by applying the factor function again with new order of the levels.
# Create the vector for data frame.

data <- c("East", "West", "East", "North", "North", "East",
          "West", "West", "West", "East", "North")

# Create the factors
factor_data <- factor(data)
print(factor_data)

# Apply the factor function with required order of the level.

new_order_data <- factor(factor_data,levels = c("East","West","North"))
print(new_order_data)

# # When we execute the above code, it produces the following result ???
# 
# [1] East  West  East  North North East  West  West  West  East  North
# Levels: East North West
# [1] East  West  East  North North East  West  West  West  East  North
# Levels: East West North

# # Generating Factor Levels
# We can generate factor levels by using the gl() function. It takes two integers as input 
# which indicates how many levels and how many times each level.

# # Syntax
# gl(n, k, labels)
# Following is the description of the parameters used ???
# n is a integer giving the number of levels.
# k is a integer giving the number of replications.
# labels is a vector of labels for the resulting factor levels.

# as an example
v <- gl(3, 4, labels = c("Chitedze", "Liempe", "Domboshava"))
print(v)

# # When we execute the above code, it produces the following result ???
# 
# [1] Chitedze   Chitedze   Chitedze   Chitedze   Liempe     Liempe     Liempe     Liempe
# [9] Domboshava Domboshava Domboshava Domboshava
# Levels: Chitedze Liempe Domboshava

# You can subset factors in a similar way that you subset vectors. 
# As usual, [ ] is the key! However, R has some interesting behavior 
# when you want to remove a factor level from your analysis.
# subsetting the first name from the above example

v[1]              # select the first name still returns all the levels
v[1, drop=TRUE]   # using the drop=TRUE only gives the desired factor & levels


v[1:7]              # select the first 2 names but still returns all the levels
v[1:7, drop=TRUE]   # using the drop=TRUE only gives the desired factor & levels


# You can subset factors in a similar way that you subset vectors. 
#As usual, [ ] is the key! However, R has some interesting behavior when you want to remove a factor level from your analysis. 

# Factor Subsetting
(x = factor(c("BS", "MS", "PhD", "MS")))

x[1:2]

# R selects the BS MS at the first and second position, but left the MS level behind. A better plan would have been to tell R to drop the PhD level entirely. To do that, add drop = TRUE

x[1:2, drop=TRUE]

###########################################################################
###########################################################################
# 5. Dealing with large data sets: Explore by yourselves!
###########################################################################
###########################################################################
# # Handling large dataset in R, especially CSV data
# R has a few packages for big data support: 
# bigmemory and ff; 
# and also some uses of parallelism to accomplish the same goal using Hadoop and MapReduce;
# For datasets with size in the range 10GB, bigmemory and ff handle themselves well;
# For larger datasets, use Hadoop.

#### Big Data
## Big Data" is a catch phrase for any dataset or data application 
## that does not fit into available RAM on one system.
## R has several a few packages for big data support like:
#   1 bigmemory
#   2 ff

#### Large Datasets
## R reads data into RAM all at once, if using the usual read.table function. 
## Objects in R live in memory entirely. 
## Keeping unnecessary data in RAM will cause R to choke eventually.
### Specifically,
#   1 on most systems it is not possible to use more than 2GB of memory.
#   2 the range of indexes that can be used is limited due to 
#     lack of a 64 bit integer data type in R and R64.
#   3 on 32 bit systems, the maximum amount of virtual memory space is 
#     limited to between 2 and 4GB.
#   4 relying on virtual memory will cause the system to grind to a halt "thrashing"

### There are two major solutions in R:
##   1 bigmemory: It is ideal for problems involving the analysis in 
##    R of manageable subsets of the data, 
##    or when an analysis is conducted mostly in C++.
##    " It's part of the \big" family, some of which we will discuss.
##  2 ff: file-based access to datasets that cannot fit in memory.
##  3 can also use databases which provide fast read/write access for piecemeal analysis.

#### The big family consists of several packages for performing tasks
### on large datasets.
##   1 bigmemory is our focus.
##   2 biganalytics provides analysis routines on big.matrix such as GLM and bigkmeans.
##   3 synchronicity adds Boost mutex functionality to R.
##   4 bigtabulate adds table and split-like support for 
##      R matrices and big.matrix memory efficiently.
##   5 bigalgebra provides BLAS and LAPACK linear algebra routines 
##      for native R matrices and big.matrix.
##   6 bigvideo provides video camera streaming via OpenCV.

#########################################################################
# Session 1
# install.packages(c("bigmemory", "biganalytics"))

library(bigmemory)
library(biganalytics)

options(bigmemory.typecast.warning=FALSE)

A <- big.matrix(5000, 5000, type="char", init=0)

# Fill the matrix by randomly picking 20% of the positions for a 1.
x <- sample(1:5000,size=5000,replace=TRUE)
x
y <- sample(1:5000,size=5000,replace=TRUE)
y

for(i in 1:5000) {
  A[x[i],y[i]] <- 1
}

# Get the location in RAM of the pointer to A.
desc <- describe(A)
desc

# Write it to disk.
dput(desc, file="A.desc")
sums <- colsum(A, 1:20)
sums

##########################################################################
# Session 2
library(bigmemory)
library(biganalytics)

# Read the pointer from disk.
desc <- dget("A.desc")
desc

# Attach to the pointer in RAM.
A <- attach.big.matrix(desc)
A

# Check our results.
sums <- colsum(A, 1:20)
sums

##########################################################################