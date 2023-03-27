
# Basic R training.

#########################################################################################

# PREAMBLE: OBJECTS AND OBJECT NAMES IN R

# Data are objects in R
# Data objects are assigned names
# We manipulate named objects in R

##############################################################################
##############################################################################
#
# 1.  Data in R, types and structures.
#
##############################################################################
##############################################################################
#
# The objective of this R script is to provide information on the topic under
# consideration, along with examples and exercises.  You should be able to
# work through it in R or R studio.  This particular script requires no 
# additional functions or packages to be loaded.
#
# The specific objective of the material in this script is to introduce
# you to the principal data types and data structures in R.  By the end
# you should have a better understanding of the R scripts, 
# and should be better-placed to start developing 
# and editing scripts yourself.  The particular topics we shall cover are
#
#   1. Basic data types in R: numeric, character, and logical
#   2. Data structures:  vectors.  
#   3. Factors
#   4. Data structures: Matrices, making matrices, matrix operations
#   5. Data structures: Dataframes
#   6. Data structures: Lists 
# 
##########################################################################################



#
# The principal data types in R are numeric, character, factor and logical.  
# There are others, but these are the main ones.
#
# A datum of type numeric is a numerical value, such as a soil pH value.
# A datum of type character is a string of characters, such as the name of 
# an experiment.
# A datum of type factor is the label for a set of treatments or categories
# which we might use in an analysis of variance.  
# A datum of type logical takes values TRUE or FALSE
#
#
##############################################################################
#
# 1.1 NUMERIC DATA
#

# Scalar:

# Before going into details of data types, we introduce the simplest data
# type in R, the scalar.  A scalar is a single value of some variable,
# e.g. the value 42, or the name "Bert".
#
#
# We can make a simple scalar value using the <- "assign" arrow in R
# Assignment is simply the association of a name (mass) with an object (the value 1000).

	mass<-1000

# This is a numeric scalar value.
#
# We can then use the print command to see the value of the scalar

	print(mass)

# We can then do simple mathematical manipulations with a numeric scalar value. 
# For example, the following will convert the mass, if this is in grams, to
# kg

	mass_kg<-mass/1000

# so you can see that the back-slash / denotes division of the value to the left by
# the value to the right.  As in most computer languages * denotes multiply, so
# to convert mass to milligrams we would do the following

	mass_mg<-mass*1000

#
#
# EXERCISE 1.1 Create a scalar value which is your age at last birthday then, assuming that
# there are 52 weeks in each year, calculate your age in weeks and assign this value
# to a scalar called my_age


# Here are some simple operations which we can perform on a numeric scalar in R using
# inbuilt functions
#
# Square root
# 
	sroot_mass<-sqrt(mass)
#
# Logarithm (natural or Napierian logarithm to base e)
#
	log_mass<-log(mass)
#
# Raising to a power
	mass_squared<-mass^2
#
################################################################################	
	
# ERROR ALERT: Run the script below with one parenthesis removed and see what happens
	sroot_mass<-sqrt(mass)
	#
################################################################################
	
#
# EXERCISE 1.3:  Operator precedence
#
# (i)  You may recall from school that, when faced with a string of operations, e.g
#
#  10-4*(2+1)
#
# you do the operations in a certain order, completing calculations inside brackets first.
# R follows a standard order of precedence in operations.  Using the rules from school
# work out the correct answer to the expression above, and then assign its value in R
# to a new numeric scalar, and see if you were correct.  
#
# (ii)  Work out the order of precedence in R for ^, -, +, * by examining the values 
# assigned to numeric scalars by the following

	1+2*3
	3*2^2
	3*2^3-1

###########################################################################################
#
# 1.2  CHARACTER DATA
#
# A character scalar is just a string of characters, for example, the name of a treatment
# 
	tname<-"Mulched"
#
# Note that the character string in the assignment is put inside double quotes.  See
# what happens if you run the command above, but with the quotes removed.
#
# It can be useful to use character scalars, as they can appear in commands such as 
# those for data plots, and you can assign the value once, so, for example, the following
# commands would create a pane of three graphs based on two variables, (here on some 
# random data).  In each plot command ("plot" and "hist") "xlab= " specifies the name on
# the x-axis label, and similarly "ylab=" for the y-axis.  One can put the label name in 
# directly here, xlab="pH", for example, but we can also put a character scalar here 
# which has been given a value elswhere e.g.:

# Assign variable names

	xname<-"pH"
	yname<-"SOC"

	par(mfrow=c(2,2))
	x<-runif(100,4,8)
	y<-rnorm(100,4,0.5)
	hist(x,xlab=xname,main="Histogram")
	hist(y,xlab=yname,main="Histogram")
	plot(x,y,xlab=xname,ylab=yname,pch=16,main="Scatterplot")

# If you were using the script to produce such plots from various variables you can see 
# how using a character scalar saves you from having to type the same variable name into
# the function for each plot.  Changing it once at the top ensures that you get the
# correct name in each case.
#


# Using "paste" to combine character variables.

# Imagine that I had a character scalar that denotes the block to which an experimental
# plot belongs in a RCBD experiment, and another one that denotes the treatment:

	block<-"Block1"
	treatment<-"CA"

# I can make a plot name by combining these two using paste.  The "sep" term allows
# me to specify the separator between the two scalars:

	plot_lab<-paste(block,treatment,sep="_")

	print(plot_lab)


	
# EXERCISE 1.6 (by yourself) Create a character scalar that includes your name and then use paste to join this
# with your age in weeks(as computed in the section on numeric data types).
#
#
###########################################################################################
#
# 1.3 LOGICAL DATA
#
# A logical scalar takes the value TRUE or FALSE.  An R command which states some relation 
# between two variables will have a logical value.  For example, let us create two numeric
# scalars

	three<-3
	five<-5

# Now the R statement (three<five) will take the value TRUE, because the value of "three" is 
# less than the value of "five", so

	three_lt_five<-(three<five)

	print(three_lt_five)

# the command below will show that it is not the case that three<three ....

	three_lt_three<-(three<three)
	print(three_lt_three)

# .... but <= (less than or equal to) gives us a different outcome ...

	three_le_three<-(three<=three)
	print(three_le_three)

# Some other useful "logical connectives" are == for "equal to" and != for "not equal to"
# and, of course, > for "greater than" and >= for "greater than or equal to".
#
# Note that == is used for "equal to", a single = will allocate the value of the scalar on the 
# left to that on the right.
	

	

# EXERCISE 1.7 satisfy yourself that the set of connectives described above behave as they do,
# using the scalars three and five above, and others of your creation.
#

	
# Logical variables can be the subject of logical functions, notably "if  .. then"
# Consider the example below

	soil_pH<-4

	if(soil_pH<5.5) {management_option<-"Lime"}else{management_option<-"No_lime"}

# In the script above if soil_pH is less than 5.5 then the scalar management_option will be given the
# character value "Lime", otherwise it will be given the value "No_lime".	
	
	
	time<-12.00
	if(time<12) {learning<-"continue"}else{learning_option<-"take a break"}

# In the script above if time is less than 12.00 then the scalar learning option will be given the
	# character value "continue", otherwise it will be given the value "take a break".

# A logical variable can be defined on the basis of more than one logical condition, this can 
# be done using the conditionals  & for "and", | for "or" (&& and || are sometimes applied to vectors 
# of logical variables.  Here is an example.  We define three numeric scalars sand, silt and clay as
# the percent by mass of sand, silt and clay-sized particles in soil.

	sand<-10
	silt<-20
	clay<-70

# first, check that the values are consistent

	consistent_particle_size<-((sand+silt+clay)==100)
	print(consistent_particle_size)

# The USDA soil texture class Clay contains soils with more than 40% clay AND less than 40% silt
# and less than 45% sand,  so we can determine whether or not our soil belongs to class clay as 
# follows

	is_clay<-(clay>40)&(silt<=40)&(sand<=45)
	print(is_clay)
	


#
# EXERCISE 1.8
#
# (i) try the commands out with some different (consistent) particle size values.
# (ii) modify the commands above so as to compute a logical variable is_silty_clay.  In the
# USDA texture triangle a soil is silty clay if the clay content is greater than 40% AND the
# silt content is greater than 40%: E.G.:
	

######################################################################################################
######################################################################################################
#
######################################################################################################
#
# B.  DATA STRUCTURES:  
#
######################################################################################################
#
#
#  A data structure in R is an R object which holds one or more data objects, a data object will be a
#  a data type, such as we have encountered in section 1 (numeric, character, etc).  In this script we introduce vectors, 
#  factors, matrices, data frames and lists. The examples and exercises should help you to understand
#  better how R holds and manages data.

######################################################################################################
#
#  1.4  VECTORS
#
# A vector is a series of values of a variable (e.g. pH value measurements
# from a sensor).  The easiest way to form a vector of values in R is with the "combine" function c().
# An example of a vector of numeric values (pH readings) is shown below:
#
	pH_values<-c(5.6, 5.5, 5.0, 5.7, 5.4, 5.3, 6.0, 6.7, 6.5, 6.4, 6.2, 6.3)
#
#
# We can count the number of items in a vector with the length() function:

	length(pH_values)

# Each item in a vector can be referenced by its index (i.e. its position in the sequence of values), 
# and we can pull out a particular item using the square brackets after the vector name.  For example,
# the 7th item in pH_values can be accessed like this

	pH_values[7]

# 
# A vector is a "homogeneous" data structure.  You could make a vector of logical values or of
# numeric values or of character values, but not a vector which has a mixture. 
# See what happens if you try:

	mixture<-c(5.2, TRUE, "CA")
	print(mixture)
#
# everything is turned into a character, and so appears in quotation marks.  You could not perform
# arithmetic on the first object in the vector, as you can see if you try:

sqrt(mixture[1])

#
# The operations that we have applied above to scalar data types can be applied to vectors.  So, for
# example, see what happens with the R command

	print(2*pH_values)
  
# the operation is applied to every element in the vector, and the output is a vector.
# 
# There are other functions that we can apply to vectors for example

	mean_pH<-mean(pH_values) 
	print(mean_pH)
#
########## ########### ########### ########### ############ ########## ############ ###########
	
# EXERCISE 1.9  Explore what the functions sum() and median() do using the pH_values vector

##############################################################################################
	
# A conditional operation applied to a vector will produce a vector of logical values

	pH_le_6<-pH_values<=6.0
	print(pH_le_6)

#
# The which() function, applied to a vector, will extract the index values for all values in the 
# vector which meet certain conditions, e.g.:
#
	index_lt_6<-which(pH_values<6)
	print(index_lt_6)
#
# we can use this vector of index values to extract the pH values which meet the condition into
# a new vector:
#
	small_pH_values<-pH_values[index_lt_6]
	print(small_pH_values)
	
# 

# EXERCISE 1.11    Extract into a new vector the values of pH which exceed 5.6	
#

# If we use the multiplication operation on two vectors, a and b, which are the same length...

	a<-c(1,2,3,4,5)
	b<-c(1,10,100,1000,10000)
	c<-a*b
	print(c)

# ... then the output, c, is a vector of the same length as a or b, where c[i]= a[i]*b[i]
# If a and b were of different length you would get an error message.
# Note that this is not a standard product of vectors from matrix algebra.


#######################################################################################################
#
#
# There are various useful commands for creating vectors, rep() is one of the best.  
# Assume that our pH values are drawn from six soil samples which are sands (S) 
# and six which are clay loam (CL).  We can make a vector of character values which 
# corresponds to the pH data with the following command:

soil_type<-rep(c("S","CL"),each=6)
print(soil_type)

# and you could then extract the index of the clay loams:

	clay_loam_index<-which(soil_type=="CL")

# .... and then extract the pH values for the clay loam soils

	pH_clay_loam<-pH_values[clay_loam_index]
#
#

	
# EXERCISE 1.13  Compare the output of the rep commands below, in order to work out what "times" and 
# "each" are doing

	rep(c("CA","Conv"),times=6)
	rep(c("CA","Conv"),each=6)
	rep(c("CA","Conv"),each=2,times=3)
	rep(c("CA","Conv"),times=1:2)

# Use rep to produce a vector of treatment labels for an experiment in which conservation
# agriculture (CA) and conventional (Conv) treatments are applied, with six reps of each in
# 6 blocks.  Produce a corresponding vector of block labels.

	

######################################################################################################
#
#  1.5  FACTORS
#
#
# In section 1.2 we introduced the idea of vectors of character variables as treatment 
# labels.  However, in order to be most useful, such a vector needs to be turned into a factor.
# A factor is a variable which is not a continuous number, or is not treated as one.  It is a label
# for some variable controlled at different levels in an experiment, so a factor might be Nitrogen 
# (application rate) with levels 0, 50 and 100 kg/ha, or it might be Variety with levels
# AB_123, AB_234, CD_120, CD_130.  When we first set up a vector of N rates these
# could be numeric, the Varieties will be made character variables because they contain letters.  In
# both cases, however, we are likely to want the variable to be turned into a factor for use in an
# analysis of variance, for example.  The example below sets up two vectors with levels for these factors
# in an experiment with the two factors in factorial combination (giving 12 treatments) and with four 
# replicates of each treatment in randomized blocks
#
	Nitrogen<-rep(c(0,50,100),16)
	Variety<-rep(c("AB_123", "AB_234", "CD_120", "CD_130"),each=3,times=4)
	Blocks<-rep(c("Block_1","Block_2","Block_3","Block_4"),each=12)

# The cbind command below will print out these vectors as columns in a matrix of character values
# it is easy to see how the twelve combinations of the two factor levels are structured in each block

	cbind(Blocks,Variety,Nitrogen)

#
# Now use the factor() function to create factor data structures.
#

	Nitrogen_factor<-factor(Nitrogen)
	print(Nitrogen_factor)

	Variety_factor<-factor(Variety)
	print(Variety_factor)

	Block_factor<-factor(Blocks)
	print(Block_factor)

# Note that a factor is actually a vector, but with an associated list of levels, always presented in 
# alpha-numeric order.  These are used by R functions such as lm() which does linear modelling, such as
# the analysis of variance.  We shall see how factors can be used in the later section on data frames.


	
# EXERCISE 1.14  An experiment has been set up to examine the combined effects of zero till vs conventional
# tillage with intercropping vs monocropping.  The experiment is designed with randomized complete blocks, 
# each of five blocks contains one replicate of each treatment.  Following the example above, set up the 
# three factors required to represent one season of this experiment.
	

#
######################################################################################################
#
# 1.6 MATRICES.
#
#
#
# A matrix is a rectangular array of values, arranged in rows and columns.  A vector
# is therefore a type of matrix with just one column.  We can create a matrix in R in
# one of two main ways.  The first is the matrix command:
#
  M1<-matrix(1:6,3,2)
 	print(M1)
#
# The first term in the command is a vector of numbers, 1 to 6, the second is the number 
# rows in the matrix and the third is the number of columns.  Note that we always refer 
# to an entry in the matrix by the ROW first and the COLUMN second.  As an exercise, look
# at the effect of swapping 3 and 2 round in the command above.  Note also that the command
# enters the terms down the first column then down the second, so that the first column of
# M1 goes 1,2,3 and the second 4, 5, 6.
#
# The second way to make a matrix in R is to "bind" some vectors together, which works only
# if they are the same length.  If we start with 2 vectors:
#
	a<-c(2,4,5)
	b<-c(4,7,10)

# .. then we can make a 3 x 2 matrix (remember that means 3 rows and 2 columns) with the
# cbind() command (for binding the vectors as columns of the matrix):

	M2<-cbind(a,b)
	print(M2) 
#
# .. and we could make a 2 x 3 matrix with the rbind() command (for binding the vectors as 
# rows of the matrix):

	M3<-rbind(a,b)
	print(M3) 
#
# Note that the vector names become column or row names of the matrix.  You can look at these
# names with the colnames or rownames command, and also use these to change the names:

colnames(M2)

colnames(M2)<-c("Column_1","Column_2")
colnames(M2)

# We can refer to a particular cell of a matrix as follows
#
	M2[2,2]
# 
# .. and we can refer to a particular column of the matrix as follows
	
	M2[,2]
# 
# For example, to find the sum of the first column of M2 ..

	Sum_Col_1<-sum(M2[,1])
	print(Sum_Col_1)
#
#

# EXERCISE 1.15  Take the vectors below (from an exercise above) which contain gravimetric
# water content and bulk density of ten soils.  Combine these two as columns in a matrix
# and then change the column names to "Gravimetric_water" and "Bulk_density".
	
	gwc<-c(0.4,0.5,0.3,0.2,0.5,0.6,0.3,0.2,0.4,0.3)
	rho<-c(1.05,1.42,1.50,1.65,1.44,0.90,1.35,1.36,1.10,1.43)
	gr<-cbind(gwc,rho)
#
# When that is done, extract the mean value of bulk density from the matrix.


######################################################################################################
#
# 1.7  DATA FRAMES.  
#
#
# A matrix or a vector is a "homogeneous" data structure, which means you can't mix data types.
# Consider the previous example, we have a vector of numeric data on soil pH:

	pH_values<-c(5.6, 5.5, 5.0, 5.7, 5.4, 5.3, 6.0, 6.7, 6.5, 6.4, 6.2, 6.3)
# 
# ... and a corresponding vector of character values on the soil texture class for
# each sample

	soil_type<-rep(c("S","CL"),each=6)
#
# .. and then attempt to put the character and the numeric vector into a matrix:

	combined_data<-cbind(soil_type,pH_values)

# print the outcome..

	print(combined_data)
#
# We can can see that the pH values are now character, a string "5.6" which could not
# be used in calculations.

# Let us see what happens if we convert soil type to a factor first..

	soil_type<-factor(soil_type)
	combined_data<-cbind(soil_type,pH_values)

# print the outcome..

	print(combined_data)

# ... the factor values have been converted to numeric ones (label 1 goes to the factor whose 
# original name came first in alphabetical order ("CL").  

# This is why matrices, while important for many applications, are not the most basic data structure
# in R.  The data frame serves this purpose, which is why we will generally use commands such as
# read.table or read.csv to read data from external files into an R data frame
#
#  We can turn our two vectors into a data frame as follows:
#
	combined.df<-data.frame(soil_type=soil_type,pH_values=pH_values,stringsAsFactors = TRUE)

# The option "stringsAsFactors = TRUE" tells R that soil_type should be treated as a factor.

# We can reference the data object inside the data frame using the dollar notation, combined.df$soil_type
# so we can confirm that soil_type is indeed a factor in the data frame as follows.

print(combined.df$soil_type)

# We now have a data frame in which the soil_type variable is a factor and the pH_values are
# numeric: just what we need.
#
	print(combined.df)
#
# This allows us to do some interesting things.  For example, make a boxplot of pH within each
# level of the factor
#
	boxplot(pH_values~soil_type,data=combined.df)
#
# ... or to extract the mean value of pH for each level of the factor
#
	by(combined.df$pH_values,combined.df$soil_type,mean)
#
# ... or to do a t-test to compare the mean pH values
#
	t.test(pH_values~soil_type,data=combined.df)
#
# Note that the terms in brackets for t.test and boxplot are identical.  They are (i) a
# formula saying express the first variable in terms of different levels of the factor and
# (ii) a data statement pointing R to the data frame where it will find the variables with
# the specified names.
	

#
# EXERCISE 1.16 	Find the mean of all the pH values in combined.df
#	Find the median values of pH in the two textural classes
	



# Three useful R commands, which can be applied to data frames and are particularly useful when
# examining one created by reading in data are names(), head() and nrow().  Apply these to 
# combined.df in order to work out what they do.
#
	names(combined.df)
	head(combined.df)
	nrow(combined.df)
#
# The colnames() command can be used either to show or to set the column names (as for matrices)
#
#
	colnames(combined.df)

#
#
# EXERCISE 1.17  The two sets of values below are, respectively, the bulk density of 20 different
# topsoil samples and their textural class: sandy loam (SL) or silt loam (KL).

  1.46,1.44,1.43,1.42,0.96,1.3,1.48,1.22,1.41,1.5,1.27,1.21,1.23,1.16,1.27,1.37,1.07,1.16,1.04,1.42

  SL,SL,SL,SL,SL,SL,SL,SL,SL,SL,KL,KL,KL,KL,KL,KL,KL,KL,KL,KL
  
 
  

# Make a data frame in which textural class is a factor and bulk density is a numerical variable
# Give each column an appropriate name then (i) produce a boxplot of bulk density in the two
# texture classes (ii) compute the mean bulk density over all the data and (iii) compute the mean
# bulk density for each textural class separately then (iv) conduct a t-test to compare the two
# textural classes with respect to bulk density.
  
  

#######################################################################################################
#
# 1.8 LISTS
#
#
#
#  A list in R is a vector of data objects.  It can be a useful structure for holding outputs from
# analyses in a consistent format, and you might find that some packages you use produce lists as
# outputs.
#
#  In this example we take three vectors of values and put them together in a list:
#

list_example<-list("some_numbers"=c(1,2,3,4), 
"some_odd_numbers"=c(1,3,5,7), 
"some_even_numbers"=c(2,4,6,8))

print(list_example)

# The three components of the list are referred to as list "slices"
# You can use the name of a slice and the dollar notation to refer to a particular slice

list_example$some_even_numbers

# and because a slice is a vector you can refer to a single element of it by its index (order in the
# vector) and the square brackets notation:

list_example$some_even_numbers[3]

# Alternatively, you can refer to slices of a list by using the double square bracket with an index
# (for the first, second... slice).  For example, an alternative way to reference the vector 
# "some_even_numbers", which is the third slice, is as follows

list_example[[3]]

# .. and you can refer to an element of this vector thus

list_example[[3]][4]

# For example, one may change the value of the fourth element in the third slice as follows

list_example[[3]][4]<-10

#	
#  Home Exercise 1.18  Produce a list object based on the data in Exercise 5.2 which contains the following
#  slices. (i) all the bulk density data  (ii) the list of corresponding textural classes (iii) the 
#  mean bulk density over all the data and (iv) the mean bulk densities for the two classes.

#
