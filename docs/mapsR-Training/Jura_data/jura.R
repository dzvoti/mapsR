


#######################################################################################

# Exercise: Exploratory data analysis


# Run the scripts below and:
# Discuss the summary statistics and distributions
# Carry out the transformations and discuss the results

#####################################################################################################################

source("CEPHaStat_3.R")
source("boxcox.r")

data.df<-read.csv("Jura.csv",header=T,stringsAsFactors=T)


#  Pb analyses

# Distributions and Summary statistics

dev.new()
summa(data.df$Pb)
summaplot(data.df$Pb,"Pb /ppm")

dev.off()

dev.new()

# Remember to base transformations on residuals, after fitting the model!

mod<-lm(Pb~Rock,data=data.df)

summa(mod$residuals)
summaplot(mod$residuals,"Residual /ppm Pb")
dev.off()

##################################################################################################

# COMPARING LOG AND BOXCOX TRANSFORMATIONS

# Remember to base transformations on residuals, after fitting the model!

#log transformation

dev.new()

mod2<-lm(log(Pb)~Rock,data=data.df)
summa(mod2$residuals)
summaplot(mod2$residuals,"Residual /Bcox ppm Pb")
anova(mod2)

dev.off()

#BoxCox transformation

dev.new()

mod3<-lm(BoxCox(Pb)~Rock,data=data.df)
summa(mod3$residuals)
summaplot(mod3$residuals,"Residual /Bcox ppm Pb")

anova(mod3)
dev.off()

#####################################################################################################

# Co analyses

summa(data.df$Co)
summaplot(data.df$Co,"Co /ppm")

mod<-lm(Co~Rock,data=data.df)

summa(mod$residuals)
summaplot(mod$residuals,"Residual /ppm Co")

# What do you conclude from the Co analysis?

###########################################################################################################

