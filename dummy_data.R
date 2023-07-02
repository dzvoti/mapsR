
# Create dataframe
soil_temp<-rnorm(90,20,1)
soil_moisture<-rnorm(90,12,2)
site<-rep(c("liempe","chitedze","domboshava"),30)
data_temp <- data.frame (soil_temp,soil_moisture,site)

head(data_temp)

write.csv(data_temp, here::here("data/soil-test.csv"), row.names = FALSE)
