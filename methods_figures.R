#### Methods figures ##
## WA Perenjori Watering Experiment
# April 2022

#Sowed seeds between 12th February 2020 and 29th February 2020
# Tally germination and focal thinning 18th July?
# Neighbourhood thinning 
# Rainout shelters built in mid-August
# Surveys done on 11, 12, 13, 14 September, and 29/30th September (of things that died previously)
# Only want to plot February - October

#### Rainfall events in growing season ####
## Import BOM rainfall data from Perenjori weather station, daily for 2022
rainfalldata2020 <- read_csv("Data/BOM_2020_rainfall_data.csv")

#Calculate amount of total rainfall for growing season in 2020
rainfalldata2020 <- within(rainfalldata2020, Month[Month == 1] <- 'January')
rainfalldata2020 <- within(rainfalldata2020, Month[Month == 2] <- 'February')
rainfalldata2020 <- within(rainfalldata2020, Month[Month == 3] <- 'March')
rainfalldata2020 <- within(rainfalldata2020, Month[Month == 4] <- 'April')
rainfalldata2020 <- within(rainfalldata2020, Month[Month == 5] <- 'May')
rainfalldata2020 <- within(rainfalldata2020, Month[Month == 6] <- 'June')
rainfalldata2020 <- within(rainfalldata2020, Month[Month == 7] <- 'July')
rainfalldata2020 <- within(rainfalldata2020, Month[Month == 8] <- 'August')
rainfalldata2020 <- within(rainfalldata2020, Month[Month == 9] <- 'September')
rainfalldata2020 <- within(rainfalldata2020, Month[Month == 10] <- 'October')
rainfalldata2020 <- within(rainfalldata2020, Month[Month == 11] <- 'November')
rainfalldata2020 <- within(rainfalldata2020, Month[Month == 12] <- 'December')
#Reordering them by month order
rainfalldata2020$Month <- factor(rainfalldata2020$Month, level = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))
gs_months <- c("June", "July", "August", "September", "October")
rainfall_2020_gs <- filter(rainfalldata2020, Month %in% gs_months)
gs_2020_summary <- rainfall_2020_gs %>% summarise(total_rainfall = sum(Rainfall_mm, na.rm = TRUE))

#Plotting rainfall for 2020 growing season
#Filter to rainfall events
rainevents <- rainfalldata2020 %>% filter(Rainfall_mm > 0)
# Rained 26 times all year
rainevents_gs <- rainfall_2020_gs %>% filter(Rainfall_mm > 0)
# Rained 10 times during gs

#Plot entire year
#Need dates to be recognised as dates so they plot in order

rainevents$Date <- as.Date(rainevents$Date, "%d/%m/%Y")
ggplot(rainevents, aes(x = Date, y = Rainfall_mm))+
  geom_point()+
  geom_line()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

#Plot growing season
rainevents_gs$Date <- as.Date(rainevents_gs$Date, "%d/%m/%Y")
#Make `Period over which rainfall was measured (days)` a factor to plot it
#Some are 3 and 4 day long rain events, others are 3. Can I colour-code by that, create a legend?
rainevents_gs$`Period over which rainfall was measured (days)` <- as.factor(rainevents_gs$`Period over which rainfall was measured (days)`)

dev.off()
pdf("Output/Figures/rainfall_events_growing_season.pdf")
ggplot(rainevents_gs, aes(x = Date, y = Rainfall_mm, colour = `Period over which rainfall was measured (days)`))+
  geom_col(aes(fill = `Period over which rainfall was measured (days)`))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_text(size = 14, face = 'bold'),
        axis.title.y = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = "bottom")
dev.off()
### Calculating long-term monthly average rainfall in Perenjori ####
# Monthly average rainfall for 1918 - 2022
longtermrainfalldata <- read_csv("Data/perenjori_longterm_monthly_rainfall.csv")

#Renaming months from numbers to names
longtermrainfalldata <- within(longtermrainfalldata, Month[Month == 1] <- 'January')
longtermrainfalldata <- within(longtermrainfalldata, Month[Month == 2] <- 'February')
longtermrainfalldata <- within(longtermrainfalldata, Month[Month == 3] <- 'March')
longtermrainfalldata <- within(longtermrainfalldata, Month[Month == 4] <- 'April')
longtermrainfalldata <- within(longtermrainfalldata, Month[Month == 5] <- 'May')
longtermrainfalldata <- within(longtermrainfalldata, Month[Month == 6] <- 'June')
longtermrainfalldata <- within(longtermrainfalldata, Month[Month == 7] <- 'July')
longtermrainfalldata <- within(longtermrainfalldata, Month[Month == 8] <- 'August')
longtermrainfalldata <- within(longtermrainfalldata, Month[Month == 9] <- 'September')
longtermrainfalldata <- within(longtermrainfalldata, Month[Month == 10] <- 'October')
longtermrainfalldata <- within(longtermrainfalldata, Month[Month == 11] <- 'November')
longtermrainfalldata <- within(longtermrainfalldata, Month[Month == 12] <- 'December')

#Summarising averages by month 
rainfallsummary <- longtermrainfalldata %>% group_by(Month) %>% summarise(mean_rainfall = mean(monthly_total_precipitation_mm),
                                                               n = n(),
                                                               sd = sd(monthly_total_precipitation_mm), 
                                                               se_rainfall = sd/sqrt(n))
#Reordering them by month order
rainfallsummary$Month <- factor(rainfallsummary$Month, level = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))

#Plotting long-term average montly rainfall
ggplot(rainfallsummary, aes (x = Month, y = mean_rainfall))+
  geom_point()+
  geom_errorbar(aes(ymax = mean_rainfall + se_rainfall, ymin = mean_rainfall - se_rainfall))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

# Calculating longterm average total rainfall and se only for growing season, June - Oct
gs_months <- c("June", "July", "August", "September", "October")
gs_longterm_rainfall <- filter(longtermrainfalldata, Month %in% gs_months)
gs_summary <- gs_longterm_rainfall %>% group_by(Year) %>% summarise(total_rainfall = sum(monthly_total_precipitation_mm),
                                                                     n = n(),
                                                                     sd = sd(monthly_total_precipitation_mm), 
                                                                     se_rainfall = sd/sqrt(n))
gs_summary2 <- gs_summary %>% summarise(mean_total_rainfall = mean(total_rainfall),
                                         n = n(),
                                         sd = sd(total_rainfall), 
                                         se_rainfall = sd/sqrt(n))





