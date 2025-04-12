rm(list = ls())

setwd("C:/Users/DELL/Documents/R/projects/MP1_CHLA_NPP_ER")
data<-read.csv("NPP_ER.csv")


install.packages("dplyr")
install.packages("scales")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("patchwork")
install.packages("gridExtra")
install.packages("RColorBrewer")
install.packages("lubridate")
install.packages("lme4")
install.packages("mgcv")
install.packages("emmeans")
install.packages("reshape2")
install.packages("zoo")
install.packages("viridis")
install.packages("akima")
install.packages("car")
install.packages("readr")
install.packages("gtable")
install.packages("grid")
library(reshape2)
library(emmeans)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(RColorBrewer)
library(lubridate)
library(lme4)
library(mgcv)
library(zoo)
library(viridis)
library(akima)
library(car)
library(readr)
library(gtable)
library(grid)
library(cowplot)

#Convert the date column to Date format
data$date <- as.Date(data$date, format="%m/%d/%Y")

#Change Tank_ID based on conditions for all dates after 3/8/23 - Error in tank labeling corrected
data <- data %>%
  mutate(Tank_ID = ifelse(date > as.Date("2023-03-08") & Tank_ID == 4, 21,
                          ifelse(date > as.Date("2023-03-08") & Tank_ID == 21, 4,
                                 Tank_ID)))

#Check for this correction
view(data)

#Add unique position to each row based on tank.id
data <- data %>%
  mutate(Tank_ID = case_when(
    Tank_ID == "1" ~ "A0",
    Tank_ID == "2" ~ "C3",
    Tank_ID == "3" ~ "B6",
    Tank_ID == "4" ~ "C9",
    Tank_ID == "5" ~ "A3",
    Tank_ID == "6" ~ "B2",
    Tank_ID == "7" ~ "C5",
    Tank_ID == "8" ~ "B3",
    Tank_ID == "9" ~ "A1",
    Tank_ID == "10" ~ "B9",
    Tank_ID == "11" ~ "B5",
    Tank_ID == "12" ~ "C0",
    Tank_ID == "13" ~ "A9",
    Tank_ID == "14" ~ "C6",
    Tank_ID == "15" ~ "B0",
    Tank_ID == "16" ~ "B8",
    Tank_ID == "17" ~ "A6",
    Tank_ID == "18" ~ "C2",
    Tank_ID == "19" ~ "A8",
    Tank_ID == "20" ~ "B7",
    Tank_ID == "21" ~ "C7",
    Tank_ID == "22" ~ "C8",
    Tank_ID == "23" ~ "A7",
    Tank_ID == "24" ~ "A5",
    Tank_ID == "25" ~ "C1",
    Tank_ID == "26" ~ "B4",
    Tank_ID == "27" ~ "A4",
    Tank_ID == "28" ~ "B1",
    Tank_ID == "29" ~ "A2",
    Tank_ID == "30" ~ "C4",
  ))

#Sort the data by Tank_ID, date, and then AM_PM
data <- arrange(data, Tank_ID, date, AM_PM)

#Separate AM and PM data
am_data <- filter(data, AM_PM == "AM")
pm_data <- filter(data, AM_PM == "PM")

#Merge AM and PM data on Tank_ID and date for NPP calculation
npp_data <- merge(am_data, pm_data, by=c("Tank_ID", "date"), suffixes=c("_AM", "_PM"))

#Calculate NPP as the difference in DO_perc from AM to PM
npp_data$NPP <- npp_data$DO_perc_PM - npp_data$DO_perc_AM

#For ER calculation, shift PM data to align with the next day's AM data
pm_data$date <- pm_data$date + 1

#Merge shifted PM data with AM data on Tank_ID and date for ER calculation
er_data <- merge(pm_data, am_data, by=c("Tank_ID", "date"), suffixes=c("_PM_prev", "_AM"))

#Calculate ER as the difference in DO_perc from PM to the next AM
er_data$ER <- er_data$DO_perc_AM - er_data$DO_perc_PM_prev

#Select relevant columns for NPP and ER dataframes
npp_data <- npp_data[c("Tank_ID", "date", "NPP")]
er_data <- er_data[c("Tank_ID", "date", "ER")]

#Display first few rows of NPP and ER data
head(npp_data)
head(er_data)

#Prepare npp_data
npp_data <- npp_data %>%
  mutate(
    plastic_type = substr(as.character(Tank_ID), 1, 1),
    plastic_conc = as.numeric(substr(as.character(Tank_ID), 2, nchar(as.character(Tank_ID)))),
    days_since_start = as.numeric(date - min(date, na.rm = TRUE))
  )

#Prepare er_data
er_data <- er_data %>%
  mutate(
    plastic_type = substr(as.character(Tank_ID), 1, 1),
    plastic_conc = as.numeric(substr(as.character(Tank_ID), 2, nchar(as.character(Tank_ID)))),
    days_since_start = as.numeric(date - min(date, na.rm = TRUE))
  )


#Heatmap for plastic type A - NPP
ggplot(subset(npp_data, plastic_type == "A"), 
       aes(x = factor(days_since_start), y = plastic_conc, fill = NPP)) +
  geom_tile() +
  scale_fill_gradient(low = "gray", high = "red") +
  labs(title = "NPP for Plastic Type A", x = "Days Since Start", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Heatmap for plastic type B - NPP
ggplot(subset(npp_data, plastic_type == "B"), 
       aes(x = factor(days_since_start), y = plastic_conc, fill = NPP)) +
  geom_tile() +
  scale_fill_gradient(low = "gray", high = "red") +
  labs(title = "NPP for Plastic Type B", x = "Days Since Start", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Heatmap for plastic type C - NPP
ggplot(subset(npp_data, plastic_type == "C"), 
       aes(x = factor(days_since_start), y = plastic_conc, fill = NPP)) +
  geom_tile() +
  scale_fill_gradient(low = "gray", high = "red") +
  labs(title = "NPP for Plastic Type C", x = "Days Since Start", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Heatmap for plastic type A - ER
ggplot(subset(er_data, plastic_type == "A"), 
       aes(x = factor(days_since_start), y = plastic_conc, fill = ER)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "gray") +
  labs(title = "ER for Plastic Type A", x = "Days Since Start", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Heatmap for plastic type B - ER
ggplot(subset(er_data, plastic_type == "B"), 
       aes(x = factor(days_since_start), y = plastic_conc, fill = ER)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "gray") +
  labs(title = "ER for Plastic Type B", x = "Days Since Start", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Heatmap for plastic type C - ER
ggplot(subset(er_data, plastic_type == "C"), 
       aes(x = factor(days_since_start), y = plastic_conc, fill = ER)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "gray") +
  labs(title = "ER for Plastic Type C", x = "Days Since Start", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#NPP over time for each Tank_ID - check for removing sample read errors
ggplot(npp_data, aes(x = date, y = NPP, color = Tank_ID)) +
  geom_point() +
  labs(title = "NPP Over Time", x = "Date", y = "NPP (% DO)") +
  theme_minimal()

#ER over time for each Tank_ID - check for removing sample read errors
ggplot(er_data, aes(x = date, y = ER, color = Tank_ID)) +
  geom_point() +
  labs(title = "Ecosystem Respiration (ER) Over Time", 
       x = "Date", 
       y = "ER (% DO)") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_viridis_d(name = "Tank ID")


#Identify rows where NPP < 0 (to be removed)
npp_removed <- npp_data %>% 
  filter(NPP < 0)

#Count the number of removed samples
n_npp_removed <- nrow(npp_removed)
cat("Number of NPP samples removed due to impossible values (NPP < 0):", n_npp_removed, "\n")

#Print the removed samples
if (n_npp_removed > 0) {
  print("Removed NPP samples:")
  print(npp_removed)
}


#Identify rows where ER > 0 (to be removed)
er_removed <- er_data %>% 
  filter(ER > 0)

#Count the number of removed samples
n_er_removed <- nrow(er_removed)
cat("Number of ER samples removed due to impossible values (ER > 0):", n_er_removed, "\n")

#Print the removed samples
if (n_er_removed > 0) {
  print("Removed ER samples:")
  print(er_removed)
}


#Remove NPP sampling errors
npp_data <- npp_data %>%
  mutate(NPP = ifelse(NPP < 0, NA, NPP))

#Remove ER sampling errors
er_data <- er_data %>%
  mutate(ER = ifelse(ER > 0, NA, ER))
  
#Run all the heatmaps again
#Plastic type A - NPP
ggplot(subset(npp_data, plastic_type == "A"), 
       aes(x = factor(days_since_start), y = plastic_conc, fill = NPP)) +
  geom_tile() +
  scale_fill_gradient(low = "gray", high = "red", limits = c(0, 100), na.value = NA) +
  labs(title = "NPP for Plastic Type A", x = "Days Since Start", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Plastic type B - NPP
ggplot(subset(npp_data, plastic_type == "B"), 
       aes(x = factor(days_since_start), y = plastic_conc, fill = NPP)) +
  geom_tile() +
  scale_fill_gradient(low = "gray", high = "red", limits = c(0, 100), na.value = NA) +
  labs(title = "NPP for Plastic Type B", x = "Days Since Start", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Plastic type C - NPP
ggplot(subset(npp_data, plastic_type == "C"), 
       aes(x = factor(days_since_start), y = plastic_conc, fill = NPP)) +
  geom_tile() +
  scale_fill_gradient(low = "gray", high = "red", limits = c(0, 100), na.value = NA) +
  labs(title = "NPP for Plastic Type C", x = "Days Since Start", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
#Plastic type A - ER
ggplot(subset(er_data, plastic_type == "A"), 
       aes(x = factor(days_since_start), y = plastic_conc, fill = ER)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "gray", limits = c(-100, 0), na.value = NA) +
  labs(title = "ER for Plastic Type A", x = "Days Since Start", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Plastic type B - ER
ggplot(subset(er_data, plastic_type == "B"), 
       aes(x = factor(days_since_start), y = plastic_conc, fill = ER)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "gray", limits = c(-100, 0), na.value = NA) +
  labs(title = "ER for Plastic Type B", x = "Days Since Start", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Plastic type C - ER
ggplot(subset(er_data, plastic_type == "C"), 
       aes(x = factor(days_since_start), y = plastic_conc, fill = ER)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "gray", limits = c(-100, 0), na.value = NA) +
  labs(title = "ER for Plastic Type C", x = "Days Since Start", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



  
##Try to get contour lines
Aselected_letter_group <- "A"
Bselected_letter_group <- "B"
Cselected_letter_group <- "C"

#Filter the data for the selected letter group
Aselected_data <- filter(npp_data, plastic_type == Aselected_letter_group)
Bselected_data <- filter(npp_data, plastic_type == Bselected_letter_group)
Cselected_data <- filter(npp_data, plastic_type == Cselected_letter_group)

#For NPP
#Calculate the mean and standard deviation for the selected group
###For A###
Asummary_data_contour <- Aselected_data %>%
  group_by(date, plastic_conc) %>%
  summarise(
    Mean = mean(NPP),
    SD = sd(NPP),
    .groups = 'drop'
  )

###For B###
Bsummary_data_contour <- Bselected_data %>%
  group_by(date, plastic_conc) %>%
  summarise(
    Mean = mean(NPP),
    SD = sd(NPP),
    .groups = 'drop'
  )

###For C###
Csummary_data_contour <- Cselected_data %>%
  group_by(date, plastic_conc) %>%
  summarise(
    Mean = mean(NPP),
    SD = sd(NPP),
    .groups = 'drop'
  )

#Format the date without the year and store it in a new variable
Asummary_data_contour$date_without_year <- format(as.Date(Asummary_data_contour$date, format = "%m/%d/%Y"), "%m/%d")

Bsummary_data_contour$date_without_year <- format(as.Date(Bsummary_data_contour$date, format = "%m/%d/%Y"), "%m/%d")

Csummary_data_contour$date_without_year <- format(as.Date(Csummary_data_contour$date, format = "%m/%d/%Y"), "%m/%d")



#Change plastic_conc values to actual concentrations
###For A###
Asummary_data_contour <- Asummary_data_contour %>%
  mutate(plastic_conc = case_when(
    plastic_conc == 0 ~ 0.000,
    plastic_conc == 1 ~ 0.004,
    plastic_conc == 2 ~ 0.008,
    plastic_conc == 3 ~ 0.013,
    plastic_conc == 4 ~ 0.023,
    plastic_conc == 5 ~ 0.041,
    plastic_conc == 6 ~ 0.072,
    plastic_conc == 7 ~ 0.126,
    plastic_conc == 8 ~ 0.220,
    plastic_conc == 9 ~ 0.385,
    TRUE ~ plastic_conc  
  ))

###For B###
Bsummary_data_contour <- Bsummary_data_contour %>%
  mutate(plastic_conc = case_when(
    plastic_conc == 0 ~ 0.000,
    plastic_conc == 1 ~ 0.004,
    plastic_conc == 2 ~ 0.008,
    plastic_conc == 3 ~ 0.013,
    plastic_conc == 4 ~ 0.023,
    plastic_conc == 5 ~ 0.041,
    plastic_conc == 6 ~ 0.072,
    plastic_conc == 7 ~ 0.126,
    plastic_conc == 8 ~ 0.220,
    plastic_conc == 9 ~ 0.385,
    TRUE ~ plastic_conc  
  ))

###For C###
Csummary_data_contour <- Csummary_data_contour %>%
  mutate(plastic_conc = case_when(
    plastic_conc == 0 ~ 0.000,
    plastic_conc == 1 ~ 0.004,
    plastic_conc == 2 ~ 0.008,
    plastic_conc == 3 ~ 0.013,
    plastic_conc == 4 ~ 0.023,
    plastic_conc == 5 ~ 0.041,
    plastic_conc == 6 ~ 0.072,
    plastic_conc == 7 ~ 0.126,
    plastic_conc == 8 ~ 0.220,
    plastic_conc == 9 ~ 0.385,
    TRUE ~ plastic_conc  
  ))


#Convert date to a Date object and then to a numeric scale
###For A###
Asummary_data_contour$date <- as.Date(Asummary_data_contour$date, format = "%m/%d/%Y")
start_date <- min(Asummary_data_contour$date)
Asummary_data_contour$date_num <- as.numeric(Asummary_data_contour$date - start_date)

Asummary_data_contour$plastic_conc <- as.numeric(Asummary_data_contour$plastic_conc)

###For B###
Bsummary_data_contour$date <- as.Date(Bsummary_data_contour$date, format = "%m/%d/%Y")
start_date <- min(Bsummary_data_contour$date)
Bsummary_data_contour$date_num <- as.numeric(Bsummary_data_contour$date - start_date)

Bsummary_data_contour$plastic_conc <- as.numeric(Bsummary_data_contour$plastic_conc)


###For C###
Csummary_data_contour$date <- as.Date(Csummary_data_contour$date, format = "%m/%d/%Y")
start_date <- min(Csummary_data_contour$date)
Csummary_data_contour$date_num <- as.numeric(Csummary_data_contour$date - start_date)

Csummary_data_contour$plastic_conc <- as.numeric(Csummary_data_contour$plastic_conc)

#Get date without year
###For A###
Asummary_data_contour <- select(Asummary_data_contour, -date_without_year)
Asummary_data_contour <- select(Asummary_data_contour, -date)

###For B###
Bsummary_data_contour <- select(Bsummary_data_contour, -date_without_year)
Bsummary_data_contour <- select(Bsummary_data_contour, -date)

###For C
Csummary_data_contour <- select(Csummary_data_contour, -date_without_year)
Csummary_data_contour <- select(Csummary_data_contour, -date)

###For A###
#All combinations of date_num and plastic_conc
all_combinations <- expand.grid(
  date_num = unique(Asummary_data_contour$date_num),
  plastic_conc = unique(Asummary_data_contour$plastic_conc)
)

#Find missing combinations
Amissing_combinations <- anti_join(all_combinations, Asummary_data_contour, by = c("date_num", "plastic_conc"))

#Interpolate mean and SD for missing combinations
if (nrow(Amissing_combinations) > 0) {
  interpolated_values <- Amissing_combinations %>%
    mutate(
      Mean = NA, 
      SD = NA     
    ) %>%
    group_by(plastic_conc) %>%
    mutate(
      Mean = mean(Asummary_data_contour$Mean, na.rm = TRUE),
      SD = sd(Asummary_data_contour$Mean, na.rm = TRUE)
    )
  
  #Add interpolated values to the dataframe
  Asummary_data_contour <- bind_rows(Asummary_data_contour, interpolated_values)
} else {
  cat("No missing combinations found in the dataframe.\n")
}

#Add a pseudocount of 0.001 to treatment concentrations
Asummary_data_contour$plastic_conc_pseudocount <- Asummary_data_contour$plastic_conc + 0.009

#Asummary_data_contour
ggplot(Asummary_data_contour, aes(x = date_num, y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  theme_minimal() +
  labs(title = "NPP for Plastic Type A - Tiles Only")

#Correcting the geom_contour usage so that y is on a log axis
###For A###
heatmap_contour_plot_A_log <- ggplot(Asummary_data_contour, aes(x = date_num, y = plastic_conc_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean), width = 2.9, height = 0.10) +  # Adjust width and height parameters to minimize gaps
  geom_contour(color = "black") +  # This adds the contour lines based on `z = Mean`
  scale_fill_viridis_c(option = "plasma", direction = -1) + # Use 'plasma' palette and reverse gradient
  scale_x_continuous(name = "Day Number", breaks = seq(min(Asummary_data_contour$date_num), max(Asummary_data_contour$date_num), by = 5)) + # Add x-axis tick marks with a step of 5
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +  # Set y-axis to log scale
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  labs(fill = "NPP", title = paste("Elastollan"))

#Print the plot
print(heatmap_contour_plot_A_log)



#Day clusters for averaging - repeat measures
day_clusters <- list(
  c(0, 1, 2),
  c(6, 7, 8),
  c(16, 17, 18),
  c(26, 28, 29),
  c(36, 37, 38),
  c(56, 57, 58),
  c(76, 77),
  c(97, 98, 99)
)

#Recalcualte averages for each cluster
averaged_data <- do.call(rbind, lapply(day_clusters, function(cluster) {
  subset <- Asummary_data_contour[Asummary_data_contour$date_num %in% cluster,]
  
  #Aggregate by plastic_conc_pseudocount to calcualte mean
  averaged_subset <- aggregate(Mean ~ plastic_conc_pseudocount, data = subset, FUN = mean, na.rm = TRUE)
  
  #Calculate SD within each plastic_conc_pseudocount group for mean
  sd_subset <- aggregate(Mean ~ plastic_conc_pseudocount, data = subset, FUN = sd, na.rm = TRUE)
  
  #Add SD values to the averaged data frame
  averaged_subset$SD <- sd_subset$Mean
  
  #Assign an averaged date_num for identification
  averaged_subset$date_num <- rep(mean(cluster), nrow(averaged_subset))
  
  return(averaged_subset)
}))

#Convert into a data frame 
averaged_data <- data.frame(averaged_data)

#date_num as a factor for plotting
averaged_data$date_num <- factor(averaged_data$date_num)


#Heatmap for Group 1
g1_plot <- ggplot(subset(Asummary_data_contour, date_num %in% c(0, 6, 16, 26, 36, 56, 76, 97)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 1", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for Group 2
g2_plot <- ggplot(subset(Bsummary_data_contour, date_num %in% c(1, 7, 17, 28, 37, 57, 77, 98)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 2", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for Group 3
g3_plot <- ggplot(subset(Csummary_data_contour, date_num %in% c(2, 8, 18, 29, 38, 58, 76, 99)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 3", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for averaged data
avg_plot <- ggplot(averaged_data, aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Averaged Data", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Display plots
g1_plot
g2_plot
g3_plot
avg_plot

#Make sure 'date_num' is numeric in 'averaged_data'
averaged_data$date_num <- as.numeric(as.character(averaged_data$date_num))

#Interpolation step
interp <- with(averaged_data, {
  interp(x = date_num, y = log10(plastic_conc_pseudocount), z = Mean,
         xo = seq(min(date_num), max(date_num), length = 100),  
         yo = seq(min(log10(plastic_conc_pseudocount)), max(log10(plastic_conc_pseudocount)), length = 100),
         linear = FALSE)
})

#Transform the interpolated data back into a data frame for plotting
interp_df <- expand.grid(date_num = interp$x, plastic_conc_pseudocount = interp$y)
interp_df$Mean <- as.vector(interp$z)

#Plot the heatmap with interpolated data and contour lines
avg_plot_interp <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean)) +  
  geom_contour(bins = 10, color = "black") +  #Add contour lines
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100)) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "NPP Elastollan", x = "Averaged Day Number", y = "Plastic Concentration Pseudocount", fill = "NPP") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

#Print
print(avg_plot_interp)

#Plotting the heatmap with interpolated data without contour lines
avg_plot_interp_no_contours_ANPP <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() + 
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100)) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    expand = c(0, 0)
  ) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "Elastollan", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "NPP") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_ANPP)


#Make the edges transparent
avg_plot_interp_no_contours_ANPP <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile(na.rm = TRUE) +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100), na.value = NA) +  # Set NA values to be transparent
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    expand = c(0, 0)
  ) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "Elastollan", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "NPP") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_ANPP)

#Trim the edges so it is just within our data range
#Determine the limits based on the range of your data
x_limits <- range(interp_df$date_num, na.rm = TRUE)
y_limits <- range(10^interp_df$plastic_conc_pseudocount, na.rm = TRUE)

#Plot adjusted limits
avg_plot_interp_no_contours_ANPP <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100), na.value = NA) +  # Set NA values to be transparent
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    limits = x_limits,  
    expand = c(0, 0)
  ) +
  scale_y_log10(
    name = "Log Plastic Concentration (g/L)",
    limits = y_limits,  
    expand = c(0, 0)
  ) +
  theme_minimal() +
  labs(title = "Elastollan", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "NPP") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_ANPP)


#Format for final figure
# Define the original concentrations for Plot C with custom format for y-axis
original_concentrations <- c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385)

#Add pseudocount to the concentrations for plotting
pseudocount_concentrations <- original_concentrations + 0.009

#Set x_limits to cover the data range
x_limits <- range(interp_df$date_num)  

#Use geom_raster() 
avg_plot_interp_no_contours_ANPP <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  
  
  #Set the fill color using viridis
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100), na.value = NA, guide = "none") +  
  
  #Customize the x-axise
  scale_x_continuous(
    breaks = seq(10, 90, by = 10),  
    limits = x_limits,  
    expand = c(0, 0)  
  ) +
  
  #Customize the y-axis with pseudocounts
  scale_y_log10(
    name = "Log plastic concentration (g/L)",
    breaks = pseudocount_concentrations,  
    labels = sprintf("%.3f", original_concentrations),  
    limits = c(min(pseudocount_concentrations), max(pseudocount_concentrations)),  # Set y-axis limits
    expand = c(0, 0)  
  ) +
  
  #Apply coord_cartesian without major adjustments
  coord_cartesian(xlim = c(0, max(x_limits))) +  
  
  #Apply the minimal theme
  theme_minimal() +
  labs(title = "Elastollan", x = "Day Number", y = "Plastic Concentration Pseudocount") +
  
  #Customize the theme to show tick marks
  theme(
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    plot.title = element_text(size = 16, hjust = 0.5),  
    
    #Add visible tick marks
    axis.ticks = element_line(color = "black"), 
    axis.ticks.length = unit(0.2, "cm"),  
    
    #Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    panel.border = element_blank()  
  )

#Print
print(avg_plot_interp_no_contours_ANPP)
####Removed rows come from interpolation outside of measured values



###For B Now###
#Generate all combinations of date_num and plastic_conc
all_combinations <- expand.grid(
  date_num = unique(Bsummary_data_contour$date_num),
  plastic_conc = unique(Bsummary_data_contour$plastic_conc)
)

#Find missing combinations
Bmissing_combinations <- anti_join(all_combinations, Bsummary_data_contour, by = c("date_num", "plastic_conc"))

#Interpolate mean and SD for missing combinations
if (nrow(Bmissing_combinations) > 0) {
  interpolated_values <- Bmissing_combinations %>%
    mutate(
      Mean = NA, 
      SD = NA     
    ) %>%
    group_by(plastic_conc) %>%
    mutate(
      Mean = mean(Bsummary_data_contour$Mean, na.rm = TRUE),
      SD = sd(Bsummary_data_contour$Mean, na.rm = TRUE)
    )
  
  #Add interpolated values to the dataframe
  Bsummary_data_contour <- bind_rows(Bsummary_data_contour, interpolated_values)
} else {
  cat("No missing combinations found in the dataframe.\n")
}

#Add a pseudocount of 0.001 to treatment concentrations
Bsummary_data_contour$plastic_conc_pseudocount <- Bsummary_data_contour$plastic_conc + 0.009

#Bsummary_data_contour
ggplot(Bsummary_data_contour, aes(x = date_num, y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  theme_minimal() +
  labs(title = "NPP for Plastic Type A - Tiles Only")

#Correcting the geom_contour usage
###For B###
heatmap_contour_plot_B_log <- ggplot(Bsummary_data_contour, aes(x = date_num, y = plastic_conc_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean), width = 2.9, height = 0.10) +  
  geom_contour(color = "black") +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100)) + 
  scale_x_continuous(name = "Day Number", breaks = seq(min(Bsummary_data_contour$date_num), max(Bsummary_data_contour$date_num), by = 5)) + # Add x-axis tick marks with a step of 5
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  labs(fill = "NPP", title = paste("TPU181"))

#Print
print(heatmap_contour_plot_B_log)

#Day clusters for averaging
day_clusters <- list(
  c(0, 1, 2),
  c(6, 7, 8),
  c(16, 17, 18),
  c(26, 28, 29),
  c(36, 37, 38),
  c(56, 57, 58),
  c(76, 77),
  c(97, 98, 99)
)

#Recalculate averages for each cluster
averaged_data <- do.call(rbind, lapply(day_clusters, function(cluster) {
  subset <- Bsummary_data_contour[Bsummary_data_contour$date_num %in% cluster,]
  
  #Aggregate by plastic_conc_pseudocount to calculate mean
  averaged_subset <- aggregate(Mean ~ plastic_conc_pseudocount, data = subset, FUN = mean, na.rm = TRUE)
  
  #Calculate SD within each plastic_conc_pseudocount group for mean
  sd_subset <- aggregate(Mean ~ plastic_conc_pseudocount, data = subset, FUN = sd, na.rm = TRUE)
  
  #Add SD values to the averaged data frame
  averaged_subset$SD <- sd_subset$Mean
  
  #Assign an averaged date_num for identification
  averaged_subset$date_num <- rep(mean(cluster), nrow(averaged_subset))
  
  return(averaged_subset)
}))

#Convert the result into a data frame
averaged_data <- data.frame(averaged_data)

#date_num as a factor
averaged_data$date_num <- factor(averaged_data$date_num)


#Heatmap for Group 1
g1_plot <- ggplot(subset(Asummary_data_contour, date_num %in% c(0, 6, 16, 26, 36, 56, 76, 97)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 1", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for Group 2
g2_plot <- ggplot(subset(Bsummary_data_contour, date_num %in% c(1, 7, 17, 28, 37, 57, 77, 98)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 2", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for Group 3
g3_plot <- ggplot(subset(Csummary_data_contour, date_num %in% c(2, 8, 18, 29, 38, 58, 76, 99)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 3", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for averaged data
avg_plot <- ggplot(averaged_data, aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Averaged Data", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Display plots
g1_plot
g2_plot
g3_plot
avg_plot

#Make sure date_num is numeric in averaged_data
averaged_data$date_num <- as.numeric(as.character(averaged_data$date_num))

#Perform interpolation
interp <- with(averaged_data, {
  interp(x = date_num, y = log10(plastic_conc_pseudocount), z = Mean,
         xo = seq(min(date_num), max(date_num), length = 100),  # Directly use the numeric range of date_num
         yo = seq(min(log10(plastic_conc_pseudocount)), max(log10(plastic_conc_pseudocount)), length = 100),
         linear = FALSE)
})

#Transform the interpolated data back into a data frame
interp_df <- expand.grid(date_num = interp$x, plastic_conc_pseudocount = interp$y)
interp_df$Mean <- as.vector(interp$z)

#Plotting the heatmap with interpolated data and contour lines
avg_plot_interp <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean)) +  
  geom_contour(bins = 10, color = "black") +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100)) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "TPU 181", x = "Averaged Day Number", y = "Plastic Concentration Pseudocount", fill = "NPP") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

#Print
print(avg_plot_interp)


#Plotting the heatmap with interpolated data without contour lines
avg_plot_interp_no_contours_BNPP <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100)) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    expand = c(0, 0)
  ) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "TPU 181", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "NPP") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_BNPP)


#Make the edges transparent
avg_plot_interp_no_contours_BNPP <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile(na.rm = TRUE) +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100), na.value = NA) +  # Set NA values to be transparent
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    expand = c(0, 0)
  ) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "TPU 181", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "NPP") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_BNPP)

#Trim the edges so it is just within data range
# Determine the limits based on the range of data
x_limits <- range(interp_df$date_num, na.rm = TRUE)
y_limits <- range(10^interp_df$plastic_conc_pseudocount, na.rm = TRUE)

#Plot with adjusted limits
avg_plot_interp_no_contours_BNPP <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100), na.value = NA) +  
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    limits = x_limits,  
    expand = c(0, 0)
  ) +
  scale_y_log10(
    name = "Log Plastic Concentration (g/L)",
    limits = y_limits,  
    expand = c(0, 0)
  ) +
  theme_minimal() +
  labs(title = "TPU 181", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "NPP") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_BNPP)


#Update format for final figures
#Define the original concentrations for plot B with custom format for y-axis
original_concentrations_B <- c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385)

#Add pseudocount to the concentrations
pseudocount_concentrations_B <- original_concentrations_B + 0.009

#Set x_limits to cover the data range
x_limits <- range(interp_df$date_num)  

#Use geom_raster()
avg_plot_interp_no_contours_BNPP <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  
  
  #Set the fill color using viridis
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100), na.value = NA, guide = "none") +  # Remove legend
  
  #Customize the x-axis
  scale_x_continuous(
    breaks = seq(10, 90, by = 10),  
    limits = x_limits,  
    expand = c(0, 0)  
  ) +
  
  #Customize the y-axis with pseudocounts
  scale_y_log10(
    name = "Log plastic concentration (g/L)",
    breaks = pseudocount_concentrations_B,  
    labels = sprintf("%.3f", original_concentrations_B),  
    limits = c(min(pseudocount_concentrations_B), max(pseudocount_concentrations_B)),  
    expand = c(0, 0)  
  ) +
  
  #Apply coord_cartesian
  coord_cartesian(xlim = c(0, max(x_limits))) +  # Ensure we start at 0 with no shifting
  
  #Apply the minimal themees
  theme_minimal() +
  labs(title = "TPU 181", x = "Day Number", y = "Plastic Concentration Pseudocount") +
  
  #Customize
  theme(
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    plot.title = element_text(size = 16, hjust = 0.5),  
    
    #Add visible tick marks
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    
    #Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    panel.border = element_blank()  
  )

#Print
print(avg_plot_interp_no_contours_BNPP)
####Removed rows come from interpolation outside of measured values


###Now for C###
#Generate all combinations of date_num and plastic_conc
all_combinations <- expand.grid(
  date_num = unique(Csummary_data_contour$date_num),
  plastic_conc = unique(Csummary_data_contour$plastic_conc)
)

#Find missing combinations
Cmissing_combinations <- anti_join(all_combinations, Csummary_data_contour, by = c("date_num", "plastic_conc"))

#Interpolate mean and SD for missing combinations
if (nrow(Cmissing_combinations) > 0) {
  interpolated_values <- Cmissing_combinations %>%
    mutate(
      Mean = NA, 
      SD = NA     
    ) %>%
    group_by(plastic_conc) %>%
    mutate(
      Mean = mean(Csummary_data_contour$Mean, na.rm = TRUE),
      SD = sd(Csummary_data_contour$Mean, na.rm = TRUE)
    )
  
  #Add interpolated values to the dataframe
  Csummary_data_contour <- bind_rows(Csummary_data_contour, interpolated_values)
} else {
  cat("No missing combinations found in the dataframe.\n")
}

#Add a pseudocount of 0.001 to treatment concentrations
Csummary_data_contour$plastic_conc_pseudocount <- Csummary_data_contour$plastic_conc + 0.009

#Csummary_data_contour
ggplot(Csummary_data_contour, aes(x = date_num, y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  theme_minimal() +
  labs(title = "NPP for Plastic Type A - Tiles Only")

#Correcting the geom_contour
###For C
heatmap_contour_plot_C_log <- ggplot(Csummary_data_contour, aes(x = date_num, y = plastic_conc_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean), width = 2.9, height = 0.10) +  
  geom_contour(color = "black") +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100)) + 
  scale_x_continuous(name = "Day Number", breaks = seq(min(Csummary_data_contour$date_num), max(Csummary_data_contour$date_num), by = 5)) + # Add x-axis tick marks with a step of 5
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),  
  ) +
  labs(fill = "NPP", title = paste("TPU FC2.1"))

#Print 
print(heatmap_contour_plot_C_log)

#Day clusters for averaging
day_clusters <- list(
  c(0, 1, 2),
  c(6, 7, 8),
  c(16, 17, 18),
  c(26, 28, 29),
  c(36, 37, 38),
  c(56, 57, 58),
  c(76, 77),
  c(97, 98, 99)
)

#Recalculate averages for each cluster
averaged_data <- do.call(rbind, lapply(day_clusters, function(cluster) {
  subset <- Csummary_data_contour[Csummary_data_contour$date_num %in% cluster,]
  
  #Aggregate by plastic_conc_pseudocount to calculate mean
  averaged_subset <- aggregate(Mean ~ plastic_conc_pseudocount, data = subset, FUN = mean, na.rm = TRUE)
  
  #Caculate SD within each plastic_conc_pseudocount group for mean
  sd_subset <- aggregate(Mean ~ plastic_conc_pseudocount, data = subset, FUN = sd, na.rm = TRUE)
  
  #Add SD values to the averaged data frame
  averaged_subset$SD <- sd_subset$Mean
  
  #Assign an averaged date_num for identification
  averaged_subset$date_num <- rep(mean(cluster), nrow(averaged_subset))
  
  return(averaged_subset)
}))

#Convert the result into a data frame
averaged_data <- data.frame(averaged_data)

#date_num as a factor for plotting
averaged_data$date_num <- factor(averaged_data$date_num)


#Heatmap for Group 1
g1_plot <- ggplot(subset(Csummary_data_contour, date_num %in% c(0, 6, 16, 26, 36, 56, 76, 97)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 1", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for Group 2
g2_plot <- ggplot(subset(Csummary_data_contour, date_num %in% c(1, 7, 17, 28, 37, 57, 77, 98)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 2", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for Group 3
g3_plot <- ggplot(subset(Csummary_data_contour, date_num %in% c(2, 8, 18, 29, 38, 58, 76, 99)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 3", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for averaged data
avg_plot <- ggplot(averaged_data, aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Averaged Data", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Display plots
g1_plot
g2_plot
g3_plot
avg_plot

#Make sure date_num is numeric in averaged_data
averaged_data$date_num <- as.numeric(as.character(averaged_data$date_num))

#Interpolate
interp <- with(averaged_data, {
  interp(x = date_num, y = log10(plastic_conc_pseudocount), z = Mean,
         xo = seq(min(date_num), max(date_num), length = 100),  # Directly use the numeric range of date_num
         yo = seq(min(log10(plastic_conc_pseudocount)), max(log10(plastic_conc_pseudocount)), length = 100),
         linear = FALSE)
})

#Transform the interpolated data back into a data frame
interp_df <- expand.grid(date_num = interp$x, plastic_conc_pseudocount = interp$y)
interp_df$Mean <- as.vector(interp$z)

#Plotting the heatmap with interpolated data and contour lines
avg_plot_interp <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean)) +  
  geom_contour(bins = 10, color = "black") +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100)) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "TPU FC2.1", x = "Averaged Day Number", y = "Plastic Concentration Pseudocount", fill = "NPP") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

#Print
print(avg_plot_interp)


#Plotting the heatmap with interpolated data without contour lines
avg_plot_interp_no_contours_CNPP <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100)) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    expand = c(0, 0)
  ) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "TPU FC2.1", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "NPP") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_CNPP)

#Make the edges transparent
avg_plot_interp_no_contours_CNPP <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile(na.rm = TRUE) +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100), na.value = NA) +  # Set NA values to be transparent
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    expand = c(0, 0)
  ) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "TPU FC2.1", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "NPP") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_CNPP)

#Trim he edges so it is just within our data range
#Determine the limits based on the data range
x_limits <- range(interp_df$date_num, na.rm = TRUE)
y_limits <- range(10^interp_df$plastic_conc_pseudocount, na.rm = TRUE)

#Plot
avg_plot_interp_no_contours_CNPP <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100), na.value = NA) +  # Set NA values to be transparent
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    limits = x_limits,  
    expand = c(0, 0)
  ) +
  scale_y_log10(
    name = "Log Plastic Concentration (g/L)",
    limits = y_limits,  
    expand = c(0, 0)
  ) +
  theme_minimal() +
  labs(title = "TPU FC2.1", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "NPP") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_CNPP)



#Update for final fomrat
#Define the original concentrations for Plot C with custom format
original_concentrations_C <- c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385)

#Add pseudocount to the concentrations
pseudocount_concentrations_C <- original_concentrations_C + 0.009

#Set x_limits
x_limits <- range(interp_df$date_num)  # Use the actual range of your data

#Use geom_raster()
avg_plot_interp_no_contours_CNPP <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  
  
  #Set the fill color using viridis
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(0, 100), na.value = NA, guide = "none") +  # Remove legend
  
  #Customize the x axis
  scale_x_continuous(
    breaks = seq(10, 90, by = 10),  
    limits = x_limits,  
    expand = c(0, 0)  
  ) +
  
  #Customize the y-axis with pseudocounts
  scale_y_log10(
    name = "Log plastic concentration (g/L)",
    breaks = pseudocount_concentrations_C,  
    labels = sprintf("%.3f", original_concentrations_C),  
    limits = c(min(pseudocount_concentrations_C), max(pseudocount_concentrations_C)),  # Set y-axis limits
    expand = c(0, 0)  
  ) +
  
  #Apply coord_cartesian without major adjustments to avoid shifting issues
  coord_cartesian(xlim = c(0, max(x_limits))) +  # Ensure we start at 0 with no shifting
  
  #Apply the minimal theme
  theme_minimal() +
  labs(title = "TPU FC2.1", x = "Day Number", y = "Plastic Concentration Pseudocount") +
  
  #Customize
  theme(
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    plot.title = element_text(size = 16, hjust = 0.5),  
    
    #Add visible tick marks
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    
    #Remove grid lines for major ticks
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    panel.border = element_blank()  
  )

#Print
print(avg_plot_interp_no_contours_CNPP)
####Removed rows come from interpolation outside of measured values

#Now combine the 3 plots
#Remove x-axis title from ANPP and CNPP, and y-axis labels/titles/tickmarks from BNPP and CNPP
avg_plot_interp_no_contours_ANPP <- avg_plot_interp_no_contours_ANPP + 
  theme(axis.title.x = element_blank())  

avg_plot_interp_no_contours_BNPP <- avg_plot_interp_no_contours_BNPP + 
  theme(axis.title.y = element_blank(),  
        axis.text.y = element_blank(),   
        axis.ticks.y = element_blank())  

avg_plot_interp_no_contours_CNPP <- avg_plot_interp_no_contours_CNPP + 
  theme(axis.title.y = element_blank(),  
        axis.text.y = element_blank(),   
        axis.ticks.y = element_blank(),  
        axis.title.x = element_blank())  

#Combine plots
combined_plot <- avg_plot_interp_no_contours_ANPP + 
  avg_plot_interp_no_contours_BNPP + 
  avg_plot_interp_no_contours_CNPP + 
  plot_layout(ncol = 3)

#Add overall title
final_combined_plot <- combined_plot + 
  plot_annotation(title = "Net primary production", 
                  theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")))

#Display combined plot
print(final_combined_plot)


#Format update for the plot
#Modify individual plot ANPP
avg_plot_interp_no_contours_ANPP <- avg_plot_interp_no_contours_ANPP + 
  theme(
    axis.title.y = element_text(size = 14, face = "bold"),  
    axis.title.x = element_blank()  
  )

#Modify individual plot BNPP
avg_plot_interp_no_contours_BNPP <- avg_plot_interp_no_contours_BNPP + 
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),  
    axis.title.y = element_blank(),  
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank()   
  )

#Modify individual plot CNPP
avg_plot_interp_no_contours_CNPP <- avg_plot_interp_no_contours_CNPP + 
  theme(
    axis.title.y = element_blank(),  
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank(),  
    axis.title.x = element_blank()   
  )

#Combine plots
combined_plot <- avg_plot_interp_no_contours_ANPP + 
  avg_plot_interp_no_contours_BNPP + 
  avg_plot_interp_no_contours_CNPP + 
  plot_layout(ncol = 3)

#Add overall title
final_combined_plot <- combined_plot

#Display final combined plot ### NPP data for MANUSCRIPT [Figure S]
print(final_combined_plot)
####Removed rows come from interpolation outside of measured values



######NOW RUN THE GAMS###########
#Start the GAMs analysis
#npp_data
#Convert date to a numeric scale (number of days since the start of the study)
npp_data$date <- as.Date(npp_data$date, format = "%m/%d/%Y")
start_date <- min(npp_data$date)
npp_data$date_num <- as.numeric(npp_data$date - start_date)

#plastic_conc is a factor and plastic_conc is numeric
npp_data$plastic_type <- as.factor(npp_data$plastic_type)
npp_data$plastic_conc <- as.numeric(npp_data$plastic_conc)

#Change date_num values to actual concentrations
npp_data <- npp_data %>%
  mutate(plastic_conc = case_when(
    plastic_conc == 0 ~ 0.000,
    plastic_conc == 1 ~ 0.004,
    plastic_conc == 2 ~ 0.008,
    plastic_conc == 3 ~ 0.013,
    plastic_conc == 4 ~ 0.023,
    plastic_conc == 5 ~ 0.041,
    plastic_conc == 6 ~ 0.072,
    plastic_conc == 7 ~ 0.126,
    plastic_conc == 8 ~ 0.220,
    plastic_conc == 9 ~ 0.385,
    TRUE ~ plastic_conc  
  ))

#Check data distribution to meet assumptions
#Histogram of response variable
hist(npp_data$NPP, main="Histogram of NPP", xlab="NPP")
#Density plot for smoother
plot(density(na.omit(npp_data$NPP)), main="Density Plot of NPP")
#Plastic type
ggplot(npp_data, aes(x=plastic_type, y=NPP)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess", color="blue") +
  theme_minimal() +
  labs(title="Relationship between plastic type and NPP", x="Predictor", y="NPP")
#Day number
ggplot(npp_data, aes(x=date_num, y=NPP)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess", color="blue") +
  theme_minimal() +
  labs(title="Relationship between date_num and NPP", x="Predictor", y="NPP")
#Plastic concentration
ggplot(npp_data, aes(x=plastic_conc, y=NPP)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess", color="blue") +
  theme_minimal() +
  labs(title="Relationship between plastic concentration and NPP", x="Predictor", y="NPP")
#Check auto correlations
acf(npp_data$date_num)
table(npp_data$plastic_type)

#Fix the non-normal data issue by setting family distribution (Not this one)
gam_model5 <- gam(NPP ~ s(plastic_conc, by = plastic_type) + s(date_num, by = plastic_type) + plastic_type, family = Gamma(link = "log"), data = npp_data)
summary(gam_model5)
plot(gam_model5)
check_result5 <- gam.check(gam_model5)

#Fix the non-normal data issue by transforming the data (This one is best)
gam_model_transformed <- gam(I(log(NPP + 1)) ~ s(plastic_conc, by = plastic_type) + s(date_num, by = plastic_type) + plastic_type, data = npp_data)
summary(gam_model_transformed)
check_result_transformed <- gam.check(gam_model_transformed)

#From the gam.check we can increase the k (Not this One)
gam_model_adjusted <- gam(I(log(NPP + 1)) ~ s(plastic_conc, by = plastic_type, k = 10) + 
                            s(date_num, by = plastic_type, k = 10) + plastic_type, 
                          data = npp_data)
summary(gam_model_adjusted)
check_result_adjusted <- gam.check(gam_model_adjusted)

#Fit the GAM model (Not this one)
gam_model6 <- gam(NPP ~ s(plastic_conc, by = plastic_type) + s(date_num, by = plastic_type) + plastic_type, data = npp_data)
summary(gam_model6)
plot(gam_model6)
check_result6 <- gam.check(gam_model6)


#Compare models
summary(gam_model_transformed)
gam.check(gam_model_transformed)

summary(gam_model_adjusted)
gam.check(gam_model_adjusted)

summary(gam_model6)
gam.check(gam_model6)

#Compare AIC values
AIC(gam_model_transformed, gam_model_adjusted, gam_model6)

#Visualize smooths
par(mfrow = c(1, 2)) 
plot(gam_model_transformed, pages = 1)
plot(gam_model_adjusted, pages = 1)

#transformed data fits best and is appropriate, build more complex models

NPPgam_model1 <- gam(I(log(NPP + 1)) ~ 
                    s(date_num) + 
                    s(plastic_conc), 
                  data = npp_data)


NPPgam_model2 <- gam(I(log(NPP + 1)) ~ 
                    s(date_num) + 
                    s(plastic_conc) + 
                    plastic_type, 
                  data = npp_data)


NPPgam_model3 <- gam(I(log(NPP + 1)) ~ 
                    s(date_num, by = plastic_type) + 
                    s(plastic_conc, by = plastic_type), 
                  data = npp_data)


NPPgam_model4 <- gam(I(log(NPP + 1)) ~ 
                    s(date_num, by = plastic_type) + 
                    s(plastic_conc, by = plastic_type) + 
                    plastic_type, 
                  data = npp_data)

#Fit models of increasing complexity
NPPgam_model1 <- gam(I(log(NPP + 1)) ~ s(date_num) + s(plastic_conc), data = npp_data)
NPPgam_model2 <- gam(I(log(NPP + 1)) ~ s(date_num) + s(plastic_conc) + plastic_type, data = npp_data)
NPPgam_model3 <- gam(I(log(NPP + 1)) ~ s(date_num, by = plastic_type) + s(plastic_conc, by = plastic_type), data = npp_data)
NPPgam_model4 <- gam(I(log(NPP + 1)) ~ s(date_num, by = plastic_type) + s(plastic_conc, by = plastic_type) + plastic_type, data = npp_data)

#Compare models using AIC
NPPmodel_comparison <- AIC(NPPgam_model1, NPPgam_model2, NPPgam_model3, NPPgam_model4)
print(NPPmodel_comparison) ###NPP data used for MANUSCRIPT [TABLE S1]

#Check diagnostics for the best model
NPPbest_model <- NPPgam_model4
summary(NPPbest_model) ###NPP data used for MANUSCRIPT [Table S2]
gam.check(NPPbest_model)
plot(NPPbest_model, pages = 1)

#plot_difference
#Comparison between plastic types: Group A and Group B
npp_diff_AB_date <- plot_difference(
  NPPbest_model,
  series = date_num,
  difference = list(plastic_type = c("A", "B"))
) +
  labs(x = "Day Number", y = "log NPP (Difference smooth)") +
  ggtitle("Difference between Elastollan and TPU 181") +
  theme_minimal()


#Comparison between plastic types - Group A and Group C
npp_diff_AC_date <- plot_difference(
  NPPbest_model,
  series = date_num,
  difference = list(plastic_type = c("A", "C"))
) +
  labs(x = "Day Number", y = "log NPP (Difference smooth)") +
  ggtitle("Difference between Elastollan and TPU FC2.1") +
  theme_minimal()


#Comparison between plastic types - Group B and Group C
npp_diff_BC_date <- plot_difference(
  NPPbest_model,
  series = date_num,
  difference = list(plastic_type = c("B", "C"))
) +
  labs(x = "Day Number", y = "log NPP (Difference smooth)") +
  ggtitle("Difference between TPU 181 and TPU FC2.1") +
  theme_minimal()



#Comparison between plastic types - Group A and Group B
npp_diff_AB_conc <- plot_difference(
  NPPbest_model,
  series = plastic_conc,
  difference = list(plastic_type = c("A", "B"))
) +
  labs(x = "Plastic Concentration", y = "log NPP (Difference smooth)") +
  ggtitle("Difference between Elastollan and TPU 181") +
  theme_minimal()


#Comparison between plastic type - Group A and Group C
npp_diff_AC_conc <- plot_difference(
  NPPbest_model,
  series = plastic_conc,
  difference = list(plastic_type = c("A", "C"))
) +
  labs(x = "Plastic Concentration", y = "log NPP (Difference smooth)") +
  ggtitle("Difference between Elastollan and TPU FC2.1") +
  theme_minimal()


#Comparison between plastic types - Group B and Group C
npp_diff_BC_conc <- plot_difference(
  NPPbest_model,
  series = plastic_conc,
  difference = list(plastic_type = c("B", "C"))
) +
  labs(x = "Plastic Concentration", y = "log NPP (Difference smooth)") +
  ggtitle("Difference between TPU181 and TPU FC2.1") +
  theme_minimal()


print(npp_diff_AB_date)
print(npp_diff_AC_date)
print(npp_diff_BC_date)
print(npp_diff_AB_conc)
print(npp_diff_AC_conc)
print(npp_diff_BC_conc)
### NPP data used for MANUSCRIPT [Figure S10]



##Now for ER data
##Try to get contour lines
Aselected_letter_group <- "A"
Bselected_letter_group <- "B"
Cselected_letter_group <- "C"

#Filter the data for the selected letter group
Aselected_data <- filter(er_data, plastic_type == Aselected_letter_group)
Bselected_data <- filter(er_data, plastic_type == Bselected_letter_group)
Cselected_data <- filter(er_data, plastic_type == Cselected_letter_group)

### Calculate the mean and standard deviation for ER ###
##Plastic type A
Asummary_data_contour_ER <- Aselected_data %>%
  group_by(date, plastic_conc) %>%
  summarise(
    Mean = mean(ER),
    SD = sd(ER),
    .groups = 'drop'
  )

##Plastic type B
Bsummary_data_contour_ER <- Bselected_data %>%
  group_by(date, plastic_conc) %>%
  summarise(
    Mean = mean(ER),
    SD = sd(ER),
    .groups = 'drop'
  )

##Plastic type C
Csummary_data_contour_ER <- Cselected_data %>%
  group_by(date, plastic_conc) %>%
  summarise(
    Mean = mean(ER),
    SD = sd(ER),
    .groups = 'drop'
  )


#Format date without year and store it in a new variable
Asummary_data_contour_ER$date_without_year <- format(as.Date(Asummary_data_contour_ER$date, format = "%m/%d/%Y"), "%m/%d")

Bsummary_data_contour_ER$date_without_year <- format(as.Date(Bsummary_data_contour_ER$date, format = "%m/%d/%Y"), "%m/%d")

Csummary_data_contour_ER$date_without_year <- format(as.Date(Csummary_data_contour_ER$date, format = "%m/%d/%Y"), "%m/%d")



#Change plastic_conc values to actual concentrations
###For A###
Asummary_data_contour_ER <- Asummary_data_contour_ER %>%
  mutate(plastic_conc = case_when(
    plastic_conc == 0 ~ 0.000,
    plastic_conc == 1 ~ 0.004,
    plastic_conc == 2 ~ 0.008,
    plastic_conc == 3 ~ 0.013,
    plastic_conc == 4 ~ 0.023,
    plastic_conc == 5 ~ 0.041,
    plastic_conc == 6 ~ 0.072,
    plastic_conc == 7 ~ 0.126,
    plastic_conc == 8 ~ 0.220,
    plastic_conc == 9 ~ 0.385,
    TRUE ~ plastic_conc 
  ))

###For B###
Bsummary_data_contour_ER <- Bsummary_data_contour_ER %>%
  mutate(plastic_conc = case_when(
    plastic_conc == 0 ~ 0.000,
    plastic_conc == 1 ~ 0.004,
    plastic_conc == 2 ~ 0.008,
    plastic_conc == 3 ~ 0.013,
    plastic_conc == 4 ~ 0.023,
    plastic_conc == 5 ~ 0.041,
    plastic_conc == 6 ~ 0.072,
    plastic_conc == 7 ~ 0.126,
    plastic_conc == 8 ~ 0.220,
    plastic_conc == 9 ~ 0.385,
    TRUE ~ plastic_conc  
  ))

###For C###
Csummary_data_contour_ER <- Csummary_data_contour_ER %>%
  mutate(plastic_conc = case_when(
    plastic_conc == 0 ~ 0.000,
    plastic_conc == 1 ~ 0.004,
    plastic_conc == 2 ~ 0.008,
    plastic_conc == 3 ~ 0.013,
    plastic_conc == 4 ~ 0.023,
    plastic_conc == 5 ~ 0.041,
    plastic_conc == 6 ~ 0.072,
    plastic_conc == 7 ~ 0.126,
    plastic_conc == 8 ~ 0.220,
    plastic_conc == 9 ~ 0.385,
    TRUE ~ plastic_conc  
  ))


#Convert date to a date object and then to a numeric scale
###For A###
Asummary_data_contour_ER$date <- as.Date(Asummary_data_contour_ER$date, format = "%m/%d/%Y")
start_date <- min(Asummary_data_contour_ER$date)
Asummary_data_contour_ER$date_num <- as.numeric(Asummary_data_contour_ER$date - start_date)

Asummary_data_contour_ER$plastic_conc <- as.numeric(Asummary_data_contour_ER$plastic_conc)

###For B###
Bsummary_data_contour_ER$date <- as.Date(Bsummary_data_contour_ER$date, format = "%m/%d/%Y")
start_date <- min(Bsummary_data_contour_ER$date)
Bsummary_data_contour_ER$date_num <- as.numeric(Bsummary_data_contour_ER$date - start_date)

Bsummary_data_contour_ER$plastic_conc <- as.numeric(Bsummary_data_contour_ER$plastic_conc)


###For C###
Csummary_data_contour_ER$date <- as.Date(Csummary_data_contour_ER$date, format = "%m/%d/%Y")
start_date <- min(Csummary_data_contour_ER$date)
Csummary_data_contour_ER$date_num <- as.numeric(Csummary_data_contour_ER$date - start_date)

Csummary_data_contour_ER$plastic_conc <- as.numeric(Csummary_data_contour_ER$plastic_conc)

#Get date without year
###For A###
Asummary_data_contour_ER <- select(Asummary_data_contour_ER, -date_without_year)
Asummary_data_contour_ER <- select(Asummary_data_contour_ER, -date)

###For B###
Bsummary_data_contour_ER <- select(Bsummary_data_contour_ER, -date_without_year)
Bsummary_data_contour_ER <- select(Bsummary_data_contour_ER, -date)

###For C###
Csummary_data_contour_ER <- select(Csummary_data_contour_ER, -date_without_year)
Csummary_data_contour_ER <- select(Csummary_data_contour_ER, -date)

##For A drop row 78 as it is a true and confirmed error in the sensor prob reading##
#Remove row 78 from Asummary_data_contour_ER
Asummary_data_contour_ER <- Asummary_data_contour_ER[-78, ]

###For A###
#Generate all combinations of date_num and plastic_conc
all_combinations <- expand.grid(
  date_num = unique(Asummary_data_contour_ER$date_num),
  plastic_conc = unique(Asummary_data_contour_ER$plastic_conc)
)

#Find missing combinations
Amissing_combinations <- anti_join(all_combinations, Asummary_data_contour_ER, by = c("date_num", "plastic_conc"))

#Interpolate
if (nrow(Amissing_combinations) > 0) {
  interpolated_values <- Amissing_combinations %>%
    mutate(
      Mean = NA, 
      SD = NA     
    ) %>%
    group_by(plastic_conc) %>%
    mutate(
      Mean = mean(Asummary_data_contour_ER$Mean, na.rm = TRUE),
      SD = sd(Asummary_data_contour_ER$Mean, na.rm = TRUE)
    )
  
  #Add interpolated values to the dataframe
  Asummary_data_contour_ER <- bind_rows(Asummary_data_contour_ER, interpolated_values)
} else {
  cat("No missing combinations found in the dataframe.\n")
}

#Add a pseudocount of 0.001 to treatment concentrations
Asummary_data_contour_ER$plastic_conc_pseudocount <- Asummary_data_contour_ER$plastic_conc + 0.009

#Asummary_data_contour
ggplot(Asummary_data_contour_ER, aes(x = date_num, y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  theme_minimal() +
  labs(title = "ER for Plastic Type A - Tiles Only")

#Correcting the geom_contour usage so that y is on a log axis
###For A###
heatmap_contour_plot_A_log <- ggplot(Asummary_data_contour_ER, aes(x = date_num, y = plastic_conc_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean), width = 2.9, height = 0.10) +  # Adjust width and height parameters to minimize gaps
  geom_contour(color = "black") +  # This adds the contour lines based on `z = Mean`
  scale_fill_viridis_c(option = "plasma", direction = -1) + # Use 'plasma' palette and reverse gradient
  scale_x_continuous(name = "Day Number", breaks = seq(min(Asummary_data_contour_ER$date_num), max(Asummary_data_contour_ER$date_num), by = 5)) + # Add x-axis tick marks with a step of 5
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +  # Set y-axis to log scale
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  labs(fill = "ER", title = paste("Elastollan"))

#Print 
print(heatmap_contour_plot_A_log)



#Day clusters for averaging
day_clusters <- list(
  c(0, 1),
  c(6, 7),
  c(16, 17),
  c(26, 28),
  c(36, 37),
  c(56, 57),
  c(76),
  c(97, 98)
)

#Recalcualte averages for each cluster
averaged_data <- do.call(rbind, lapply(day_clusters, function(cluster) {
  subset <- Asummary_data_contour_ER[Asummary_data_contour_ER$date_num %in% cluster,]
  
  #Aggregate by plastic_conc_pseudocount to calculate mean
  averaged_subset <- aggregate(Mean ~ plastic_conc_pseudocount, data = subset, FUN = mean, na.rm = TRUE)
  
  #Compute SD within each plastic_conc_pseudocount group for mean
  sd_subset <- aggregate(Mean ~ plastic_conc_pseudocount, data = subset, FUN = sd, na.rm = TRUE)
  
  #Add SD values to the averaged data frame
  averaged_subset$SD <- sd_subset$Mean
  
  #Assign an averaged date_num for identification
  averaged_subset$date_num <- rep(mean(cluster), nrow(averaged_subset))
  
  return(averaged_subset)
}))

#Convert the result into a data frame 
averaged_data <- data.frame(averaged_data)

#date_num as a factor
averaged_data$date_num <- factor(averaged_data$date_num)


#Heatmap for Group 1
g1_plot <- ggplot(subset(Asummary_data_contour_ER, date_num %in% c(0, 6, 16, 26, 36, 56, 76, 97)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 1 ER", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for Group 2
g2_plot <- ggplot(subset(Asummary_data_contour_ER, date_num %in% c(1, 7, 17, 28, 37, 57, 76, 98)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 2 ER", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for averaged data
avg_plot <- ggplot(averaged_data, aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Averaged Data ER", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Display plots
g1_plot
g2_plot
avg_plot

#Make sure date_num is numeric in averaged_data
averaged_data$date_num <- as.numeric(as.character(averaged_data$date_num))

#Interpolate
interp <- with(averaged_data, {
  interp(x = date_num, y = log10(plastic_conc_pseudocount), z = Mean,
         xo = seq(min(date_num), max(date_num), length = 100),  # Directly use the numeric range of date_num
         yo = seq(min(log10(plastic_conc_pseudocount)), max(log10(plastic_conc_pseudocount)), length = 100),
         linear = FALSE)
})

#Transform the interpolated data back into a data frame
interp_df <- expand.grid(date_num = interp$x, plastic_conc_pseudocount = interp$y)
interp_df$Mean <- as.vector(interp$z)

#Plotting the heatmap with interpolated data
avg_plot_interp <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean)) +  # Continuous heatmap representation
  geom_contour(bins = 10, color = "black") +  # Add contour lines
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0)) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "ER Elastollan", x = "Averaged Day Number", y = "Plastic Concentration Pseudocount", fill = "ER") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

#Print 
print(avg_plot_interp)


#Plotting the heatmap with interpolated data without contour lines
avg_plot_interp_no_contours_AER <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  # Continuous heatmap representation
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0)) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  # Manually specify breaks starting at 0, ending at 90, at intervals of 10
    expand = c(0, 0)
  ) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "Elastollan", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "ER") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_AER)


#Update
avg_plot_interp_no_contours_AER <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0), na.value = NA) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    expand = c(0, 0)
  ) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "Elastollan", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "ER") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_AER)

#Determine the limits based on the range of data
x_limits <- range(interp_df$date_num, na.rm = TRUE)
y_limits <- range(10^interp_df$plastic_conc_pseudocount, na.rm = TRUE)

#Adjust the plot
avg_plot_interp_no_contours_AER <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0), na.value = NA) +  # Set NA values to be transparent
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    limits = x_limits,  
    expand = c(0, 0)
  ) +
  scale_y_log10(
    name = "Log Plastic Concentration (g/L)",
    limits = y_limits,  
    expand = c(0, 0)
  ) +
  theme_minimal() +
  labs(title = "Elastollan", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "ER") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_AER)

#Format changes to figure
#Define the original concentrations for plot A
original_concentrations_A <- c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385)

#Add pseudocount to the concentrations
pseudocount_concentrations_A <- original_concentrations_A + 0.009

#Set x_limits to cover the data range
x_limits <- range(interp_df$date_num)  

#Use geom_raster()
avg_plot_interp_no_contours_AER <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  
  
  #Set the fill color using viridis 
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0), na.value = NA, guide = "none") +  # Remove legend
  
  #Customize the x-axis
  scale_x_continuous(
    breaks = seq(10, 90, by = 10),  
    limits = x_limits,  
    expand = c(0, 0)  
  ) +
  
  #Customize the y-axis with pseudocounts
  scale_y_log10(
    name = "Log plastic concentration (g/L)",
    breaks = pseudocount_concentrations_A,  
    labels = sprintf("%.3f", original_concentrations_A),  
    limits = c(min(pseudocount_concentrations_A), max(pseudocount_concentrations_A)),  # Set y-axis limits
    expand = c(0, 0)  
  ) +
  
  #Apply coord_cartesian
  coord_cartesian(xlim = c(0, max(x_limits))) +  
  
  #Apply minimal themes
  theme_minimal() +
  labs(title = "Elastollan", x = "Day Number", y = "Plastic Concentration Pseudocount") +
  
  #Customizes
  theme(
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    plot.title = element_text(size = 16, hjust = 0.5),  
    
    #Add tick marks
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    
    #Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    panel.border = element_blank()  
  )

#Print
print(avg_plot_interp_no_contours_AER)




###For ER B###
###For B###
# Generate all combinations of date_num and plastic_conc
all_combinations <- expand.grid(
  date_num = unique(Bsummary_data_contour_ER$date_num),
  plastic_conc = unique(Bsummary_data_contour_ER$plastic_conc)
)

#Find missing combinations
Bmissing_combinations <- anti_join(all_combinations, Bsummary_data_contour_ER, by = c("date_num", "plastic_conc"))

#Interpolate mean and SD for missing combinations
if (nrow(Bmissing_combinations) > 0) {
  interpolated_values <- Bmissing_combinations %>%
    mutate(
      Mean = NA, 
      SD = NA    
    ) %>%
    group_by(plastic_conc) %>%
    mutate(
      Mean = mean(Bsummary_data_contour_ER$Mean, na.rm = TRUE),
      SD = sd(Bsummary_data_contour_ER$Mean, na.rm = TRUE)
    )
  
  #Add interpolated values to the dataframe
  Bsummary_data_contour_ER <- bind_rows(Bsummary_data_contour_ER, interpolated_values)
} else {
  cat("No missing combinations found in the dataframe.\n")
}

#Add a pseudocount of 0.001 to treatment concentrations
Bsummary_data_contour_ER$plastic_conc_pseudocount <- Bsummary_data_contour_ER$plastic_conc + 0.009

#Asummary_data_contour
ggplot(Bsummary_data_contour_ER, aes(x = date_num, y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  theme_minimal() +
  labs(title = "ER for Plastic Type B - Tiles Only")

#Correcting the geom_contour usage so that y is on a log axis
###For B###
heatmap_contour_plot_B_log <- ggplot(Bsummary_data_contour_ER, aes(x = date_num, y = plastic_conc_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean), width = 2.9, height = 0.10) +  
  geom_contour(color = "black") +  
  scale_fill_viridis_c(option = "plasma", direction = -1) + 
  scale_x_continuous(name = "Day Number", breaks = seq(min(Bsummary_data_contour_ER$date_num), max(Bsummary_data_contour_ER$date_num), by = 5)) + 
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  labs(fill = "ER", title = paste("TPU 181"))

#Print
print(heatmap_contour_plot_B_log)



#Day clusters for averaging
day_clusters <- list(
  c(0, 1),
  c(6, 7),
  c(16, 17),
  c(26, 28),
  c(36, 37),
  c(56, 57),
  c(76),
  c(97, 98)
)

#Recalculate averages for each cluster
averaged_data <- do.call(rbind, lapply(day_clusters, function(cluster) {
  subset <- Bsummary_data_contour_ER[Bsummary_data_contour_ER$date_num %in% cluster,]
  
  #Aggregate by plastic_conc_pseudocount to calculate mean
  averaged_subset <- aggregate(Mean ~ plastic_conc_pseudocount, data = subset, FUN = mean, na.rm = TRUE)
  
  #Calculate SD within each plastic_conc_pseudocount group
  sd_subset <- aggregate(Mean ~ plastic_conc_pseudocount, data = subset, FUN = sd, na.rm = TRUE)
  
  #Add SD values to the averaged data frame
  averaged_subset$SD <- sd_subset$Mean
  
  #Assign an averaged date_num for identification 
  averaged_subset$date_num <- rep(mean(cluster), nrow(averaged_subset))
  
  return(averaged_subset)
}))

#Convert the result into a data frame
averaged_data <- data.frame(averaged_data)

#Treating date_num as a factor
averaged_data$date_num <- factor(averaged_data$date_num)


#Heatmap for Group 1
g1_plot <- ggplot(subset(Bsummary_data_contour_ER, date_num %in% c(0, 6, 16, 26, 36, 56, 76, 97)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 1 ER", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for Group 2
g2_plot <- ggplot(subset(Bsummary_data_contour_ER, date_num %in% c(1, 7, 17, 28, 37, 57, 76, 98)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 2 ER", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for Averaged Data
avg_plot <- ggplot(averaged_data, aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Averaged Data ER", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Display plots
g1_plot
g2_plot
avg_plot

#Make 'date_num' is numeric in 'averaged_data'
averaged_data$date_num <- as.numeric(as.character(averaged_data$date_num))

#Interpolate
interp <- with(averaged_data, {
  interp(x = date_num, y = log10(plastic_conc_pseudocount), z = Mean,
         xo = seq(min(date_num), max(date_num), length = 100),  
         yo = seq(min(log10(plastic_conc_pseudocount)), max(log10(plastic_conc_pseudocount)), length = 100),
         linear = FALSE)
})

#Transform the interpolated data back into a data frame
interp_df <- expand.grid(date_num = interp$x, plastic_conc_pseudocount = interp$y)
interp_df$Mean <- as.vector(interp$z)

#Plot the heatmap with interpolated data
avg_plot_interp <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean)) +  
  geom_contour(bins = 10, color = "black") +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0)) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "TPU 181", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "ER") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

#Print 
print(avg_plot_interp)


#Plot the heatmap with interpolated data
avg_plot_interp_no_contours_BER <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0)) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    expand = c(0, 0)
  ) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "TPU 181", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "ER") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_BER)


#Make gray boxes transparent
avg_plot_interp_no_contours_BER <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0), na.value = NA) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    expand = c(0, 0)
  ) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "TPU 181", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "ER") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_BER)

#Determine the limits based on the range of data
x_limits <- range(interp_df$date_num, na.rm = TRUE)
y_limits <- range(10^interp_df$plastic_conc_pseudocount, na.rm = TRUE)

#Adjust plot
avg_plot_interp_no_contours_BER <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0), na.value = NA) +  # Set NA values to be transparent
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    limits = x_limits,  
    expand = c(0, 0)
  ) +
  scale_y_log10(
    name = "Log Plastic Concentration (g/L)",
    limits = y_limits,  
    expand = c(0, 0)
  ) +
  theme_minimal() +
  labs(title = "TPU 181", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "ER") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_BER)

#Figure format changes
#Define the original concentrations for Plot B
original_concentrations_B <- c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385)

#Add pseudocount to the concentrations
pseudocount_concentrations_B <- original_concentrations_B + 0.009

#Set x_limits to cover the data range
x_limits <- range(interp_df$date_num)  # Use the actual range of your data

#geom_raster()
avg_plot_interp_no_contours_BER <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_raster(interpolate = TRUE) + 
  
  #Set the fill color using viridis 
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0), na.value = NA, guide = "none") +  # Remove legend
  
  #Customize the x-axis 
  scale_x_continuous(
    breaks = seq(10, 90, by = 10), 
    limits = x_limits,  
    expand = c(0, 0)  
  ) +
  
  #Customize the y-axis with pseudocounts
  scale_y_log10(
    name = "Log plastic concentration (g/L)",
    breaks = pseudocount_concentrations_B,  
    labels = sprintf("%.3f", original_concentrations_B),  
    limits = c(min(pseudocount_concentrations_B), max(pseudocount_concentrations_B)),  
    expand = c(0, 0)  
  ) +
  
  #Apply coord_cartesian
  coord_cartesian(xlim = c(0, max(x_limits))) +  # Ensure we start at 0 with no shifting
  
  #Apply the minimal theme
  theme_minimal() +
  labs(title = "TPU 181", x = "Day Number", y = "Plastic Concentration Pseudocount") +
  
  #Customize the theme 
  theme(
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    plot.title = element_text(size = 16, hjust = 0.5),  
    
    #Add visible tick marks (excluding 0)
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    
    #Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    panel.border = element_blank()  
  )

#Print
print(avg_plot_interp_no_contours_BER)






####For C####
###For ER C###
#Generate all combinations of date_num and plastic_conc
all_combinations <- expand.grid(
  date_num = unique(Csummary_data_contour_ER$date_num),
  plastic_conc = unique(Csummary_data_contour_ER$plastic_conc)
)

#Find missing combinations
Cmissing_combinations <- anti_join(all_combinations, Csummary_data_contour_ER, by = c("date_num", "plastic_conc"))

#Interpolate mean and SD for missing combinations
if (nrow(Cmissing_combinations) > 0) {
  interpolated_values <- Cmissing_combinations %>%
    mutate(
      Mean = NA, 
      SD = NA     
    ) %>%
    group_by(plastic_conc) %>%
    mutate(
      Mean = mean(Csummary_data_contour_ER$Mean, na.rm = TRUE),
      SD = sd(Csummary_data_contour_ER$Mean, na.rm = TRUE)
    )
  
  #Add interpolated values to dataframe
  Csummary_data_contour_ER <- bind_rows(Csummary_data_contour_ER, interpolated_values)
} else {
  cat("No missing combinations found in the dataframe.\n")
}

#Add a pseudocount of 0.001 to treatment concentrations
Csummary_data_contour_ER$plastic_conc_pseudocount <- Csummary_data_contour_ER$plastic_conc + 0.009

#Csummary_data_contour
ggplot(Csummary_data_contour_ER, aes(x = date_num, y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  theme_minimal() +
  labs(title = "ER for Plastic Type A - Tiles Only")

#Correcting the geom_contour usage
###For C###
heatmap_contour_plot_C_log <- ggplot(Csummary_data_contour_ER, aes(x = date_num, y = plastic_conc_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean), width = 2.9, height = 0.10) +  
  geom_contour(color = "black") +  
  scale_fill_viridis_c(option = "plasma", direction = -1) + 
  scale_x_continuous(name = "Day Number", breaks = seq(min(Csummary_data_contour_ER$date_num), max(Csummary_data_contour_ER$date_num), by = 5)) + # Add x-axis tick marks with a step of 5
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  labs(fill = "ER", title = paste("TPU FC2.1"))

#Print 
print(heatmap_contour_plot_C_log)



#Day clusters for averaging
day_clusters <- list(
  c(0, 1),
  c(6, 7),
  c(16, 17),
  c(26, 28),
  c(36, 37),
  c(56, 57),
  c(76),
  c(97, 98)
)

#Recalculate averages for each cluster
averaged_data <- do.call(rbind, lapply(day_clusters, function(cluster) {
  subset <- Csummary_data_contour_ER[Csummary_data_contour_ER$date_num %in% cluster,]
  
  #Aggregate by plastic_conc_pseudocount to calculate mean
  averaged_subset <- aggregate(Mean ~ plastic_conc_pseudocount, data = subset, FUN = mean, na.rm = TRUE)
  
  #Compute SD within each plastic_conc_pseudocount group
  sd_subset <- aggregate(Mean ~ plastic_conc_pseudocount, data = subset, FUN = sd, na.rm = TRUE)
  
  #Add SD values to the averaged data frame
  averaged_subset$SD <- sd_subset$Mean
  
  #Assign an averaged date_num
  averaged_subset$date_num <- rep(mean(cluster), nrow(averaged_subset))
  
  return(averaged_subset)
}))

#Convert the result into a data frame
averaged_data <- data.frame(averaged_data)

#date_num as a factor
averaged_data$date_num <- factor(averaged_data$date_num)


#Heatmap for Group 1
g1_plot <- ggplot(subset(Csummary_data_contour_ER, date_num %in% c(0, 6, 16, 26, 36, 56, 76, 97)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 1 ER", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for Group 2
g2_plot <- ggplot(subset(Csummary_data_contour_ER, date_num %in% c(1, 7, 17, 28, 37, 57, 76, 98)),
                  aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Group 2 ER", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Heatmap for Averaged Data
avg_plot <- ggplot(averaged_data, aes(x = factor(date_num), y = plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap for Averaged Data ER", x = "Day Number", y = "Plastic Concentration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Display plots
g1_plot
g2_plot
avg_plot

#Make sure date_num is numeric in averaged_data
averaged_data$date_num <- as.numeric(as.character(averaged_data$date_num))

#Interpolate
interp <- with(averaged_data, {
  interp(x = date_num, y = log10(plastic_conc_pseudocount), z = Mean,
         xo = seq(min(date_num), max(date_num), length = 100),  
         yo = seq(min(log10(plastic_conc_pseudocount)), max(log10(plastic_conc_pseudocount)), length = 100),
         linear = FALSE)
})

#Transform the interpolated data back into a data frame
interp_df <- expand.grid(date_num = interp$x, plastic_conc_pseudocount = interp$y)
interp_df$Mean <- as.vector(interp$z)

#Plotting the heatmap with interpolated data
avg_plot_interp <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean)) +  
  geom_contour(bins = 10, color = "black") +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0)) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "ER TPU FC2.1", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "ER") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

#Print
print(avg_plot_interp)


#Plotting the heatmap with interpolated data
avg_plot_interp_no_contours_CER <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0)) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    expand = c(0, 0)
  ) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "TPU FC2.1", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "ER") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold")
  )

print(avg_plot_interp_no_contours_CER)


#Make gray boxes transparent
avg_plot_interp_no_contours_CER <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0), na.value = NA) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    expand = c(0, 0)
  ) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +
  theme_minimal() +
  labs(title = "TPU FC2.1", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "ER") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_CER)

#Determine the limits based on the range of data
x_limits <- range(interp_df$date_num, na.rm = TRUE)
y_limits <- range(10^interp_df$plastic_conc_pseudocount, na.rm = TRUE)

#Adjust the plot
avg_plot_interp_no_contours_CER <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0), na.value = NA) +  # Set NA values to be transparent
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),  
    limits = x_limits,  
    expand = c(0, 0)
  ) +
  scale_y_log10(
    name = "Log Plastic Concentration (g/L)",
    limits = y_limits,  
    expand = c(0, 0)
  ) +
  theme_minimal() +
  labs(title = "TPU FC2.1", x = "Day Number", y = "Plastic Concentration Pseudocount", fill = "ER") +
  theme(axis.text.x = element_text(angle = 90), 
        plot.title = element_text(size = 16, hjust = 0.5))

print(avg_plot_interp_no_contours_CER)

#Format changes
#Define the original concentrations for Plot C
original_concentrations_C <- c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385)

#Add pseudocount to the concentrations
pseudocount_concentrations_C <- original_concentrations_C + 0.009

#Set x_limits to cover the data range
x_limits <- range(interp_df$date_num)  

#geom_raster()
avg_plot_interp_no_contours_CER <- ggplot(interp_df, aes(x = date_num, y = 10^plastic_conc_pseudocount, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  
  
  #Set the fill color using viridis
  scale_fill_viridis_c(option = "plasma", direction = -1, limits = c(-100, 0), na.value = NA, guide = "none") +  # Remove legend
  
  #Customize the x-axis
  scale_x_continuous(
    breaks = seq(10, 90, by = 10),  
    limits = x_limits,  
    expand = c(0, 0)  
  ) +
  
  #Customize the y-axis with pseudocounts 
  scale_y_log10(
    name = "Log plastic concentration (g/L)",
    breaks = pseudocount_concentrations_C,  
    labels = sprintf("%.3f", original_concentrations_C),  
    limits = c(min(pseudocount_concentrations_C), max(pseudocount_concentrations_C)),  # Set y-axis limits
    expand = c(0, 0)  
  ) +
  
  #Apply coord_cartesian
  coord_cartesian(xlim = c(0, max(x_limits))) +  
  
  #Apply minimal theme
  theme_minimal() +
  labs(title = "TPU FC2.1", x = "Day Number", y = "Plastic Concentration Pseudocount") +
  
  #Customize
  theme(
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    plot.title = element_text(size = 16, hjust = 0.5),  
    
    # Add visible tick marks (excluding 0)
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),  
    
    #Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),  
    panel.border = element_blank()  
  )

#Print the updated heatmap plot for ER Plot C with three decimal places on the y-axis
print(avg_plot_interp_no_contours_CER)


#Now produce the combined figure
#Remove x-axis titles from AER and CER
# Adjust AER plot
avg_plot_interp_no_contours_AER <- avg_plot_interp_no_contours_AER + 
  theme(axis.title.x = element_blank())  

#Adjust BER plot 
avg_plot_interp_no_contours_BER <- avg_plot_interp_no_contours_BER + 
  theme(axis.title.y = element_blank(),  
        axis.text.y = element_blank(),   
        axis.ticks.y = element_blank())  

#Adjust CER plot
avg_plot_interp_no_contours_CER <- avg_plot_interp_no_contours_CER + 
  theme(axis.title.y = element_blank(),  
        axis.text.y = element_blank(),   
        axis.ticks.y = element_blank(),  
        axis.title.x = element_blank())  

#Combine ER plots
combined_ER_plot <- avg_plot_interp_no_contours_AER + 
  avg_plot_interp_no_contours_BER + 
  avg_plot_interp_no_contours_CER + 
  plot_layout(ncol = 3)

#Add overall title
final_combined_ER_plot <- combined_ER_plot + 
  plot_annotation(title = "Ecosystem Respiration", 
                  theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")))

#Display final combined ER plot
print(final_combined_ER_plot)


#Adjust AER plot
avg_plot_interp_no_contours_AER <- avg_plot_interp_no_contours_AER + 
  theme(
    axis.title.x = element_blank(),       
    axis.title.y = element_text(size = 14, face = "bold")  
  )

#Adjust BER plot
avg_plot_interp_no_contours_BER <- avg_plot_interp_no_contours_BER + 
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),  
    axis.title.y = element_blank(),       
    axis.text.y = element_blank(),        
    axis.ticks.y = element_blank()        
  )

#Adjust CER plot
avg_plot_interp_no_contours_CER <- avg_plot_interp_no_contours_CER + 
  theme(
    axis.title.y = element_blank(),       
    axis.text.y = element_blank(),        
    axis.ticks.y = element_blank(),       
    axis.title.x = element_blank()        
  )

#Combine ER plots
combined_ER_plot <- avg_plot_interp_no_contours_AER + 
  avg_plot_interp_no_contours_BER + 
  avg_plot_interp_no_contours_CER + 
  plot_layout(ncol = 3)

#Add an overall title
final_combined_ER_plot <- combined_ER_plot

#Display final combined ER plot
print(final_combined_ER_plot) ### ER data used for MANUSCRIPT [FIGURE S2]







####GAMs for ER data########
#Convert date to numeric scale (number of days since the start of the study)
er_data$date <- as.Date(er_data$date, format = "%m/%d/%Y")
start_date <- min(er_data$date)
er_data$date_num <- as.numeric(er_data$date - start_date)

#Make plastic_conc a factor and plastic_conc numeric
er_data$plastic_type <- as.factor(er_data$plastic_type)
er_data$plastic_conc <- as.numeric(er_data$plastic_conc)

#Change date_num values to actual concentrations
er_data <- er_data %>%
  mutate(plastic_conc = case_when(
    plastic_conc == 0 ~ 0.000,
    plastic_conc == 1 ~ 0.004,
    plastic_conc == 2 ~ 0.008,
    plastic_conc == 3 ~ 0.013,
    plastic_conc == 4 ~ 0.023,
    plastic_conc == 5 ~ 0.041,
    plastic_conc == 6 ~ 0.072,
    plastic_conc == 7 ~ 0.126,
    plastic_conc == 8 ~ 0.220,
    plastic_conc == 9 ~ 0.385,
    TRUE ~ plastic_conc  # Default case to leave values unchanged
  ))

#Check data distribution to meet assumptions
#Histogram of response variable
hist(er_data$ER, main="Histogram of ER", xlab="ER")
#Density plot for smoother visualization
plot(density(na.omit(er_data$ER)), main="Density Plot of ER")
#Plastic type
ggplot(er_data, aes(x=plastic_type, y=ER)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess", color="blue") +
  theme_minimal() +
  labs(title="Relationship between plastic type and ER", x="Predictor", y="ER")
#Day number
ggplot(er_data, aes(x=date_num, y=ER)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess", color="blue") +
  theme_minimal() +
  labs(title="Relationship between date_num and ER", x="Predictor", y="ER")
#Plastic concentration
ggplot(er_data, aes(x=plastic_conc, y=ER)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess", color="blue") +
  theme_minimal() +
  labs(title="Relationship between plastic concentration and ER", x="Predictor", y="ER")
#Check auto correlations
acf(er_data$date_num)
table(er_data$plastic_type)

#Fix the non-normal data issue by setting family distribution (Not this one)
gam_model5 <- gam(ER ~ s(plastic_conc, by = plastic_type) + s(date_num, by = plastic_type) + plastic_type, family = Gamma(link = "log"), data = er_data)
summary(gam_model5)
plot(gam_model5)
check_result5 <- gam.check(gam_model5)

#Fix the non-normal data issue by transforming the data (This one is best)
gam_model_transformed <- gam(I(log(ER + 1)) ~ s(plastic_conc, by = plastic_type) + s(date_num, by = plastic_type) + plastic_type, data = er_data)
summary(gam_model_transformed)
check_result_transformed <- gam.check(gam_model_transformed)

#From gam.check we should increase the k (Not this One)
gam_model_adjusted <- gam(I(log(ER + 1)) ~ s(plastic_conc, by = plastic_type, k = 10) + 
                            s(date_num, by = plastic_type, k = 10) + plastic_type, 
                          data = er_data)
summary(gam_model_adjusted)
check_result_adjusted <- gam.check(gam_model_adjusted)

#Fit the GAM model (Not this one)
gam_model6 <- gam(ER ~ s(plastic_conc, by = plastic_type) + s(date_num, by = plastic_type) + plastic_type, data = er_data)
summary(gam_model6)
plot(gam_model6)
#Check GAM model
check_result6 <- gam.check(gam_model6)
plot(check_result6)

#Compare models
summary(gam_model_transformed)
gam.check(gam_model_transformed)

summary(gam_model_adjusted)
gam.check(gam_model_adjusted)

summary(gam_model6)
gam.check(gam_model6)

summary(gam_model5)
gam.check(gam_model5)

#Compare AIC values
AIC(gam_model_transformed, gam_model_adjusted, gam_model6, gam_model5)

#Visualize smooths
par(mfrow = c(1, 2)) 
plot(gam_model_transformed, pages = 1)
plot(gam_model_adjusted, pages = 1)


#New column for absolute value of ER
er_data <- er_data %>% mutate(ER_abs = abs(ER))

#More complex models
ERgam_model1 <- gam(I(log(ER_abs + 1)) ~ 
                      s(date_num) + 
                      s(plastic_conc), 
                    data = er_data)

ERgam_model2 <- gam(I(log(ER_abs + 1)) ~ 
                      s(date_num) + 
                      s(plastic_conc) + 
                      plastic_type, 
                    data = er_data)

ERgam_model3 <- gam(I(log(ER_abs + 1)) ~ 
                      s(date_num, by = plastic_type) + 
                      s(plastic_conc, by = plastic_type), 
                    data = er_data)

ERgam_model4 <- gam(I(log(ER_abs + 1)) ~ 
                      s(date_num, by = plastic_type) + 
                      s(plastic_conc, by = plastic_type) + 
                      plastic_type, 
                    data = er_data)

#Compare models with AIC
ERmodel_comparison <- AIC(ERgam_model1, ERgam_model2, ERgam_model3, ERgam_model4)
print(ERmodel_comparison) ### ER data used for MANUSCRIPT [Table S1]

#Check diagnostics for best model
ERbest_model <- ERgam_model4
summary(ERbest_model) ### ER data used for MANUSCRIPT [Table S2]
gam.check(ERbest_model)
plot(ERbest_model, pages = 1)

#plot_difference
#Group A and Group B
er_diff_AB_date <- plot_difference(
  ERbest_model,
  series = date_num,
  difference = list(plastic_type = c("A", "B"))
) +
  labs(x = "Day Number", y = "log ER (Difference smooth)") +
  ggtitle("Difference between Elastollan and TPU 181") +
  theme_minimal()

#Group A and Group C
er_diff_AC_date <- plot_difference(
  ERbest_model,
  series = date_num,
  difference = list(plastic_type = c("A", "C"))
) +
  labs(x = "Day Number", y = "log ER (Difference smooth)") +
  ggtitle("Difference between Elastollan and TPU FC2.1") +
  theme_minimal()

#Group B and Group C
er_diff_BC_date <- plot_difference(
  ERbest_model,
  series = date_num,
  difference = list(plastic_type = c("B", "C"))
) +
  labs(x = "Day Number", y = "log ER (Difference smooth)") +
  ggtitle("Difference between TPU 181 and TPU FC2.1") +
  theme_minimal()

#Group A and Group B
er_diff_AB_conc <- plot_difference(
  ERbest_model,
  series = plastic_conc,
  difference = list(plastic_type = c("A", "B"))
) +
  labs(x = "Plastic Concentration", y = "log ER (Difference smooth)") +
  ggtitle("Difference between Elastollan and TPU 181") +
  theme_minimal()

#Group A and Group C
er_diff_AC_conc <- plot_difference(
  ERbest_model,
  series = plastic_conc,
  difference = list(plastic_type = c("A", "C"))
) +
  labs(x = "Plastic Concentration", y = "log ER (Difference smooth)") +
  ggtitle("Difference between Elastollan and TPU FC2.1") +
  theme_minimal()

# Comparison between plastic types: Group B vs Group C
er_diff_BC_conc <- plot_difference(
  ERbest_model,
  series = plastic_conc,
  difference = list(plastic_type = c("B", "C"))
) +
  labs(x = "Plastic Concentration", y = "log ER (Difference smooth)") +
  ggtitle("Difference between TPU 181 and TPU FC2.1") +
  theme_minimal()

#Print
print(er_diff_AB_date)
print(er_diff_AC_date)
print(er_diff_BC_date)
print(er_diff_AB_conc)
print(er_diff_AC_conc)
print(er_diff_BC_conc)
### ER data used for MANUSCRIPT [Figure S11]
