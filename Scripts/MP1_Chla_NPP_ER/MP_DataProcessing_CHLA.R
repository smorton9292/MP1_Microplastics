#This R script runs the chla analysis and figure creation


rm(list = ls())

setwd("C:/Users/DELL/Documents/R/projects/MP1_CHLA_NPP_ER")
data<-read.csv("MP_Chla.csv")

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
install.packages("cowplot")
install.packages("tidymv")

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
library(cowplot)
library(tidymv)

#Now gather the data to long format
data_long <- data %>% 
  gather(key = "replicate", value = "CHLA_UGL", -c(1:3))

#Filter out 4 true and verified machine error sample reads
data_long <- data_long %>% 
  filter(!(row_number() %in% c(87, 185, 785, 1385)))

summary_data <- data_long %>%
  group_by(trt_code) %>%
  summarize(
    Mean = mean(CHLA_UGL),
    Median = median(CHLA_UGL),
    SD = sd(CHLA_UGL),
    Min = min(CHLA_UGL),
    Max = max(CHLA_UGL)
  )

#Check normality - hist
ggplot(data_long, aes(CHLA_UGL)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  labs(title="Histogram for CHLA_UGL")

#QQ-plot
qqnorm(data_long$CHLA_UGL, main = "Q-Q plot for CHLA_MGM3")
qqline(data_long$CHLA_UGL)

#check date format
str(data_long)

#Format dates to as.Date
data_long$date <- as.Date(data_long$date, format = "%m/%d/%Y")

#Plot of each channel over time
ggplot(data_long, aes(x = as.factor(date), y = CHLA_UGL)) +
  geom_boxplot() +
  facet_wrap(~ trt_code) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Date", y = "CHLA_UGL", title = "Boxplot for each Mesocosm")

#Each y axis scaled to 500 (Not Necessary to run)
ggplot(data_long, aes(x = as.factor(date), y = CHLA_UGL)) +
  geom_boxplot() +
  facet_wrap(~ trt_code) +
  ylim(0, 500) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Date", y = "CHLA_UGL", title = "Boxplot for each Mesocosm")

#Plot each plastic type on one figure with shades depicting concentrations
#Split the trt_code into letter and number parts
data_long <- data_long %>%
  mutate(
    letter_group = substr(trt_code, 1, 1),
    number_designation = as.numeric(substr(trt_code, 2, 2))
  )

#Calculate the mean and standard error for each group
summary_data <- data_long %>%
  group_by(letter_group, number_designation, date) %>%
  summarize(
    Mean = mean(CHLA_UGL),
    SE = sd(CHLA_UGL) / sqrt(n()),
    .groups = 'drop'
  )

#Scatter plot
ggplot(summary_data, aes(x = as.factor(date), y = Mean, color = as.factor(number_designation))) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  facet_wrap(~ letter_group) +
  ylim(0, 500) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Date", y = "CHLA_UGL", title = "Scatter Plot with Error Bars for each Mesocosm")

#Improve figure
#Function to create a color gradient
get_colors <- function(base_color, n) {
  seq_color <- colorRampPalette(c("grey", base_color))(n)
  return(seq_color)
}

#Define the unique letter groups
letter_groups <- unique(data_long$letter_group)

#Define the desired order for the letter groups
letter_groups <- c("A", "B", "C")

#Define the base color for each letter group
colors <- c("darkred", "darkblue", "darkgreen")

#Create a list to store the plots
plots <- list()

max(data_long$CHLA_UGL)
max(summary_data$Mean)

#Loop through the letter groups and create a plot for each
for (i in seq_along(letter_groups)) {
  letter_data <- filter(data_long, letter_group == letter_groups[i])
  
  #Calculate the mean and standard deviation for each group
  summary_data <- letter_data %>%
    group_by(date, number_designation) %>%
    summarise(
      Mean = mean(CHLA_UGL),
      SD = sd(CHLA_UGL),
      .groups = 'drop' 
    )
  
  #Format the date without the year and store it in a new variable
  summary_data$date_without_year <- format(as.Date(summary_data$date, format = "%m/%d/%Y"), "%m/%d")
  
  #Add a numeric date variable for trend line
  summary_data$date_num <- as.numeric(as.factor(summary_data$date))
  
  #Define the color palette using the custom function
  palette <- get_colors(colors[i], length(unique(letter_data$number_designation)))
  
  plot <- ggplot(summary_data, aes(x = as.factor(date_without_year), y = Mean, color = as.factor(number_designation))) +
    geom_point() +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
    geom_smooth(aes(x = date_num), se = FALSE) +
    ylim(0, 500) +
    scale_color_manual(values = palette, name = "Concentration") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.title = element_text(size = 9), 
      legend.text = element_text(size = 7),   
      legend.box.margin = margin(45, 0, 0, 0)
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = "Time", y = "CHLA_UGL", title = paste("Type", letter_groups[i])) +
    guides(color = guide_legend(ncol = 2))
  
  #Add the plot to the list
  plots[[i]] <- plot
}

#Display the plots in a grid
do.call(grid.arrange, plots)
#Error does not effect actual data - just in how the smooth was added to figure. Only exploratory OK to move forward.


#Set up unique titles for final plots
custom_titles <- c("Elastollan", "TPU 181", "TPU 2.1")  

#Starting heatmap plotting
#Loop through the letter groups and create a plot for each
for (i in seq_along(letter_groups)) {
  letter_data <- filter(data_long, letter_group == letter_groups[i])
  
  #Calculate the mean and standard deviation for each group
  summary_data <- letter_data %>%
    group_by(date, number_designation) %>%
    summarise(
      Mean = mean(CHLA_UGL),
      SD = sd(CHLA_UGL),
      .groups = 'drop' 
    )
  
  #Filter the data
  summary_data_filtered <- filter(summary_data, Mean <= 500)
  
  #Format the date without the yeare
  summary_data_filtered$date_without_year <- format(as.Date(summary_data_filtered$date, format = "%m/%d/%Y"), "%m/%d")
  
  #Add a numeric date variable for trendline
  summary_data_filtered$date_num <- as.numeric(as.factor(summary_data_filtered$date))
  
  #Convert date to a numeric scale (number of days since the start of the study)
  summary_data_filtered$date <- as.Date(summary_data_filtered$date, format = "%m/%d/%Y")
  start_date <- min(summary_data_filtered$date)
  summary_data_filtered$date <- as.numeric(summary_data_filtered$date - start_date)
  
  #Change number_designation values to actual concentrations

  summary_data_filtered <- summary_data_filtered %>%
    mutate(number_designation = case_when(
      number_designation == 0 ~ "0.000",
      number_designation == 1 ~ "0.004",
      number_designation == 2 ~ "0.008",
      number_designation == 3 ~ "0.013",
      number_designation == 4 ~ "0.023",
      number_designation == 5 ~ "0.041",
      number_designation == 6 ~ "0.072",
      number_designation == 7 ~ "0.126",
      number_designation == 8 ~ "0.220",
      number_designation == 9 ~ "0.385",
      TRUE ~ as.character(number_designation)  
    ))
  
  #Make it categorical again
  summary_data_filtered$number_designation <- as.factor(summary_data_filtered$number_designation)
  
  
  str(summary_data_filtered$number_designation)
  
  #Heatmap Plot
  heatmap_plot <- ggplot(summary_data_filtered, aes(x = as.factor(date), y = as.factor(number_designation), fill = Mean)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = colors[i]) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(angle = 0),
      plot.title = element_text(size = 16, hjust = 0.5)
    ) +
    labs(x = "Day", y = "Plastic Concentration", fill = "CHLA_UGL", title = custom_titles[i])
  
  #Add the heatmap plot to the list
  plots[[i]] <- heatmap_plot
}

#Display the plots in a grid
do.call(grid.arrange, plots)

str(data_long)



##Try to get contour lines##
Aselected_letter_group <- "A"
Bselected_letter_group <- "B"
Cselected_letter_group <- "C"

#Filter the data for the selected letter group
Aselected_data <- filter(data_long, letter_group == Aselected_letter_group)
Bselected_data <- filter(data_long, letter_group == Bselected_letter_group)
Cselected_data <- filter(data_long, letter_group == Cselected_letter_group)

#Calculate the mean and standard deviation for the selected group
###For A###
Asummary_data_contour <- Aselected_data %>%
  group_by(date, number_designation) %>%
  summarise(
    Mean = mean(CHLA_UGL),
    SD = sd(CHLA_UGL),
    .groups = 'drop'
  )

###For B###
Bsummary_data_contour <- Bselected_data %>%
  group_by(date, number_designation) %>%
  summarise(
    Mean = mean(CHLA_UGL),
    SD = sd(CHLA_UGL),
    .groups = 'drop'
  )

###For C###
Csummary_data_contour <- Cselected_data %>%
  group_by(date, number_designation) %>%
  summarise(
    Mean = mean(CHLA_UGL),
    SD = sd(CHLA_UGL),
    .groups = 'drop'
  )

#Format the date without the year and store it in a new variable
Asummary_data_contour$date_without_year <- format(as.Date(Asummary_data_contour$date, format = "%m/%d/%Y"), "%m/%d")

Bsummary_data_contour$date_without_year <- format(as.Date(Bsummary_data_contour$date, format = "%m/%d/%Y"), "%m/%d")

Csummary_data_contour$date_without_year <- format(as.Date(Csummary_data_contour$date, format = "%m/%d/%Y"), "%m/%d")



#Change number_designation values to actual concentrations
###For A###
Asummary_data_contour <- Asummary_data_contour %>%
  mutate(number_designation = case_when(
    number_designation == 0 ~ 0.000,
    number_designation == 1 ~ 0.004,
    number_designation == 2 ~ 0.008,
    number_designation == 3 ~ 0.013,
    number_designation == 4 ~ 0.023,
    number_designation == 5 ~ 0.041,
    number_designation == 6 ~ 0.072,
    number_designation == 7 ~ 0.126,
    number_designation == 8 ~ 0.220,
    number_designation == 9 ~ 0.385,
    TRUE ~ number_designation  
  ))

#Change number_designation values to actual concentrations
###For B###
Bsummary_data_contour <- Bsummary_data_contour %>%
  mutate(number_designation = case_when(
    number_designation == 0 ~ 0.000,
    number_designation == 1 ~ 0.004,
    number_designation == 2 ~ 0.008,
    number_designation == 3 ~ 0.013,
    number_designation == 4 ~ 0.023,
    number_designation == 5 ~ 0.041,
    number_designation == 6 ~ 0.072,
    number_designation == 7 ~ 0.126,
    number_designation == 8 ~ 0.220,
    number_designation == 9 ~ 0.385,
    TRUE ~ number_designation  
  ))

#Change number_designation values to actual concentrations
###For C###
Csummary_data_contour <- Csummary_data_contour %>%
  mutate(number_designation = case_when(
    number_designation == 0 ~ 0.000,
    number_designation == 1 ~ 0.004,
    number_designation == 2 ~ 0.008,
    number_designation == 3 ~ 0.013,
    number_designation == 4 ~ 0.023,
    number_designation == 5 ~ 0.041,
    number_designation == 6 ~ 0.072,
    number_designation == 7 ~ 0.126,
    number_designation == 8 ~ 0.220,
    number_designation == 9 ~ 0.385,
    TRUE ~ number_designation  
  ))

# Convert date to a Date object and then to a numeric scale
###For A###
Asummary_data_contour$date <- as.Date(Asummary_data_contour$date, format = "%m/%d/%Y")
start_date <- min(Asummary_data_contour$date)
Asummary_data_contour$date_num <- as.numeric(Asummary_data_contour$date - start_date)

Asummary_data_contour$number_designation <- as.numeric(Asummary_data_contour$number_designation)

###For B###
Bsummary_data_contour$date <- as.Date(Bsummary_data_contour$date, format = "%m/%d/%Y")
start_date <- min(Bsummary_data_contour$date)
Bsummary_data_contour$date_num <- as.numeric(Bsummary_data_contour$date - start_date)

Bsummary_data_contour$number_designation <- as.numeric(Bsummary_data_contour$number_designation)


###For C###
Csummary_data_contour$date <- as.Date(Csummary_data_contour$date, format = "%m/%d/%Y")
start_date <- min(Csummary_data_contour$date)
Csummary_data_contour$date_num <- as.numeric(Csummary_data_contour$date - start_date)

Csummary_data_contour$number_designation <- as.numeric(Csummary_data_contour$number_designation)

#Get date without year
###For A###
Asummary_data_contour <- dplyr::select(Asummary_data_contour, -date_without_year)
Asummary_data_contour <- dplyr::select(Asummary_data_contour, -date)

###For B###
Bsummary_data_contour <- dplyr::select(Bsummary_data_contour, -date_without_year)
Bsummary_data_contour <- dplyr::select(Bsummary_data_contour, -date)

###For C###
Csummary_data_contour <- dplyr::select(Csummary_data_contour, -date_without_year)
Csummary_data_contour <- dplyr::select(Csummary_data_contour, -date)

###FOr A###
#Generate all combinations of date_num and number_designation
all_combinations <- expand.grid(
  date_num = unique(Asummary_data_contour$date_num),
  number_designation = unique(Asummary_data_contour$number_designation)
)

#Find missing combinations
Amissing_combinations <- anti_join(all_combinations, Asummary_data_contour, by = c("date_num", "number_designation"))


#Ensure the data is combined and sorted
Afull_data <- bind_rows(Asummary_data_contour, Amissing_combinations) %>%
  arrange(number_designation, date_num)

#Interpolate Mean and SD for missing combinations - next closest data point
Afull_data <- Afull_data %>%
  group_by(number_designation) %>%
  mutate(
    Mean = ifelse(is.na(Mean), na.approx(Mean, na.rm = FALSE), Mean),
    SD = ifelse(is.na(SD), na.approx(SD, na.rm = FALSE), SD)
  )

#Splitting back the data
Asummary_data_contour <- filter(Afull_data, !is.na(Mean) & !is.na(SD))
Amissing_combinations_interpolated <- filter(Afull_data, is.na(Mean) | is.na(SD))

#Check and append the interpolated values if missing combinations were filled
if (nrow(Amissing_combinations_interpolated) > 0) {
  #This should now have interpolated values where possible
  Asummary_data_contour <- bind_rows(Asummary_data_contour, Amissing_combinations_interpolated)
} else {
  cat("No missing combinations found or interpolatable in the dataframe.\n")
}

#Add a pseudocount of 0.009 to treatment concentrations
Asummary_data_contour$number_designation_pseudocount <- Asummary_data_contour$number_designation + 0.009

#Correcting the geom_contour usage so that y is on a log axis
###For A###
heatmap_contour_plot_A_log <- ggplot(Asummary_data_contour, aes(x = date_num, y = number_designation_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean), width = 2.9, height = 0.10) +  #Adjust width and height parameters to minimize gaps
  geom_contour(color = "black") +  #This adds the contour lines based on `z = Mean`
  scale_fill_viridis_c(option = "plasma", direction = -1) + #Use 'plasma' palette and reverse gradient
  scale_x_continuous(name = "Day Number", breaks = seq(min(Asummary_data_contour$date_num), max(Asummary_data_contour$date_num), by = 5)) + #Add x-axis tick marks with a step of 5
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +  #Set y-axis to log scale
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  labs(fill = "Chlα µg/L", title = paste("Elastollan"))

#Print the plot
print(heatmap_contour_plot_A_log)

#Interpolate data
interp_dataA <- with(Asummary_data_contour, {
  interp(x = date_num, y = log10(number_designation_pseudocount), z = Mean, 
         xo = seq(min(date_num), max(date_num), length = 100),  #Adjust granularity as needed
         yo = seq(min(log10(number_designation_pseudocount)), max(log10(number_designation_pseudocount)), length = 100), 
         linear = FALSE)  #Set to TRUE for linear interpolation, FALSE for natural spline
})

#Convert interpolated data back to a data frame for ggplot
interp_dfA <- expand.grid(date_num = interp_dataA$x, log_concentration = interp_dataA$y)
interp_dfA$Mean <- as.vector(interp_dataA$z)


#Create the continuous heatmap with geom_raster()
heatmap_contour_plot_A_log_lines <- ggplot(interp_dfA, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(
    name = "Day Number", 
    breaks = seq(min(interp_dfA$date_num), max(interp_dfA$date_num), by = 10),  #Adjust tick marks
    expand = c(0, 0),  #No expansion
    limits = range(interp_dfA$date_num)  #Explicitly set x-axis limits
  ) +
  scale_y_log10(
    name = "Log Plastic Concentration (g/L)", 
    expand = c(0, 0),  #No expansion
    limits = c(min(10^interp_dfA$log_concentration, na.rm = TRUE), max(10^interp_dfA$log_concentration, na.rm = TRUE))  # Explicitly set y-axis limits
  ) +
  stat_contour(aes(z = Mean), color = "black", bins = 5) +  #Retain smoother contour lines with fewer bins
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines")  
  ) +
  labs(fill = "Chlα µg/L", title = "Elastollan")

#Print the adjusted heatmap plot with smoother contour lines
print(heatmap_contour_plot_A_log_lines)

#No Contours
heatmap_contour_plot_A_log <- ggplot(interp_dfA, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(name = "Day Number", 
                     breaks = seq(10, 90, by = 10),  
                     expand = c(0, 0),  #No expansion
                     limits = range(interp_dfA$date_num)) +  #Explicitly set x-axis limits
  scale_y_log10(name = "Log Plastic Concentration (g/L)", 
                expand = c(0, 0),  #No expansion
                limits = c(min(10^interp_dfA$log_concentration, na.rm = TRUE), max(10^interp_dfA$log_concentration, na.rm = TRUE))) +  # Explicitly set y-axis limits
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines")  
  ) +
  labs(fill = "Chlα µg/L", title = "Elastollan")

#Print the adjusted heatmap plot
print(heatmap_contour_plot_A_log)

#Add a vector for plastic concentrations corresponding to the number designations
plastic_concentrations <- c(0.009, 0.013, 0.017, 0.022, 0.032, 0.050, 0.081, 0.135, 0.229, 0.394)


#No Contours
heatmap_contour_plot_A_log <- ggplot(interp_dfA, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(name = "Day Number", 
                     breaks = seq(10, 90, by = 10),  
                     expand = c(0, 0),  #No expansion
                     limits = range(interp_dfA$date_num)) +  s
  scale_y_log10(name = "Log Plastic Concentration (g/L)", 
                expand = c(0, 0),  #No expansion
                limits = c(min(10^interp_dfA$log_concentration, na.rm = TRUE), max(10^interp_dfA$log_concentration, na.rm = TRUE)),
                sec.axis = sec_axis(~., 
                                    breaks = c(0.0093, 0.013, 0.017, 0.022, 0.032, 0.050, 0.081, 0.135, 0.229, 0.385),  # Your plastic concentrations
                                    labels = NULL,  #No labels
                                    name = NULL)) +  #Add tick marks to the secondary y-axis
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    axis.text.y.right = element_blank(),  #Hide the labels on the right y-axis
    axis.ticks.y.right = element_line(color = "black", size = 1.5),  #Show tick marks on the right y-axis
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines")  
  ) +
  labs(fill = "Chlα µg/L", title = "Elastollan")

#Print the adjusted heatmap plot
print(heatmap_contour_plot_A_log)




#Try to fix the y axis with hard coding the treatment concentrations
heatmap_contour_plot_A_log <- ggplot(interp_dfA, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(name = "Day Number", 
                     breaks = seq(10, 90, by = 10),  
                     expand = c(0, 0),  #No expansion
                     limits = range(interp_dfA$date_num)) +  
  scale_y_log10(name = "Plastic Concentration (g/L)", 
                breaks = c(0.0093, 0.013, 0.017, 0.022, 0.032, 0.050, 0.081, 0.135, 0.229, 0.385),  # Hard-coded plastic concentrations
                labels = c("0.0093", "0.013", "0.017", "0.022", "0.032", "0.050", "0.081", "0.135", "0.229", "0.385"),  # Custom labels
                expand = c(0, 0),  #No expansion
                limits = c(0.0093, 0.385)) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines")  
  ) +
  labs(fill = "Chlα µg/L", title = "Elastollan")

#Print the adjusted heatmap plot
print(heatmap_contour_plot_A_log)






#Original plastic concentrations
original_concentrations <- c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385)

#Add pseudocount to the concentrations for plotting
pseudocount_concentrations <- original_concentrations + 0.009

heatmap_contour_plot_A_log_y <- ggplot(interp_dfA, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(name = "Day Number", 
                     breaks = seq(10, 90, by = 10),  
                     expand = c(0, 0),  #No expansion
                     limits = range(interp_dfA$date_num)) +  
  scale_y_log10(name = "Plastic Concentration (g/L)", 
                breaks = pseudocount_concentrations,  #Set breaks with pseudocount
                labels = original_concentrations,  #Set labels to original concentrations
                expand = c(0, 0),  # No expansion
                limits = c(min(pseudocount_concentrations), max(pseudocount_concentrations))) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  #Center-align tick labels
    axis.text.y = element_text(angle = 0),
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.ticks = element_line(color = "black"),   
    axis.ticks.length = unit(0.15, "cm"),  #Set the length of the tick marks
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines")  
  ) +
  guides(fill = "none") +  #Remove the legend
  labs(fill = "Chlα µg/L", title = "Elastollan")

#Print the adjusted heatmap plot
print(heatmap_contour_plot_A_log_y)

#Final Adjustments
#Increase the tickmarks
heatmap_contour_plot_A_log_y <- ggplot(interp_dfA, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  
  # Customize x-axis
  scale_x_continuous(name = "Day Number", 
                     breaks = seq(10, 90, by = 10),  #Adjust tick marks
                     expand = c(0, 0),  #No expansion
                     limits = range(interp_dfA$date_num)) +  
  
  #Customize y-axis with three decimal places
  scale_y_log10(name = "Log plastic concentration (g/L)", 
                breaks = pseudocount_concentrations,  #Set breaks with pseudocount
                labels = sprintf("%.3f", original_concentrations),  #Set labels to three decimal places
                expand = c(0, 0),  # No expansion
                limits = c(min(pseudocount_concentrations), max(pseudocount_concentrations))) +  
  
  #Customize
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 0, vjust = 0.5, hjust = 0.5),  
    axis.text.y = element_text(size = 14, angle = 0),  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.15, "cm"),  
    plot.title = element_text(size = 16, hjust = 0.5),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines")  
  ) +
  
  guides(fill = "none") +  
  labs(fill = "Chlα µg/L", title = "Elastollan")

#Print the adjusted heatmap plot
print(heatmap_contour_plot_A_log_y)
####Error is associated with interpolation beyond the range of observed values. Ok to move forward. Collected data represented.  





###FOr B###
#Generate all combinations of date_num and number_designation
all_combinations <- expand.grid(
  date_num = unique(Bsummary_data_contour$date_num),
  number_designation = unique(Bsummary_data_contour$number_designation)
)

#Find missing combinations
Bmissing_combinations <- anti_join(all_combinations, Bsummary_data_contour, by = c("date_num", "number_designation"))
#THERE WERE NO MISSING COMBINATIONS FOR B



###For B###
#Add a pseudocount of 0.001 to treatment concentrations
Bsummary_data_contour$number_designation_pseudocount <- Bsummary_data_contour$number_designation + 0.009


###For B###
heatmap_contour_plot_B_log <- ggplot(Bsummary_data_contour, aes(x = date_num, y = number_designation_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean), width = 2.9, height = 0.10) +  
  geom_contour(color = "black") +  #This adds the contour lines based on `z = Mean`
  scale_fill_viridis_c(option = "plasma", direction = -1) + 
  scale_x_continuous(name = "Day Number", breaks = seq(min(Bsummary_data_contour$date_num), max(Bsummary_data_contour$date_num), by = 5)) +
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +  # Set y-axis to log scale
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  labs(fill = "Chlα µg/L", title = paste("TPU 181"))

#Print the plot for B
print(heatmap_contour_plot_B_log)


#Interpolation for B
interp_dataB <- with(Bsummary_data_contour, {
  interp(
    x = date_num, 
    y = log10(number_designation_pseudocount), 
    z = Mean, 
    xo = seq(min(date_num), max(date_num), length = 100),  
    yo = seq(min(log10(number_designation_pseudocount)), max(log10(number_designation_pseudocount)), length = 100),
    linear = FALSE  
  )
})

#Convert interpolated data back to a data frame for ggplot
interp_dfB <- expand.grid(date_num = interp_dataB$x, log_concentration = interp_dataB$y)
interp_dfB$Mean <- as.vector(interp_dataB$z)

#For B Heatmap#
heatmap_contour_plot_B_log_lines <- ggplot(interp_dfB, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(
    name = "Day Number", 
    breaks = seq(min(interp_dfB$date_num), max(interp_dfB$date_num), by = 10),  
    expand = c(0, 0),  #No expansion
    limits = range(interp_dfB$date_num)  
  ) +
  scale_y_log10(
    name = "Log Plastic Concentration (g/L)", 
    expand = c(0, 0),  #No expansion
    limits = c(min(10^interp_dfB$log_concentration, na.rm = TRUE), max(10^interp_dfB$log_concentration, na.rm = TRUE))  # Explicitly set y-axis limits
  ) +
  stat_contour(aes(z = Mean), color = "black", bins = 5) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),  
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines")  
  ) +
  labs(fill = "Chlα µg/L", title = "TPU 181")  

#Print the adjusted heatmap plot for B with smoother contour lines
print(heatmap_contour_plot_B_log_lines)

#Without Contours
heatmap_contour_plot_B_log <- ggplot(interp_dfB, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(
    name = "Day Number", 
    breaks = seq(10, 90, by = 10),  
    expand = c(0, 0),  # No expansion
    limits = range(interp_dfB$date_num)  
  ) +
  scale_y_log10(
    name = "Log Plastic Concentration (g/L)", 
    expand = c(0, 0),  # No expansion
    limits = c(min(10^interp_dfB$log_concentration, na.rm = TRUE), max(10^interp_dfB$log_concentration, na.rm = TRUE))  
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),  
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines")  
  ) +
  labs(fill = "Chlα µg/L", title = "TPU 181")  

#Print the adjusted heatmap plot for B with smoother contour lines
print(heatmap_contour_plot_B_log)


#Make the formating changes for B
#Original plastic concentrations for Plot B
original_concentrations_B <- c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385)

#Add pseudocount to the concentrations for plotting
pseudocount_concentrations_B <- original_concentrations_B + 0.009

heatmap_contour_plot_B_log_y <- ggplot(interp_dfB, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(
    name = "Day Number", 
    breaks = seq(10, 90, by = 10),  
    expand = c(0, 0),  #No expansion
    limits = range(interp_dfB$date_num)  
  ) +
  scale_y_log10(
    name = "Plastic Concentration (g/L)", 
    breaks = pseudocount_concentrations_B,  # et breaks with pseudocount
    labels = original_concentrations_B,  
    expand = c(0, 0),  #No expansion
    limits = c(min(pseudocount_concentrations_B), max(pseudocount_concentrations_B))  
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(angle = 0),
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.15, "cm"),  
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines")   
  ) +
  guides(fill = "none") +  
  labs(fill = "Chlα µg/L", title = "TPU 181")  

#Print the adjusted heatmap plot for B
print(heatmap_contour_plot_B_log_y)


#Final format changes
heatmap_contour_plot_B_log_y <- ggplot(interp_dfB, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  
  # Customize x-axis
  scale_x_continuous(
    name = "Day Number", 
    breaks = seq(10, 90, by = 10),  
    expand = c(0, 0),  #No expansion
    limits = range(interp_dfB$date_num)  
  ) +
  
  # Customize y-axis
  scale_y_log10(
    name = "Plastic Concentration (g/L)", 
    breaks = pseudocount_concentrations_B,  #Set breaks with pseudocount
    labels = original_concentrations_B,  
    expand = c(0, 0),  #No expansion
    limits = c(min(pseudocount_concentrations_B), max(pseudocount_concentrations_B))  
  ) +
  
  # Customize the theme
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 0, vjust = 0.5, hjust = 0.5),  
    axis.text.y = element_text(size = 14, angle = 0),  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.15, "cm"),  
    plot.title = element_text(size = 16, hjust = 0.5),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    plot.margin = unit(c(1, 1, 1, 1), "lines")  
  ) +
  
  guides(fill = "none") +  
  labs(fill = "Chlα µg/L", title = "TPU 181") 

#Print the adjusted heatmap plot for B
print(heatmap_contour_plot_B_log_y)

#Remove y axis labels and title
heatmap_contour_plot_B_log_y <- ggplot(interp_dfB, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  
  #Customize x-axis
  scale_x_continuous(
    name = "Day Number", 
    breaks = seq(10, 90, by = 10),  
    expand = c(0, 0),  #No expansion
    limits = range(interp_dfB$date_num)  
  ) +
  
  #Customize y-axis
  scale_y_log10(
    breaks = pseudocount_concentrations_B,  #Set breaks with pseudocount
    labels = original_concentrations_B,  
    expand = c(0, 0),  #No expansion
    limits = c(min(pseudocount_concentrations_B), max(pseudocount_concentrations_B))  
  ) +
  
  # Customize the theme
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 0, vjust = 0.5, hjust = 0.5),  
    axis.text.y = element_blank(),  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_blank(),  
    axis.ticks.y = element_blank(),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.15, "cm"),  
    plot.title = element_text(size = 16, hjust = 0.5),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines")  
  ) +
  
  guides(fill = "none") +  
  labs(fill = "Chlα µg/L", title = "TPU 181")  

#Print the adjusted heatmap plot for B
print(heatmap_contour_plot_B_log_y)
####Error is associated with interpolation beyond the range of observed values. Ok to move forward. Collected data represented.  



###For C###
#Generate all combinations of date_num and number_designation
all_combinations <- expand.grid(
  date_num = unique(Csummary_data_contour$date_num),
  number_designation = unique(Csummary_data_contour$number_designation)
)

#Find missing combinations
Cmissing_combinations <- anti_join(all_combinations, Csummary_data_contour, by = c("date_num", "number_designation"))

#Add a pseudocount of 0.001 to treatment concentrations
Csummary_data_contour$number_designation_pseudocount <- Csummary_data_contour$number_designation + 0.009

#For C#
heatmap_contour_plot_C_log <- ggplot(Csummary_data_contour, aes(x = date_num, y = number_designation_pseudocount, z = Mean)) +
  geom_tile(aes(fill = Mean), width = 2.9, height = 0.10) +  
  geom_contour(color = "black") +  #This adds the contour lines based on `z = Mean`
  scale_fill_viridis_c(option = "plasma", direction = -1) + 
  scale_x_continuous(name = "Day Number", breaks = seq(min(Csummary_data_contour$date_num), max(Csummary_data_contour$date_num), by = 5)) + 
  scale_y_log10(name = "Log Plastic Concentration (g/L)") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  labs(fill = "Chlα µg/L", title = paste("TPU FC2.1"))

#Print the plot for C
print(heatmap_contour_plot_C_log)

#Interpolation for C
interp_dataC <- with(Csummary_data_contour, {
  interp(
    x = date_num, 
    y = log10(number_designation_pseudocount), 
    z = Mean, 
    xo = seq(min(date_num), max(date_num), length = 100),
    yo = seq(min(log10(number_designation_pseudocount)), max(log10(number_designation_pseudocount)), length = 100),
    linear = FALSE
  )
})

#Convert interpolated data back to a data frame for ggplot
interp_dfC <- expand.grid(date_num = interp_dataC$x, log_concentration = interp_dataC$y)
interp_dfC$Mean <- as.vector(interp_dataC$z)

#For C#
heatmap_contour_plot_C_log_lines <- ggplot(interp_dfC, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(
    name = "Day Number", 
    breaks = seq(min(interp_dfC$date_num), max(interp_dfC$date_num), by = 10),
    expand = c(0, 0),
    limits = range(interp_dfC$date_num)
  ) +
  scale_y_log10(
    name = "Log Plastic Concentration (g/L)", 
    expand = c(0, 0),
    limits = c(min(10^interp_dfC$log_concentration, na.rm = TRUE), max(10^interp_dfC$log_concentration, na.rm = TRUE))
  ) +
  stat_contour(aes(z = Mean), color = "black", bins = 5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  ) +
  labs(fill = "Chlα µg/L", title = "TPU FC2.1") 

#Print the adjusted heatmap plot for C with smoother contour lines
print(heatmap_contour_plot_C_log_lines)

#For C No Contour
#For C#
heatmap_contour_plot_C_log <- ggplot(interp_dfC, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(
    name = "Day Number", 
    breaks = seq(10, 90, by = 10),
    expand = c(0, 0),
    limits = range(interp_dfC$date_num)
  ) +
  scale_y_log10(
    name = "Log Plastic Concentration (g/L)", 
    expand = c(0, 0),
    limits = c(min(10^interp_dfC$log_concentration, na.rm = TRUE), max(10^interp_dfC$log_concentration, na.rm = TRUE))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),
    axis.text.y = element_text(angle = 0),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  ) +
  labs(fill = "Chlα µg/L", title = "TPU FC2.1")  

#Print the adjusted heatmap plot for C with smoother contour lines
print(heatmap_contour_plot_C_log)



#Original plastic concentrations for Plot C
original_concentrations_C <- c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385)

#Add pseudocount to the concentrations for plotting
pseudocount_concentrations_C <- original_concentrations_C + 0.009

heatmap_contour_plot_C_log_y <- ggplot(interp_dfC, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(
    name = "Day Number", 
    breaks = seq(10, 90, by = 10),  
    expand = c(0, 0),
    limits = range(interp_dfC$date_num)  
  ) +
  scale_y_log10(
    name = "Plastic Concentration (g/L)", 
    breaks = pseudocount_concentrations_C,  #Set breaks with pseudocount
    labels = original_concentrations_C,  
    expand = c(0, 0),
    limits = c(min(pseudocount_concentrations_C), max(pseudocount_concentrations_C))  # Set y-axis limits
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  
    axis.text.y = element_text(angle = 0),
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.15, "cm"),  
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines"),  
    legend.title = element_text(face = "bold", size = 10)  
  ) +
  labs(fill = "Chlα µg/L", title = "TPU FC2.1")  

#Print the adjusted heatmap plot for C
print(heatmap_contour_plot_C_log_y)


#Remove legend
heatmap_contour_plot_C_log_y <- ggplot(interp_dfC, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  scale_x_continuous(
    name = "Day Number", 
    breaks = seq(10, 90, by = 10),  
    expand = c(0, 0),
    limits = range(interp_dfC$date_num)  
  ) +
  scale_y_log10(
    name = "Plastic Concentration (g/L)", 
    breaks = pseudocount_concentrations_C,  
    labels = original_concentrations_C,  
    expand = c(0, 0),
    limits = c(min(pseudocount_concentrations_C), max(pseudocount_concentrations_C))  
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  
    axis.text.y = element_text(angle = 0),
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.15, "cm"),  
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines")  
  ) +
  guides(fill = "none") +  
  labs(title = "TPU FC2.1")  

#Print the adjusted heatmap plot for C without the legend
print(heatmap_contour_plot_C_log_y)


#Format changes for final combined plot
heatmap_contour_plot_C_log_y <- ggplot(interp_dfC, aes(x = date_num, y = 10^log_concentration, fill = Mean)) +
  geom_raster(interpolate = TRUE) +  #Continuous heatmap
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  
  #Customize x-axis
  scale_x_continuous(
    name = "Day Number", 
    breaks = seq(10, 90, by = 10),  
    expand = c(0, 0),
    limits = range(interp_dfC$date_num) 
  ) +
  
  #Customize y-axis
  scale_y_log10(
    breaks = pseudocount_concentrations_C,  #Set breaks with pseudocount
    labels = original_concentrations_C,  
    expand = c(0, 0),
    limits = c(min(pseudocount_concentrations_C), max(pseudocount_concentrations_C))  
  ) +
  
  # Customize the theme
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 0, vjust = 0.5, hjust = 0.5),  
    axis.text.y = element_blank(),  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_blank(),  
    axis.ticks.y = element_blank(),  
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.15, "cm"),  
    plot.title = element_text(size = 16, hjust = 0.5),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.margin = unit(c(1, 1, 1, 1), "lines")  
  ) +
  
  guides(fill = "none") +  
  labs(title = "TPU FC2.1")  

#Print the adjusted heatmap plot for C without the legend
print(heatmap_contour_plot_C_log_y)
####Error is associated with interpolation beyond the range of observed values. Ok to move forward. Collected data represented.  


#Create the combined plot
combined_heatmap_plot <- heatmap_contour_plot_A_log_y + 
  heatmap_contour_plot_B_log_y + 
  heatmap_contour_plot_C_log_y +
  plot_layout(ncol = 3)  

#Display the combined plot
print(combined_heatmap_plot)


#Fix the x axis for combined plot
heatmap_contour_plot_A_log_y <- heatmap_contour_plot_A_log_y + 
  theme(axis.title.x = element_blank())  

#No changes needed for plot B

#Update heatmap plot C
heatmap_contour_plot_C_log_y <- heatmap_contour_plot_C_log_y + 
  theme(axis.title.x = element_blank()) 

#Combine plots horizontally
combined_heatmap_plot <- heatmap_contour_plot_A_log_y + 
  heatmap_contour_plot_B_log_y + 
  heatmap_contour_plot_C_log_y +
  plot_layout(ncol = 3)  

#Display the combined plot 
print(combined_heatmap_plot) ### MANUSCRIPT [Fig. 1]





#Start the GAMs analysis
#Convert date to a numeric scale (number of days since the start of the study)
data_long$date <- as.Date(data_long$date, format = "%m/%d/%Y")
start_date <- min(data_long$date)
data_long$date_num <- as.numeric(data_long$date - start_date)

#Change letter_group to a factor and number_designation to numeric
data_long$letter_group <- as.factor(data_long$letter_group)
data_long$number_designation <- as.numeric(data_long$number_designation)

#Change date_num values to actual concentrations
data_long <- data_long %>%
  mutate(number_designation = case_when(
    number_designation == 0 ~ 0.000,
    number_designation == 1 ~ 0.004,
    number_designation == 2 ~ 0.008,
    number_designation == 3 ~ 0.013,
    number_designation == 4 ~ 0.023,
    number_designation == 5 ~ 0.041,
    number_designation == 6 ~ 0.072,
    number_designation == 7 ~ 0.126,
    number_designation == 8 ~ 0.220,
    number_designation == 9 ~ 0.385,
    TRUE ~ number_designation  
  ))

#check data distribution to meat assumptions
#Histogram
hist(data_long$CHLA_UGL, main="Histogram of CHLA_UGL", xlab="CHLA_UGL")
#Density plot
plot(density(na.omit(data_long$CHLA_UGL)), main="Density Plot of CHLA_UGL")
#Plastic type
ggplot(data_long, aes(x=letter_group, y=CHLA_UGL)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess", color="blue") +
  theme_minimal() +
  labs(title="Relationship between plastic type and CHLA_UGL", x="Predictor", y="CHLA_UGL")
#Day number
ggplot(data_long, aes(x=date_num, y=CHLA_UGL)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess", color="blue") +
  theme_minimal() +
  labs(title="Relationship between date_num and CHLA_UGL", x="Predictor", y="CHLA_UGL")
#Plastic concentration
ggplot(data_long, aes(x=number_designation, y=CHLA_UGL)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="loess", color="blue") +
  theme_minimal() +
  labs(title="Relationship between plastic concentration and CHLA_UGL", x="Predictor", y="CHLA_UGL")
#Check auto correlations
acf(data_long$date_num)

table(data_long$letter_group)


#Model 1 - Smooth terms for date_num and number_designation
model1_gamma <- gam(CHLA_UGL ~ 
                      s(date_num) + 
                      s(number_designation), 
                    family = Gamma(link = "log"), 
                    data = data_long)

#Model 2 - Add letter_group as a parametric effect
model2_gamma <- gam(CHLA_UGL ~ 
                      s(date_num) + 
                      s(number_designation) + 
                      letter_group, 
                    family = Gamma(link = "log"), 
                    data = data_long)

#Model 3 - Smooth terms vary by letter_group
model3_gamma <- gam(CHLA_UGL ~ 
                      s(date_num, by = letter_group) + 
                      s(number_designation, by = letter_group) + 
                      letter_group, 
                    family = Gamma(link = "log"), 
                    data = data_long)

#Model 3.5 - manual log of plastic concentration, Smooth terms vary by letter_group
model3.5_gamma <- gam(CHLA_UGL ~ 
                      s(date_num, by = letter_group) + 
                      s(log(number_designation + 0.001), by = letter_group) +  
                      letter_group, 
                    family = Gamma(link = "log"), 
                    data = data_long)

#Model 4 - Full model, allowing smooth terms to vary by letter_group
#And adding interaction with letter_group
model4_gamma <- gam(CHLA_UGL ~ 
                      s(date_num) +                         #Global smooth for date_num
                      s(date_num, by = letter_group) +      #Group-specific smooth for date_num
                      s(number_designation, by = letter_group) + #Group-specific smooth for number_designation
                      letter_group,                        #Parametric effect of letter_group
                    family = Gamma(link = "log"), 
                    data = data_long)

#Summaries and diagnostics for each model
models <- list(model1_gamma, model2_gamma, model3_gamma, model3.5_gamma, model4_gamma)
model_names <- c("Model 1", "Model 2", "Model 3", "Model 3.5", "Model 4")

for (i in seq_along(models)) {
  cat("Summary for", model_names[i], "\n")
  print(summary(models[[i]]))
  cat("\nDiagnostics for", model_names[i], "\n")
  gam.check(models[[i]])
}

#Compare models with AIC and degrees of freedom
aic_comparison_gamma <- data.frame(
  Model = model_names,
  AIC = sapply(models, AIC),
  df = sapply(models, function(m) m$df.residual)
)

#Print the AIC comparison
print(aic_comparison_gamma)

#ANOVA between models
anova_results <- anova(models[[1]], models[[2]], models[[3]], models[[4]], models[[5]])
print(anova_results)

#Plot the best model
best_model <- models[[which.min(aic_comparison_gamma$AIC)]]
plot(best_model, pages = 1)
summary(best_model)

#Add pairwise comparisons for each letter group

#Fitted vs. residuals
plot(fitted(model3_gamma), residuals(model3_gamma), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
#QQ plots
qqnorm(residuals(model3_gamma))
qqline(residuals(model3_gamma), col = "red")
#Auto correlations
acf(residuals(model3_gamma))
#Visualize smooth terms
plot(model3_gamma, pages = 1, all.terms = TRUE)

#Update the generate_predictions function to use model3_gamma
generate_predictions <- function(data) {
  preds <- predict(model3_gamma, newdata = data, type = "response", se = TRUE)
  data$fit <- preds$fit
  data$upper <- preds$fit + 1.96 * preds$se
  data$lower <- preds$fit - 1.96 * preds$se
  return(data)
}

#Sequences for predictors
num_designation_seq <- seq(from = min(data_long$number_designation), 
                           to = max(data_long$number_designation), 
                           length.out = 100)
date_num_seq <- seq(from = min(data_long$date_num), 
                    to = max(data_long$date_num), 
                    length.out = 100)

#Grid for number_designation predictions
prediction_data <- expand.grid(number_designation = num_designation_seq,
                               date_num = mean(data_long$date_num),  
                               letter_group = factor(c("A", "B", "C")))

#Grid for date_num predictions
prediction_data_date <- expand.grid(number_designation = mean(data_long$number_designation),  
                                    date_num = date_num_seq,
                                    letter_group = factor(c("A", "B", "C")))

#Apply function prediction data
prediction_data <- generate_predictions(prediction_data)
prediction_data_date <- generate_predictions(prediction_data_date)

#Plot for plastic concentration - Group A
p1 <- ggplot(prediction_data[prediction_data$letter_group == "A", ], aes(x = number_designation, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(title = "Plastic concentration: Group A", x = "Plastic concentration", y = "Fitted Value") +
  theme_minimal()

#Plot for plastic Concentration - Group B
p2 <- ggplot(prediction_data[prediction_data$letter_group == "B", ], aes(x = number_designation, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(title = "Plastic Concentration: Group B", x = "Plastic Concentration", y = "Fitted Value") +
  theme_minimal()

#Plot for plastic concentration - Group C
p3 <- ggplot(prediction_data[prediction_data$letter_group == "C", ], aes(x = number_designation, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(title = "Plastic Concentration: Group C", x = "Plastic Concentration", y = "Fitted Value") +
  theme_minimal()

#Plot for Time - Group A
p4 <- ggplot(prediction_data_date[prediction_data_date$letter_group == "A", ], aes(x = date_num, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(title = "Time: Group A", x = "Time", y = "Fitted Value") +
  theme_minimal()

#Plot for Time - Group B
p5 <- ggplot(prediction_data_date[prediction_data_date$letter_group == "B", ], aes(x = date_num, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(title = "Time: Group B", x = "Time", y = "Fitted Value") +
  theme_minimal()

#Plot for Time - Group C
p6 <- ggplot(prediction_data_date[prediction_data_date$letter_group == "C", ], aes(x = date_num, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(title = "Time: Group C", x = "Time", y = "Fitted Value") +
  theme_minimal()

#Arrange the plots
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)

#Post-hoc test for comparisons of B and C for model3_gamma
emm <- emmeans(model3_gamma, ~ letter_group)

#Pairwise comparisons
pairwise_comp <- pairs(emm)
print(pairwise_comp)

#Comparison specifically between B and C
BC_comparison <- contrast(emm, method = "revpairwise", levels = c("B", "C"))
print(BC_comparison)





#Model 1 - Smooth terms for date_num and number_designation
model1_gamma <- gam(CHLA_UGL ~ 
                      s(date_num) + 
                      s(number_designation), 
                    family = Gamma(link = "log"), 
                    data = data_long)

#Model 2-  Add letter_group as a parametric effect
model2_gamma <- gam(CHLA_UGL ~ 
                      s(date_num) + 
                      s(number_designation) + 
                      letter_group, 
                    family = Gamma(link = "log"), 
                    data = data_long)

#Model 3 - Smooth terms vary by letter_group
model3_gamma <- gam(CHLA_UGL ~ 
                      s(date_num, by = letter_group) + 
                      s(number_designation, by = letter_group) + 
                      letter_group, 
                    family = Gamma(link = "log"), 
                    data = data_long)

#Model 4 - Full model, allowing smooth terms to vary by letter_group
#Adding interaction with letter_group
model4_gamma <- gam(CHLA_UGL ~ 
                      s(date_num) +                         #Global smooth for date_num
                      s(date_num, by = letter_group) +      #Group-specific smooth for date_num
                      s(number_designation, by = letter_group) + #Group-specific smooth for number_designation
                      letter_group,                        #Parametric effect of letter_group
                    family = Gamma(link = "log"), 
                    data = data_long)

#Summaries and diagnostics for each model
models <- list(model1_gamma, model2_gamma, model3_gamma, model4_gamma)
model_names <- c("Model 1", "Model 2", "Model 3", "Model 4")

for (i in seq_along(models)) {
  cat("Summary for", model_names[i], "\n")
  print(summary(models[[i]]))
  cat("\nDiagnostics for", model_names[i], "\n")
  gam.check(models[[i]])
}

#Compare models with AIC and degrees of freedom
aic_comparison_gamma <- data.frame(
  Model = model_names,
  AIC = sapply(models, AIC),
  df = sapply(models, function(m) m$df.residual)
)

#Print AIC comparison
print(aic_comparison_gamma) ####USed for Chla MANUSCRIPT [Extended Data Table 1]###

#ANOVA between models
anova_results <- anova(models[[1]], models[[2]], models[[3]], models[[4]])
print(anova_results)

#Plot best model
best_model <- models[[which.min(aic_comparison_gamma$AIC)]]
plot(best_model, pages = 1)
summary(best_model) ####Used for Chla MANUSCRIPT [Extended Data Table 2]###
#anova.gam(best_model)

#Add pairwise comparisons for each letter group
#Get estimated marginal means (EMMs) for plastic type
plastic_emm <- emmeans(best_model, ~ letter_group)

#Pairwise comparisons with Tukey correction
pairwise_results <- contrast(plastic_emm, method = "pairwise", adjust = "tukey")

#View the pairwise comparisons
summary(pairwise_results)
###############################################

#Fitted vs. residuals
plot(fitted(model3_gamma), residuals(model3_gamma), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
#QQ plots
qqnorm(residuals(model3_gamma))
qqline(residuals(model3_gamma), col = "red")
#Auto correlations
acf(residuals(model3_gamma))
#Visualize smooth terms
plot(model3_gamma, pages = 1, all.terms = TRUE)

#Update the generate_predictions function to use model3_gamma
generate_predictions <- function(data) {
  preds <- predict(model3_gamma, newdata = data, type = "response", se = TRUE)
  data$fit <- preds$fit
  data$upper <- preds$fit + 1.96 * preds$se
  data$lower <- preds$fit - 1.96 * preds$se
  return(data)
}

#Sequences for predictors
num_designation_seq <- seq(from = min(data_long$number_designation), 
                           to = max(data_long$number_designation), 
                           length.out = 100)
date_num_seq <- seq(from = min(data_long$date_num), 
                    to = max(data_long$date_num), 
                    length.out = 100)

#Grid for number_designation predictions
prediction_data <- expand.grid(number_designation = num_designation_seq,
                               date_num = mean(data_long$date_num),  # Placeholder, will vary later
                               letter_group = factor(c("A", "B", "C")))

#Grid for date_num predictions
prediction_data_date <- expand.grid(number_designation = mean(data_long$number_designation),  # Placeholder
                                    date_num = date_num_seq,
                                    letter_group = factor(c("A", "B", "C")))

#Apply function to your prediction data
prediction_data <- generate_predictions(prediction_data)
prediction_data_date <- generate_predictions(prediction_data_date)

#Plot for plastic concentration - Group A
p1 <- ggplot(prediction_data[prediction_data$letter_group == "A", ], aes(x = number_designation, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(title = "Plastic concentration: Group A", x = "Plastic concentration", y = "Fitted Value") +
  theme_minimal()

# Plot for plastic concentration - Group B
p2 <- ggplot(prediction_data[prediction_data$letter_group == "B", ], aes(x = number_designation, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(title = "Plastic Concentration: Group B", x = "Plastic Concentration", y = "Fitted Value") +
  theme_minimal()

#Plot for plastic concentration - Group C
p3 <- ggplot(prediction_data[prediction_data$letter_group == "C", ], aes(x = number_designation, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(title = "Plastic Concentration: Group C", x = "Plastic Concentration", y = "Fitted Value") +
  theme_minimal()

#Plot for time - Group A
p4 <- ggplot(prediction_data_date[prediction_data_date$letter_group == "A", ], aes(x = date_num, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(title = "Time: Group A", x = "Time", y = "Fitted Value") +
  theme_minimal()

#Plot for time - Group B
p5 <- ggplot(prediction_data_date[prediction_data_date$letter_group == "B", ], aes(x = date_num, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(title = "Time: Group B", x = "Time", y = "Fitted Value") +
  theme_minimal()

#Plot for time - Group C
p6 <- ggplot(prediction_data_date[prediction_data_date$letter_group == "C", ], aes(x = date_num, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "blue") +
  labs(title = "Time: Group C", x = "Time", y = "Fitted Value") +
  theme_minimal()

#Arrange plots
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)

#Post-hoc test for comparisons of B and C for model3_gamma
emm <- emmeans(model3_gamma, ~ letter_group)

#Pairwise comparisons
pairwise_comp <- pairs(emm)
print(pairwise_comp)

#Comparison between B and C
BC_comparison <- contrast(emm, method = "revpairwise", levels = c("B", "C"))
print(BC_comparison)






#Group A and Group B
chla_diff_AB_date <- plot_difference(
  model3_gamma,
  series = date_num,
  difference = list(letter_group = c("A", "B"))
) + 
  labs(x = "Day Number", y = "Chlorophyll-a (Difference smooth)") +
  ggtitle("Difference between Elastollan and TPU 181") +
  theme_minimal()

#Group A and Group C
chla_diff_AC_date <- plot_difference(
  model3_gamma,
  series = date_num,
  difference = list(letter_group = c("A", "C"))
) + 
  labs(x = "Day Number", y = "Chlorophyll-a (Difference smooth)") +
  ggtitle("Difference between Elastollan and TPU FC2.1") +
  theme_minimal()

#Group B and Group C
chla_diff_BC_date <- plot_difference(
  model3_gamma,
  series = date_num,
  difference = list(letter_group = c("B", "C"))
) + 
  labs(x = "Day Number", y = "Chlorophyll-a (Difference smooth)") +
  ggtitle("Difference between TPU181 and TPU FC2.1") +
  theme_minimal()

#Group A and Group B
chla_diff_AB_num <- plot_difference(
  model3_gamma,
  series = number_designation,
  difference = list(letter_group = c("A", "B"))
) + 
  labs(x = "Plastic Concentration", y = "Chlorophyll-a (Difference smooth)") +
  ggtitle("Difference between Elastollan and TPU181") +
  theme_minimal()

#Group A and Group C
chla_diff_AC_num <- plot_difference(
  model3_gamma,
  series = number_designation,
  difference = list(letter_group = c("A", "C"))
) + 
  labs(x = "Plastic Concentration", y = "Chlorophyll-a (Difference smooth)") +
  ggtitle("Difference between Elastollan and TPU FC2.1") +
  theme_minimal()

#Group B and Group C
chla_diff_BC_num <- plot_difference(
  model3_gamma,
  series = number_designation,
  difference = list(letter_group = c("B", "C"))
) + 
  labs(x = "Plastic Concentration", y = "Chlorophyll-a (Difference smooth)") +
  ggtitle("Difference between TPU 181 and TPU FC2.1") +
  theme_minimal()

####Used for MANUSCRIPT [Supplementary Fig. 9]###
print(chla_diff_AB_date)
print(chla_diff_AC_date)
print(chla_diff_BC_date)
print(chla_diff_AB_num)
print(chla_diff_AC_num)
print(chla_diff_BC_num)
