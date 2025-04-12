#Clear if need
rm(list = ls())

setwd("C:/Users/DELL/Documents/R/projects/MP1_Zoop")
data<-read.csv("MP Zooplankton Counts - zoop_counts.csv")

#Install additional packages
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("tidymv")
install.packages("patchwork")
install.packages("gridExtra")
install.packages("pals")
install.packages("RColorBrewer")
install.packages("lubridate")
install.packages("lme4")
install.packages("vegan")
install.packages("mgcv")
install.packages("RVAideMemoire")
install.packages("car")
install.packages("viridisLite")
install.packages("mgcv")
install.packages("margins")
install.packages("MASS")
install.packages("ggeffects")
install.packages("glmmTMB")
install.packages("kableExtra")
install.packages("knitr")
install.packages("emmeans")
install.packages("patchwork")

update.packages(ask = FALSE)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidymv)
library(pals)
library(patchwork)
library(gridExtra)
library(RColorBrewer)
library(lubridate)
library(lme4)
library(vegan)
library(mgcv)
library(RVAideMemoire)
library(car)
library(viridisLite)



#Sort data in data frame w/pipe format
data%>%arrange(data$sample.date)

data%>%arrange(data$tank.id)

#Reset data for date and tank.ID
data<-data%>%arrange(data$sample.date,data$tank.id)


#Mutate the subsample for tank.id
data<- data%>% mutate(tank.id = case_when (
  tank.id == "Sub1 C2" ~ "C2",
  tank.id == "Sub2 C2" ~ "C2",
  tank.id == "Sub1 C3" ~ "C3",
  tank.id == "Sub2 C3" ~ "C3",
  TRUE ~ tank.id
))

#Add unique position to each row based on tank.id
data <- data %>%
  mutate(tank.number = case_when(
    tank.id == "A0" ~ "1",
    tank.id == "C3" ~ "2",
    tank.id == "B6" ~ "3",
    tank.id == "C9" ~ "4",
    tank.id == "A3" ~ "5",
    tank.id == "B2" ~ "6",
    tank.id == "C5" ~ "7",
    tank.id == "B3" ~ "8",
    tank.id == "A1" ~ "9",
    tank.id == "B9" ~ "10",
    tank.id == "B5" ~ "11",
    tank.id == "C0" ~ "12",
    tank.id == "A9" ~ "13",
    tank.id == "C6" ~ "14",
    tank.id == "B0" ~ "15",
    tank.id == "B8" ~ "16",
    tank.id == "A6" ~ "17",
    tank.id == "C2" ~ "18",
    tank.id == "A8" ~ "19",
    tank.id == "B7" ~ "20",
    tank.id == "C7" ~ "21",
    tank.id == "C8" ~ "22",
    tank.id == "A7" ~ "23",
    tank.id == "A5" ~ "24",
    tank.id == "C1" ~ "25",
    tank.id == "B4" ~ "26",
    tank.id == "A4" ~ "27",
    tank.id == "B1" ~ "28",
    tank.id == "A2" ~ "29",
    tank.id == "C4" ~ "30",
  ))

#Getting rid of weird columns (cleaning)
data<-subset(data, select = -c(photo.analyzed, X,X.1,X.2,X.3))

#Copy column taxa code and counts + update counts based on the prop.counted
data$new.taxa = data$taxa.code
data$new.count = data$count
data <- data %>% mutate(new.count = new.count / prop.counted)

#Column for plastic type
data$plastic.type=data$tank.id


data<- data %>% mutate(plastic.type = case_when(
  plastic.type == "A0" ~ "No Plastic",
  plastic.type == "A1" ~ "Petroleum Based TPU",
  plastic.type == "A2" ~ "Petroleum Based TPU",
  plastic.type == "A3" ~ "Petroleum Based TPU",
  plastic.type == "A4" ~ "Petroleum Based TPU",
  plastic.type == "A5" ~ "Petroleum Based TPU",
  plastic.type == "A6" ~ "Petroleum Based TPU",
  plastic.type == "A7" ~ "Petroleum Based TPU",
  plastic.type == "A8" ~ "Petroleum Based TPU", 
  plastic.type == "A9" ~ "Petroleum Based TPU",
  plastic.type == "B0" ~ "No Plastic",
  plastic.type == "B1" ~ "Algal Based TPU",
  plastic.type =="B2" ~ "Algal Based TPU",
  plastic.type == "B3" ~ "Algal Based TPU",
  plastic.type == "B4" ~ "Algal Based TPU",
  plastic.type == "B5" ~ "Algal Based TPU",
  plastic.type == "B6" ~ "Algal Based TPU",
  plastic.type == "B7" ~ "Algal Based TPU",
  plastic.type == "B8" ~ "Algal Based TPU", 
  plastic.type == "B9" ~ "Agal Based TPU",
  plastic.type =="C0" ~ "No Plastic",
  plastic.type == "C1" ~ "Biodegradable TPU",
  plastic.type == "C2" ~ "Biodegradable TPU",
  plastic.type == "C3" ~ "Biodegradable TPU",
  plastic.type == "C4" ~ "Biodegradable TPU",
  plastic.type == "C5" ~ "Biodegradable TPU",
  plastic.type == "C6" ~ "Biodegradable TPU",
  plastic.type == "C7" ~ "Biodegradable TPU",
  plastic.type == "C8" ~ "Biodegradable TPU", 
  plastic.type == "C9" ~ "Biodegradable TPU",
))



#Sorting the taxa 
unique(data$taxa.code)
unique(data$new.taxa)

data <- data %>% mutate(new.taxa = case_when (
  taxa.code == "chy" ~ "cladoceran", 
  taxa.code == "bos" ~ "cladoceran",
  taxa.code == "eur" ~ "cladoceran",
  taxa.code == "unc" ~ "cladoceran",
  taxa.code == "dap" ~ "daphnia",
  taxa.code == "ost" ~ "ostracoda",
  taxa.code == "ostra" ~ "ostracoda",
  taxa.code == "hydra" ~ "mite",
  taxa.code == "cyc" ~ "large copepod",
  taxa.code == "euy" ~ "large copepod",
  taxa.code == "cal" ~ "large copepod",
  taxa.code == "erg" ~ "large copepod",
  taxa.code == "nau" ~ "nauplii",
  taxa.code == "Un7" ~ "rotifer",
  taxa.code == "plm" ~ "rotifer",
  taxa.code == "ker" ~ "rotifer",
  taxa.code == "ephi" ~ "ephippia",
  taxa.code == "mon" ~ "rotifer",
  taxa.code == "zor" ~ "insect",
  taxa.code == "in1" ~ "insect",
  taxa.code == "uni" ~ "insect",
  taxa.code == "i16" ~ "insect",
  taxa.code == "uix" ~ "insect",
  taxa.code == "i17" ~ "insect",
  taxa.code == "in5" ~ "insect",
  taxa.code == "mqa" ~ "insect",
  taxa.code == "dyt" ~ "insect",
  taxa.code == "spi" ~ "insect",
  taxa.code == "uniz" ~ "insect",
  taxa.code == "uniy" ~ "insect",
  taxa.code == "ant" ~ "insect",
  taxa.code == "in11" ~ "insect",
  taxa.code == "mql" ~ "insect",
  taxa.code == "chi" ~ "insect",
  taxa.code == "chir" ~ "insect",
  taxa.code == "mos" ~ "insect",
  taxa.code == "thy" ~ "insect",
  taxa.code == "thrip" ~ "insect",
  taxa.code == "thrips" ~ "insect",
  taxa.code == "mgl" ~ "insect",
  taxa.code == "i15" ~ "insect",
  taxa.code == "in2" ~ "insect",
  taxa.code == "bzn" ~ "insect",
  taxa.code == "bzr" ~ "insect",
  taxa.code == "i12" ~ "insect",
  taxa.code == "un6" ~ "unknown",
  taxa.code == "unkown" ~ "unknown",
  taxa.code == "unk 1" ~ "unknown",
  taxa.code == "unk 2" ~ "unknown",
  taxa.code == "unk 3" ~ "unknown",
  taxa.code == "u13" ~ "unknown",
  taxa.code == "un12" ~ "unknown",
  taxa.code == "un2" ~ "unknown",
  taxa.code == "u14" ~ "unknown",
  taxa.code == "un0" ~ "unknown",
  taxa.code == "uno" ~ "unknown",
  taxa.code == "unk 4" ~ "unknown",
  taxa.code == "unk" ~ "unknown",
  taxa.code == "un13" ~ "unknown",
  taxa.code == "un7" ~ "unknown",
  TRUE ~ taxa.code 
))

print(unique(data$new.taxa))


#Now get the unique dates and sort them
unique_dates <- sort(unique(data$sample.date))
print(unique_dates)

#Add sample.period based on sample.date
data <- data %>%
  mutate(sample.period = case_when(
    sample.date == as.Date("2023-03-07") ~ "0",
    sample.date == as.Date("2023-03-13") ~ "6",
    sample.date == as.Date("2023-03-23") ~ "16",
    sample.date == as.Date("2023-04-12") ~ "36",
    sample.date == as.Date("2023-05-02") ~ "56",
    sample.date == as.Date("2023-06-12") ~ "97",
    TRUE ~ NA_character_  #Assign NA to any date that doesn't match
  ))




#Subset data and check how many entries are in each sample.period
sample1 <- data %>%
  filter(sample.period == "0")

#Convert to numeric and then get unique sorted values for tank.number in sample1
unique_tank_numbers_sample1 <- sort(as.numeric(unique(sample1$tank.number)))
print(unique_tank_numbers_sample1)

sample4 <- data %>%
  filter(sample.period == "36")
unique_tank_numbers_sample4<-sort(as.numeric(unique(sample4$tank.number)))
print(unique_tank_numbers_sample4)

#Convert to numeric and then get unique sorted values for tank.number in sample4
unique_tank_numbers_sample4 <- sort(as.numeric(unique(sample4$tank.number)))
print(unique_tank_numbers_sample4)

#sampling 5/2/23 for tank numbers
sample5 <- data%>%
  filter(sample.period == "56")
unique_tank_numbers_sample5 <- sort(as.numeric(unique(sample5$tank.number)))
print(unique_tank_numbers_sample5)


#Convert to numeric and then get unique sorted values for tank.number in sample6
sample6 <- data %>%
  filter(sample.period == "97")
unique_tank_numbers_sample6 <- sort(as.numeric(unique(sample6$tank.number)))
print(unique_tank_numbers_sample6)

#Summarize the count by sample_period, channel_id, and taxa_code
data_wide <- data %>%
  group_by(sample.period, tank.id, new.taxa, plastic.type) %>%
  summarise(count = sum(new.count)) %>%
  pivot_wider(names_from = new.taxa, values_from = count, values_fill = list(count = 0))


data_long <- data_wide %>%
  pivot_longer(
    cols = -c(sample.period, tank.id, plastic.type),
    names_to = "Taxa",
    values_to = "Count"
  )

#If you want to only display specific taxa codes (removing ephippia)
print(unique(data_long$Taxa))

desired_values <- c("cladoceran","daphnia","nauplii","unknown","insect",
                    "large copepod", "rotifer","mite", "ostracoda")
data_long <- data_long %>%
  filter( Taxa%in% desired_values)


#Generate the color palette after filtering the data
unique_entries <- unique(data$new.taxa)

#color palette
color_palette <- pals::watlington(length(unique_entries))



#Convert sample.period to factor if it's not already
data_long$sample.period <- factor(data_long$sample.period, levels = unique(data_long$sample.period))

data_long$sample.period

#Filter the data for sample periods 1, 4, 5
filtered_data <- data_long %>%
  filter(sample.period %in% c("0", "36", "56", "97"))

#ggplot

c.plot<-ggplot(filtered_data, aes(x = sample.period, y = (log((Count)/4))^1.0001, fill = Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual( values = color_palette,
                     unique_entries,
                     guide = guide_legend( title = "Taxa")) +
  facet_wrap(~ tank.id, nrow=3) +
  labs(title= "Zooplankton Concentrations Over 97 Days", y = "log(Number of zooplankton per Liter)", x = "Days of Plastic Exposure") +
  
  theme_minimal()
c.plot
#Warning is zero counts

#Do the same figure but give manual colors
#Define your color palette with 9 color
color_palette_manual <- c("midnightblue", "green4", "lightgreen", "darkorchid", "magenta", "pink4", "orange2", "red2", "yellow1")

#Plot
c.plot.color <- ggplot(filtered_data, aes(x = sample.period, y = (log((Count)/4))^1.0001, fill = Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette_manual,
                    guide = guide_legend(title = "Taxa")) +
  facet_wrap(~ tank.id, nrow = 3) +
  labs(title = "Zooplankton Concentrations Over 97 Days", 
       y = "log(Number of zooplankton per Liter)", 
       x = "Days of Plastic Exposure") +
  theme_minimal()

c.plot.color
#Warning is zero counts

####Try it with the plasma color palette####

#Define color palette using the plasma color scale from viridisLite
color_palette_plasma <- plasma(9)

#Plot
c.plot.plasma <- ggplot(filtered_data, aes(x = sample.period, y = (log((Count)/4))^1.0001, fill = Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette_plasma,
                    guide = guide_legend(title = "Taxa")) +
  facet_wrap(~ tank.id, nrow = 3) +
  labs(title = "Zooplankton concentrations over 97 day exsposure period", 
       y = "log(Zooplankton abundance [individuals/L])", 
       x = "Days of plastic exposure") +
  theme_minimal()

c.plot.plasma

################Try to manually change all plot aesthetics######################

#Define color palette using the plasma color scale from viridisLite
color_palette_plasma <- plasma(9)

#Define the labels for the legend components for custom names
legend_labels <- c("Cladoceran", "Daphnia", "Insecta", "Large Copepod", "Mite", "Nauplii", "Ostracoda", "Rotifera", "Unknown")

#Plot
c.plot.plasma.custom <- ggplot(filtered_data, aes(x = sample.period, y = (log((Count)/4))^1.0001, fill = Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette_plasma,
                    labels = legend_labels,
                    guide = guide_legend(title = "Taxon")) +
  facet_wrap(~ tank.id, nrow = 3) +
  labs(title = "Zooplankton concentrations [Individuals/L] over 97 days", 
       y = "Abundance - log(Zooplankton per liter)", 
       x = "Time - (Days of plastic exsposure)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"), # Center and bold title
    axis.title.x = element_text(size = 13), # Custom x-axis title
    axis.title.y = element_text(size = 13), # Custom y-axis title
    legend.title = element_text(face = "bold", size = 10), # Custom legend title
    legend.text = element_text(size = 9) # Custom legend text size
  )

c.plot.plasma.custom

####Keep editing format####
#Define color palette using the plasma color scale from viridisLite
color_palette_plasma <- plasma(9)

#Plot
c.plot.plasma <- ggplot(filtered_data, aes(x = sample.period, y = (log((Count)/4))^1.0001, fill = Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette_plasma,
                    guide = guide_legend(title = "Taxa")) +
  facet_wrap(~ tank.id, nrow = 3) +
  labs(title = "Zooplankton abundance with 97 day exsposure", 
       y = "log(Zooplnakton abundance [individuals/L])", 
       x = "Days of plastic exposure") +
  theme_minimal() +
  theme(
    legend.title = element_text(face = "bold"),  # Bold legend title
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Box around each plot
    panel.grid = element_blank()  # Remove gridlines
  )

#Display the plot
c.plot.plasma


####Combine formats####
#Define color palette using the plasma color scale from viridisLite
color_palette_plasma <- plasma(9)

#Define the labels for the legend components if you need custom names
legend_labels <- c("Cladoceran", "Daphnia", "Insecta", "Large Copepod", "Mite", "Nauplii", "Ostracoda", "Rotifera", "Unknown")

#Plot
c.plot.plasma.custom <- ggplot(filtered_data, aes(x = sample.period, y = (log((Count)/4))^1.0001, fill = Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette_plasma,
                    labels = legend_labels,
                    guide = guide_legend(title = "Taxon")) +
  facet_wrap(~ tank.id, nrow = 3) +
  labs( 
       y = "log(Zooplankton abundance [individuals/L])", 
       x = "Day Number") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 13, face = "bold"),  # Custom x-axis title
    axis.title.y = element_text(size = 13, face = "bold"),  # Custom y-axis title
    legend.title = element_text(face = "bold", size = 10),  # Bold and custom size legend title
    legend.text = element_text(size = 9),  # Custom legend text size
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Box around each plot
    panel.grid = element_blank()  # Remove gridlines
  )

#Display the plot
c.plot.plasma.custom

#One more
#Define color palette with 9 color
color_palette_manual <- c("deepskyblue3", "green4", "lightgreen", "darkorchid", "magenta", "midnightblue", "orange2", "red2", "firebrick")

#Define the labels for the legend components for custom names
legend_labels <- c("Cladoceran", "Daphnia", "Insecta", "Large Copepod", "Mite", "Nauplii", "Ostracoda", "Rotifera", "Unknown")

#Plot
c.plot.new.custom <- ggplot(filtered_data, aes(x = sample.period, y = (log((Count)/4))^1.0001, fill = Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette_manual,
                    labels = legend_labels,
                    guide = guide_legend(title = "Taxon")) +
  facet_wrap(~ tank.id, nrow = 3) +
  labs( 
    y = "log(Zooplankton abundance [individuals/L])", 
    x = "Day Number") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 13, face = "bold"),  # Custom x-axis title
    axis.title.y = element_text(size = 13, face = "bold"),  # Custom y-axis title
    legend.title = element_text(face = "bold", size = 10),  # Bold and custom size legend title
    legend.text = element_text(size = 9),  # Custom legend text size
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Box around each plot
    panel.grid = element_blank()  # Remove gridlines
  )

#Display the plot
c.plot.new.custom ### Zoop data for MANUSCRIPT [Figure S4]
#Warning is zero counts

#taxa counts
#time point 1
abundance_counts <- filtered_data %>%
  filter(sample.period=="0")%>%
  group_by(tank.id, Taxa) %>%
  summarise(Count = sum(Count)) %>%
  mutate(Count = sum(Count))

abundance_counts<-abundance_counts %>%
  filter(Taxa =="cladoceran") 

abundance_counts<- abundance_counts %>%
  mutate( Taxa= case_when(
    Taxa== "cladoceran" ~ "organism")) %>%
  mutate(sample.period = "1")




#time point 4
abundance_counts1 <- filtered_data %>%
  filter(sample.period=="36")%>%
  group_by(tank.id, Taxa) %>%
  summarise(Count = sum(Count)) %>%
  mutate(Count = sum(Count))

abundance_counts1<-abundance_counts1 %>%
  filter(Taxa =="cladoceran") 

abundance_counts1 <- abundance_counts1 %>%
  mutate( Taxa= case_when(
    Taxa== "cladoceran" ~ "organism"))%>%
  mutate(sample.period = "4")

#time point 5
abundance_counts2 <- filtered_data %>%
  filter(sample.period=="56")%>%
  group_by(tank.id, Taxa) %>%
  summarise(Count = sum(Count)) %>%
  mutate(Count = sum(Count))

abundance_counts2<-abundance_counts2 %>%
  filter(Taxa =="cladoceran") 

abundance_counts2 <- abundance_counts2 %>%
  mutate( Taxa= case_when(
    Taxa== "cladoceran" ~ "organism"))%>%
  mutate(sample.period = "5")

#time point 6
abundance_counts3 <- filtered_data %>%
  filter(sample.period=="97")%>%
  group_by(tank.id, Taxa) %>%
  summarise(Count = sum(Count)) %>%
  mutate(Count = sum(Count))

abundance_counts3<-abundance_counts3 %>%
  filter(Taxa =="cladoceran") 

abundance_counts3 <- abundance_counts3 %>%
  mutate( Taxa= case_when(
    Taxa== "cladoceran" ~ "organism"))%>%
  mutate(sample.period = "6")


############Library additional packages for next set of code############
library(knitr)
library(kableExtra)
library(glmmTMB)
library(ggeffects)
library(margins)
library(mgcv)
library(MASS)
library(emmeans)
library(patchwork)

#combine the abundances
Total_abundance0 <- bind_rows(abundance_counts, abundance_counts1, abundance_counts2, abundance_counts3)
Total_abundance0$logcount<-log(Total_abundance0$Count)
Total_abundance0$count.per.l<-(Total_abundance0$Count)/4
Total_abundance <- subset(Total_abundance0, !grepl("0", tank.id))
Control_abundance <-subset(Total_abundance0, grepl("0", tank.id))
Control_abundance$logcount <- log(Control_abundance$Count)


#GOING TO MAKE THE COPPIES OF CONTROL ABUNDANCES
#linear including the negative controls
Control_abundance$treatment<- 0
Control_abundance$plastic_type<- Control_abundance$tank.id
Control_abundanceE<- Control_abundance
Control_abundanceE$tank.id<-"A0"
Control_abundance1<-Control_abundance
Control_abundance1$tank.id<-"B0"
Control_abundance2<-Control_abundance 
Control_abundance2$tank.id<-"C0"
Control_abundance<-rbind(Control_abundance1, Control_abundance2, Control_abundanceE)
Total_abundance0<- rbind(Total_abundance,Control_abundance)
#Giving each each tank the indv control tank for plastic.type

#Run this
#subsample for sample periods
Total_abundance1 <-subset(rbind(Total_abundance,Control_abundance), grepl("1", sample.period))
hist(Total_abundance1$Count)


#run this for the 4 periods
Total_abundance4 <-subset(rbind(Total_abundance,Control_abundance), grepl("4", sample.period))
hist(Total_abundance4$Count)




Total_abundance5 <-subset(rbind(Total_abundance,Control_abundance), grepl("5", sample.period))
hist(Total_abundance5$Count)


Total_abundance6 <-subset(rbind(Total_abundance,Control_abundance), grepl("6", sample.period))
hist(Total_abundance6$Count)







Total_abundance0$plastic_type<-Total_abundance0$tank.id
Total_abundance0$concentration <-Total_abundance0$tank.id

Total_abundance0<- Total_abundance0 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Total_abundance0<- Total_abundance0 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Total_abundance0$log.con<-Total_abundance0$tank.id

Total_abundance0<- Total_abundance0 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))



#LOG trans above
Total_abundance0$log.con<-as.numeric(Total_abundance0$log.con)
Total_abundance0$log.con<-log(Total_abundance0$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Total_abundance0$log.con)

hist(Total_abundance0$log.con)

########Total_abundance objects which are Total_abundance1, Total_abundance4, Total_abundance5, Total_abundance6

Total_abundance1$plastic_type<-Total_abundance1$tank.id
Total_abundance1$concentration <-Total_abundance1$tank.id

Total_abundance1<- Total_abundance1 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Total_abundance1<- Total_abundance1 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Total_abundance1$log.con<-Total_abundance1$tank.id

Total_abundance1<- Total_abundance1 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))



#LOG trans above
Total_abundance1$log.con<-as.numeric(Total_abundance1$log.con)
Total_abundance1$log.con<-log(Total_abundance1$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Total_abundance1$log.con)

hist(Total_abundance1$log.con)
#######################

Total_abundance4$plastic_type<-Total_abundance4$tank.id
Total_abundance4$concentration <-Total_abundance4$tank.id

Total_abundance4<- Total_abundance4 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Total_abundance4<- Total_abundance4 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Total_abundance4$log.con<-Total_abundance4$tank.id

Total_abundance4<- Total_abundance4 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))



#LOG TRANS ABOVE!!!!
Total_abundance4$log.con<-as.numeric(Total_abundance4$log.con)
Total_abundance4$log.con<-log(Total_abundance4$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Total_abundance4$log.con)

hist(Total_abundance4$log.con)

####################################
Total_abundance5$plastic_type<-Total_abundance5$tank.id
Total_abundance5$concentration <-Total_abundance5$tank.id

Total_abundance5<- Total_abundance5 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Total_abundance5<- Total_abundance5 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Total_abundance5$log.con<-Total_abundance5$tank.id

Total_abundance5<- Total_abundance5 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))



#LOG trans above
Total_abundance5$log.con<-as.numeric(Total_abundance5$log.con)
Total_abundance5$log.con<-log(Total_abundance5$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Total_abundance5$log.con)

hist(Total_abundance5$log.con)

###################################

Total_abundance6$plastic_type<-Total_abundance6$tank.id
Total_abundance6$concentration <-Total_abundance6$tank.id

Total_abundance6<- Total_abundance6 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Total_abundance6<- Total_abundance6 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Total_abundance6$log.con<-Total_abundance6$tank.id

Total_abundance6<- Total_abundance6 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))



#LOG trans above
Total_abundance6$log.con<-as.numeric(Total_abundance6$log.con)
Total_abundance6$log.con<-log(Total_abundance6$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Total_abundance6$log.con)

hist(Total_abundance6$log.con)


#factoring 
Total_abundance0$plastic_type<-as.factor(Total_abundance0$plastic_type)
Total_abundance0$concentration<- as.numeric(Total_abundance0$concentration)
Total_abundance0$log.con<-as.numeric(Total_abundance0$log.con)
#Total_abundance$concentration<- as.numeric(Total_abundance$concentration)

Total_abundance1$plastic_type<-as.factor(Total_abundance1$plastic_type)
Total_abundance1$concentration<- as.numeric(Total_abundance1$concentration)
Total_abundance1$log.con<- as.numeric(Total_abundance1$log.con)


Total_abundance4$plastic_type<-as.factor(Total_abundance4$plastic_type)
Total_abundance4$concentration<- as.numeric(Total_abundance4$concentration)
Total_abundance4$log.con<- as.numeric(Total_abundance4$log.con)


Total_abundance5$plastic_type<-as.factor(Total_abundance5$plastic_type)
Total_abundance5$concentration<- as.numeric(Total_abundance5$concentration)
Total_abundance5$log.con<- as.numeric(Total_abundance5$log.con)


Total_abundance6$plastic_type<-as.factor(Total_abundance6$plastic_type)
Total_abundance6$concentration<- as.numeric(Total_abundance6$concentration)
Total_abundance6$log.con<- as.numeric(Total_abundance6$log.con)






#distribution

hist(Total_abundance$Count, main= "Histogram of Zooplankton Counts From All Samples", xlab="Number of Zooplankton", ylab="Frequency")
hist(Total_abundance1$Count, main= "Histogram of Zooplankton Count From Day 0", xlab="Number of Zooplankton", ylab="Frequency")
hist(Total_abundance4$Count, main= "Histogram of Zooplankton Count From Day 36", xlab="Number of Zooplankton", ylab="Frequency")
hist(Total_abundance5$Count, main= "Histogram of Zooplankton Count", xlab="Number of Zooplankton", ylab="Frequency")
hist(Total_abundance6$Count, main= "Histogram of Zooplankton Count", xlab="Number of Zooplankton", ylab="Frequency")

#Model selection
#TIME POINT 0
test.run.t1<- gam(count.per.l~ plastic_type+s(log.con, by=plastic_type), data=Total_abundance1, method= 'REML', family = Gamma(link='log')) 

test2.run.t1<- gam(count.per.l~ s(log.con, by=plastic_type), data=Total_abundance1, method= 'REML', family = Gamma(link='log'))

test3.run.t1<- gam(count.per.l~ s(log.con), data=Total_abundance1, method= 'REML', Gamma(link='log')) 

test.run4.t1<- gam(count.per.l~ plastic_type+s(log.con), data=Total_abundance1, method= 'REML', family = Gamma(link='log')) 


AIC(test.run.t1,test2.run.t1,test3.run.t1,test.run4.t1) #Abundance data used for MANUSCRIPT [Table S3] this is time point 0
plot(test.run.t1, all.terms = TRUE, page=1)
plot(test2.run.t1, all.terms = TRUE, page=1) 

summary(test.run.t1)
summary(test2.run.t1)
summary(test3.run.t1) 
summary(test.run4.t1)

#----------------------------------------------------

#TIME POINT 4
test.run.t4<- gam(count.per.l~ plastic_type+s(concentration, by=plastic_type), data=Total_abundance4, method= 'REML', family = Gamma(link='log')) #FIRST IS THE BEST 
test2.run.t4<- gam(count.per.l~ s(concentration, by=plastic_type, k=10), data=Total_abundance4, method= 'REML', family = Gamma(link='log')) 
test3.run.t4<- gam(count.per.l~ s(concentration), data=Total_abundance4, method= 'REML', family = Gamma(link='log')) 
test4.run.t4<- gam(count.per.l~ plastic_type+s(concentration), data=Total_abundance4, method= 'REML', family = Gamma(link='log'))
summary(test.run.t4)

AIC(test.run.t4,test2.run.t4,test3.run.t4,test4.run.t4) #Abundance data used for MANUSCRIPT [Table S3]
plot(test.run.t4, all.terms = TRUE, page=1)
plot(test4.run.t4, all.terms = TRUE, page=1)

summary(test.run.t4) 
summary(test2.run.t4)
summary(test3.run.t4)
summary(test4.run.t4) 
#----------------------------------------------------

#TIME POINT 5
test.run.t5<- gam(count.per.l~ plastic_type+s(concentration, by=plastic_type), data=Total_abundance5, method= 'REML', family = Gamma(link='log')) #best model
test2.run.t5<- gam(count.per.l~ s(concentration, by=plastic_type), data=Total_abundance5, method= 'REML', family = Gamma(link='log')) 
test3.run.t5<- gam(count.per.l~ s(concentration), data=Total_abundance5, method= 'REML', family = Gamma(link='log')) 
test4.run.t5<- gam(count.per.l~ plastic_type+s(concentration), data=Total_abundance5, method= 'REML', family = Gamma(link='log')) 

AIC(test.run.t5,test2.run.t5,test3.run.t5,test4.run.t5) #Abundance data used for MANUSCRIPT [Table S3]
plot(test.run.t5, all.terms = TRUE, page=1)
plot(test2.run.t5, all.terms = TRUE, page=1)
plot(test3.run.t5, all.terms = TRUE, page=1)

summary(test.run.t5) 
summary(test2.run.t5)
summary(test3.run.t5)
summary(test4.run.t5)
#----------------------------------------------------

#TIME POINT 6
test.run.t6<- gam(count.per.l~ plastic_type+s(log.con, by=plastic_type), data=Total_abundance6, method= 'REML', family = Gamma(link='log')) 
test2.run.t6<- gam(count.per.l~ s(concentration, by=plastic_type), data=Total_abundance6, method= 'REML', family = Gamma(link='log')) 
test3.run.t6<- gam(count.per.l~ s(concentration), data=Total_abundance6, method= 'REML', family = Gamma(link='log')) 
test4.run.t6<- gam(count.per.l~ plastic_type+s(log.con), data=Total_abundance6, method= 'REML', family = Gamma(link='log')) 


AIC(test.run.t6,test2.run.t6,test3.run.t6,test4.run.t6) #Abundance data used for MANUSCRIPT [Table S3]
plot(test.run.t6, all.terms = TRUE, page=1)
plot(test2.run.t6, all.terms = TRUE, page=1)
plot(test3.run.t6, all.terms = TRUE, page=1)


summary(test.run.t6) 
summary(test2.run.t6)
summary(test3.run.t6)
summary(test4.run.t6)

#GAMs with lowest AIC
#CONGREGATE THE MODELS THAT WE ARE GOING TO USE AND FURTHER TEST
##########gam.t1 this is for time 0
gam.t1<- gam(count.per.l~ s(log.con), data=Total_abundance1, method= 'REML', family = Gamma(link='log'))#1
#gam.t1.0<- gam(count.per.l~ plastic_type+s(log.con, by=plastic_type), data=Total_abundance1, method= 'REML', family = Gamma(link='log'))#1


gam.t4<- gam(count.per.l~ plastic_type+s(log.con, by=plastic_type), data=Total_abundance4, method= 'REML', family = Gamma(link='log')) #4
gam.t5<- gam(count.per.l~ plastic_type+s(log.con, by=plastic_type), data=Total_abundance5, method= 'REML', family = Gamma(link='log')) #5
gam.t6<- gam(count.per.l~ plastic_type+s(log.con, by=plastic_type, k=6), data=Total_abundance6, method= 'REML', family = Gamma(link='log')) #6



hist(Total_abundance0$count.per.l)
#hist(Total_abundance0$logcount.per.l)

hist(Total_abundance1$count.per.l)
#hist(Total_abundance1$logcount.per.l)
hist(Total_abundance4$count.per.l)
#hist(Total_abundance4$logcount.per.l)
hist(Total_abundance5$count.per.l)
#hist(Total_abundance5$logcount.per.l)
hist(Total_abundance6$count.per.l)
#hist(Total_abundance6$logcount.per.l)


#Looking at summary and ANOVA.GAM and plot
AIC(gam.t1,gam.t4, gam.t5, gam.t6)

print(gam.t1)
print(gam.t4)
print(gam.t5)
print(gam.t6)

summary(gam.t1) ### Abundance data used for MANUSCRIPT [Table S4]
summary(gam.t4) ### Abundance data used for MANUSCRIPT [Table S4]
summary(gam.t5) ### Abundance data used for MANUSCRIPT [Table S4]
summary(gam.t6) ### Abundance data used for MANUSCRIPT [Table S4]

anova.gam(gam.t1)
anova.gam(gam.t4)
anova.gam(gam.t5)
anova.gam(gam.t6)

plot(gam.t1, all.terms = TRUE, page=1)
plot(gam.t4, all.terms = TRUE, page=1)
plot(gam.t5, all.terms = TRUE, page=1)
plot(gam.t6, all.terms = TRUE, page=1)

concurvity(gam.t1)#LOW CONCURVITY VALUES ARE DESIRABLE
concurvity(gam.t4)#LOW CONCURVITY VALUES ARE DESIRABLE
concurvity(gam.t5)#LOW CONCURVITY VALUES ARE DESIRABLE
concurvity(gam.t6)#LOW CONCURVITY VALUES ARE DESIRABLE

#FINAL GAM CHECK
gam.check(gam.t1)
gam.check(gam.t4)
gam.check(gam.t5)
gam.check(gam.t6)

#setting concentrations to look at (every concentration)
log_con_levels <- c(-0.9493306,-1.5050779,-2.0557250,-2.6036902,-3.1465552, -3.6888795,-4.1997051, -4.6051702, -5.1159958, -6.2146081)

# Calculate pairwise comparisons for each model
###Best pairwise comparisons for data
pw1 <- emmeans(gam.t1, ~ log.con, at = list(log.con = log_con_levels))
pw4 <- emmeans(gam.t4, ~ plastic_type | log.con, at = list(log.con = log_con_levels))
pw5 <- emmeans(gam.t5, ~ plastic_type | log.con, at = list(log.con = log_con_levels))
pw6 <- emmeans(gam.t6, ~ plastic_type | log.con, at = list(log.con = log_con_levels))

#Perform pairwise comparisons
pw1c <- pairs(pw1)
pw4c <- pairs(pw4)
pw5c <- pairs(pw5)
pw6c <- pairs(pw6)

#table to copy and paste into excel
table1<-kable(pw1c)
table4<-kable(pw4c)
table5<-kable(pw5c)
table6<-kable(pw6c)

view(pw1c)
view(pw4c)
view(pw5c)
view(pw6c)

#######################################################
###############Adding Summary Table Creation###########
#######################################################

#Define log concentrations and matching real concentrations
log_con_levels_sumtable <- c(-6.2146081, -5.1159958, -4.6051702, -4.1997051, -3.6888795, 
                    -3.1465552, -2.6036902, -2.0557250, -1.5050779, -0.9493306)

plastic_conc_sumtable <- c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385)

#Create a lookup table
log_conc_map <- data.frame(
  log_con = log_con_levels_sumtable,
  plastic_concentration = plastic_conc_sumtable
)

#Summarize each pairwise comparison object
summary_pw1 <- as.data.frame(summary(pw1c))  #Day 0 (before plastics)
summary_pw4 <- as.data.frame(summary(pw4c))  #Day 36
summary_pw5 <- as.data.frame(summary(pw5c))  #Day 56
summary_pw6 <- as.data.frame(summary(pw6c))  #Day 97

#Add sampling day to each summary
summary_pw1$sampling_day <- 0
summary_pw4$sampling_day <- 36
summary_pw5$sampling_day <- 56
summary_pw6$sampling_day <- 97

#Combine all summaries into one dataframe
full_summarytable <- bind_rows(summary_pw1, summary_pw4, summary_pw5, summary_pw6)

#Merge with real plastic concentration values
full_summarytable <- full_summarytable %>%
  left_join(log_conc_map, by = c("log.con" = "log_con"))

colnames(full_summarytable)

#Create the final clean summary table
finished_summarytable <- full_summarytable %>%
  mutate(
    Comparison = contrast,
    Estimate_SE = paste0(round(estimate, 3), " ± ", round(SE, 3)),
    p_value = round(p.value, 4),
    Significant = ifelse(p.value < 0.05, "Yes", "No")
  ) %>%
  dplyr::select(   # <--- IMPORTANT: Force dplyr::select
    sampling_day,
    plastic_concentration,
    Comparison,
    Estimate_SE,
    t.ratio,
    p_value,
    Significant
  ) %>%
  arrange(sampling_day, plastic_concentration)

#View the final table
kable(finished_summarytable)
view(finished_summarytable)

#Export to CSV
#write.csv(finished_summarytable, "zooplankton_pairwise_comparisons_final.csv", row.names = FALSE)

#Clean p-values
finished_summarytable <- finished_summarytable %>%
  mutate(
    p_value_clean = ifelse(p_value < 0.0001, "<0.0001", sprintf("%.4f", p_value)),
    Estimate_SE_clean = gsub("Â", "", Estimate_SE) #Remove character
  )

#Filter to show only significant comparisons
significant_summarytable <- finished_summarytable %>%
  filter(Significant == "Yes")

#Select columns and rename
final_summarytable <- significant_summarytable %>%
  dplyr::select(          #Force dplyr::select
    `Sampling Day` = sampling_day,
    `Plastic Concentration (g/L)` = plastic_concentration,
    `Comparison` = Comparison,
    `Estimate ± SE` = Estimate_SE_clean,
    `t-ratio` = t.ratio,
    `p-value` = p_value_clean,
    `Significant Difference` = Significant
  ) %>%
  arrange(`Sampling Day`, `Plastic Concentration (g/L)`)

#Print the final polished table
kable(final_summarytable, caption = "Supplementary Table X | Significant pairwise comparisons of zooplankton biomass across plastic treatments and concentrations.") %>%
  kable_styling(full_width = FALSE, position = "center")

#write.csv(final_summarytable, "zooplankton_pairwise_comparisons_final_clean.csv", row.names = FALSE)


####################################################
####################################################

fitted_values <- fitted(gam.t5)
residuals <- residuals(gam.t5, type = "response")

#Plot these transformed values manually to inspect
plot(fitted_values, residuals)

#CONTINUE-----------CONTINUE----------------CONTINUE 

#TIME POINT 4
t4.Evs21<-plot_difference(
  gam.t4,
  series = log.con,
  difference = list(plastic_type = c("Elastollan", "TPU 2.1")))+
  labs(x="log(Concentration of Plastic)", y="Abundance (Difference smooth)")+
  ggtitle("Differnce between Elastollan and TPU FC2.1, Day 36")+coord_cartesian(ylim = c(-4.6,2.8))

t4.Evs181<-plot_difference(
  gam.t4,
  series = log.con,
  difference = list(plastic_type = c("Elastollan", "TPU 181")))+
  labs(x="log(Concentration of Plastic)", y="Abundance (Difference smooth)")+
  ggtitle("Differnce between Elastollan and TPU 181, Day 36")+coord_cartesian(ylim = c(-4.6,2.8))

t4.21vs181<-plot_difference(
  gam.t4,
  series = log.con,
  difference = list(plastic_type = c("TPU 2.1", "TPU 181")))+
  labs(x="log(Concentration of Plastic)", y="Abundance (Difference smooth)")+
  ggtitle("Differnce between TPU FC2.1 and TPU 181, Day 36")+coord_cartesian(ylim = c(-4.6,2.8))

t4.Evs21 ### Abundance data used for MANUSCRIPT [Table S13]
t4.Evs181 ### Abundance data used for MANUSCRIPT [Table S13]
t4.21vs181 ### Abundance data used for MANUSCRIPT [Table S13]

#---------------------------
#TIME POINT 5
t5.Evs21<-plot_difference(
  gam.t5,
  series = log.con,
  difference = list(plastic_type = c("Elastollan", "TPU 2.1")))+
  labs(x="log(Concentration of Plastic)", y="Abundance (Difference smooth)")+
  ggtitle("Differnce between Elastollan and TPU FC2.1, Day 56")+coord_cartesian(ylim = c(-4.6,2.8))

t5.Evs181<-plot_difference(
  gam.t5,
  series = log.con,
  difference = list(plastic_type = c("Elastollan", "TPU 181")))+
  labs(x="log(Concentration of Plastic)", y="Abundance (Difference smooth)")+
  ggtitle("Differnce between Elastollan and TPU 181, 56 days exposure")+coord_cartesian(ylim = c(-4.6,2.8))

t5.21vs181<-plot_difference(
  gam.t5,
  series = log.con,
  difference = list(plastic_type = c("TPU 2.1", "TPU 181")))+
  labs(x="log(Concentration of Plastic)", y="Abundance (Difference smooth)")+
  ggtitle("Differnce between TPU FC2.1 and TPU 181, Day 56")+coord_cartesian(ylim = c(-4.6,2.8))

t5.Evs21 ### Abundance data used for MANUSCRIPT [Table S13]
t5.Evs181 ### Abundance data used for MANUSCRIPT [Table S13]
t5.21vs181 ### Abundance data used for MANUSCRIPT [Table S13]


#---------------------------
#TIME POINT 6
t6.Evs21<-plot_difference(
  gam.t6,
  series = log.con,
  difference = list(plastic_type = c("Elastollan", "TPU 2.1")))+
  labs(x="log(Concentration of Plastic)", y="Abundance (Difference smooth)")+
  ggtitle("Differnce between Elastollan and TPU FC2.1, Day 97")+coord_cartesian(ylim = c(-4.6,2.8))

t6.Evs181<-plot_difference(
  gam.t6,
  series = log.con,
  difference = list(plastic_type = c("Elastollan", "TPU 181")))+
  labs(x="log(Concentration of Plastic)", y="Abundance (Difference smooth)")+
  ggtitle("Differnce between Elastollan and TPU 181, Day 97")+coord_cartesian(ylim = c(-4.6,2.8))

t6.21vs181<-plot_difference(
  gam.t6,
  series = log.con,
  difference = list(plastic_type = c("TPU 2.1", "TPU 181")))+
  labs(x="log(Concentration of Plastic)", y="Abundance (Difference smooth)")+
  ggtitle("Differnce between TPU FC2.1 and TPU 181, Day 97")+coord_cartesian(ylim = c(-4.6,2.8))

t6.Evs21 ### Abundance data used for MANUSCRIPT [Table S13]
t6.Evs181 ### Abundance data used for MANUSCRIPT [Table S13]
t6.21vs181 ### Abundance data used for MANUSCRIPT [Table S13]


#STOP-------------STOP--------------STOP#


T1.mod.plot<-
  plot_smooths(
    model = gam.t1,
    series = log.con
    
  )+ 
  geom_point(data=Total_abundance1, 
             aes(x=log.con, y=log(count.per.l), color=plastic_type))+
  labs(x="log(Concentration of Plastic)", y="log(# of Zooplankton)")+
  ggtitle("Concentration of Plastic v.s. Zoop. Abundance, 0 days of exposure")+coord_cartesian(ylim = c(2,8.1))
#Manually set the maximum value for the y-axis


T1.mod.plot


##############Add customization###############
T1.mod.plot <- plot_smooths(
  model = gam.t1,
  series = log.con
) + 
  geom_point(data = Total_abundance1, 
             aes(x = log.con, y = log(count.per.l), color = plastic_type), 
             size = 4) +  # Adjust the size as needed
  scale_color_manual(
    values = c("darkorchid4", "deeppink3", "darkorange"),
    labels = c("Elastollan", "TPU 181", "TPU FC2.1"),
    guide = guide_legend(title = "Plastic type")
  ) +
  labs(x = "log(Plastic concentration [g/L])", y = "log(Zooplankton abundance [individuals/L])") +
  ggtitle("Zooplankton abundance by plastic concentration - 0 days of exposure") +
  coord_cartesian(ylim = c(2, 8.1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
    axis.title.x = element_text(size = 12), # Custom x-axis title
    axis.title.y = element_text(size = 12), # Custom y-axis title
    legend.title = element_text(face = "bold", size = 10), # Custom legend title
    legend.text = element_text(size = 9) # Custom legend text size
  )

#Display the plot
T1.mod.plot

#Edit some of the customization
T1.mod.plot.custom <- plot_smooths(
  model = gam.t1,
  series = log.con
) + 
  geom_point(data = Total_abundance1, 
             aes(x = log.con, y = log(count.per.l), color = plastic_type), 
             size = 3, position = position_dodge(width = 0.2)) +  # Adjust width as needed
  scale_color_manual(
    values = c("darkorchid4", "deeppink3", "darkorange"),
    labels = c("Elastollan", "TPU 181", "TPU FC2.1"),
    guide = guide_legend(title = "Plastic type")
  ) +
  labs(x = "log(Plastic concentration [g/L])", y = "log(Zooplankton abundance [individuals/L])") +
  ggtitle("Zooplankton abundance by plastic concentration - 0 days of exposure") +
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
    axis.title.x = element_text(size = 12), # Custom x-axis title
    axis.title.y = element_text(size = 12), # Custom y-axis title
    legend.title = element_text(face = "bold", size = 10), # Custom legend title
    legend.text = element_text(size = 9) # Custom legend text size
  )

#Display the plot
T1.mod.plot.custom

#Edit the customization
T1.mod.plot.custom <- plot_smooths(
  model = gam.t1,
  series = log.con
) + 
  #Without jitter
  geom_point(data = subset(Total_abundance1, log.con == min(log.con)),
             aes(x = log.con, y = log(count.per.l), color = plastic_type),
             size = 3) +  # No position adjustment here
  
  #With jitter
  geom_point(data = subset(Total_abundance1, log.con != min(log.con)),
             aes(x = log.con, y = log(count.per.l), color = plastic_type),
             size = 3, position = position_dodge(width = 0.2)) +  # Jitter applied to non-minimal concentration
  
  #Customize the color scale
  scale_color_manual(
    values = c("darkorchid4", "deeppink3", "darkorange"),
    labels = c("Elastollan", "TPU 181", "TPU FC2.1"),
    guide = guide_legend(title = "Plastic type")
  ) +
  
  #Axis labels and title
  labs(x = "log(Plastic concentration [g/L])", 
       y = "log(Zooplankton abundance [individuals/L])") +
  ggtitle("Zooplankton abundance by plastic concentration - 0 days of exposure") +
  
  #Set the plot limits
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75)) +
  
  #Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12),                 # Custom x-axis title
    axis.title.y = element_text(size = 12),                 # Custom y-axis title
    legend.title = element_text(face = "bold", size = 10),  # Custom legend title
    legend.text = element_text(size = 9)                    # Custom legend text size
  )

#Display the updated plot
T1.mod.plot.custom



# Same plot but now points with plastic_conc == 0.000 will be black and added to legend
T1.mod.plot.custom <- plot_smooths(
  model = gam.t1,
  series = log.con
) + 
  #Without jitter
  geom_point(data = subset(Total_abundance1, log.con == min(log.con)),
             aes(x = log.con, y = log(count.per.l), color = plastic_type),
             size = 3) +  # No position adjustment here
  
  #Points where plastic_conc == 0.000 (no plastic) colored black
  geom_point(data = subset(Total_abundance1, concentration == 0.000),
             aes(x = log.con, y = log(count.per.l)),
             color = "black",  # Black color for no plastic points
             size = 3) +
  
  #With jitter
  geom_point(data = subset(Total_abundance1, log.con != min(log.con) & concentration != 0.000),
             aes(x = log.con, y = log(count.per.l), color = plastic_type),
             size = 3, position = position_dodge(width = 0.2)) +  # Jitter applied to non-minimal concentration
  
  #Customize the color scale
  scale_color_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic"),  # Add 'No Plastic' label
    guide = guide_legend(title = "Plastic type")
  ) +
  
  #Axis labels and title
  labs(x = "log(Plastic concentration [g/L])", 
       y = "log(Zooplankton abundance [individuals/L])") +
  ggtitle("0 Days of exposure") +
  
  #Set the plot limits
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75)) +
  
  #Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12),                 # Custom x-axis title
    axis.title.y = element_text(size = 12),                 # Custom y-axis title
    legend.title = element_text(face = "bold", size = 10),  # Custom legend title
    legend.text = element_text(size = 9)                    # Custom legend text size
  )

#Display the updated plot
T1.mod.plot.custom




#Making final figure with a few changes.
T1.mod.plot.custom.final <- plot_smooths(
  model = gam.t1,
  series = log.con
) + 
  #Without jitter
  geom_point(data = subset(Total_abundance1, log.con == min(log.con)),
             aes(x = log.con, y = log(count.per.l), color = plastic_type),
             size = 3) +  # No position adjustment here
  
  #Points where plastic_conc == 0.000 (no plastic) colored black
  geom_point(data = subset(Total_abundance1, concentration == 0.000),
             aes(x = log.con, y = log(count.per.l)),
             color = "black",  # Black color for no plastic points
             size = 3) +
  
  #With jitter
  geom_point(data = subset(Total_abundance1, log.con != min(log.con) & concentration != 0.000),
             aes(x = log.con, y = log(count.per.l), color = plastic_type),
             size = 3, position = position_dodge(width = 0.2)) +  # Jitter applied to non-minimal concentration
  
  #Customize the color scale
  scale_color_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic"),  # Add 'No Plastic' label
    guide = guide_legend(title = "Plastic type")
  ) +
  
  #Axis labels and title
  labs(x = "log(Plastic concentration [g/L])", 
       y = "log(Zooplankton abundance [individuals/L])") +
  ggtitle("0 Days of exposure") +
  
  #Set the plot limits
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75)) +
  
  #Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12),                 # Custom x-axis title
    axis.title.y = element_text(size = 14),                 # Increase size of y-axis title
    axis.text.y = element_text(size = 14),                  # Increase size of y-axis tick labels
    axis.text.x = element_text(size = 14),                  # Optional: you can also adjust the x-axis tick labels
    legend.position = "none",                               # Remove legend
    panel.grid = element_blank(),                           # Remove gridlines
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Add black box around figure
    axis.line = element_line(colour = "black")              # Ensure axes are visible
  )

#Display the updated plot
T1.mod.plot.custom.final ### Abundance data used for MANUSCRIPT [Figure S3]









#TIME POINT 4
T4.mod.plot<-
  plot_smooths(
    model =gam.t4,
    series = log.con,
    comparison = plastic_type
    
  )+ 
  geom_point(data=Total_abundance4, 
             aes(x=log.con, y=log(count.per.l), color=plastic_type))+
  labs(x="log(Concentration of Plastic)", y="log(# of Zooplankton)")+
  ggtitle("Concentration of Plastic v.s. Zoop. Abundance, 36 days of exposure")+coord_cartesian(ylim = c(2,8.1))

T4.mod.plot

#Change the colors and point size
# TIME POINT 4
T4.mod.plot <- plot_smooths(
  model = gam.t4,
  series = log.con,
  comparison = plastic_type
) + 
  geom_point(data = Total_abundance4, 
             aes(x = log.con, y = log(count.per.l), color = plastic_type), size = 3) +
  scale_color_manual(values = c("darkorchid4", "deeppink3", "darkorange")) +
  scale_fill_manual(values = c("darkorchid4", "deeppink3", "darkorange")) +
  labs(x = "log(Concentration of Plastic)", y = "log(# of Zooplankton)") +
  ggtitle("Concentration of Plastic vs. Zoop. Abundance, 36 days of exposure") +
  coord_cartesian(ylim = c(2, 8.1))

T4.mod.plot

# TIME POINT 4
T4.mod.plot <- plot_smooths(
  model = gam.t4,
  series = log.con,
  comparison = plastic_type
) + 
  geom_point(data = Total_abundance4, 
             aes(x = log.con, y = log(count.per.l), color = plastic_type), size = 3) +
  scale_color_manual(values = c("darkorchid4", "deeppink3", "darkorange")) +
  scale_fill_manual(values = c("darkorchid4", "deeppink3", "darkorange")) +
  labs(x = "log(Concentration of plastic)", y = "log(# of zooplankton)") +
  ggtitle("Concentration of Plastic vs. Zoop. Abundance, 36 days of exposure") +
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75))  # Adjust limits as necessary

T4.mod.plot

#Update the second plot to include black points for "No Plastic"
T4.mod.plot <- plot_smooths(
  model = gam.t4,
  series = log.con,
  comparison = plastic_type
) + 
  #Points for other plastic concentrations
  geom_point(data = subset(Total_abundance4, concentration != 0.000),
             aes(x = log.con, y = log(count.per.l), color = plastic_type),
             size = 3) +
  
  #Points where plastic_conc == 0.000 (no plastic) colored black
  geom_point(data = subset(Total_abundance4, concentration == 0.000),
             aes(x = log.con, y = log(count.per.l)),
             color = "black",  # Black color for no plastic points
             size = 3) +
  
  #Customize the color scale to include black for "No Plastic"
  scale_color_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic"),  # Add 'No Plastic' label
    guide = guide_legend(title = "Plastic type")
  ) +
  
  #Customize the fill scale for consistency
  scale_fill_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic")   # Add 'No Plastic' label
  ) +
  
  #Axis labels and title
  labs(x = "log(Concentration of plastic)", y = "log(# of zooplankton)") +
  ggtitle("Concentration of Plastic vs. Zoop. Abundance, 56 days of exposure") +
  
  #Set the plot limits
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75)) +
  
  #Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12),                 # Custom x-axis title
    axis.title.y = element_text(size = 12),                 # Custom y-axis title
    legend.title = element_text(face = "bold", size = 10),  # Custom legend title
    legend.text = element_text(size = 9)                    # Custom legend text size
  )

#Display the updated plot
T4.mod.plot

#Now make adjustments for final plot
T4.mod.plot.custom.final <- plot_smooths(
  model = gam.t4,
  series = log.con,
  comparison = plastic_type
) + 
  #Points for other plastic concentrations
  geom_point(data = subset(Total_abundance4, concentration != 0.000),
             aes(x = log.con, y = log(count.per.l), color = plastic_type),
             size = 3) +
  
  #Points where plastic_conc == 0.000 (no plastic) colored black
  geom_point(data = subset(Total_abundance4, concentration == 0.000),
             aes(x = log.con, y = log(count.per.l)),
             color = "black",  # Black color for no plastic points
             size = 3) +
  
  #Customize the color scale to include black for "No Plastic"
  scale_color_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic"),  # Add 'No Plastic' label
    guide = guide_legend(title = "Plastic type")
  ) +
  
  #Customize the fill scale for consistency
  scale_fill_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic")   # Add 'No Plastic' label
  ) +
  
  #Axis labels and title
  labs(x = "log(Concentration of plastic)", y = NULL) +  # Remove y-axis label
  ggtitle("Day 36") +
  
  #Set the plot limits
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75)) +
  
  #Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12),                 # Custom x-axis title
    axis.title.y = element_blank(),                         # Remove y-axis label
    axis.text.y = element_blank(),                          # Remove y-axis tick marks
    axis.ticks.y = element_blank(),                         # Remove y-axis tick lines
    axis.text.x = element_text(size = 14),                  # Optional: adjust x-axis tick label size
    legend.position = "none",                               # Remove legend
    panel.grid = element_blank(),                           # Remove gridlines
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Add black box around figure
    axis.line = element_line(colour = "black")              # Ensure axes are visible
  )

#Display the updated plot
T4.mod.plot.custom.final ### Abundance data used for MANUSCRIPT [Figure S3]









#TIME POINT 5
T5.mod.plot<-
  plot_smooths(
    model =gam.t5,
    series = log.con,
    comparison = plastic_type
    
  )+ 
  geom_point(data=Total_abundance5, 
             aes(x=log.con, y=log(count.per.l), color=plastic_type))+
  labs(x="log(Concentration of Plastic)", y="log(# of Zooplankton)")+
  ggtitle("Concentration of Plastic v.s. Zoop. Abundance, 56 days of exposure")+coord_cartesian(ylim = c(2,8.1))

T5.mod.plot


# TIME POINT 5
T5.mod.plot <- plot_smooths(
  model = gam.t5,
  series = log.con,
  comparison = plastic_type
) + 
  geom_point(data = Total_abundance5, 
             aes(x = log.con, y = log(count.per.l), color = plastic_type), size = 3) +
  scale_color_manual(values = c("darkorchid4", "deeppink3", "darkorange")) +
  scale_fill_manual(values = c("darkorchid4", "deeppink3", "darkorange")) +
  labs(x = "log(Concentration of plastic)", y = "log(# of zooplankton)") +
  ggtitle("Concentration of Plastic vs. Zoop. Abundance, 56 days of exposure") +
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75))

T5.mod.plot


#Update the plot for time point 5 to include black points for "No Plastic"
T5.mod.plot <- plot_smooths(
  model = gam.t5,
  series = log.con,
  comparison = plastic_type
) + 
  #Points for other plastic concentrations
  geom_point(data = subset(Total_abundance5, concentration != 0.000),
             aes(x = log.con, y = log(count.per.l), color = plastic_type),
             size = 3) +
  
  #Points where plastic_conc == 0.000 (no plastic) colored black
  geom_point(data = subset(Total_abundance5, concentration == 0.000),
             aes(x = log.con, y = log(count.per.l)),
             color = "black",  # Black color for no plastic points
             size = 3) +
  
  #Customize the color scale to include black for "No Plastic"
  scale_color_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic"),  # Add 'No Plastic' label
    guide = guide_legend(title = "Plastic type")
  ) +
  
  #Customize the fill scale for consistency
  scale_fill_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic")   # Add 'No Plastic' label
  ) +
  
  #Axis labels and title
  labs(x = "log(Concentration of plastic)", y = "log(# of zooplankton)") +
  ggtitle("Concentration of Plastic vs. Zoop. Abundance, 56 days of exposure") +
  
  #Set the plot limits
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75)) +
  
  #Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12),                 # Custom x-axis title
    axis.title.y = element_text(size = 12),                 # Custom y-axis title
    legend.title = element_text(face = "bold", size = 10),  # Custom legend title
    legend.text = element_text(size = 9)                    # Custom legend text size
  )

#Display the updated plot
T5.mod.plot


#Now make the final figure
T5.mod.plot.custom.final <- plot_smooths(
  model = gam.t5,
  series = log.con,
  comparison = plastic_type
) + 
  #Points for other plastic concentrations
  geom_point(data = subset(Total_abundance5, concentration != 0.000),
             aes(x = log.con, y = log(count.per.l), color = plastic_type),
             size = 3) +
  
  #Points where plastic_conc == 0.000 (no plastic) colored black
  geom_point(data = subset(Total_abundance5, concentration == 0.000),
             aes(x = log.con, y = log(count.per.l)),
             color = "black",  # Black color for no plastic points
             size = 3) +
  
  #Customize the color scale to include black for "No Plastic"
  scale_color_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic"),  # Add 'No Plastic' label
    guide = guide_legend(title = "Plastic type")
  ) +
  
  #Customize the fill scale for consistency
  scale_fill_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic")   # Add 'No Plastic' label
  ) +
  
  #Axis labels and title
  labs(x = "log(Concentration of plastic)", y = NULL) +  # Remove y-axis label
  ggtitle("Day 56") +
  
  #Set the plot limits
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75)) +
  
  #Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12),                 # Custom x-axis title
    axis.title.y = element_blank(),                         # Remove y-axis label
    axis.text.y = element_blank(),                          # Remove y-axis tick marks
    axis.ticks.y = element_blank(),                         # Remove y-axis tick lines
    axis.text.x = element_text(size = 14),                  # Optional: adjust x-axis tick label size
    legend.position = "none",                               # Remove legend
    panel.grid = element_blank(),                           # Remove gridlines
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Add black box around figure
    axis.line = element_line(colour = "black")              # Ensure axes are visible
  )

#Display the updated plot
T5.mod.plot.custom.final ### Abundance data used for MANUSCRIPT [Figure S3]









#TIME POINT 6
T6.mod.plot<-
  plot_smooths(
    model = gam.t6,
    series = log.con,
    comparison = plastic_type,
    ci_z=2.576
    
  )+ 
  geom_point(data=Total_abundance6, 
             aes(x=log.con, y=log(count.per.l), color=plastic_type))+
  labs(x="log(Concentration of Plastic)", y="log(# of Zooplankton)")+
  ggtitle("Concentration of Plastic v.s. Zoop. Abundance, 97 days of exposure")+coord_cartesian(ylim = c(2,8.1))

T6.mod.plot
summary(test.run.t6)

# TIME POINT 6
T6.mod.plot <- plot_smooths(
  model = gam.t6,
  series = log.con,
  comparison = plastic_type,
  ci_z = 2.576
) + 
  geom_point(data = Total_abundance6, 
             aes(x = log.con, y = log(count.per.l), color = plastic_type), size = 3) +
  scale_color_manual(values = c("darkorchid4", "deeppink3", "darkorange")) +
  scale_fill_manual(values = c("darkorchid4", "deeppink3", "darkorange")) +
  labs(x = "log(Concentration of plastic)", y = "log(# of zooplankton)") +
  ggtitle("Concentration of Plastic vs. Zoop. Abundance, 97 days of exposure") +
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75))

T6.mod.plot


#Update the plot for Time Point 6 to include black points for "No Plastic"
T6.mod.plot <- plot_smooths(
  model = gam.t6,
  series = log.con,
  comparison = plastic_type,
  ci_z = 2.576  #Keep this confidence interval setting
) + 
  #Points for other plastic concentrations
  geom_point(data = subset(Total_abundance6, concentration != 0.000),
             aes(x = log.con, y = log(count.per.l), color = plastic_type),
             size = 3) +
  
  #Points where plastic_conc == 0.000 (no plastic) colored black
  geom_point(data = subset(Total_abundance6, concentration == 0.000),
             aes(x = log.con, y = log(count.per.l)),
             color = "black",  # Black color for no plastic points
             size = 3) +
  
  #Customize the color scale to include black for "No Plastic"
  scale_color_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic"),  # Add 'No Plastic' label
    guide = guide_legend(title = "Plastic type")
  ) +
  
  #Customize the fill scale for consistency
  scale_fill_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic")   # Add 'No Plastic' label
  ) +
  
  #Axis labels and title
  labs(x = "log(Concentration of plastic)", y = "log(# of zooplankton)") +
  ggtitle("Concentration of Plastic vs. Zoop. Abundance, 97 days of exposure") +
  
  #Set the plot limits
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75)) +
  
  #Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12),                 # Custom x-axis title
    axis.title.y = element_text(size = 12),                 # Custom y-axis title
    legend.title = element_text(face = "bold", size = 10),  # Custom legend title
    legend.text = element_text(size = 9)                    # Custom legend text size
  )

#Display the updated plot
T6.mod.plot



#Now make final plot
T6.mod.plot.custom.final <- plot_smooths(
  model = gam.t6,
  series = log.con,
  comparison = plastic_type,
  ci_z = 2.576  # Keep this confidence interval setting
) + 
  #Points for other plastic concentrations
  geom_point(data = subset(Total_abundance6, concentration != 0.000),
             aes(x = log.con, y = log(count.per.l), color = plastic_type),
             size = 3) +
  
  #Points where plastic_conc == 0.000 (no plastic) colored black
  geom_point(data = subset(Total_abundance6, concentration == 0.000),
             aes(x = log.con, y = log(count.per.l)),
             color = "black",  # Black color for no plastic points
             size = 3) +
  
  #Customize the color scale to include black for "No Plastic"
  scale_color_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic"),  # Add 'No Plastic' label
    guide = guide_legend(title = "Plastic type")
  ) +
  
  #Customize the fill scale for consistency
  scale_fill_manual(
    values = c("darkorchid4", "deeppink3", "darkorange", "black"),  # Add black for 'No Plastic'
    labels = c("Elastollan", "TPU 181", "TPU FC2.1", "No Plastic")   # Add 'No Plastic' label
  ) +
  
  #Axis labels and title
  labs(x = "log(Concentration of plastic [g/L])", y = NULL) +  # Remove y-axis label
  ggtitle("Day 97") +
  
  #Set the plot limits
  coord_cartesian(ylim = c(-0.5, 7), xlim = c(-6.25, -0.75)) +
  
  #Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12),                 # Custom x-axis title
    axis.title.y = element_blank(),                         # Remove y-axis label
    axis.text.y = element_blank(),                          # Remove y-axis tick marks
    axis.ticks.y = element_blank(),                         # Remove y-axis tick lines
    axis.text.x = element_text(size = 14),                  # Optional: adjust x-axis tick label size
    legend.position = "none",                               # Remove legend
    panel.grid = element_blank(),                           # Remove gridlines
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Add black box around figure
    axis.line = element_line(colour = "black")              # Ensure axes are visible
  )

#Display the updated plot
T6.mod.plot.custom.final ### Abundance data used for MANUSCRIPT [Figure S3]



#Combine the plots
library(patchwork)


#Add overall title and remove and center x axis label
#Remove x-axis labels from the individual plots, keep the tick marks, and bold titles and axis labels
T1.mod.plot.custom.final <- T1.mod.plot.custom.final + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),  # Bold and add space to y-axis title
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Centered and bold for individual figure titles
        plot.margin = margin(0, 0, 0, 0))

T4.mod.plot.custom.final <- T4.mod.plot.custom.final + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),  # Bold and add space to y-axis title
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Centered and bold for individual figure titles
        plot.margin = margin(0, 0, 0, 0))

T5.mod.plot.custom.final <- T5.mod.plot.custom.final + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),  # Bold and add space to y-axis title
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Centered and bold for individual figure titles
        plot.margin = margin(0, 0, 0, 0))

T6.mod.plot.custom.final <- T6.mod.plot.custom.final + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),  # Bold and add space to y-axis title
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Centered and bold for individual figure titles
        plot.margin = margin(0, 0, 0, 0))

#Combine the plots horizontally
combined_plot <- (T1.mod.plot.custom.final + 
                    T4.mod.plot.custom.final + 
                    T5.mod.plot.custom.final + 
                    T6.mod.plot.custom.final) +
  plot_layout(ncol = 4)

#Create an empty plot to serve as a placeholder for the single x-axis label
x_label_plot <- ggplot() + 
  theme_void() +   #Remove all elements from this empty plot
  labs(x = "log(Concentration of plastic [g/L])") +  # Add the single x-axis label here
  theme(axis.title.x = element_text(size = 14, face = "bold", hjust = 0.5))  # Bold and center the x-axis label

#Combine the plots with the title and x-axis label
final_plot <- (combined_plot / x_label_plot) + 
  plot_layout(heights = c(10, 0.5)) +  #Further reduce the space by adjusting heights
  
  #This part controls ONLY the overall title (Zooplankton abundance)
  plot_annotation(title = "Zooplankton abundance") &  #Add overall title
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  #Center and keep the overall title bold

#Display the final combined plot
print(final_plot)





#BIOMASS CODE
#Library any additional packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidymv)
library(patchwork)
library(purrr)
library(mgcv)

#Read in length data
measurements0<-read.csv("MP Zoop Measurements - Zoop_Lengths.csv")


measurements0 <- subset(measurements0, select = -c(who.took.measurment, hours, who.counted, number.measured, count.date, measured, sample.period))
measurements0$days.of.exposure<-measurements0$sample.date
measurements0 <- measurements0 %>%
  mutate(days.of.exposure = case_when(
    sample.date == as.Date("2023-03-07") ~ "0",
    sample.date == as.Date("2023-03-13") ~ "6",
    sample.date == as.Date("2023-03-23") ~ "16",
    sample.date == as.Date("2023-04-12") ~ "36",
    sample.date == as.Date("2023-05-02") ~ "56",
    sample.date == as.Date("2023-06-12") ~ "97",
    TRUE ~ NA_character_  #Assign NA to any date that doesn't match
  ))

measurements0 <- measurements0 %>%
  mutate(tank.id = str_replace(tank.id, "A0 ", "A0"))

unique(measurements0$sample.date)
measurements0 <- measurements0[!is.na(measurements0$prop.counted)]


measurements0<- subset(measurements0, !is.na(length..mm.))
unique(measurements0$tank.id)
measurements0 <-subset(measurements0, !grepl("2023-04-12", sample.date))

#################################

#GOING TO LOOK AT DISTRIBUTION OF MEASURMENT OF TAXA FOR DATES
unique(measurements0$taxa.code)
measurements0$length..mm.<-as.numeric(measurements0$length..mm.)


#CHY
measurements.chy<-measurements0 %>% filter(days.of.exposure %in%c("97","0"), taxa.code=="chy")
measurements.chy <- na.omit(measurements.chy)
hist(measurements.chy$length..mm.)

measurements97.chy<-measurements.chy %>% filter(days.of.exposure %in% c("97"))
measurements97.chy$length..mm.<-as.numeric(measurements97.chy$length..mm.)
hist(measurements97.chy$length..mm.)

measurements0.chy<-measurements.chy %>% filter(days.of.exposure %in% c("0"))
measurements0.chy$length..mm.<-as.numeric(measurements0.chy$length..mm.)
hist(measurements0.chy$length..mm.)

#DAPHNIA
measurements.dap<-measurements0 %>% filter(days.of.exposure %in% c("97","0"), taxa.code=="dap")
measurements.dap <-subset(measurements.dap, !grepl("C6", tank.id))
measurements.dap <- na.omit(measurements.dap)
hist(measurements.dap$length..mm.)

measurements97.dap<-measurements.dap %>% filter(days.of.exposure %in% c("97"))
measurements97.dap$length..mm.<-as.numeric(measurements97.dap$length..mm.)
hist(measurements97.dap$length..mm.)

measurements0.dap<-measurements.dap %>% filter(days.of.exposure %in% c("0"))
measurements0.dap$length..mm.<-as.numeric(measurements0.dap$length..mm.)
hist(measurements0.dap$length..mm.)

#HYDRA
measurements.hydra<-measurements0 %>% filter(days.of.exposure %in% c("97","0"), taxa.code=="hydra")
measurements.hydra <- na.omit(measurements.hydra)
measurements97.hydra<-measurements.hydra %>% filter(days.of.exposure %in% c("97"))
measurements97.hydra$length..mm.<-as.numeric(measurements97.hydra$length..mm.)
hist(measurements97.hydra$length..mm.)

measurements0.hydra<-measurements.hydra %>% filter(days.of.exposure %in% c("0"))
measurements0.hydra$length..mm.<-as.numeric(measurements0.hydra$length..mm.)
hist(measurements0.hydra$length..mm.)

#OST
measurements.ost<-measurements0 %>% filter(days.of.exposure %in% c("97","0"), taxa.code=="ost")
measurements.ost <- na.omit(measurements.ost)
hist(measurements.ost$length..mm.)

measurements97.ost<-measurements.ost %>% filter(days.of.exposure %in% c("97"))
measurements97.ost$length..mm.<-as.numeric(measurements97.ost$length..mm.)
hist(measurements97.ost$length..mm.)

measurements0.ost<-measurements.ost %>% filter(days.of.exposure %in% c("0"))
measurements0.ost$length..mm.<-as.numeric(measurements0.ost$length..mm.)
hist(measurements0.ost$length..mm.)

#MON
measurements.mon<-measurements0 %>% filter(days.of.exposure %in% c("97","0"), taxa.code=="mon")
measurements.mon <- na.omit(measurements.mon)
hist(measurements.mon$length..mm.)

measurements97.mon<-measurements.mon %>% filter(days.of.exposure %in% c("97"))
measurements97.mon$length..mm.<-as.numeric(measurements97.mon$length..mm.)
hist(measurements97.mon$length..mm.)

measurements0.mon<-measurements.mon %>% filter(days.of.exposure %in% c("0"))
measurements0.mon$length..mm.<-as.numeric(measurements0.mon$length..mm.)
hist(measurements0.mon$length..mm.)#no monostyla at time 0

#NAU
measurements.nau<-measurements0 %>% filter(days.of.exposure %in% c("97","0"), taxa.code=="nau")
measurements.nau <- na.omit(measurements.nau)
hist(measurements.nau$length..mm.)

measurements97.nau<-measurements.nau %>% filter(days.of.exposure %in% c("97"))
measurements97.nau$length..mm.<-as.numeric(measurements97.nau$length..mm.)
hist(measurements97.nau$length..mm.)

measurements0.nau<-measurements.nau %>% filter(days.of.exposure %in% c("0"))
measurements0.nau$length..mm.<-as.numeric(measurements0.nau$length..mm.)
hist(measurements0.nau$length..mm.)

#ERG
measurements.erg<-measurements0 %>% filter(days.of.exposure %in% c("97","0"), taxa.code=="erg")
measurements.erg <- na.omit(measurements.erg)
hist(measurements.erg$length..mm.)


measurements97.erg<-measurements.erg %>% filter(days.of.exposure %in% c("97"))
measurements97.erg$length..mm.<-as.numeric(measurements97.erg$length..mm.)
hist(measurements97.erg$length..mm.)

measurements0.erg<-measurements.erg %>% filter(days.of.exposure %in% c("0"))
measurements0.erg$length..mm.<-as.numeric(measurements0.erg$length..mm.)
hist(measurements0.erg$length..mm.)

#CYC
measurements.cyc<-measurements0 %>% filter(days.of.exposure %in% c("97","0"), taxa.code=="cyc")
measurements.cyc <- na.omit(measurements.cyc)
hist(measurements.cyc$length..mm.)

measurements97.cyc<-measurements.cyc %>% filter(days.of.exposure %in% c("97"))
measurements97.cyc$length..mm.<-as.numeric(measurements97.cyc$length..mm.)
hist(measurements97.cyc$length..mm.)

measurements0.cyc<-measurements.cyc %>% filter(days.of.exposure %in% c("0"))
measurements0.cyc$length..mm.<-as.numeric(measurements0.cyc$length..mm.)
hist(measurements0.cyc$length..mm.)

#CAL
measurements.cal<-measurements0 %>% filter(days.of.exposure %in% c("97","0"), taxa.code=="cal")
measurements.cal <- na.omit(measurements.cal)
hist(measurements.cal$length..mm.)

measurements97.cal<-measurements.cal %>% filter(days.of.exposure %in% c("97"))
measurements97.cal$length..mm.<-as.numeric(measurements97.cal$length..mm.)
hist(measurements97.cal$length..mm.)

measurements0.cal<-measurements.cal %>% filter(days.of.exposure %in% c("0"))
measurements0.cal$length..mm.<-as.numeric(measurements0.cal$length..mm.)
hist(measurements0.cal$length..mm.)

#large copapods

measurements.cop<-measurements0 %>% filter(days.of.exposure %in% c("97","0"), taxa.code %in% c("cal", "erg", "cyc"))
measurements.cop <- na.omit(measurements.cal)
hist(measurements.cop$length..mm.)



#Remove true and verified outlier measurement in DAP
measurements.dap <- subset(measurements.dap, !grepl("A6", tank.id) & days.of.exposure=="97" )


#GOING TO CALCULATE THE MEDIAN FROM ACROSS ALL TANKS ACROSS ALL TIMES 
#GOING TO USE ^ TO REDUCE VARIATION IN SYSTEMATIC ERROR
median.chy<-median(measurements.chy$length..mm.)
median.dap<- median(measurements.dap$length..mm.)
median.hydra<-median(measurements.hydra$length..mm.)
median.ost<- median(measurements.ost$length..mm.)
median.mon<-median(measurements.mon$length..mm.)
median.nau<-median(measurements.nau$length..mm.)
median.erg<-median(measurements.erg$length..mm.)
median.cyc<-median(measurements.cyc$length..mm.)
median.cal<-median(measurements.cal$length..mm.)
median.cop<- median(measurements.cop$length..mm.)

#use the mean
mean.chy<-mean(measurements.chy$length..mm.)
mean.dap<- mean(measurements.dap$length..mm.)
mean.hydra<-mean(measurements.hydra$length..mm.)
mean.ost<- mean(measurements.ost$length..mm.)
mean.mon<-mean(measurements.mon$length..mm.)
mean.nau<-mean(measurements.nau$length..mm.)
mean.erg<-mean(measurements.erg$length..mm.)
mean.cyc<-mean(measurements.cyc$length..mm.)
mean.cal<-mean(measurements.cal$length..mm.)
mean.cop<- mean(measurements.cop$length..mm.)


#going to apply median measurements to all taxa over all time periods
all.time.measurements<-filtered_data
#removing non-primary taxa (insects, hydra, unknown)
desired_values <- c("cladoceran","daphnia","nauplii",
                    "large copepod", "rotifer", "ostracoda")

all.time.measurements <- all.time.measurements %>%
  filter( Taxa%in% desired_values)
all.time.measurements<- all.time.measurements[, -which(names(all.time.measurements) == "plastic.type")]


Total_measurements0 <- subset(all.time.measurements, !grepl("0", tank.id))
Control_measurements1 <- subset(all.time.measurements, grepl("0", tank.id) & sample.period == 0)
Control_measurements4 <- subset(all.time.measurements, grepl("0", tank.id) & sample.period == 36)
Control_measurements5 <- subset(all.time.measurements, grepl("0", tank.id) & sample.period == 56)
Control_measurements6 <- subset(all.time.measurements, grepl("0", tank.id) & sample.period == 97)



Total_measurements0 <- subset(Total_measurements0, Count != 0)
Control_measurements1 <- subset(Control_measurements1, Count != 0)
Control_measurements4 <- subset(Control_measurements4, Count != 0)
Control_measurements5<- subset(Control_measurements5, Count != 0)
Control_measurements6<- subset(Control_measurements6, Count != 0) 

#APPLYING MEAN MEASURMENTS TO ALL TAXA AND MERGE so we work with measurments0 and controls from now on
unique(Total_measurements0$Taxa)

Total_measurements0$mean.length<-Total_measurements0$Taxa
Total_measurements0<- Total_measurements0 %>% mutate(
  mean.length = case_when(
    mean.length == "cladoceran" ~ mean.chy,
    mean.length == "daphnia" ~ mean.dap,
    mean.length == "ostracoda" ~ mean.ost,
    mean.length == "rotifer" ~ mean.mon,
    mean.length == "nauplii" ~ mean.nau,
    mean.length == "large copepod" ~ mean.cop
  ))

Control_measurements1$mean.length <-Control_measurements1$Taxa
Control_measurements1<- Control_measurements1 %>% mutate(
  mean.length = case_when(
    mean.length == "cladoceran" ~ mean.chy,
    mean.length == "daphnia" ~ mean.dap,
    mean.length == "ostracoda" ~ mean.ost,
    mean.length == "rotifer" ~ mean.mon,
    mean.length == "nauplii" ~ mean.nau,
    mean.length == "large copepod" ~ mean.cop
  ))
Control_measurements4$mean.length<-Control_measurements4$Taxa
Control_measurements4<- Control_measurements4 %>% mutate(
  mean.length = case_when(
    mean.length == "cladoceran" ~ mean.chy,
    mean.length == "daphnia" ~ mean.dap,
    mean.length == "ostracoda" ~ mean.ost,
    mean.length == "rotifer" ~ mean.mon,
    mean.length == "nauplii" ~ mean.nau,
    mean.length == "large copepod" ~ mean.cop
  ))
Control_measurements5$mean.length<-Control_measurements5$Taxa
Control_measurements5<- Control_measurements5 %>% mutate(
  mean.length = case_when(
    mean.length == "cladoceran" ~ mean.chy,
    mean.length == "daphnia" ~ mean.dap,
    mean.length == "ostracoda" ~ mean.ost,
    mean.length == "rotifer" ~ mean.mon,
    mean.length == "nauplii" ~ mean.nau,
    mean.length == "large copepod" ~ mean.cop
  ))
Control_measurements6$mean.length<-Control_measurements6$Taxa
Control_measurements6<- Control_measurements6 %>% mutate(
  mean.length = case_when(
    mean.length == "cladoceran" ~ mean.chy,
    mean.length == "daphnia" ~ mean.dap,
    mean.length == "ostracoda" ~ mean.ost,
    mean.length == "rotifer" ~ mean.mon,
    mean.length == "nauplii" ~ mean.nau,
    mean.length == "large copepod" ~ mean.cop
  ))



#MEANS HAVE BEEN APPLIED, NOW UPDATE PLASTIC TYPE AND CONCENTRATION

#subsample for sample period
Total_measurements1 <-subset(Total_measurements0, grepl("0", sample.period))
Total_measurements4 <-subset(Total_measurements0, grepl("36", sample.period))
Total_measurements5 <-subset(Total_measurements0, grepl("56", sample.period))
Total_measurements6 <-subset(Total_measurements0, grepl("97", sample.period))


#master mutate code below#

Control_measurements1$plastic_type<-Control_measurements1$tank.id
Control_measurements1$concentration <-Control_measurements1$tank.id

Control_measurements1<- Control_measurements1 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Control_measurements1<- Control_measurements1 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Control_measurements1$log.con<-Control_measurements1$tank.id

Control_measurements1<- Control_measurements1 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))

#LOG trans above
Control_measurements1$log.con<-as.numeric(Control_measurements1$log.con)
Control_measurements1$log.con<-log(Control_measurements1$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Control_measurements1$log.con)

Control_measurements4$plastic_type<-Control_measurements4$tank.id
Control_measurements4$concentration <-Control_measurements4$tank.id

Control_measurements4<- Control_measurements4 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Control_measurements4<- Control_measurements4 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Control_measurements4$log.con<-Control_measurements4$tank.id

Control_measurements4<- Control_measurements4 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))

#LOG trans above
Control_measurements4$log.con<-as.numeric(Control_measurements4$log.con)
Control_measurements4$log.con<-log(Control_measurements4$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Control_measurements4$log.con)

Control_measurements5$plastic_type<-Control_measurements5$tank.id
Control_measurements5$concentration <-Control_measurements5$tank.id

Control_measurements5<- Control_measurements5 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Control_measurements5<- Control_measurements5 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Control_measurements5$log.con<-Control_measurements5$tank.id

Control_measurements5<- Control_measurements5 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))

#LOG trans above
Control_measurements5$log.con<-as.numeric(Control_measurements5$log.con)
Control_measurements5$log.con<-log(Control_measurements5$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Control_measurements5$log.con)

Control_measurements6$plastic_type<-Control_measurements6$tank.id
Control_measurements6$concentration <-Control_measurements6$tank.id

Control_measurements6<- Control_measurements6 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Control_measurements6<- Control_measurements6 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Control_measurements6$log.con<-Control_measurements6$tank.id

Control_measurements6<- Control_measurements6 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))

#LOG trans above
Control_measurements6$log.con<-as.numeric(Control_measurements6$log.con)
Control_measurements6$log.con<-log(Control_measurements6$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Total_measurements6$log.con)

Total_measurements1$plastic_type<-Total_measurements1$tank.id
Total_measurements1$concentration <-Total_measurements1$tank.id

Total_measurements1<- Total_measurements1 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Total_measurements1<- Total_measurements1 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Total_measurements1$log.con<-Total_measurements1$tank.id

Total_measurements1<- Total_measurements1 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))

#LOG trans above
Total_measurements1$log.con<-as.numeric(Total_measurements1$log.con)
Total_measurements1$log.con<-log(Total_measurements1$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Total_measurements1$log.con)

Total_measurements4$plastic_type<-Total_measurements4$tank.id
Total_measurements4$concentration <-Total_measurements4$tank.id

Total_measurements4<- Total_measurements4 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Total_measurements4<- Total_measurements4 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Total_measurements4$log.con<-Total_measurements4$tank.id

Total_measurements4<- Total_measurements4 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))

#LOG tran above
Total_measurements4$log.con<-as.numeric(Total_measurements4$log.con)
Total_measurements4$log.con<-log(Total_measurements4$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Total_measurements4$log.con)

Total_measurements5$plastic_type<-Total_measurements5$tank.id
Total_measurements5$concentration <-Total_measurements5$tank.id

Total_measurements5<- Total_measurements5 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Total_measurements5<- Total_measurements5 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Total_measurements5$log.con<-Total_measurements5$tank.id

Total_measurements5<- Total_measurements5 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))

#LOG trans above
Total_measurements5$log.con<-as.numeric(Total_measurements5$log.con)
Total_measurements5$log.con<-log(Total_measurements5$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Total_measurements5$log.con)



Total_measurements6$plastic_type<-Total_measurements6$tank.id
Total_measurements6$concentration <-Total_measurements6$tank.id

Total_measurements6<- Total_measurements6 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Total_measurements6<- Total_measurements6 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))



#LOG TRANS CONCENTRATION
Total_measurements6$log.con<-Total_measurements6$tank.id

Total_measurements6<- Total_measurements6 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))

#LOG trans above
Total_measurements6$log.con<-as.numeric(Total_measurements6$log.con)
Total_measurements6$log.con<-log(Total_measurements6$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10





unique(Control_measurements1$log.con)
unique(Control_measurements4$log.con)
unique(Control_measurements5$log.con)
unique(Control_measurements6$log.con)
unique(Total_measurements1$log.con)
unique(Total_measurements4$log.con)
unique(Total_measurements5$log.con)
unique(Total_measurements6$log.con)

Control1.E<-Control_measurements1
Control1.E <- Control1.E %>% 
  mutate(plastic_type = "Elastollan")

Control1.181<-Control_measurements1
Control1.181 <- Control1.181 %>% 
  mutate(plastic_type = "TPU 181")

Control1.2.1<-Control_measurements1
Control1.2.1 <- Control1.2.1 %>% 
  mutate(plastic_type = "TPU 2.1")

Control4.E<-Control_measurements4
Control4.E <- Control4.E %>% 
  mutate(plastic_type = "Elastollan")

Control4.181<-Control_measurements4
Control4.181 <- Control4.181 %>% 
  mutate(plastic_type = "TPU 181")

Control4.2.1<-Control_measurements4
Control4.2.1 <- Control4.2.1 %>% 
  mutate(plastic_type = "TPU 2.1")

Control5.E<-Control_measurements5
Control5.E <- Control5.E %>% 
  mutate(plastic_type = "Elastollan")

Control5.181<-Control_measurements5
Control5.181 <- Control5.181 %>% 
  mutate(plastic_type = "TPU 181")

Control5.2.1<-Control_measurements5
Control5.2.1 <- Control5.2.1 %>% 
  mutate(plastic_type = "TPU 2.1")

Control6.E<-Control_measurements6
Control6.E <- Control6.E %>% 
  mutate(plastic_type = "Elastollan")

Control6.181<-Control_measurements6
Control6.181 <- Control6.181 %>% 
  mutate(plastic_type = "TPU 181")

Control6.2.1<-Control_measurements6
Control6.2.1 <- Control6.2.1 %>% 
  mutate(plastic_type = "TPU 2.1")


Control_measurements1<- bind_rows(Control1.181,Control1.2.1, Control1.E)
Control_measurements4<- bind_rows(Control4.181,Control4.2.1, Control4.E)
Control_measurements5<- bind_rows(Control5.181,Control5.2.1, Control5.E)
Control_measurements6<- bind_rows(Control6.181,Control6.2.1, Control6.E)


Total_measurements6$plastic_type<-Total_measurements6$tank.id
Total_measurements6$concentration <-Total_measurements6$tank.id

Total_measurements6<- Total_measurements6 %>% mutate(plastic_type = case_when(
  plastic_type == "A0" ~ "Elastollan",
  plastic_type == "B0" ~ "TPU 181",
  plastic_type == "C0" ~ "TPU 2.1",
  plastic_type == "A1" ~ "Elastollan",
  plastic_type == "A2" ~ "Elastollan",
  plastic_type == "A3" ~ "Elastollan",
  plastic_type == "A4" ~ "Elastollan",
  plastic_type == "A5" ~ "Elastollan",
  plastic_type == "A6" ~ "Elastollan",
  plastic_type == "A7" ~ "Elastollan",
  plastic_type == "A8" ~ "Elastollan", 
  plastic_type == "A9" ~ "Elastollan",
  plastic_type == "B1" ~ "TPU 181",
  plastic_type =="B2" ~ "TPU 181",
  plastic_type == "B3" ~ "TPU 181",
  plastic_type == "B4" ~ "TPU 181",
  plastic_type == "B5" ~ "TPU 181",
  plastic_type == "B6" ~ "TPU 181",
  plastic_type == "B7" ~ "TPU 181",
  plastic_type == "B8" ~ "TPU 181", 
  plastic_type == "B9" ~ "TPU 181",
  plastic_type == "C1" ~ "TPU 2.1",
  plastic_type == "C2" ~ "TPU 2.1",
  plastic_type == "C3" ~ "TPU 2.1",
  plastic_type == "C4" ~ "TPU 2.1",
  plastic_type == "C5" ~ "TPU 2.1",
  plastic_type == "C6" ~ "TPU 2.1",
  plastic_type == "C7" ~ "TPU 2.1",
  plastic_type == "C8" ~ "TPU 2.1", 
  plastic_type == "C9" ~ "TPU 2.1",
))



Total_measurements6<- Total_measurements6 %>% mutate(concentration = case_when(
  concentration == "A0" ~ "0.000",
  concentration == "A1" ~ "0.004",
  concentration == "A2" ~ "0.008",
  concentration == "A3" ~ "0.013",
  concentration == "A4" ~ "0.023",
  concentration == "A5" ~ "0.041",
  concentration == "A6" ~ "0.072",
  concentration == "A7" ~ "0.126",
  concentration == "A8" ~ "0.220", 
  concentration == "A9" ~ "0.385",
  concentration == "B0" ~ "0.000",
  concentration == "B1" ~ "0.004",
  concentration == "B2" ~ "0.008",
  concentration == "B3" ~ "0.013",
  concentration == "B4" ~ "0.023",
  concentration == "B5" ~ "0.041",
  concentration == "B6" ~ "0.072",
  concentration == "B7" ~ "0.126",
  concentration == "B8" ~ "0.220", 
  concentration == "B9" ~ "0.385",
  concentration == "C0" ~ "0.000",
  concentration == "C1" ~ "0.004",
  concentration == "C2" ~ "0.008",
  concentration == "C3" ~ "0.013",
  concentration == "C4" ~ "0.023",
  concentration == "C5" ~ "0.041",
  concentration == "C6" ~ "0.072",
  concentration == "C7" ~ "0.126",
  concentration == "C8" ~ "0.220", 
  concentration == "C9" ~ "0.385",
))

#LOG TRANS CONCENTRATION
Total_measurements6$log.con<-Total_measurements6$tank.id

Total_measurements6<- Total_measurements6 %>% mutate(log.con = case_when(
  log.con == "A0" ~ "0",    
  log.con == "A1" ~ "0.004",
  log.con == "A2" ~ "0.008",
  log.con== "A3"  ~ "0.013",
  log.con == "A4" ~ "0.023",
  log.con == "A5" ~ "0.041",
  log.con == "A6" ~ "0.072",
  log.con == "A7" ~ "0.126",
  log.con == "A8" ~ "0.220", 
  log.con == "A9" ~ "0.385",
  log.con == "B0" ~ "0",
  log.con == "B1" ~ "0.004",
  log.con == "B2" ~ "0.008",
  log.con == "B3" ~ "0.013",
  log.con == "B4" ~ "0.023",
  log.con == "B5" ~ "0.041",
  log.con == "B6" ~ "0.072",
  log.con == "B7" ~ "0.126",
  log.con == "B8" ~ "0.220", 
  log.con == "B9" ~ "0.385",
  log.con == "C0" ~ "0",
  log.con == "C1" ~ "0.004",
  log.con == "C2" ~ "0.008",
  log.con == "C3" ~ "0.013",
  log.con == "C4" ~ "0.023",
  log.con == "C5" ~ "0.041",
  log.con == "C6" ~ "0.072",
  log.con == "C7" ~ "0.126",
  log.con == "C8" ~ "0.220", 
  log.con == "C9" ~ "0.385",
))

#LOG trans above
Total_measurements6$log.con<-as.numeric(Total_measurements6$log.con)
Total_measurements6$log.con<-log(Total_measurements6$log.con+.002)#THIS IS TAKING THE LN!!! IF YOU WANT LOG10 YOU NEED TO SPECIFY LOG10
unique(Total_measurements6$log.con)

hist(Total_measurements6$log.con)

#PRODUCING BIOMASS DATA FIRST NEED TO GET INTO THE CORRECT UNITS...COUNT IS PER L (NEED TO DIVIDE BY 4)
#Count.per.L
Total_measurements1$Count.per.L <-(Total_measurements1$Count)/4
Total_measurements4$Count.per.L <-(Total_measurements4$Count)/4
Total_measurements5$Count.per.L <-(Total_measurements5$Count)/4
Total_measurements6$Count.per.L <-(Total_measurements6$Count)/4
Control_measurements1$Count.per.L <- (Control_measurements1$Count)/4
Control_measurements4$Count.per.L <- (Control_measurements4$Count)/4
Control_measurements5$Count.per.L <- (Control_measurements5$Count)/4
Control_measurements6$Count.per.L <- (Control_measurements6$Count)/4

#########Biomass estimation starting with weights and then biomass in next column
Total_measurements1$weights<-Total_measurements1$Taxa
Total_measurements1 <- Total_measurements1 %>% mutate(weights = case_when(
  weights == "cladoceran" ~ exp(3.93 * log(mean.length) + 4.493),
  weights == "large copepod" ~ exp(2.40 * log(mean.length) + 1.953),
  weights == "daphnia" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "nauplii" ~ exp(2.40* log(mean.length)+1.95 ),
  weights == "ostracoda" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "rotifer" ~ 0.22*(mean.length^3),
  TRUE ~ as.numeric(weights)
))


Total_measurements4$weights<-Total_measurements4$Taxa
Total_measurements4 <- Total_measurements4 %>% mutate(weights = case_when(
  weights == "cladoceran" ~ exp(3.93 * log(mean.length) + 4.493),
  weights == "large copepod" ~ exp(2.40 * log(mean.length) + 1.953),
  weights == "daphnia" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "nauplii" ~ exp(2.40* log(mean.length)+1.95 ),
  weights == "ostracoda" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "rotifer" ~ 0.22*(mean.length^3),
  TRUE ~ as.numeric(weights)
))

Total_measurements5$weights<-Total_measurements5$Taxa
Total_measurements5 <- Total_measurements5 %>% mutate(weights = case_when(
  weights == "cladoceran" ~ exp(3.93 * log(mean.length) + 4.493),
  weights == "large copepod" ~ exp(2.40 * log(mean.length) + 1.953),
  weights == "daphnia" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "nauplii" ~ exp(2.40* log(mean.length)+1.95 ),
  weights == "ostracoda" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "rotifer" ~ 0.22*(mean.length^3),
  TRUE ~ as.numeric(weights)
))

Total_measurements6$weights<-Total_measurements6$Taxa
Total_measurements6 <- Total_measurements6 %>% mutate(weights = case_when(
  weights == "cladoceran" ~ exp(3.93 * log(mean.length) + 4.493),
  weights == "large copepod" ~ exp(2.40 * log(mean.length) + 1.953),
  weights == "daphnia" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "nauplii" ~ exp(2.40* log(mean.length)+1.95 ),
  weights == "ostracoda" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "rotifer" ~ 0.22*(mean.length^3),
  TRUE ~ as.numeric(weights)
))

Control_measurements1$weights<-Control_measurements1$Taxa
Control_measurements1 <- Control_measurements1 %>% mutate(weights = case_when(
  weights == "cladoceran" ~ exp(3.93 * log(mean.length) + 4.493),
  weights == "large copepod" ~ exp(2.40 * log(mean.length) + 1.953),
  weights == "daphnia" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "nauplii" ~ exp(2.40* log(mean.length)+1.95 ),
  weights == "ostracoda" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "rotifer" ~ 0.22*(mean.length^3),
  TRUE ~ as.numeric(weights)
))

Control_measurements4$weights<-Control_measurements4$Taxa
Control_measurements4 <- Control_measurements4 %>% mutate(weights = case_when(
  weights == "cladoceran" ~ exp(3.93 * log(mean.length) + 4.493),
  weights == "large copepod" ~ exp(2.40 * log(mean.length) + 1.953),
  weights == "daphnia" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "nauplii" ~ exp(2.40* log(mean.length)+1.95 ),
  weights == "ostracoda" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "rotifer" ~ 0.22*(mean.length^3),
  TRUE ~ as.numeric(weights)
))

Control_measurements5$weights<-Control_measurements5$Taxa
Control_measurements5 <- Control_measurements5 %>% mutate(weights = case_when(
  weights == "cladoceran" ~ exp(3.93 * log(mean.length) + 4.493),
  weights == "large copepod" ~ exp(2.40 * log(mean.length) + 1.953),
  weights == "daphnia" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "nauplii" ~ exp(2.40* log(mean.length)+1.95 ),
  weights == "ostracoda" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "rotifer" ~ 0.22*(mean.length^3),
  TRUE ~ as.numeric(weights)
))

Control_measurements6$weights<-Control_measurements6$Taxa
Control_measurements6 <- Control_measurements6 %>% mutate(weights = case_when(
  weights == "cladoceran" ~ exp(3.93 * log(mean.length) + 4.493),
  weights == "large copepod" ~ exp(2.40 * log(mean.length) + 1.953),
  weights == "daphnia" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "nauplii" ~ exp(2.40* log(mean.length)+1.95 ),
  weights == "ostracoda" ~ exp(2.83 * log(mean.length) + 1.468),
  weights == "rotifer" ~ 0.22*(mean.length^3),
  TRUE ~ as.numeric(weights)
))

Total_measurements1$biomass<-Total_measurements1$weights*Total_measurements1$Count.per.L
Total_measurements4$biomass<-Total_measurements4$weights*Total_measurements4$Count.per.L
Total_measurements5$biomass<-Total_measurements5$weights*Total_measurements5$Count.per.L
Total_measurements6$biomass<-Total_measurements6$weights*Total_measurements6$Count.per.L
Control_measurements1$biomass<-Control_measurements1$weights*Control_measurements1$Count.per.L
Control_measurements4$biomass<-Control_measurements4$weights*Control_measurements4$Count.per.L
Control_measurements5$biomass<-Control_measurements5$weights*Control_measurements5$Count.per.L
Control_measurements6$biomass<-Control_measurements6$weights*Control_measurements6$Count.per.L

#Sum the biomass

Total_measurements1 <- subset(Total_measurements1, select = c(tank.id, plastic_type, concentration, log.con, biomass))
Total_measurements4 <- subset(Total_measurements4, select = c(tank.id, plastic_type, concentration, log.con, biomass))
Total_measurements5 <- subset(Total_measurements5, select = c(tank.id, plastic_type, concentration, log.con, biomass))
Total_measurements6 <- subset(Total_measurements6, select = c(tank.id, plastic_type, concentration, log.con, biomass))
Control_measurements1 <- subset(Control_measurements1, select = c(tank.id, plastic_type, concentration, log.con, biomass))
Control_measurements4 <- subset(Control_measurements4, select = c(tank.id, plastic_type, concentration, log.con, biomass))
Control_measurements5 <- subset(Control_measurements5, select = c(tank.id, plastic_type, concentration, log.con, biomass))
Control_measurements6 <- subset(Control_measurements6, select = c(tank.id, plastic_type, concentration, log.con, biomass))

Total_measurements1 <- Total_measurements1 %>%
  group_by(tank.id) %>%
  mutate(biomass = sum(biomass))
Total_measurements4 <- Total_measurements4 %>%
  group_by(tank.id) %>%
  mutate(biomass = sum(biomass))
Total_measurements5 <- Total_measurements5 %>%
  group_by(tank.id) %>%
  mutate(biomass = sum(biomass))
Total_measurements6 <- Total_measurements6 %>%
  group_by(tank.id) %>%
  mutate(biomass = sum(biomass))
Control_measurements1 <- Control_measurements1 %>%
  group_by(tank.id) %>%
  mutate(biomass = sum(biomass))
Control_measurements4 <- Control_measurements4 %>%
  group_by(tank.id) %>%
  mutate(biomass = sum(biomass))
Control_measurements5 <- Control_measurements5 %>%
  group_by(tank.id) %>%
  mutate(biomass = sum(biomass))
Control_measurements6 <- Control_measurements6 %>%
  group_by(tank.id) %>%
  mutate(biomass = sum(biomass))

Total_measurements1 <- distinct(Total_measurements1)
Total_measurements4 <- distinct(Total_measurements4)
Total_measurements5 <- distinct(Total_measurements5)
Total_measurements6 <- distinct(Total_measurements6)
Control_measurements1 <- distinct(Control_measurements1)
Control_measurements4 <- distinct(Control_measurements4)
Control_measurements5 <- distinct(Control_measurements5)
Control_measurements6 <- distinct(Control_measurements6)

#A9 in total_measurements4 has biomass of 0ug/L, no row, need to introduce fabricated row
A9.4<- data.frame(
  tank.id="A9",
  plastic_type="Elastollan",
  concentration = 0.385,
  log.con= -0.9493306,
  biomass = 0.0000001
)
str(A9.4)

summary(A9.4)
summary(Total_measurements4)
Total_measurements4$concentration<- as.numeric(Total_measurements4$concentration)

str(Total_measurements4)
Total_measurements4 <- rbind(Total_measurements4, A9.4)
Total_measurements4 <- subset(Total_measurements4, biomass != 0)

hist(Total_measurements1$biomass)
hist(Total_measurements4$biomass)
hist(Total_measurements5$biomass)
hist(Total_measurements6$biomass)

hist(Control_measurements1$biomass)
hist(Control_measurements4$biomass)
hist(Control_measurements5$biomass)
hist(Control_measurements6$biomass)



#Recombining controls with the measurements

Total_measurements1<-rbind(Total_measurements1,Control_measurements1)

str(Control_measurements4)
str(Total_measurements4)

Total_measurements4$plastic_type<-as.factor(Total_measurements4$plastic_type)
Control_measurements4$plastic_type<-as.factor(Control_measurements4$plastic_type)
Control_measurements4$concentration<-as.numeric(Control_measurements4$concentration)

str(Control_measurements4)
str(Total_measurements4)

Total_measurements4<-rbind(Total_measurements4,Control_measurements4)
Total_measurements5<-rbind(Total_measurements5,Control_measurements5)
Total_measurements6<-rbind(Total_measurements6,Control_measurements6)


#################################################################################################
###LOOKING AT DISTRIBUTION### GAMs are starting here #biomass is already in ug/L
#################################################################################################

hist(Total_measurements1$biomass)
hist(Total_measurements4$biomass)
hist(Total_measurements5$biomass)
hist(Total_measurements6$biomass)


hist(log(Total_measurements1$biomass))
hist(log(Total_measurements4$biomass))
hist(log(Total_measurements5$biomass))
hist(log(Total_measurements6$biomass))
#Close to gaussian, check 5


hist(log10(Total_measurements1$biomass))
hist(log10(Total_measurements4$biomass))
hist(log10(Total_measurements5$biomass))
hist(log10(Total_measurements6$biomass))

#converting everything to factors
Total_measurements1$plastic_type<- as.factor(Total_measurements1$plastic_type)
Total_measurements1$log.con<- as.numeric(Total_measurements1$log.con)
Total_measurements1$biomass <- as.numeric(Total_measurements1$biomass)

Total_measurements4$plastic_type<- as.factor(Total_measurements4$plastic_type)
Total_measurements4$log.con<- as.numeric(Total_measurements4$log.con)
Total_measurements4$biomass <- as.numeric(Total_measurements4$biomass)

Total_measurements5$plastic_type<- as.factor(Total_measurements5$plastic_type)
Total_measurements5$log.con<- as.numeric(Total_measurements5$log.con)
Total_measurements5$biomass <- as.numeric(Total_measurements5$biomass)

Total_measurements6$plastic_type<- as.factor(Total_measurements6$plastic_type)
Total_measurements6$log.con<- as.numeric(Total_measurements6$log.con)
Total_measurements6$biomass <- as.numeric(Total_measurements6$biomass)


#START OF GAMS

measurement1.gam1<-gam(biomass~plastic_type+s(log.con, by=plastic_type, k=7), data=Total_measurements1, method='REML', family= Gamma(link="log"))
measurement1.gam2<-gam(biomass~s(log.con, by=plastic_type, k=7), data=Total_measurements1, method='REML', family= Gamma(link="log"))
measurement1.gam3<-gam(biomass~plastic_type+s(log.con, k=7), data=Total_measurements1, method='REML', family= Gamma(link="log"))
measurement1.gam4<-gam(biomass~ s(log.con, k=7), data=Total_measurements1, method='REML', family= Gamma(link="log"))
AIC(measurement1.gam1,measurement1.gam2,measurement1.gam3,measurement1.gam4) ### Biomass data used for MANUSCRIPT [Table S3]


measurement4.gam1<-gam(biomass~plastic_type+s(log.con, by=plastic_type, k=7), data=Total_measurements4, method='REML', family= Gamma(link="log"))
measurement4.gam2<-gam(biomass~s(log.con, by=plastic_type, k=7), data=Total_measurements4, method='REML', family= Gamma(link="log"))
measurement4.gam3<-gam(biomass~plastic_type+s(log.con, k=7), data=Total_measurements4, method='REML', family= Gamma(link="log"))
measurement4.gam4<-gam(biomass~ s(log.con, k=7), data=Total_measurements4, method='REML', family= Gamma(link="log"))
AIC(measurement4.gam1,measurement4.gam2,measurement4.gam3,measurement4.gam4) ### Biomass data used for MANUSCRIPT [Table S3]

####Check log.con###
#Model 1 -Plastic type as a parametric term + smooth of raw concentration per plastic type
measurement4.gam1_raw <- gam(biomass ~ plastic_type + s(concentration, by = plastic_type, k = 7), 
                             data = Total_measurements4, method = 'REML', 
                             family = Gamma(link = "log"))

#Model 2 - Smooth of raw concentration per plastic type (no parametric plastic type effect)
measurement4.gam2_raw <- gam(biomass ~ s(concentration, by = plastic_type, k = 7), 
                             data = Total_measurements4, method = 'REML', 
                             family = Gamma(link = "log"))

#Model 3 -Plastic type as a parametric term + smooth of raw concentration (not by plastic type)
measurement4.gam3_raw <- gam(biomass ~ plastic_type + s(concentration, k = 7), 
                             data = Total_measurements4, method = 'REML', 
                             family = Gamma(link = "log"))

#Model 4 - Smooth of raw concentration only (no plastic type effect)
measurement4.gam4_raw <- gam(biomass ~ s(concentration, k = 7), 
                             data = Total_measurements4, method = 'REML', 
                             family = Gamma(link = "log"))

AIC(measurement4.gam1, measurement4.gam2, measurement4.gam3, measurement4.gam4, 
    measurement4.gam1_raw, measurement4.gam2_raw, measurement4.gam3_raw, measurement4.gam4_raw)

#Check log concentration models
###AIC for non log concentrations is much higher, keep using log.con#########

measurement5.gam1<-gam(biomass~plastic_type+s(log.con, by=plastic_type, k=7), data=Total_measurements5, method='REML', family= Gamma(link="log"))
measurement5.gam2<-gam(biomass~s(log.con, by=plastic_type, k=7), data=Total_measurements5, method='REML', family= Gamma(link="log"))
measurement5.gam3<-gam(biomass~plastic_type+s(log.con, k=7), data=Total_measurements5, method='REML', family= Gamma(link="log"))
measurement5.gam4<-gam(biomass~ s(log.con, k=7), data=Total_measurements5, method='REML', family= Gamma(link="log"))
AIC(measurement5.gam1,measurement5.gam2,measurement5.gam3,measurement5.gam4) ### Biomass data used for MANUSCRIPT [Table S3]

measurement6.gam1<-gam(biomass~plastic_type+s(log.con, by=plastic_type, k=7), data=Total_measurements6, method='REML', family= Gamma(link="log"))
measurement6.gam2<-gam(biomass~s(log.con, by=plastic_type, k=7), data=Total_measurements6, method='REML', family= Gamma(link="log"))
measurement6.gam3<-gam(biomass~plastic_type+s(log.con, k=7), data=Total_measurements6, method='REML', family= Gamma(link="log"))
measurement6.gam4<-gam(biomass~ s(log.con, k=7), data=Total_measurements6, method='REML', family= Gamma(link="log"))
AIC(measurement6.gam1,measurement6.gam2,measurement6.gam3,measurement6.gam4) ### Biomass data used for MANUSCRIPT [Table S13]

measurement1.gam.check<-gam.check(measurement1.gam4) ### Biomass data used for MANUSCRIPT [Table S4]
measurement4.gam.check<-gam.check(measurement4.gam1) ### Biomass data used for MANUSCRIPT [Table S4]
measurement5.gam.check<-gam.check(measurement5.gam1) ### Biomass data used for MANUSCRIPT [Table S4]
measurement6.gam.check<-gam.check(measurement6.gam1) ### Biomass data used for MANUSCRIPT [Table S4]

anova.gam(measurement1.gam4)
anova.gam(measurement4.gam1)
anova.gam(measurement5.gam1)
anova.gam(measurement6.gam1)

plot(measurement1.gam4, all.terms = TRUE, page=1)
plot(measurement4.gam1, all.terms = TRUE, page=1)
plot(measurement5.gam1, all.terms = TRUE, page=1)
plot(measurement6.gam1, all.terms = TRUE, page=1)


#Setting concentrations to look at every concentration
log_con_levels <- c(-0.9493306,-1.5050779,-2.0557250,-2.6036902,-3.1465552, -3.6888795,-4.1997051, -4.6051702, -5.1159958, -6.2146081)

#Calculate pairwise comparisons for each model
pw1.biomass <- emmeans(measurement1.gam4, ~ log.con, at = list(log.con = log_con_levels))
pw4.biomass <- emmeans(measurement4.gam1, ~ plastic_type | log.con, at = list(log.con = log_con_levels))
pw5.biomass <- emmeans(measurement5.gam1, ~ plastic_type | log.con, at = list(log.con = log_con_levels))
pw6.biomass <- emmeans(measurement6.gam1, ~ plastic_type | log.con, at = list(log.con = log_con_levels))

#Perform pairwise comparisons
pw1c.biomass <- pairs(pw1.biomass)
pw4c.biomass <- pairs(pw4.biomass)
pw5c.biomass <- pairs(pw5.biomass)
pw6c.biomass <- pairs(pw6.biomass)

#Table to copy and paste into excel
table1.biomass<-kable(pw1c.biomass)
table4.biomass<-kable(pw4c.biomass)
table5.biomass<-kable(pw5c.biomass)
table6.biomass<-kable(pw6c.biomass)

view(pw1c.biomass)
view(pw4c.biomass)
view(pw5c.biomass)
view(pw6c.biomass)

#Get pairwise comparisons table for tiempoint1
pw1c.biomass <- as.data.frame(pw1c.biomass)
writeClipboard(
  paste(
    paste(names(pw1c.biomass), collapse = "\t"),  # Column headers
    paste(apply(pw1c.biomass, 1, paste, collapse = "\t"), collapse = "\n"),  # Data
    sep = "\n"
  )
)

#For timepoint4
pw4c.biomass <- as.data.frame(pw4c.biomass)
writeClipboard(
  paste(
    paste(names(pw4c.biomass), collapse = "\t"),  # Column headers
    paste(apply(pw4c.biomass, 1, paste, collapse = "\t"), collapse = "\n"),  # Data
    sep = "\n"
  )
)

#For timepoint5
pw5c.biomass <- as.data.frame(pw5c.biomass)
writeClipboard(
  paste(
    paste(names(pw5c.biomass), collapse = "\t"),  # Column headers
    paste(apply(pw5c.biomass, 1, paste, collapse = "\t"), collapse = "\n"),  # Data
    sep = "\n"
  )
)

#For timepoint6
pw6c.biomass <- as.data.frame(pw6c.biomass)
writeClipboard(
  paste(
    paste(names(pw6c.biomass), collapse = "\t"),  # Column headers
    paste(apply(pw6c.biomass, 1, paste, collapse = "\t"), collapse = "\n"),  # Data
    sep = "\n"
  )
)
#########plot.difference##############################
#No T1 because no plastic types added to ponds at T1#

#TIME POINT 4
biomass.t4.Evs21<-plot_difference(
  measurement4.gam1,
  series = log.con,
  difference = list(plastic_type = c("Elastollan", "TPU 2.1")))+
  labs(x="log(Concentration of Plastic)", y="Biomass (Difference smooth)")+
  ggtitle("Differnce between Elastollan and TPU FC2.1, Day 36")+coord_cartesian(ylim = c(-4.6,2.8))

biomass.t4.Evs181<-plot_difference(
  measurement4.gam1,
  series = log.con,
  difference = list(plastic_type = c("Elastollan", "TPU 181")))+
  labs(x="log(Concentration of Plastic)", y="Biomass (Difference smooth)")+
  ggtitle("Differnce between Elastollan and TPU 181, Day 36")+coord_cartesian(ylim = c(-4.6,2.8))

biomass.t4.21vs181<-plot_difference(
  measurement4.gam1,
  series = log.con,
  difference = list(plastic_type = c("TPU 2.1", "TPU 181")))+
  labs(x="log(Concentration of Plastic)", y="Biomass (Difference smooth)")+
  ggtitle("Differnce between TPU FC2.1 and TPU 181, Day 36")+coord_cartesian(ylim = c(-4.6,2.8))

biomass.t4.Evs21 ### Biomass data used for MANUSCRIPT [FIGURE S12]
biomass.t4.Evs181 ### Biomass data used for MANUSCRIPT [FIGURE S12]
biomass.t4.21vs181 ### Biomass data used for MANUSCRIPT [FIGURE S12]

#---------------------------
#TIME POINT 5
biomass.t5.Evs21<-plot_difference(
  measurement5.gam1,
  series = log.con,
  difference = list(plastic_type = c("Elastollan", "TPU 2.1")))+
  labs(x="log(Concentration of Plastic)", y="Biomass (Difference smooth)")+
  ggtitle("Differnce between Elastollan and TPU FC2.1, Day 56")+coord_cartesian(ylim = c(-4.6,2.8))

biomass.t5.Evs181<-plot_difference(
  measurement5.gam1,
  series = log.con,
  difference = list(plastic_type = c("Elastollan", "TPU 181")))+
  labs(x="log(Concentration of Plastic)", y="Biomass (Difference smooth)")+
  ggtitle("Differnce between Elastollan and TPU 181, Day 56")+coord_cartesian(ylim = c(-4.6,2.8))

biomass.t5.21vs181<-plot_difference(
  measurement5.gam1,
  series = log.con,
  difference = list(plastic_type = c("TPU 2.1", "TPU 181")))+
  labs(x="log(Concentration of Plastic)", y="Biomass (Difference smooth)")+
  ggtitle("Differnce between TPU FC2.1 and TPU 181, Day 56")+coord_cartesian(ylim = c(-4.6,2.8))

biomass.t5.Evs21 ### Biomass data used for MANUSCRIPT [FIGURE S12]
biomass.t5.Evs181 ### Biomass data used for MANUSCRIPT [FIGURE S12]
biomass.t5.21vs181 ### Biomass data used for MANUSCRIPT [FIGURE S12]
#---------------------------
#TIME POINT 6
biomass.t6.Evs21<-plot_difference(
  measurement6.gam1,
  series = log.con,
  difference = list(plastic_type = c("Elastollan", "TPU 2.1")))+
  labs(x="log(Concentration of Plastic)", y="Biomass (Difference smooth)")+
  ggtitle("Differnce between Elastollan and TPU FC2.1, Day 97")+coord_cartesian(ylim = c(-4.6,2.8))

biomass.t6.Evs181<-plot_difference(
  measurement6.gam1,
  series = log.con,
  difference = list(plastic_type = c("Elastollan", "TPU 181")))+
  labs(x="log(Concentration of Plastic)", y="Biomass (Difference smooth)")+
  ggtitle("Differnce between Elastollan and TPU 181, Day 97")+coord_cartesian(ylim = c(-4.6,2.8))

biomass.t6.21vs181<-plot_difference(
  measurement6.gam1,
  series = log.con,
  difference = list(plastic_type = c("TPU 2.1", "TPU 181")))+
  labs(x="log(Concentration of Plastic)", y="Biomass (Difference smooth)")+
  ggtitle("Differnce between TPU FC2.1 and TPU 181, Day 97")+coord_cartesian(ylim = c(-4.6,2.8))

biomass.t6.Evs21 ### Biomass data used for MANUSCRIPT [FIGURE S12]
biomass.t6.Evs181 ### Biomass data used for MANUSCRIPT [FIGURE S12]
biomass.t6.21vs181 ### Biomass data used for MANUSCRIPT [FIGURE S12]



#################Figure creation######################


Total_measurements1 <- Total_measurements1 %>%
  mutate(plastic_type = case_when(
    tank.id %in% c("A0", "B0", "C0") ~ "No plastic",
    TRUE ~ plastic_type  # Keep the original value for all other cases
  ))

Total_measurements1 <- Total_measurements1 %>%
  mutate(plastic_type = case_when(
    tank.id %in% c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9") ~ "TPU FC2.1",
    TRUE ~ plastic_type  # Keep the original value for all other cases
  ))

Total_measurements1 <- Total_measurements1 %>%
  mutate(plastic_type = factor(plastic_type, levels = c("No plastic", "Elastollan", "TPU 181", "TPU FC2.1")))

plotplot1<-
  plot_smooths(
    model = measurement1.gam4,
    series = log.con,
    
  )+ 
  geom_point(data=Total_measurements1, 
             aes(x=log.con, y=log(biomass), color=plastic_type))+
  labs(x="ln(Concentration of Plastic)", y="ln(Zooplankton biomass[ug/L])")+
  ggtitle("Concentration of Plastic v.s. Zoop. Biomass, 0 days of exposure")+coord_cartesian(ylim = c(-2,8.1))

plot(plotplot1)

#Change colors
plotplot1_color <- plot_smooths(
  model = measurement1.gam4,
  series = log.con
) + 
  geom_point(data = Total_measurements1, 
             aes(x = log.con, y = log(biomass), color = plastic_type), 
             size = 4) +  # Adjust the size as needed
  scale_color_manual(
    values = c("No plastic" = "black", 
               "Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU FC2.1" = "darkorange"), # Explicit color mapping
    guide = guide_legend(title = "Plastic type")
  ) +
  labs(x = "log(Plastic concentration [g/L])", y = "log(Zooplankton biomass [ug/L])") +
  ggtitle("0 Days of exposure") +
  coord_cartesian(ylim = c(-2, 8.1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
    axis.title.x = element_text(size = 12), # Custom x-axis title
    axis.title.y = element_text(size = 12), # Custom y-axis title
    legend.title = element_text(face = "bold", size = 10), # Custom legend title
    legend.text = element_text(size = 9) # Custom legend text size
  )

#Display the plot
plotplot1_color




#Jitter the two overlayed points for A6 and C6
#Define which specific points (e.g., based on tank.id or other conditions) need jitter
overlapping_points <- Total_measurements1 %>% filter(tank.id %in% c("A6", "C6"))

plotplot1_color <- plot_smooths(
  model = measurement1.gam4,
  series = log.con
) + 
  #Plot all points without jitter
  geom_point(data = Total_measurements1 %>% filter(!tank.id %in% c("A6", "C6")), 
             aes(x = log.con, y = log(biomass), color = plastic_type), 
             size = 4) +  # Points without jitter
  
  #Apply jitter only to the specific overlapping points
  geom_jitter(data = overlapping_points, 
              aes(x = log.con, y = log(biomass), color = plastic_type), 
              size = 4, width = 0.1, height = 0) +  #Jitter applied only to specific overlapping points
  
  #Color mapping
  scale_color_manual(
    values = c("No plastic" = "black", 
               "Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU FC2.1" = "darkorange"),  # Explicit color mapping
    guide = guide_legend(title = "Plastic type")
  ) +
  
  labs(x = "log(Plastic concentration [g/L])", y = "log(Zooplankton biomass [ug/L])") +
  ggtitle("0 Days of exposure") +
  coord_cartesian(ylim = c(-2, 8.1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(size = 12),  # Custom x-axis title
    axis.title.y = element_text(size = 12),  # Custom y-axis title
    legend.title = element_text(face = "bold", size = 10),  # Custom legend title
    legend.text = element_text(size = 9)  # Custom legend text size
  )

#Display the plot
plotplot1_color



#Apply jitter
plotplot1_color <- plot_smooths(
  model = measurement1.gam4,
  series = log.con
) + 
  #Apply jitter except 'No plastic'
  geom_point(data = Total_measurements1 %>% filter(plastic_type == "No plastic"), 
             aes(x = log.con, y = log(biomass), color = plastic_type), 
             size = 4) +  #No jitter for 'No plastic'
  geom_jitter(data = Total_measurements1 %>% filter(plastic_type != "No plastic"), 
              aes(x = log.con, y = log(biomass), color = plastic_type), 
              size = 4, width = 0.1, height = 0) +  #Apply jitter for other types
  scale_color_manual(
    values = c("No plastic" = "black", 
               "Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU FC2.1" = "darkorange"), #Explicit color mapping
    guide = guide_legend(title = "Plastic type")
  ) +
  labs(x = "log(Plastic concentration [g/L])", y = "log(Zooplankton biomass [ug/L])") +
  ggtitle("0 Days of exposure") +
  coord_cartesian(ylim = c(-2, 8.1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

#Display the plot
plotplot1_color



#Change the format for final plot
plotplot1_T1_mod.custom.final <- plot_smooths(
  model = measurement1.gam4,
  series = log.con
) + 
  #Apply jitter for all points except 'No plastic'
  geom_point(data = Total_measurements1 %>% filter(plastic_type == "No plastic"), 
             aes(x = log.con, y = log(biomass), color = plastic_type), 
             size = 4) +  #No jitter for 'No plastic'
  geom_jitter(data = Total_measurements1 %>% filter(plastic_type != "No plastic"), 
              aes(x = log.con, y = log(biomass), color = plastic_type), 
              size = 4, width = 0.1, height = 0) +  #Apply jitter for other types
  
  #Explicit color mapping
  scale_color_manual(
    values = c("No plastic" = "black", 
               "Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU FC2.1" = "darkorange")
  ) +
  
  labs(x = "log(Plastic concentration [g/L])", y = "log(Zooplankton biomass [ug/L])") +
  ggtitle("0 Days of exposure") +
  coord_cartesian(ylim = c(-2, 8.1)) +
  
  # Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.title.x = element_text(size = 14, face = "bold"),  
    axis.title.y = element_text(size = 14, face = "bold"),  
    axis.text.x = element_text(size = 14),                  
    axis.text.y = element_text(size = 14),                  
    panel.grid = element_blank(),                       
    legend.position = "none",                               
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    axis.line = element_line(colour = "black")              
  )

#Display the plot
plotplot1_T1_mod.custom.final ### Biomass data used for MANUSCRIPT [FIGURE 2]




#Next plot
plotplot2<-
  plot_smooths(
    model = measurement4.gam1,
    series = log.con,
    comparison = plastic_type
    
  )+ 
  geom_point(data=Total_measurements4, 
             aes(x=log.con, y=log(biomass), color=plastic_type))+
  labs(x="ln(Concentration of Plastic)", y="ln(Zooplankton biomass[ug/L])")+
  ggtitle("Concentration of Plastic v.s. Zoop. Biomass, 36 days of exposure")+coord_cartesian(ylim = c(-2,8.1))

plot(plotplot2)

#Black control points
plotplot2_color <- plot_smooths(
  model = measurement4.gam1,
  series = log.con,
  comparison = plastic_type
) + 
  #Plot black points for tanks A0, B0, or C0 (the control group)
  geom_point(data = Total_measurements4 %>% filter(tank.id %in% c("A0", "B0", "C0")), 
             aes(x = log.con, y = log(biomass)), 
             color = "black", size = 4) +  # Black points for control tanks
  
  #Plot points for other tanks (non-control)
  geom_point(data = Total_measurements4 %>% filter(!tank.id %in% c("A0", "B0", "C0")), 
             aes(x = log.con, y = log(biomass), color = plastic_type), 
             size = 4) +  # Color-coded points for others
  
  #Define colors for the points and smoothers
  scale_color_manual(
    values = c("Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU 2.1" = "darkorange"),  #Colors for other plastic types
    labels = c("Elastollan", "TPU 181", "TPU FC2.1"),
    guide = guide_legend(title = "Plastic type")
  ) +
  
  #Define fill colors for the smoothers (without black for the control group)
  scale_fill_manual(
    values = c("Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU 2.1" = "darkorange"),  #Same colors for smoothers
    labels = c("Elastollan", "TPU 181", "TPU FC2.1"),
    guide = guide_legend(title = "Plastic type")
  ) +
  
  labs(x = "log(Plastic concentration [g/L])", y = "log(Zooplankton biomass [ug/L])") +
  ggtitle("Day 36") +
  coord_cartesian(ylim = c(-2, 8.1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 12),  
    axis.title.y = element_text(size = 12),  
    legend.title = element_text(face = "bold", size = 10),  
    legend.text = element_text(size = 9)  
  )

#Display the plot
plotplot2_color


#Change format for final plot
plotplot2_T2_mod.custom.final <- plot_smooths(
  model = measurement4.gam1,
  series = log.con,
  comparison = plastic_type
) + 
  #Plot black points for tanks A0, B0, or C0 (the control group)
  geom_point(data = Total_measurements4 %>% filter(tank.id %in% c("A0", "B0", "C0")), 
             aes(x = log.con, y = log(biomass)), 
             color = "black", size = 4) +  # Black points for control tanks
  
  #Plot points for other tanks (non-control)
  geom_point(data = Total_measurements4 %>% filter(!tank.id %in% c("A0", "B0", "C0")), 
             aes(x = log.con, y = log(biomass), color = plastic_type), 
             size = 4) +  # Color-coded points for others
  
  #Define colors for the points and smoothers
  scale_color_manual(
    values = c("Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU 2.1" = "darkorange")  #Colors for other plastic types
  ) +
  
  #Define fill colors for the smoothers (without black for the control group)
  scale_fill_manual(
    values = c("Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU 2.1" = "darkorange")  #Same colors for smoothers
  ) +
  
  labs(x = "log(Plastic concentration [g/L])") +  
  ggtitle("Day 36") +
  coord_cartesian(ylim = c(-2, 8.1)) +
  
  # Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.title.x = element_text(size = 14, face = "bold"),  
    axis.text.x = element_text(size = 14),                  
    axis.title.y = element_blank(),                         
    axis.text.y = element_blank(),                          
    panel.grid = element_blank(),                           
    legend.position = "none",                               
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    axis.line = element_line(colour = "black")              
  )

# Display the plot
plotplot2_T2_mod.custom.final ### Biomass data used for MANUSCRIPT [FIGURE 2]




#next time point
plotplot3<-
  plot_smooths(
    model = measurement5.gam1,
    series = log.con,
    comparison = plastic_type
    
  )+ 
  geom_point(data=Total_measurements5, 
             aes(x=log.con, y=log(biomass), color=plastic_type))+
  labs(x="ln(Concentration of Plastic)", y="ln(Zooplankton biomass[ug/L])")+
  ggtitle("Concentration of Plastic v.s. Zoop. Biomass, 56 days of exposure")+coord_cartesian(ylim = c(-2,8.1))

plot(plotplot3)

#Format update
plotplot3_color <- plot_smooths(
  model = measurement5.gam1,
  series = log.con,
  comparison = plastic_type
) + 
  #Plot black points for tanks A0, B0, or C0 (the control group)
  geom_point(data = Total_measurements5 %>% filter(tank.id %in% c("A0", "B0", "C0")), 
             aes(x = log.con, y = log(biomass)), 
             color = "black", size = 4) +  #Black points for control tanks
  
  #Plot points for other tanks (non-control)
  geom_point(data = Total_measurements5 %>% filter(!tank.id %in% c("A0", "B0", "C0")), 
             aes(x = log.con, y = log(biomass), color = plastic_type), 
             size = 4) +  #Color-coded points for others
  
  #Define colors for the points and smoothers
  scale_color_manual(
    values = c("Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU 2.1" = "darkorange"),  #Colors for other plastic types
    labels = c("Elastollan", "TPU 181", "TPU FC2.1"),  
    guide = guide_legend(title = "Plastic type")
  ) +
  
  #Define fill colors for the smoothers (without black for the control group)
  scale_fill_manual(
    values = c("Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU 2.1" = "darkorange"),  #Same colors for smoothers
    labels = c("Elastollan", "TPU 181", "TPU FC2.1"),  
    guide = guide_legend(title = "Plastic type")
  ) +
  
  labs(x = "log(Plastic concentration [g/L])", y = "log(Zooplankton biomass [ug/L])") +
  ggtitle("Day 56") +
  coord_cartesian(ylim = c(-2, 8.1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.title.x = element_text(size = 12),  
    axis.title.y = element_text(size = 12),  
    legend.title = element_text(face = "bold", size = 10),  
    legend.text = element_text(size = 9)  
  )

#Display the plot
plotplot3_color


#Final plot updates
plotplot3_T3_mod.custom.final <- plot_smooths(
  model = measurement5.gam1,
  series = log.con,
  comparison = plastic_type
) + 
  #Plot black points for tanks A0, B0, or C0 (the control group)
  geom_point(data = Total_measurements5 %>% filter(tank.id %in% c("A0", "B0", "C0")), 
             aes(x = log.con, y = log(biomass)), 
             color = "black", size = 4) +  # Black points for control tanks
  
  #Plot points for other tanks (non-control)
  geom_point(data = Total_measurements5 %>% filter(!tank.id %in% c("A0", "B0", "C0")), 
             aes(x = log.con, y = log(biomass), color = plastic_type), 
             size = 4) +  #Color-coded points for others
  
  #Define colors for the points and smoothers
  scale_color_manual(
    values = c("Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU 2.1" = "darkorange")  #Colors for other plastic types
  ) +
  
  #Define fill colors for the smoothers (without black for the control group)
  scale_fill_manual(
    values = c("Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU 2.1" = "darkorange")  #Same colors for smoothers
  ) +
  
  labs(x = "log(Plastic concentration [g/L])") +  
  ggtitle("Day 56") +
  coord_cartesian(ylim = c(-2, 8.1)) +
  
  #Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.title.x = element_text(size = 14, face = "bold"),  
    axis.text.x = element_text(size = 14),                  
    axis.title.y = element_blank(),                         
    axis.text.y = element_blank(),                          
    panel.grid = element_blank(),                           
    legend.position = "none",                               
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    axis.line = element_line(colour = "black")              
  )

#Display the plot
plotplot3_T3_mod.custom.final ### Biomass data used for MANUSCRIPT [FIGURE 2]









#Next time point
plotplot4<-
  plot_smooths(
    model = measurement6.gam1,
    series = log.con,
    comparison = plastic_type
    
  )+ 
  geom_point(data=Total_measurements6, 
             aes(x=log.con, y=log(biomass), color=plastic_type))+
  labs(x="ln(Concentration of Plastic)", y="ln(Zooplankton biomass[ug/L])")+
  ggtitle("Concentration of Plastic v.s. Zoop. Biomass, 97 days of exposure")+coord_cartesian(ylim = c(-2,8.1))

plot(plotplot4)

#Change format
plotplot4_color <- plot_smooths(
  model = measurement6.gam1,
  series = log.con,
  comparison = plastic_type
) + 
  #Plot black points for tanks A0, B0, or C0 (the control group)
  geom_point(data = Total_measurements6 %>% filter(tank.id %in% c("A0", "B0", "C0")), 
             aes(x = log.con, y = log(biomass)), 
             color = "black", size = 4) +  #Black points for control tanks
  
  #Plot points for other tanks (non-control)
  geom_point(data = Total_measurements6 %>% filter(!tank.id %in% c("A0", "B0", "C0")), 
             aes(x = log.con, y = log(biomass), color = plastic_type), 
             size = 4) +  #Color-coded points for others
  
  #Define colors for the points and smoothers
  scale_color_manual(
    values = c("Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU 2.1" = "darkorange"),  #Colors for other plastic types
    labels = c("Elastollan", "TPU 181", "TPU FC2.1"),  # Update legend
    guide = guide_legend(title = "Plastic type")
  ) +
  
  #Define fill colors for the smoothers (without black for the control group)
  scale_fill_manual(
    values = c("Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU 2.1" = "darkorange"),  #Same colors for smoothers
    labels = c("Elastollan", "TPU 181", "TPU FC2.1"),  
    guide = guide_legend(title = "Plastic type")
  ) +
  
  labs(x = "log(Plastic concentration [g/L])", y = "log(Zooplankton biomass [ug/L])") +
  ggtitle("Day 97") +
  coord_cartesian(ylim = c(-2, 8.1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.title.x = element_text(size = 12),  
    axis.title.y = element_text(size = 12),  
    legend.title = element_text(face = "bold", size = 10),  
    legend.text = element_text(size = 9)  
  )

#Display the plot
plotplot4_color

#Final format updates
plotplot4_T4_mod.custom.final <- plot_smooths(
  model = measurement6.gam1,
  series = log.con,
  comparison = plastic_type
) + 
  #Plot black points for tanks A0, B0, or C0 (the control group)
  geom_point(data = Total_measurements6 %>% filter(tank.id %in% c("A0", "B0", "C0")), 
             aes(x = log.con, y = log(biomass)), 
             color = "black", size = 4) +  #Black points for control tanks
  
  #Plot points for other tanks (non-control)
  geom_point(data = Total_measurements6 %>% filter(!tank.id %in% c("A0", "B0", "C0")), 
             aes(x = log.con, y = log(biomass), color = plastic_type), 
             size = 4) +  #Color-coded points for others
  
  #Define colors for the points and smoothers
  scale_color_manual(
    values = c("Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU 2.1" = "darkorange")  #Colors for other plastic types
  ) +
  
  #Define fill colors for the smoothers (without black for the control group)
  scale_fill_manual(
    values = c("Elastollan" = "darkorchid4", 
               "TPU 181" = "deeppink3", 
               "TPU 2.1" = "darkorange")
  ) +
  
  labs(x = "log(Plastic concentration [g/L])") +  
  ggtitle("Day 97") +
  coord_cartesian(ylim = c(-2, 8.1)) +
  
  #Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.title.x = element_text(size = 14, face = "bold"),  
    axis.text.x = element_text(size = 14),                  
    axis.title.y = element_blank(),                         
    axis.text.y = element_blank(),                          
    panel.grid = element_blank(),                           
    legend.position = "none",                               
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    axis.line = element_line(colour = "black")              
  )

#Display the plot
plotplot4_T4_mod.custom.final ### Biomass data used for MANUSCRIPT [FIGURE 2]



#combine plotplot1_T1_mod.custom.final, plotplot2_T2_mod.custom.final, plotplot3_T3_mod.custom.final, and plotplot4_T4_mod.custom.final
#Format them similarly to the abundance plots.
plotplot1_T1_mod.custom.final <- plotplot1_T1_mod.custom.final + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),  
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  
        plot.margin = margin(0, 0, 0, 0))

plotplot2_T2_mod.custom.final <- plotplot2_T2_mod.custom.final + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),  
        axis.text.y = element_blank(),   
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  
        plot.margin = margin(0, 0, 0, 0))

plotplot3_T3_mod.custom.final <- plotplot3_T3_mod.custom.final + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),  
        axis.text.y = element_blank(),   
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  
        plot.margin = margin(0, 0, 0, 0))

plotplot4_T4_mod.custom.final <- plotplot4_T4_mod.custom.final + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),  
        axis.text.y = element_blank(),   
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  
        plot.margin = margin(0, 0, 0, 0))

#Combine the plots horizontally
combined_plot_biomass <- (plotplot1_T1_mod.custom.final + 
                            plotplot2_T2_mod.custom.final + 
                            plotplot3_T3_mod.custom.final + 
                            plotplot4_T4_mod.custom.final) +
  plot_layout(ncol = 4)

#Create an empty plot to serve as a placeholder for the single x-axis label
x_label_plot_biomass <- ggplot() + 
  theme_void() +   #Remove all elements from this empty plot
  labs(x = "log(Concentration of plastic [g/L])") +  
  theme(axis.title.x = element_text(size = 14, face = "bold", hjust = 0.5))  

#Combine the plots with the title and x-axis label, reducing the space more
final_plot_biomass <- (combined_plot_biomass / x_label_plot_biomass) + 
  plot_layout(heights = c(10, 0.5)) +  
  
  #This part controls ONLY the overall title (Biomass over exposure periods)
  plot_annotation(title = "Zooplankton biomass") &
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  

#Display the final combined plot
print(final_plot_biomass) ### Biomass data used for MANUSCRIPT [FIGURE 2]



#Cheking other factors
plot(plotplot1)
plot(plotplot2)
plot(plotplot3)
plot(plotplot4)

plot(resid(measurement1.gam1))
plot(Total_measurements1$log.con, resid(measurement1.gam1), xlab = "log.con", ylab = "Residuals")


