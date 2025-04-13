#This R script will import the exports from my QIIME2 pipeline script so that they can be combined into a single Phyloseq object for further processing, object extraction, and statistical analysis. 
#Lines beginning with (#) were not in the final script

#Clear work space if need
rm(list = ls())

#check R version
R.version.string
#install packages
#install.packages("BiocManager")
BiocManager::install("microbiome")
install.packages("vegetarian")
#install.packages("/Users/lab/vegetarian", repos = NULL, type = "source")
BiocManager::install("phyloseq")
#BiocManager::install("phyloseq", force = TRUE)
install.packages("gdata")
install.packages("ecodist")
install.packages("ape")
install.packages("phytools")
install.packages("castor")
install.packages("doParallel")
install.packages("viridis")
install.packages("RColorBrewer")
BiocManager::install("DESeq2")
install.packages("indicspecies")
install.packages("ggrepel")
install.packages("tibble")
BiocManager::install("metagMisc")
install.packages("ggplot2")
install.packages("vegan")
install.packages("carData")
install.packages("car")
BiocManager::install("biomformat")
install.packages("dplyr")
install.packages("maps")
install.packages("Rcpp")
install.packages("tidyr")

#load packages
library(phyloseq); packageVersion("phyloseq")
library(microbiome)
library(vegetarian)
library(ggplot2); packageVersion("ggplot2")
library(gdata)
library(ecodist)
library(vegan) 
library(carData)
library(car)   
library(dplyr)
library(biomformat)
library(ape) 
library(maps)
library(phytools) 
library(Rcpp)
library(castor)
library(doParallel) 
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(metagMisc)
library(tidyr)
library(DESeq2)
library(indicspecies)
library(tibble)

#Set seed
set.seed(77)

#set working directory
setwd("C:/Users/DELL/Documents/R/projects/MP_16S")

#Read in biom table
ASV_reads <- read_biom("16S_Final_feature_table.biom") #uses biomformat package
ASV_table <- as.data.frame(as.matrix(biom_data(ASV_reads))) #shape shifting
otu_tab<-t(ASV_table) #transpose cols and rows
dim(otu_tab) #are the dimensions correct? #150 samps x 7317 ASVS 
colnames(otu_tab) #what about col names? #feature IDs
rownames(otu_tab) #what about row names? #Sample names
df.ASV<-as.data.frame(otu_tab)
df.ASV[1:10,1:10]
rownames(df.ASV) #nothing has changed
colnames(df.ASV) #nothing has changed

#Read in metadata file
metadata<-read.csv("16S_BioP_metadata.csv",header=TRUE)
dim(metadata) #150 samps x 11 columns 

#Read in tax table
taxa2<-read.csv("16S_taxonomy_Final.csv")
ph.headers<-c("ID", "Kingdom","Phylum","Class","Order","Family","Genus","species")
colnames(taxa2)<-ph.headers
rownames(taxa2) #just numbers
rownames(taxa2)<-taxa2[,1]
rownames(taxa2) #much better
dim(taxa2) #7317 x 8
taxa3<-taxa2[,2:8]
colnames(taxa3)

#Read in tree
dated.16Stree<-read.tree(file="BioP-rooted-tree-FINAL.nwk")
is.rooted(dated.16Stree) #TRUE
sample_names(dated.16Stree) #NULL
dated.16Stree$tip.label #for example "bf08d62b32cced86e829cba893bdf318" 

#Creating phyloseq object
str(taxa3) #data.frame
taxa3.m<-as.matrix(taxa3) #check matrix 
str(taxa3.m)
colnames(taxa3.m)
str(df.ASV) #data.frame 150 obs. of 7217 vars
str(metadata) #data.frame 150 obs. of 11 vars
colnames(metadata)

rownames(df.ASV)<-as.character(rownames(df.ASV))
as.character(rownames(df.ASV))==as.character(metadata[,1]) # Now TRUE, these need to match; sort your sample_name column in excel from low to high.
colnames(df.ASV) #accession numbers
rownames(metadata)<-as.character(metadata[,1])
rownames(taxa3.m)<-as.character(rownames(taxa3.m)) #accession numbers
samp.names<-as.character(metadata[,1])

#To set up sample names to match (originally marked NULL)
sample_names(df.ASV)<-samp.names
sample_names(metadata)<-samp.names
sample_names(taxa3.m)<-samp.names
sample_names(dated.16Stree)<-samp.names

############################

#Check if rownames of df.ASV match sample names
identical(rownames(df.ASV), samp.names)  # Should be TRUE
#Check if rownames of metadata match sample names
identical(rownames(metadata), samp.names)  # Should be TRUE
#Check if column names of df.ASV match the ASV IDs in taxa3
identical(colnames(df.ASV), rownames(taxa3.m))  #Should be TRUE

#Fix the taxonomy matrix to match the OTU table
taxa3.m <- taxa3.m[colnames(df.ASV), ]  #Align rows with OTU columns
#Reconfirm
identical(colnames(df.ASV), rownames(taxa3.m))  #Should return TRUE

#Check if tree labels match ASV IDs
identical(sort(dated.16Stree$tip.label), sort(colnames(df.ASV)))  #Should be TRUE

#Here is the actual phyloseq object
BioP.phylo<-phyloseq(otu_table(df.ASV, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxa3.m), phy_tree(dated.16Stree))

#Confirm
identical(sample_names(BioP.phylo), samp.names)  #Should be TRUE

###############################

#Replace 'taxa names' to short hand so that they are easier to view in R
BioP.phylo@otu_table[1:10,1:10]
dim(BioP.phylo@otu_table) #150 x 7317
rowSums(BioP.phylo@otu_table) 
colSums(BioP.phylo@otu_table) 

#Check taxa data and filter
colnames(BioP.phylo@tax_table) #"Kingdom", "Phylum", "Class", "Family", "Genus", "species"  
table(tax_table(BioP.phylo)[, "Kingdom"], exclude = NULL) 
table(tax_table(BioP.phylo)[, "Phylum"], exclude = NULL) 
table(tax_table(BioP.phylo)[, "Class"], exclude = NULL) 
table(tax_table(BioP.phylo)[, "Order"], exclude = NULL) #Reads assigned as chloroplast
table(tax_table(BioP.phylo)[, "Family"], exclude = NULL) #Reads assigned as mitochondria
table(tax_table(BioP.phylo)[, "Genus"], exclude = NULL)
table(tax_table(BioP.phylo)[, "species"], exclude = NULL)

#Filter all to Phylum - Keep Archaea
#Trim up leading or trailing spaces - Some found good, to run for safety
tax_table(BioP.phylo) <- apply(tax_table(BioP.phylo), 2, trimws)
#Keep sequences with assigned Kingdoms (removes Unassigned)
p1 <- subset_taxa(BioP.phylo, !Kingdom %in% "Unassigned")
#Remove Eukaryotic sequences
p2 <- subset_taxa(p1, !Kingdom %in% "Eukaryota")
table(tax_table(p2)[, "Kingdom"], exclude = NULL) #Check
#Remove chloroplast contaminants
p3 <- subset_taxa(p2, !Order %in% "Chloroplast")
table(tax_table(p3)[, "Order"], exclude = NULL) #Check
#Remove mitochondrial contaminants
p4 <- subset_taxa(p3, !Family %in% "Mitochondria")
table(tax_table(p4)[, "Family"], exclude = NULL) #Check
#Remove poorly classified taxa at the Phylum level
p5 <- subset_taxa(p4, !is.na(Phylum) &
                    !Phylum %in% c("Unassigned", "unknown", "Unknown", "uncultured", "NA"))
table(tax_table(p5)[, "Phylum"], exclude = NULL) #Check

###########################
#Check to make sure the chloroplasts were over abundant and a real containment
table(tax_table(p2)[, "Order"])
chloroplast_taxa <- taxa_names(subset_taxa(p2, Order == "Chloroplast"))
length(chloroplast_taxa)  #How many ASVs match
#Check if taxa are rows
taxa_are_rows(p2)
#Convert OTU table to matrix and transpose if required
otu_mat <- as(otu_table(p2), "matrix")
otu_mat <- t(otu_mat)  #Transpose to make taxa the rows
#Subset and sum chloroplast reads per sample
chloroplast_counts <- colSums(otu_mat[chloroplast_taxa, , drop = FALSE])
#Total reads per sample
total_counts <- sample_sums(p2)
#Relative abundance of chloroplast reads
chloroplast_rel_abund <- chloroplast_counts / total_counts
#View summary
summary(chloroplast_rel_abund) #Yes, mean is 13%, good to remove the chloroplasts


#Continue
p5@sam_data
rowSums(p5@otu_table)

#Assign cleaned object
BioPlastics_phylo <- p5
dim(BioPlastics_phylo@otu_table)  #150 x 6,114

#Get OTU table and convert to data frame for rarefaction step
BioPlasticsASVtable <- BioPlastics_phylo@otu_table
BioPlasticsASVtable1 <- BioPlasticsASVtable  #Naming
dim(BioPlasticsASVtable1) #150
colSums(BioPlasticsASVtable1)
min(colSums(BioPlasticsASVtable1)) #Need at least 1
BioP_ASV_table.df <- as.data.frame(BioPlasticsASVtable1)

#Find the minimum read depth
min(rowSums(BioPlasticsASVtable1)) #41,255
sort(rowSums(BioPlasticsASVtable1),decreasing=FALSE)

##############################
### Multi-rarefaction Step ###
##############################

# Install and load metagMisc if needed
#if (!requireNamespace("metagMisc", quietly = TRUE)) {
  #devtools::install_github("vmikk/metagMisc")
#}
#library(metagMisc)

# Find rarefaction depth
min_depth <- min(sample_sums(BioPlastics_phylo))  
print(min_depth) #41255
taxa_are_rows(BioPlastics_phylo) #Check if taxa are rows
if (!taxa_are_rows(BioPlastics_phylo)) { #Check
  otu_table(BioPlastics_phylo) <- t(otu_table(BioPlastics_phylo))
  taxa_are_rows(BioPlastics_phylo) <- TRUE
}
#Reconfirm
taxa_are_rows(BioPlastics_phylo)

#Run multiple rarefaction and average
set.seed(91)
averaged_phylo <- phyloseq_mult_raref_avg(
  BioPlastics_phylo,
  SampSize = min_depth,
  iter = 100,
  parallel = TRUE
)

#Create rarefied ASV-level phyloseq object
UniFrac_Phylo_Object <- phyloseq(otu_table(averaged_phylo, taxa_are_rows=FALSE),
                                 phy_tree(BioPlastics_phylo@phy_tree),
                                 tax_table(BioPlastics_phylo@tax_table),
                                 sample_data(BioPlastics_phylo@sam_data))

#Fix tree if needed
if (!ape::is.binary(phy_tree(UniFrac_Phylo_Object))) {
  phy_tree(UniFrac_Phylo_Object) <- ape::multi2di(phy_tree(UniFrac_Phylo_Object))
}

#Aggregate from rarefied ASV object - for multilevel downstream analysis
BioPlastics_family <- tax_glom(UniFrac_Phylo_Object, taxrank = "Family")
BioPlastics_phylum <- tax_glom(UniFrac_Phylo_Object, taxrank = "Phylum")

#Weighted UniFrac for ASV level
registerDoParallel(cores=4)
set.seed(88)
wdistUni <- UniFrac(UniFrac_Phylo_Object, weighted=TRUE, parallel=TRUE, fast=TRUE)

#Bray-Curtis from ASV level
set.seed(71)
bray_phylo <- transform_sample_counts(UniFrac_Phylo_Object, function(x) x / sum(x))
bray_dist <- phyloseq::distance(bray_phylo, method = "bray")

#Bray-Curtis from Family level
set.seed(76)
bray_family_phylo <- transform_sample_counts(BioPlastics_family, function(x) x / sum(x))
bray_dist_family <- phyloseq::distance(bray_family_phylo, method = "bray")

#Bray-Curtis from Phylum level
set.seed(99)
bray_phylum_phylo <- transform_sample_counts(BioPlastics_phylum, function(x) x / sum(x))
bray_dist_phylum <- phyloseq::distance(bray_phylum_phylo, method = "bray")

#Weighted UniFrac for Family and Phylum levels
set.seed(65)
wdistUni_family <- UniFrac(BioPlastics_family, weighted=TRUE, parallel=TRUE, fast=TRUE)
set.seed(66)
wdistUni_phylum <- UniFrac(BioPlastics_phylum, weighted=TRUE, parallel=TRUE, fast=TRUE)

###################################################
#Visualize the data before running further analysis
###################################################

#Transforming phyloseq data for ggplot2
physeq_df <- p5 %>%
  transform_sample_counts(function(x) x / sum(x)) %>% #Normalizes the counts
  psmelt() %>%
  filter(!is.na(Phylum)) #Makes it so only rows with defined Phyla are considered - we already filtered so gtg

#Example date order
custom_order <- c("March_7", "March_23", "April_12", "May_3", "June_9")

#Converting date to a factor with a specified order
physeq_df <- physeq_df %>%
  mutate(date = factor(date, levels = custom_order))

#Creating the plot
phyplot <- ggplot(physeq_df, aes(fill=Phylum, y=Abundance, x=date)) +
  geom_bar(stat="identity", position="fill") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample Date", title = "Phylum Abundance Over Time by Plastic Type") +
  facet_grid(.~plastic_type, space="free", scales="free") #Adjust the faceting based on metadata

print(phyplot)

#Aggregating data by phylum, plastic_type, and plastic_conc
physeq_df2 <- physeq_df %>%
  group_by(date, Phylum, plastic_type, plastic_conc) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop')

#Create stacked bar plot with flipped facet display
phyplot2 <- ggplot(physeq_df2, aes(x=date, y=Abundance, fill=Phylum)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(plastic_type ~ plastic_conc, scales="free_x", space="free") +  #Flipped facet orientation
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample Date", y = "Relative Abundance", title = "Microbial Diversity by Plastic Type and Concentration [g/L] Over Time")

#Print plot
print(phyplot2)

#Aggregating data by phylum, plastic_type, and plastic_conc
physeq_df3 <- physeq_df2 %>%
  group_by(date, Phylum, plastic_type, plastic_conc) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
  #Update the plastic_type to desired labels
  mutate(plastic_type = case_when(
    plastic_type == "TPU181" ~ "TPU 181",
    plastic_type == "TPUFC2.1" ~ "TPU FC2.1",
    TRUE ~ plastic_type
  ),
  #Make sure plastic concentration is shown with 3 decimal places
  plastic_conc = sprintf("%.3f", plastic_conc),
  #Replace the dates with corresponding day numbers
  date = case_when(
    date == "March_7" ~ 0,
    date == "March_23" ~ 16,
    date == "April_12" ~ 36,
    date == "May_3" ~ 56,
    date == "June_9" ~ 94,
    TRUE ~ as.numeric(date)  #In case there are additional dates - there are not so gtg
  ),
  #Convert date to a factor to avoid spacing issues
  date = factor(date, levels = c(0, 16, 36, 56, 94))  
  )

#Creating stacked bar plot with updated labels and plasma color scheme
phyplot3 <- ggplot(physeq_df3, aes(x = date, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(plastic_type ~ plastic_conc, scales = "free_x", space = "free") +  #Flipped facet orientation
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for date
  scale_fill_viridis(option = "plasma", discrete = TRUE) +  #Use plasma color scheme
  #Update theme to bold titles, remove grids, and add box around each facet
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
    axis.text.x = element_text(size = 9, angle = 90, hjust = 0.5),  #Center-align x-axis tick labels
    axis.text.y = element_text(size = 9),  #Increase size of y-axis tick labels
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),  #Bold and center the title
    legend.title = element_text(face = "bold", size = 10),  #Bold legend title
    legend.text = element_text(size = 8),  #Customize legend text size
    legend.key.size = unit(0.6, "lines"),  #Adjust spacing between legend items
    legend.key.height = unit(0.4, "lines"),  #Reduce height of legend keys
    panel.background = element_blank(),  #Remove grey background
    panel.grid.major = element_blank(),  #Remove major gridlines
    panel.grid.minor = element_blank(),  #Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  #Box around each plot with correct linewidth argument
  ) +
  guides(fill = guide_legend(ncol = 1)) +  #Adjust number of columns in legend
  labs(x = "Day Number", y = "Relative Abundance", title = "Microbial diversity")

#Print plot
print(phyplot3)

#Custom color palette with 42 colors
custom_col42 <- c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3",
                  "#114578","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0",
                  "#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0",
                  "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4",
                  "#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098",
                  "#F7F7C5","#784511","#A55E18","#D2781E","#E4913F","#EAAB6C","#F0C498",
                  "#781122","#A5182F","#D21E2C","#E43F5B","#EA6C81","#F098A7", "black")

#Custom color palette with 27 colors
custom_col29 <- c("#781156", "#A51876", "#D21E96", "#E43FAD", "#EA6CC0", "#F098D3",  
                  "#114578", "#185EA5", "#1E78D2", "#3F91E4", "#6CABEA", "#98C4F0",  
                  "#D2D21E", "#E4E43F", "#EAEA6C", "#F0F098",                      
                  "#A55E18", "#D2781E", "#E4913F", "#EAAB6C", "#F0C498",            
                  "#A5182F", "#D21E2C", "#E43F5B", "#EA6C81", "#F098A7", "black", "darkgreen", "green3")            

#Custom color palette with 42 colors
custom_col29_used <- c("#114578","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0",
                  "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4",
                  "#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C",
                  "#A55E18","#D2781E","#E4913F","#EAAB6C","#F0C498",
                  "#781122","#A5182F","#D21E2C","#E43F5B","#EA6C81","#F098A7", "black")

#Stacked bar plot with custom colors for 27 phyla
phyplot4 <- ggplot(physeq_df3, aes(x = date, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(plastic_type ~ plastic_conc, scales = "free_x", space = "free") +  #Flipped facet orientation
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for date
  scale_fill_manual(values = custom_col29_used) +  #Use custom color palette with 27 colors
  #Update theme to bold titles, remove grids, and add box around each facet
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5),  #Center-align x-axis tick labels
    axis.text.y = element_text(size = 7),  #Increase size of y-axis tick labels
    legend.title = element_text(face = "bold", size = 10),  #Bold legend title
    legend.text = element_text(size = 8),  #Customize legend text size
    legend.key.size = unit(0.6, "lines"),  #Adjust spacing between legend items
    legend.key.height = unit(0.4, "lines"),  #Reduce height of legend keys
    panel.background = element_blank(),  #Remove grey background
    panel.grid.major = element_blank(),  #Remove major gridlines
    panel.grid.minor = element_blank(),  #Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  #Box around each plot with correct linewidth argument
  ) +
  guides(fill = guide_legend(ncol = 1)) +  #Adjust number of columns in legend
  labs(x = "Day Number", y = "Relative Abundance")  #Removed the title

#Print plot
print(phyplot4) #MANUSCRIPT [Supplementary Fig. 5]

#Move legend
phyplot4 <- ggplot(physeq_df3, aes(x = date, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +  #Increase bar width
  facet_grid(plastic_type ~ plastic_conc, scales = "free_x", space = "free") +  # Flipped facet orientation
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for date
  scale_fill_manual(values = custom_col29_used) +  #Use custom color palette with 27 colors
  #Update theme to bold titles, remove grids, add box around each facet, and modify legend position
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
    axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),  #Center-align x-axis tick labels
    axis.text.y = element_text(size = 8),  #Increase size of y-axis tick labels
    legend.title = element_text(face = "bold", size = 10),  #Bold legend title
    legend.text = element_text(size = 8),  #Customize legend text size
    legend.key.size = unit(0.6, "lines"),  #Adjust spacing between legend items
    legend.key.height = unit(0.4, "lines"),  #Reduce height of legend keys
    legend.position = "bottom",  #Move legend below the plot
    panel.background = element_blank(),  #Remove grey background
    panel.grid.major = element_blank(),  #Remove major gridlines
    panel.grid.minor = element_blank(),  #Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  #Box around each plot with correct linewidth argument
  ) +
  guides(fill = guide_legend(ncol = 7)) +  #Adjust number of columns in legend for better fit
  labs(x = "Day Number", y = "Relative Abundance")  #Removed the title

#Print plot
print(phyplot4)


#Aggregating data by phylum, plastic_type, and plastic_conc
physeq_df4 <- physeq_df %>%
  group_by(date, Phylum, Family, plastic_type, plastic_conc) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop') %>%
  #Update the plastic_type to correct labels
  mutate(plastic_type = case_when(
    plastic_type == "TPU181" ~ "TPU 181",
    plastic_type == "TPUFC2.1" ~ "TPU FC2.1",
    TRUE ~ plastic_type
  ),
  #Plastic concentration shown with 3 decimal places
  plastic_conc = sprintf("%.3f", plastic_conc),
  #Replace the dates with corresponding day numbers
  date = case_when(
    date == "March_7" ~ 0,
    date == "March_23" ~ 16,
    date == "April_12" ~ 36,
    date == "May_3" ~ 56,
    date == "June_9" ~ 94,
    TRUE ~ as.numeric(date)
  ),
  #Convert date to a factor to avoid spacing issues
  date = factor(date, levels = c(0, 16, 36, 56, 94))  #Discrete date
  )

#Generate 300 distinct colors using colorRampPalette()
custom_col300 <- colorRampPalette(c("red", "blue", "green", "orange", "purple", "brown", "pink", "yellow", "cyan", "black"))(300)

#Create bar plot faceted by phylum, with families as fill
phyplot_family <- ggplot(physeq_df4, aes(x = date, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +  #Stack families for each date
  facet_wrap(~ Phylum, scales = "free_y") +  #Create a separate facet for each phylum
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for date
  scale_fill_manual(values = custom_col300) +  #Use custom color palette for families
  #Update theme to bold titles, remove grids, and add box around each facet
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),  #Tilt x-axis labels to reduce overlap
    axis.text.y = element_text(size = 7),  #Increase size of y-axis tick labels
    legend.position = "none",  #Remove the legend
    panel.background = element_blank(),  #Remove grey background
    panel.grid.major = element_blank(),  #Remove major gridlines
    panel.grid.minor = element_blank(),  #Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  #Box around each plot with correct linewidth argument
  ) +
  labs(x = "Day Number", y = "Relative Abundance", title = "Family-level Abundance for Each Phylum")  # Add plot title

#Print plot
print(phyplot_family)

##Look over all 27 Phyla
#Generate distinct colors for all the families (we need 300 colors)
custom_col300 <- colorRampPalette(c("red", "blue", "green", "orange", "purple", "brown", "pink", "yellow", "cyan", "black"))(300)

#Get unique values of plastic types and concentrations
plastic_types <- unique(physeq_df4$plastic_type)
plastic_concs <- unique(physeq_df4$plastic_conc)

#DONT RUN TOO LONG TO GENERERATE
#Loop through each plastic type and concentration to create individual plots
#for (ptype in plastic_types) {
  #for (pconc in plastic_concs) {
    # Subset the data for the current plastic type and concentration
    #subset_data <- physeq_df4 %>%
      #filter(plastic_type == ptype, plastic_conc == pconc)
    
    # Create the bar plot for the current plastic type and concentration
    #phyplot_family <- ggplot(subset_data, aes(x = date, y = Abundance, fill = Family)) +
      #geom_bar(stat = "identity", position = "stack", width = 0.8) +  #Stack families for each date
      #facet_wrap(~ Phylum, scales = "free_y") +  #Create a separate facet for each phylum
      #scale_x_discrete(name = "Day Number") +  #Use discrete scale for date
      #scale_fill_manual(values = custom_col294) +  #Use custom color palette with 294 colors
      # Update theme to bold titles, remove grids, and add box around each facet
      #theme(
        #axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
        #axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
        #axis.text.x = element_text(size = 7, angle = 45, hjust = 1),  #Tilt x-axis labels to reduce overlap
        #axis.text.y = element_text(size = 7),  #Increase size of y-axis tick labels
        #legend.position = "none",  #Remove the legend
        #panel.background = element_blank(),  #Remove grey background
        #panel.grid.major = element_blank(),  #Remove major gridlines
        #panel.grid.minor = element_blank(),  #Remove minor gridlines
        #panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  #Box around each plot with correct linewidth argument
      #) +
      #labs(
        #x = "Day Number",
        #y = "Relative Abundance",
        #title = paste("Family-level Abundance for Phylum\nPlastic Type:", ptype, "| Plastic Concentration:", pconc)
      #)  #Add plot title with context of plastic type and concentration
    
    #Print the plot
    #print(phyplot_family)
  #}
#}


#ALL top 3 phyla included, colors are families
#Identify the top 3 phyla for each combination of plastic_type, plastic_conc, and date
top_phyla_df <- physeq_df4 %>%
  group_by(plastic_type, plastic_conc, date, Phylum) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = 'drop') %>%
  arrange(plastic_type, plastic_conc, date, desc(TotalAbundance)) %>%
  group_by(plastic_type, plastic_conc, date) %>%
  slice_max(order_by = TotalAbundance, n = 3) %>%
  ungroup()

#Filter the original data to keep only families within the top 3 phyla for each tank and date
top_families_df <- physeq_df4 %>%
  semi_join(top_phyla_df, by = c("plastic_type", "plastic_conc", "date", "Phylum"))

#Bar plot faceted by tank (plastic type and concentration) with families as stacked bars
phyplot_top_families <- ggplot(top_families_df, aes(x = date, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  #Stack families within each phylum for each date
  facet_grid(plastic_type ~ plastic_conc) +  #Create a separate facet for each tank (plastic type and concentration)
  scale_x_discrete(name = "Sampling Date") +  #Use discrete scale for dates
  scale_fill_manual(values = colorRampPalette(c("red", "blue", "green", "orange", "purple", "yellow", "cyan", "pink", "brown", "grey"))(length(unique(top_families_df$Family)))) +  # Use custom color palette for families
  #Update theme to bold titles, remove grids, and add box around each facet
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5),  #Tilt x-axis labels to reduce overlap
    axis.text.y = element_text(size = 7),  #Increase size of y-axis tick labels
    legend.position = "none",  #Remove the legend
    panel.background = element_blank(),  #Remove grey background
    panel.grid.major = element_blank(),  #Remove major gridlines
    panel.grid.minor = element_blank(),  #Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  #Box around each plot with correct linewidth argument
  ) +
  labs(
    x = "Sampling Date",
    y = "Relative Abundance",
    title = "Family-level Abundance of Top 3 Phyla for Each Tank"
  )  #Add plot title

#Print plot
print(phyplot_top_families)

#Separate figure for each of the 3 phyla
#Identify the top 3 phyla across the entire dataset
top_phyla_overall <- physeq_df4 %>%
  group_by(Phylum) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = 'drop') %>%
  arrange(desc(TotalAbundance)) %>%
  slice_max(order_by = TotalAbundance, n = 3) %>%
  pull(Phylum)

#Create separate data subsets for each of the top 3 phyla
physeq_df_phylum1 <- physeq_df4 %>% filter(Phylum == top_phyla_overall[1])
physeq_df_phylum2 <- physeq_df4 %>% filter(Phylum == top_phyla_overall[2])
physeq_df_phylum3 <- physeq_df4 %>% filter(Phylum == top_phyla_overall[3])

#Create separate plots for each of the top 3 phyla, with legends
#Plot for phylum 1
phyplot_phylum1_Proteobacteria <- ggplot(physeq_df_phylum1, aes(x = date, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  #Stack families for each date
  facet_grid(plastic_type ~ plastic_conc) +  #Create a separate facet for each tank
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for dates
  scale_fill_manual(values = colorRampPalette(c("red", "blue", "green", "orange", "purple", "yellow", "cyan", "pink", "brown", "grey"))(length(unique(physeq_df_phylum1$Family)))) +  # Use custom color palette for families
  #Update theme to bold titles, remove grids, add box around each facet, and add legend back in
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5),  #Center-align x-axis tick labels
    axis.text.y = element_text(size = 7),  #Increase size of y-axis tick labels
    legend.position = "right",  #Move the legend below the plot
    legend.title = element_text(face = "bold", size = 10),  #Bold legend title
    legend.text = element_text(size = 6),  #Customize legend text size
    legend.key.size = unit(0.6, "lines"),  #Adjust spacing between legend items
    legend.key.height = unit(0.4, "lines"),  #Reduce height of legend keys
    panel.background = element_blank(),  #Remove grey background
    panel.grid.major = element_blank(),  #Remove major gridlines
    panel.grid.minor = element_blank(),  #Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  #Box around each plot with correct linewidth argument
  ) +
  guides(fill = guide_legend(ncol = 3)) +  #Adjust the number of columns in the legend for better fit
  labs(
    x = "Day Number",
    y = "Relative Abundance",
    title = paste("Family-level Abundance of Top Phylum:", top_phyla_overall[1])
  )  #Add plot title

#Plot for Phylum 2
phyplot_phylum2_Bacteroidota <- ggplot(physeq_df_phylum2, aes(x = date, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  #Stack families for each date
  facet_grid(plastic_type ~ plastic_conc) +  #Create a separate facet for each tank
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for dates
  scale_fill_manual(values = colorRampPalette(c("red", "blue", "green", "orange", "purple", "yellow", "cyan", "pink", "brown", "grey"))(length(unique(physeq_df_phylum2$Family)))) +  # Use custom color palette for families
  #Update theme to bold titles, remove grids, add box around each facet, and add legend back in
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5),  #Center-align x-axis tick labels
    axis.text.y = element_text(size = 7),  #Increase size of y-axis tick labels
    legend.position = "right",  #Move the legend below the plot
    legend.title = element_text(face = "bold", size = 10),  #Bold legend title
    legend.text = element_text(size = 8),  #Customize legend text size
    legend.key.size = unit(0.6, "lines"),  #Adjust spacing between legend items
    legend.key.height = unit(0.4, "lines"),  #Reduce height of legend keys
    panel.background = element_blank(),  #Remove grey background
    panel.grid.major = element_blank(),  #Remove major gridlines
    panel.grid.minor = element_blank(),  #Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  #Box around each plot with correct linewidth argument
  ) +
  guides(fill = guide_legend(ncol = 2)) +  #Adjust the number of columns in the legend for better fit
  labs(
    x = "Day Number",
    y = "Relative Abundance",
    title = paste("Family-level Abundance of Top Phylum:", top_phyla_overall[2])
  )  #Add plot title

#Plot for Phylum 3
phyplot_phylum3_Cyanobacteria <- ggplot(physeq_df_phylum3, aes(x = date, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  #Stack families for each date
  facet_grid(plastic_type ~ plastic_conc) +  #Create a separate facet for each tank
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for dates
  scale_fill_manual(values = colorRampPalette(c("red", "blue", "green", "orange", "purple", "yellow", "cyan", "pink", "brown", "grey"))(length(unique(physeq_df_phylum3$Family)))) +  # Use custom color palette for families
  #Update theme to bold titles, remove grids, add box around each facet, and add legend back in
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5),  #Center-align x-axis tick labels
    axis.text.y = element_text(size = 7),  #Increase size of y-axis tick labels
    legend.position = "right",  #Move the legend below the plot
    legend.title = element_text(face = "bold", size = 10),  #Bold legend title
    legend.text = element_text(size = 8),  #Customize legend text size
    legend.key.size = unit(0.6, "lines"),  #Adjust spacing between legend items
    legend.key.height = unit(0.4, "lines"),  #Reduce height of legend keys
    panel.background = element_blank(),  #Remove grey background
    panel.grid.major = element_blank(),  #Remove major gridlines
    panel.grid.minor = element_blank(),  #Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  #Box around each plot with correct linewidth argument
  ) +
  guides(fill = guide_legend(ncol = 1)) +  #Adjust the number of columns in the legend for better fit
  labs(
    x = "Day Number",
    y = "Relative Abundance",
    title = paste("Family-level Abundance of Top Phylum:", top_phyla_overall[3])
  )  #Add plot title

#Print the plots
print(phyplot_phylum1_Proteobacteria)
print(phyplot_phylum2_Bacteroidota)
print(phyplot_phylum3_Cyanobacteria)


###################################################
###########PERMANOVA analysis below################
###################################################

#Begin addonis2 analysis - Permutation Multivariate Analysis of Variance (PERMANOVA)
#Convert Phyloseq data to usable formats

meta_df <- as(sample_data(UniFrac_Phylo_Object), "data.frame")
#Make sure factors and numeric variables are specified correctly 
meta_df$plastic_type <- as.factor(meta_df$plastic_type)
meta_df$plastic_conc <- as.numeric(as.character(meta_df$plastic_conc))
meta_df$date <- as.factor(meta_df$date)
unique(meta_df$date)
meta_df$date <- factor(meta_df$date, levels = c("March_7", "March_23", "April_12", "May_3", "June_9"))
#cat("Checking overlap between meta_df and bray_dist...\n")
#print(all(rownames(meta_df) %in% rownames(as.matrix(bray_dist))))


#Run the global PERMANOVAs
#ASV level with Bray-Curtis
set.seed(55)
adonis2_bray_asv <- adonis2(bray_dist ~ plastic_type * plastic_conc * date + chla, data = meta_df, permutations = 1000, by = "terms")
print("ASV Level - Bray-Curtis")
print(adonis2_bray_asv)

#ASV Level with Weighted UniFrac
set.seed(66)
adonis2_unifrac_asv <- adonis2(wdistUni ~ plastic_type * plastic_conc * date + chla, data = meta_df, permutations = 1000, by = "terms")
print("ASV Level - Weighted UniFrac")
print(adonis2_unifrac_asv) #MANUSCRIPT [Extended Data Table 6]

#Family level with Bray-Curtis
set.seed(77)
adonis2_bray_family <- adonis2(bray_dist_family ~ plastic_type * plastic_conc * date + chla, data = meta_df, permutations = 1000, by = "terms")
print("Family Level - Bray-Curtis")
print(adonis2_bray_family)

#Family level with Weighted UniFrac
set.seed(88)
adonis2_unifrac_family <- adonis2(wdistUni_family ~ plastic_type * plastic_conc * date + chla, data = meta_df, permutations = 1000, by = "terms")
print("Family Level - Weighted UniFrac")
print(adonis2_unifrac_family)

#Phylum level with Bray-Curtis
set.seed(99)
adonis2_bray_phylum <- adonis2(bray_dist_phylum ~ plastic_type * plastic_conc * date + chla, data = meta_df, permutations = 1000, by = "terms")
print("Phylum Level - Bray-Curtis")
print(adonis2_bray_phylum)

#Phylum level with Weighted UniFrac
set.seed(11)
adonis2_unifrac_phylum <- adonis2(wdistUni_phylum ~ plastic_type * plastic_conc * date + chla, data = meta_df, permutations = 1000, by = "terms")
print("Phylum Level - Weighted UniFrac")
print(adonis2_unifrac_phylum)

#Compare models
#library(dplyr)
#library(tibble)
#library(tidyr)

#Organize adonis2 models into a list
models <- list(
  ASV_Bray       = adonis2_bray_asv,
  ASV_UniFrac    = adonis2_unifrac_asv,
  Family_Bray    = adonis2_bray_family,
  Family_UniFrac = adonis2_unifrac_family,
  Phylum_Bray    = adonis2_bray_phylum,
  Phylum_UniFrac = adonis2_unifrac_phylum
)

#Function to summarize adonis2 model output
summarize_adonis2 <- function(model, name) {
  aov_df <- as.data.frame(model)
  aov_df$term <- rownames(aov_df)
  
  #Total R2: exclude residual and total
  total_r2 <- sum(aov_df$R2[!(aov_df$term %in% c("Residual", "Total"))])
  
  #Sample size (n) and number of predictors (p)
  n <- sum(aov_df$Df)
  p <- sum(!(aov_df$term %in% c("Residual", "Total")))
  
  #Adjusted R2
  adj_r2 <- 1 - ((1 - total_r2) * (n - 1) / (n - p - 1))
  
  #Significant terms (p <= 0.05)
  significant_terms <- aov_df %>%
    filter(!(term %in% c("Residual", "Total")) & `Pr(>F)` <= 0.05) %>%
    pull(term) %>%
    paste(collapse = ", ")
  
  #Create summary tibble
  tibble(
    Model = name,
    Total_R2 = round(total_r2, 3),
    Adjusted_R2 = round(adj_r2, 3),
    Significant_Terms = ifelse(nchar(significant_terms) > 0, significant_terms, "None")
  )
}

#Create and print the summary table
summary_table <- bind_rows(
  lapply(names(models), function(n) summarize_adonis2(models[[n]], n))
)

print(summary_table)
view(summary_table)
#write.csv(summary_table, file = "Distance_Matrix_comparisions_table_global_PERMANOVAs.csv", row.names = FALSE)


#Convert summary_table to long format for plotting
summary_long <- summary_table %>%
  pivot_longer(cols = c(Total_R2, Adjusted_R2), names_to = "Metric", values_to = "R2")

#Order models for readability
summary_long$Model <- factor(summary_long$Model, levels = summary_table$Model)

#Create barplot to compare visually
ggplot(summary_long, aes(x = Model, y = R2, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("#4B9CD3", "#FFA500")) +
  labs(title = "Total and Adjusted R2 Across PERMANOVA Models",
       x = "Model",
       y = "R² Value",
       fill = "Metric") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#################################################
###Homogeneity of dispersion (PERMDISP) checks###
#################################################

#ASV level - Bray-Curtis
disp_bray_asv <- betadisper(bray_dist, meta_df$plastic_type)
set.seed(42)
anova(disp_bray_asv); permutest(disp_bray_asv); TukeyHSD(disp_bray_asv)

#ASV level - UniFrac
disp_unifrac_asv <- betadisper(wdistUni, meta_df$plastic_type)
set.seed(43)
anova(disp_unifrac_asv); permutest(disp_unifrac_asv); TukeyHSD(disp_unifrac_asv)

#Family level - Bray-Curtis
disp_bray_family <- betadisper(bray_dist_family, meta_df$plastic_type)
set.seed(44)
anova(disp_bray_family); permutest(disp_bray_family); TukeyHSD(disp_bray_family)

#Family level - UniFrac
disp_unifrac_family <- betadisper(wdistUni_family, meta_df$plastic_type)
set.seed(45)
anova(disp_unifrac_family); permutest(disp_unifrac_family); TukeyHSD(disp_unifrac_family)

#Phylum level - Bray-Curtis
disp_bray_phylum <- betadisper(bray_dist_phylum, meta_df$plastic_type)
set.seed(46)
anova(disp_bray_phylum); permutest(disp_bray_phylum); TukeyHSD(disp_bray_phylum)

#Phylum level - UniFrac
disp_unifrac_phylum <- betadisper(wdistUni_phylum, meta_df$plastic_type)
set.seed(47)
anova(disp_unifrac_phylum); permutest(disp_unifrac_phylum); TukeyHSD(disp_unifrac_phylum)



#OK now run for each date separately
#List of dates to run PERMANOVA for each
dates <- c("March_7", "March_23", "April_12", "May_3", "June_9")

#Convert distance matrices to matrix for subsetting
bray_asv_matrix     <- as.matrix(bray_dist)
unifrac_asv_matrix  <- as.matrix(wdistUni)
bray_fam_matrix     <- as.matrix(bray_dist_family)
unifrac_fam_matrix  <- as.matrix(wdistUni_family)
bray_phy_matrix     <- as.matrix(bray_dist_phylum)
unifrac_phy_matrix  <- as.matrix(wdistUni_phylum)


#Set up results storage
results_bray_asv     <- list()
results_unifrac_asv  <- list()
results_bray_family  <- list()
results_unifrac_family <- list()
results_bray_phylum  <- list()
results_unifrac_phylum <- list()


#Loop over each date
for (d in dates) {
  subset_df <- meta_df[meta_df$date == d, ]
  sample_ids <- rownames(subset_df)
  
  #Subset all distance matrices
  subset_bray_asv    <- as.dist(bray_asv_matrix[sample_ids, sample_ids])
  subset_unifrac_asv <- as.dist(unifrac_asv_matrix[sample_ids, sample_ids])
  subset_bray_fam    <- as.dist(bray_fam_matrix[sample_ids, sample_ids])
  subset_unifrac_fam <- as.dist(unifrac_fam_matrix[sample_ids, sample_ids])
  subset_bray_phy    <- as.dist(bray_phy_matrix[sample_ids, sample_ids])
  subset_unifrac_phy <- as.dist(unifrac_phy_matrix[sample_ids, sample_ids])
  
  #ASV level
  set.seed(31)
  res_bray_asv <- adonis2(subset_bray_asv ~ plastic_type * plastic_conc + chla, data = subset_df, permutations = 1000, by = "terms")
  set.seed(32)
  res_unifrac_asv <- adonis2(subset_unifrac_asv ~ plastic_type * plastic_conc + chla, data = subset_df, permutations = 1000, by = "terms")
  
  #Family level
  set.seed(33)
  res_bray_fam <- adonis2(subset_bray_fam ~ plastic_type * plastic_conc + chla, data = subset_df, permutations = 1000, by = "terms")
  set.seed(34)
  res_unifrac_fam <- adonis2(subset_unifrac_fam ~ plastic_type * plastic_conc + chla, data = subset_df, permutations = 1000, by = "terms")
  
  #Phylum level
  set.seed(35)
  res_bray_phy <- adonis2(subset_bray_phy ~ plastic_type * plastic_conc + chla, data = subset_df, permutations = 1000, by = "terms")
  set.seed(36)
  res_unifrac_phy <- adonis2(subset_unifrac_phy ~ plastic_type * plastic_conc + chla, data = subset_df, permutations = 1000, by = "terms")
  
  #Save the results
  results_bray_asv[[d]]       <- res_bray_asv
  results_unifrac_asv[[d]]    <- res_unifrac_asv
  results_bray_family[[d]]    <- res_bray_fam
  results_unifrac_family[[d]] <- res_unifrac_fam
  results_bray_phylum[[d]]    <- res_bray_phy
  results_unifrac_phylum[[d]] <- res_unifrac_phy
  
  #Print the results
  cat(paste0("\n--- ", d, " ---\n"))
  cat("ASV - Bray-Curtis:\n"); print(res_bray_asv)
  cat("ASV - UniFrac:\n"); print(res_unifrac_asv)
  cat("Family - Bray-Curtis:\n"); print(res_bray_fam)
  cat("Family - UniFrac:\n"); print(res_unifrac_fam)
  cat("Phylum - Bray-Curtis:\n"); print(res_bray_phy)
  cat("Phylum - UniFrac:\n"); print(res_unifrac_phy)
}

#Function to print results by taxonomic level and distance across all dates
print_results_by_type <- function(result_list, label) {
  cat(paste0("\n### ", label, " ###\n"))
  for (d in names(result_list)) {
    cat(paste0("\n--- ", d, " ---\n"))
    print(result_list[[d]])
  }
}

#Print results for each taxonomic level and distance metric
print_results_by_type(results_bray_asv,       "ASV - Bray-Curtis")
print_results_by_type(results_unifrac_asv,    "ASV - UniFrac") #MANUSCRIPT [Extended Data Table 7]
print_results_by_type(results_bray_family,    "Family - Bray-Curtis")
print_results_by_type(results_unifrac_family, "Family - UniFrac")
print_results_by_type(results_bray_phylum,    "Phylum - Bray-Curtis")
print_results_by_type(results_unifrac_phylum, "Phylum - UniFrac")

#############################################################
#############################################################
#############################################################

#Ok now Run the PCoA Code for each metric

###########PCoA code##########

#Double check metadata and phyloseq object alignment
all(rownames(meta_df) %in% sample_names(UniFrac_Phylo_Object))
all(rownames(meta_df) %in% sample_names(BioPlastics_family))
all(rownames(meta_df) %in% sample_names(BioPlastics_phylum))

#Add metadata to phyloseq objects
sample_data(UniFrac_Phylo_Object) <- sample_data(meta_df)
sample_data(BioPlastics_family)   <- sample_data(meta_df)
sample_data(BioPlastics_phylum)   <- sample_data(meta_df)

#Subset by date
asv_physeq    <- subset_samples(UniFrac_Phylo_Object, date %in% c("March_7", "March_23", "April_12", "May_3", "June_9"))
family_physeq <- subset_samples(BioPlastics_family,   date %in% c("March_7", "March_23", "April_12", "May_3", "June_9"))
phylum_physeq <- subset_samples(BioPlastics_phylum,   date %in% c("March_7", "March_23", "April_12", "May_3", "June_9"))

#Remove taxa with zero counts
asv_physeq    <- prune_taxa(taxa_sums(asv_physeq) > 0, asv_physeq)
family_physeq <- prune_taxa(taxa_sums(family_physeq) > 0, family_physeq)
phylum_physeq <- prune_taxa(taxa_sums(phylum_physeq) > 0, phylum_physeq)

#Run each ordination (PCoA)
set.seed(12)
pcoa_bray_asv    <- ordinate(physeq = asv_physeq,    method = "PCoA", distance = "bray")
set.seed(13)
pcoa_bray_family <- ordinate(physeq = family_physeq, method = "PCoA", distance = "bray")
set.seed(14)
pcoa_bray_phylum <- ordinate(physeq = phylum_physeq, method = "PCoA", distance = "bray")

set.seed(15)
pcoa_unifrac_asv    <- ordinate(physeq = asv_physeq,    method = "PCoA", distance = "wunifrac")
set.seed(16)
pcoa_unifrac_family <- ordinate(physeq = family_physeq, method = "PCoA", distance = "wunifrac")
set.seed(17)
pcoa_unifrac_phylum <- ordinate(physeq = phylum_physeq, method = "PCoA", distance = "wunifrac")

#Format date labels for better display
meta_df$date <- gsub("_", " ", meta_df$date)

#Define plot aesthetics
custom_shapes <- c("Elastollan" = 16, "TPU181" = 17, "TPUFC2.1" = 18)
custom_colors <- c("March 7" = "goldenrod2", "March 23" = "darkorchid3", 
                   "April 12" = "deeppink3", "May 3" = "darkorange3", "June 9" = "lightblue3")

#Function for plotting
plot_pcoa <- function(physeq, ordination, title) {
  sample_data(physeq) <- sample_data(meta_df)  #Sync metadata
  plot_ordination(physeq, ordination, color = "date", shape = "plastic_type", axes = c(1, 2)) +
    geom_point(size = 4, alpha = 0.8) +
    scale_color_manual(values = custom_colors) +
    scale_shape_manual(values = custom_shapes, labels = c("Elastollan", "TPU 181", "TPU FC2.1")) +
    theme_minimal() +
    labs(title = title, color = "Date", shape = "Plastic Type") +
    theme(
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      legend.box.just = "left",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )
}

#Generate plots
plot_bray_asv    <- plot_pcoa(asv_physeq,    pcoa_bray_asv,    "Bray-Curtis PCoA - ASV Level")
plot_bray_family <- plot_pcoa(family_physeq, pcoa_bray_family, "Bray-Curtis PCoA - Family Level")
plot_bray_phylum <- plot_pcoa(phylum_physeq, pcoa_bray_phylum, "Bray-Curtis PCoA - Phylum Level")

plot_unifrac_asv    <- plot_pcoa(asv_physeq,    pcoa_unifrac_asv,    "Weighted UniFrac PCoA - ASV Level")
plot_unifrac_family <- plot_pcoa(family_physeq, pcoa_unifrac_family, "Weighted UniFrac PCoA - Family Level")
plot_unifrac_phylum <- plot_pcoa(phylum_physeq, pcoa_unifrac_phylum, "Weighted UniFrac PCoA - Phylum Level")

#Print plots
print(plot_bray_asv)
print(plot_bray_family)
print(plot_bray_phylum)

print(plot_unifrac_asv)
print(plot_unifrac_family)
print(plot_unifrac_phylum)

#############################
#############################

#Sticking with UniFrac - ASV level for further visualizations

###OK try to do a three panel by plastic type
#Make sure plastic_type is a character in metadata
meta_df$plastic_type <- as.character(meta_df$plastic_type)
sample_data(UniFrac_Phylo_Object) <- sample_data(meta_df)

#Custom shape mapping for plastic_type
custom_shapes <- c("Elastollan" = 16, "TPU181" = 17, "TPUFC2.1" = 18)

#Custom color mapping for date
custom_colors <- c("March 7" = "goldenrod2", "March 23" = "darkorchid3", 
                   "April 12" = "deeppink3", "May 3" = "darkorange3", "June 9" = "lightblue3")

#Subset the phyloseq object by each plastic_type
elastollan_phylo <- subset_samples(UniFrac_Phylo_Object, plastic_type == "Elastollan")
elastollan_phylo <- prune_taxa(taxa_sums(elastollan_phylo) > 0, elastollan_phylo)

tpu181_phylo <- subset_samples(UniFrac_Phylo_Object, plastic_type == "TPU181")
tpu181_phylo <- prune_taxa(taxa_sums(tpu181_phylo) > 0, tpu181_phylo)

tpufc2_phylo <- subset_samples(UniFrac_Phylo_Object, plastic_type == "TPUFC2.1")
tpufc2_phylo <- prune_taxa(taxa_sums(tpufc2_phylo) > 0, tpufc2_phylo)



#Function to create PCoA plot for a given phyloseq object and plastic type
create_pcoa_plot <- function(phylo_obj, plastic_type_label) {
  set.seed(77)
  #Perform PCoA using weighted UniFrac distance
  pcoa_result <- ordinate(physeq = phylo_obj, method = "PCoA", distance = "wunifrac")
  
  #Create PCoA plot
  pcoa_plot <- plot_ordination(
    physeq = phylo_obj, 
    ordination = pcoa_result, 
    color = "date", 
    shape = "plastic_type",
    axes = c(1, 2)
  ) +
    geom_point(size = 4, alpha = 0.8) +
    scale_color_manual(
      values = custom_colors,
      breaks = c("March 7", "March 23", "April 12", "May 3", "June 9"),
      labels = c("March 7", "March 23", "April 12", "May 3", "June 9")
    ) +
    scale_shape_manual(
      values = custom_shapes,
      breaks = c("Elastollan", "TPU181", "TPUFC2.1"),
      labels = c("Elastollan", "TPU 181", "TPU FC2.1")
    ) +
    theme_minimal() +
    labs(
      color = "Date",
      shape = "Plastic Type",
      title = paste("PCoA for Plastic Type:", plastic_type_label)
    ) +
    theme(
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      legend.box.just = "left",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )
  
  return(pcoa_plot)
}

#Create PCoA plots for each plastic type
pcoa_plot_elastollan <- create_pcoa_plot(elastollan_phylo, "Elastollan")
pcoa_plot_tpu181 <- create_pcoa_plot(tpu181_phylo, "TPU 181")
pcoa_plot_tpufc2 <- create_pcoa_plot(tpufc2_phylo, "TPU FC2.1")

#Display the plots
print(pcoa_plot_elastollan)
print(pcoa_plot_tpu181)
print(pcoa_plot_tpufc2)


###################Continue########

#Define custom colors explicitly for correct mapping
custom_colors <- c(
  "March 7" = "goldenrod2",   #Day 0
  "March 23" = "darkorchid3",  #Day 16
  "April 12" = "deeppink3",  #Day 36
  "May 3" = "darkorange3",     #Day 56
  "June 9" = "lightblue3"     #Day 94
)


#Update the legend labels for the plots
pcoa_plot_elastollan <- pcoa_plot_elastollan +
  scale_color_manual(
    values = custom_colors,                #Make sure colors align with original dates
    breaks = names(custom_colors),         #Match to dates used in the data
    labels = c("0", "16", "36", "56", "94"), #Replace labels with day numbers
    name = "Day Number"                    #Update legend title
  )

pcoa_plot_tpu181 <- pcoa_plot_tpu181 +
  scale_color_manual(
    values = custom_colors,
    breaks = names(custom_colors),
    labels = c("0", "16", "36", "56", "94"),
    name = "Day Number"
  )

pcoa_plot_tpufc2 <- pcoa_plot_tpufc2 +
  scale_color_manual(
    values = custom_colors,
    breaks = names(custom_colors),
    labels = c("0", "16", "36", "56", "94"),
    name = "Day Number"
  )

#Print updated plots
print(pcoa_plot_elastollan)
print(pcoa_plot_tpu181)
print(pcoa_plot_tpufc2)


#Update plots
#Define manual concentration labels directly
#Assuming concentrations are in the original phyloseq metadata as plastic_conc
manual_concentration_labels <- data.frame(
  plastic_conc = c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385),
  plastic_conc_label = as.character(0:9)  #Labels from 0 to 9
)

#Function to update PCoA plot by adding annotations without modifying meta_df
update_pcoa_plot_manual_labels <- function(pcoa_plot, phylo_obj, ordination) {
  #Extract metadata for the given phyloseq object
  local_meta_df <- as(sample_data(phylo_obj), "data.frame")
  local_meta_df$SampleID <- rownames(local_meta_df)  #Add SampleID to use for merging
  
  #Extract ordination coordinates as a data frame
  ordination_df <- as.data.frame(ordination$vectors[, 1:2])
  ordination_df$SampleID <- rownames(ordination_df)  #Add SampleID to use for merging
  colnames(ordination_df)[1:2] <- c("Axis.1", "Axis.2")  #Rename columns for plotting
  
  #Merge ordination data with metadata
  combined_df <- left_join(local_meta_df, ordination_df, by = "SampleID")
  
  #Merge with manual labels for plastic concentration
  combined_df <- left_join(combined_df, manual_concentration_labels, by = "plastic_conc")
  
  #Update the existing plot with annotations (no lines)
  updated_plot <- pcoa_plot +
    geom_point(data = combined_df, aes(x = Axis.1, y = Axis.2), size = 5, alpha = 0.8) + #Reduce size of points
    geom_text(data = combined_df, aes(x = Axis.1, y = Axis.2, label = plastic_conc_label), size = 2, color = "black", vjust = 0.5, hjust = 0.5) + #Annotate with numbers 0-9
    theme(
      plot.title = element_blank(),   #Remove title
      legend.position = "none"        #Remove legends
    )
  
  return(updated_plot)
}

#Update each PCoA plot with manual labels
pcoa_plot_elastollan_updated <- update_pcoa_plot_manual_labels(pcoa_plot_elastollan, elastollan_phylo, ordinate(elastollan_phylo, method = "PCoA", distance = "wunifrac"))
pcoa_plot_tpu181_updated <- update_pcoa_plot_manual_labels(pcoa_plot_tpu181, tpu181_phylo, ordinate(tpu181_phylo, method = "PCoA", distance = "wunifrac"))
pcoa_plot_tpufc2_updated <- update_pcoa_plot_manual_labels(pcoa_plot_tpufc2, tpufc2_phylo, ordinate(tpufc2_phylo, method = "PCoA", distance = "wunifrac"))

#Print updated plots
print(pcoa_plot_elastollan_updated)
print(pcoa_plot_tpu181_updated)
print(pcoa_plot_tpufc2_updated)

update_pcoa_plot_no_legend <- function(pcoa_plot, phylo_obj, ordination) {
  #Extract metadata for the given phyloseq object
  local_meta_df <- as(sample_data(phylo_obj), "data.frame")
  local_meta_df$SampleID <- rownames(local_meta_df)  #Add SampleID to use for merging
  
  #Extract ordination coordinates as a data frame
  ordination_df <- as.data.frame(ordination$vectors[, 1:2])
  ordination_df$SampleID <- rownames(ordination_df)  #Add SampleID to use for merging
  colnames(ordination_df)[1:2] <- c("Axis.1", "Axis.2")  #Rename columns for easy plotting
  
  #Merge ordination data with metadata
  combined_df <- left_join(local_meta_df, ordination_df, by = "SampleID")
  
  #Merge with manual labels for plastic concentration
  combined_df <- left_join(combined_df, manual_concentration_labels, by = "plastic_conc")
  
  #Make sure plastic_conc_label is numeric
  combined_df$plastic_conc_label <- as.numeric(combined_df$plastic_conc_label)
  
  #Update the existing plot with transparency scaled by concentration
  updated_plot <- ggplot(combined_df, aes(x = Axis.1, y = Axis.2)) +
    geom_point(
      aes(
        size = plastic_conc_label,
        alpha = plastic_conc_label / max(plastic_conc_label, na.rm = TRUE),
        color = date,   #Retain coloring by date
        shape = plastic_type  #Retain shapes by plastic type
      ),
      show.legend = FALSE
    ) +
    geom_text(
      aes(label = plastic_conc_label),
      size = 3,
      color = "black",
      vjust = 0.5,
      hjust = 0.5
    ) +
    scale_size_continuous(
      range = c(3, 10),  #Adjust point sizes based on concentration
      guide = "none"     #Suppress size from legend
    ) +
    scale_alpha_continuous(
      range = c(0.3, 1), #Make transparency differences more dramatic
      guide = "none"     #Suppress alpha from legend
    ) +
    scale_color_manual(
      values = custom_colors,
      breaks = names(custom_colors),
      labels = c("0", "16", "36", "56", "94")  #Replace labels with day numbers
    ) +
    scale_shape_manual(
      values = custom_shapes
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",  #Remove all legends
      plot.title = element_blank(),  #Remove title
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )
  
  return(updated_plot)
}

#Rebuild each PCoA plot without any legends and adding transparency
set.seed(55)
pcoa_plot_elastollan_no_legend <- update_pcoa_plot_no_legend(
  pcoa_plot_elastollan,
  elastollan_phylo,
  ordinate(elastollan_phylo, method = "PCoA", distance = "wunifrac")
)

set.seed(56)
pcoa_plot_tpu181_no_legend <- update_pcoa_plot_no_legend(
  pcoa_plot_tpu181,
  tpu181_phylo,
  ordinate(tpu181_phylo, method = "PCoA", distance = "wunifrac")
)

set.seed(57)
pcoa_plot_tpufc2_no_legend <- update_pcoa_plot_no_legend(
  pcoa_plot_tpufc2,
  tpufc2_phylo,
  ordinate(tpufc2_phylo, method = "PCoA", distance = "wunifrac")
)

#Print updated plots
print(pcoa_plot_elastollan_no_legend)
print(pcoa_plot_tpu181_no_legend)
print(pcoa_plot_tpufc2_no_legend)

#####Put them on the same scale#####
#Manually set x and y axis limits
x_limits <- c(-0.25, 0.25)
y_limits <- c(-0.25, 0.25)

#Update each plot with fixed axis limits
pcoa_plot_elastollan_no_legend <- pcoa_plot_elastollan_no_legend +
  xlim(x_limits) +
  ylim(y_limits)

pcoa_plot_tpu181_no_legend <- pcoa_plot_tpu181_no_legend +
  xlim(x_limits) +
  ylim(y_limits)

pcoa_plot_tpufc2_no_legend <- pcoa_plot_tpufc2_no_legend +
  xlim(x_limits) +
  ylim(y_limits)

#Print updated plots
print(pcoa_plot_elastollan_no_legend) #[Supplementary Fig 6]
print(pcoa_plot_tpu181_no_legend) #[Supplementary Fig 6]
print(pcoa_plot_tpufc2_no_legend) #[Supplementary Fig 6]

#####################################################
#OK now lets do each date separately with chla vector
#Define manual concentration labels directly
unique(meta_df$date)

manual_concentration_labels <- data.frame(
  plastic_conc = c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385),
  plastic_conc_label = as.character(0:9)  #Labels from 0 to 9
)

#Custom color mapping for plastic_type
custom_colors <- c("Elastollan" = "darkorchid3", "TPU181" = "deeppink3", "TPUFC2.1" = "darkorange3")
custom_shapes <- c("Elastollan" = 16, "TPU181" = 17, "TPUFC2.1" = 18)

#Subset data by date
march_7_phylo <- subset_samples(UniFrac_Phylo_Object, date == "March 7")
march_23_phylo <- subset_samples(UniFrac_Phylo_Object, date == "March 23")
april_12_phylo <- subset_samples(UniFrac_Phylo_Object, date == "April 12")
may_3_phylo <- subset_samples(UniFrac_Phylo_Object, date == "May 3")
june_9_phylo <- subset_samples(UniFrac_Phylo_Object, date == "June 9")

#Prune taxa with zero counts
march_7_phylo <- prune_taxa(taxa_sums(march_7_phylo) > 0, march_7_phylo)
march_23_phylo <- prune_taxa(taxa_sums(march_23_phylo) > 0, march_23_phylo)
april_12_phylo <- prune_taxa(taxa_sums(april_12_phylo) > 0, april_12_phylo)
may_3_phylo <- prune_taxa(taxa_sums(may_3_phylo) > 0, may_3_phylo)
june_9_phylo <- prune_taxa(taxa_sums(june_9_phylo) > 0, june_9_phylo)

#Function to create PCoA plot with `chla` vector
create_pcoa_plot <- function(phylo_obj, date_label) {
  set.seed(10)
  # Perform PCoA using weighted UniFrac distance
  date_pcoa <- ordinate(physeq = phylo_obj, method = "PCoA", distance = "wunifrac")
  
  #Extract metadata for the given phyloseq object
  local_meta_df <- as(sample_data(phylo_obj), "data.frame")
  local_meta_df$SampleID <- rownames(local_meta_df)  # Add SampleID to use for merging
  
  #Extract ordination coordinates as a data frame from the PCoA object
  ordination_df <- as.data.frame(date_pcoa$vectors[, 1:2])
  ordination_df$SampleID <- rownames(ordination_df)  #Add SampleID to use for merging
  colnames(ordination_df)[1:2] <- c("Axis.1", "Axis.2")  #Rename columns for easy plotting
  
  #Merge ordination data with metadata for the current date
  combined_df <- left_join(local_meta_df, ordination_df, by = "SampleID")
  
  #Merge with manual labels for plastic concentration for the current date
  combined_df <- left_join(combined_df, manual_concentration_labels, by = "plastic_conc")
  
  #Fit the environmental variable (`chla`) to the ordination using `envfit`
  ord_scores <- ordination_df[, c("Axis.1", "Axis.2")]  #Extract the ordination scores for envfit
  envfit_result <- envfit(ord_scores, local_meta_df$chla, permutations = 1000)
  
  #Create PCoA plot for the current date
  pcoa_plot <- ggplot(combined_df, aes(x = Axis.1, y = Axis.2, color = plastic_type)) +
    geom_point(size = 5, alpha = 0.8) +  #Defined point size
    geom_text(aes(label = plastic_conc_label), size = 4, color = "black", vjust = 0.5, hjust = 0.5) +  # Annotations with numbers 0-9
    scale_color_manual(
      values = custom_colors,
      breaks = c("Elastollan", "TPU181", "TPUFC2.1"),
      labels = c("Elastollan", "TPU 181", "TPU FC2.1")
    ) +
    theme_minimal() +
    labs(
      title = paste("PCoA for Date:", date_label),
      color = "Plastic Type"
    ) +
    theme(
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      legend.box.just = "left",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )
  
  #Extract the `chla` vector from `envfit` results and add to the plot for the current date
  chla_arrow <- data.frame(
    Axis.1 = 0, 
    Axis.2 = 0, 
    Axis.1_end = envfit_result$vectors$arrows[1, 1] * 1.5,
    Axis.2_end = envfit_result$vectors$arrows[1, 2] * 1.5
  )
  
  pcoa_plot <- pcoa_plot +
    geom_segment(data = chla_arrow, aes(x = Axis.1, y = Axis.2, xend = Axis.1_end, yend = Axis.2_end),
                 arrow = arrow(length = unit(0.3, "cm")), color = "green4", size = 1) +
    geom_text(data = chla_arrow, aes(x = Axis.1_end, y = Axis.2_end, label = "Chla"), color = "blue", vjust = -0.5)
  
  return(pcoa_plot)
}

#Create PCoA plots for each date
pcoa_plot_march_7 <- create_pcoa_plot(march_7_phylo, "March 7")
pcoa_plot_march_23 <- create_pcoa_plot(march_23_phylo, "March 23")
pcoa_plot_april_12 <- create_pcoa_plot(april_12_phylo, "April 12")
pcoa_plot_may_3 <- create_pcoa_plot(may_3_phylo, "May 3")
pcoa_plot_june_9 <- create_pcoa_plot(june_9_phylo, "June 9")

#Print plots for each date
print(pcoa_plot_march_7)
print(pcoa_plot_march_23)
print(pcoa_plot_april_12)
print(pcoa_plot_may_3)
print(pcoa_plot_june_9)

#Update plots
#Calculate consistent axis limits for all date-specific PCoA plots
all_vectors <- list(
  ordinate(march_7_phylo, method = "PCoA", distance = "wunifrac")$vectors[, 1:2],
  ordinate(march_23_phylo, method = "PCoA", distance = "wunifrac")$vectors[, 1:2],
  ordinate(april_12_phylo, method = "PCoA", distance = "wunifrac")$vectors[, 1:2],
  ordinate(may_3_phylo, method = "PCoA", distance = "wunifrac")$vectors[, 1:2],
  ordinate(june_9_phylo, method = "PCoA", distance = "wunifrac")$vectors[, 1:2]
)

#Find the global minimum and maximum limits for the axes
x_min <- min(sapply(all_vectors, function(vec) min(vec[, 1])))
x_max <- max(sapply(all_vectors, function(vec) max(vec[, 1])))
y_min <- min(sapply(all_vectors, function(vec) min(vec[, 2])))
y_max <- max(sapply(all_vectors, function(vec) max(vec[, 2])))

#Expand axis limits slightly to accommodate the chla vector
expand_factor <- 0.1  #Expands by 10% to make sure enough space for arrows
x_range <- x_max - x_min
y_range <- y_max - y_min

axis_limits <- list(
  xlim = c(x_min - expand_factor * x_range, x_max + expand_factor * x_range),
  ylim = c(y_min - expand_factor * y_range, y_max + expand_factor * y_range)
)


#Calculate consistent axis limits for all date-specific PCoA plots
all_vectors <- list(
  ordinate(march_7_phylo, method = "PCoA", distance = "wunifrac")$vectors[, 1:2],
  ordinate(march_23_phylo, method = "PCoA", distance = "wunifrac")$vectors[, 1:2],
  ordinate(april_12_phylo, method = "PCoA", distance = "wunifrac")$vectors[, 1:2],
  ordinate(may_3_phylo, method = "PCoA", distance = "wunifrac")$vectors[, 1:2],
  ordinate(june_9_phylo, method = "PCoA", distance = "wunifrac")$vectors[, 1:2]
)

#Find the global minimum and maximum limits for the axes
x_min <- min(sapply(all_vectors, function(vec) min(vec[, 1])))
x_max <- max(sapply(all_vectors, function(vec) max(vec[, 1])))
y_min <- min(sapply(all_vectors, function(vec) min(vec[, 2])))
y_max <- max(sapply(all_vectors, function(vec) max(vec[, 2])))

#Expand axis limits slightly to accommodate the chla vector
expand_factor <- 0.1  #Expands by 10% to ensure enough space for arrows
x_range <- x_max - x_min
y_range <- y_max - y_min

axis_limits <- list(
  xlim = c(x_min - expand_factor * x_range, x_max + expand_factor * x_range),
  ylim = c(y_min - expand_factor * y_range, y_max + expand_factor * y_range)
)

#Function to create combined data frame for a given phyloseq object
create_combined_df <- function(phylo_obj) {
  set.seed(11)
  date_pcoa <- ordinate(physeq = phylo_obj, method = "PCoA", distance = "wunifrac")
  local_meta_df <- as(sample_data(phylo_obj), "data.frame")
  local_meta_df$SampleID <- rownames(local_meta_df)
  ordination_df <- as.data.frame(date_pcoa$vectors[, 1:2])
  ordination_df$SampleID <- rownames(ordination_df)
  colnames(ordination_df)[1:2] <- c("Axis.1", "Axis.2")
  combined_df <- dplyr::left_join(local_meta_df, ordination_df, by = "SampleID")
  return(combined_df)
}

#Create combined data frames for each date
combined_df_march_7 <- create_combined_df(march_7_phylo)
combined_df_march_23 <- create_combined_df(march_23_phylo)
combined_df_april_12 <- create_combined_df(april_12_phylo)
combined_df_may_3 <- create_combined_df(may_3_phylo)
combined_df_june_9 <- create_combined_df(june_9_phylo)

#Function to create `chla_arrow_data` for a given ordination and metadata
create_chla_arrow_data <- function(ord_scores, chla_values, axis_limits) {
  set.seed(12)
  envfit_result <- envfit(ord_scores, chla_values, permutations = 1000)
  arrow_data <- data.frame(
    Axis.1 = 0, 
    Axis.2 = 0, 
    Axis.1_end = envfit_result$vectors$arrows[1, 1] * 1.5,
    Axis.2_end = envfit_result$vectors$arrows[1, 2] * 1.5
  )
  while (arrow_data$Axis.1_end < axis_limits$xlim[1] || arrow_data$Axis.1_end > axis_limits$xlim[2] ||
         arrow_data$Axis.2_end < axis_limits$ylim[1] || arrow_data$Axis.2_end > axis_limits$ylim[2]) {
    arrow_data$Axis.1_end <- arrow_data$Axis.1_end * 0.9
    arrow_data$Axis.2_end <- arrow_data$Axis.2_end * 0.9
  }
  return(arrow_data)
}

#Create chla_arrow_data for each date
chla_arrow_data_march_7 <- create_chla_arrow_data(combined_df_march_7[, c("Axis.1", "Axis.2")], combined_df_march_7$chla, axis_limits)
chla_arrow_data_march_23 <- create_chla_arrow_data(combined_df_march_23[, c("Axis.1", "Axis.2")], combined_df_march_23$chla, axis_limits)
chla_arrow_data_april_12 <- create_chla_arrow_data(combined_df_april_12[, c("Axis.1", "Axis.2")], combined_df_april_12$chla, axis_limits)
chla_arrow_data_may_3 <- create_chla_arrow_data(combined_df_may_3[, c("Axis.1", "Axis.2")], combined_df_may_3$chla, axis_limits)
chla_arrow_data_june_9 <- create_chla_arrow_data(combined_df_june_9[, c("Axis.1", "Axis.2")], combined_df_june_9$chla, axis_limits)

#Consistent shape and color assignments
shape_mapping <- c("Elastollan" = 16, "TPU181" = 17, "TPUFC2.1" = 18)
color_mapping <- c("Elastollan" = "darkorchid3", "TPU181" = "deeppink3", "TPUFC2.1" = "darkorange3")

#Function to create PCoA plot for a given combined data frame and arrow data
create_pcoa_plot <- function(combined_df, chla_arrow_data, title) {
  ggplot(combined_df, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(color = plastic_type, shape = plastic_type, alpha = plastic_conc_numeric / 9), size = 5) +
    geom_text(aes(label = plastic_conc_numeric - 1), size = 4, color = "black", vjust = 0.5, hjust = 0.5) +
    scale_shape_manual(values = shape_mapping) +
    scale_color_manual(values = color_mapping) +
    scale_alpha_continuous(range = c(0.3, 1)) +  #Make sure alpha is within a visible range
    theme_minimal() +
    xlim(axis_limits$xlim) + ylim(axis_limits$ylim) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    geom_segment(data = chla_arrow_data, aes(x = Axis.1, y = Axis.2, xend = Axis.1_end, yend = Axis.2_end),
                 arrow = arrow(length = unit(0.3, "cm")), color = "green4", size = 1) +
    geom_text(data = chla_arrow_data, aes(x = Axis.1_end, y = Axis.2_end, label = "Chla"), color = "black", vjust = -0.5) +
    labs(title = title)
}

###################################
# Manually defined labels
manual_concentration_labels <- data.frame(
  plastic_conc = c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385),
  plastic_conc_numeric = 0:9  #Numeric values to be used for plotting size and alpha
)

#Function to map concentration labels to each combined_df
map_conc_labels <- function(df) {
  df <- merge(df, manual_concentration_labels, by = "plastic_conc", all.x = TRUE)
  return(df)
}

#Apply mapping to each combined_df
combined_df_march_7    <- map_conc_labels(combined_df_march_7)
combined_df_march_23   <- map_conc_labels(combined_df_march_23)
combined_df_april_12   <- map_conc_labels(combined_df_april_12)
combined_df_may_3      <- map_conc_labels(combined_df_may_3)
combined_df_june_9     <- map_conc_labels(combined_df_june_9)

###############################
#Create PCoA plots for each date
pcoa_plot_march_7 <- create_pcoa_plot(combined_df_march_7, chla_arrow_data_march_7, "March 7 PCoA Plot")
pcoa_plot_march_23 <- create_pcoa_plot(combined_df_march_23, chla_arrow_data_march_23, "March 23 PCoA Plot")
pcoa_plot_april_12 <- create_pcoa_plot(combined_df_april_12, chla_arrow_data_april_12, "April 12 PCoA Plot")
pcoa_plot_may_3 <- create_pcoa_plot(combined_df_may_3, chla_arrow_data_may_3, "May 3 PCoA Plot")
pcoa_plot_june_9 <- create_pcoa_plot(combined_df_june_9, chla_arrow_data_june_9, "June 9 PCoA Plot")

#Print updated plots for each date
print(pcoa_plot_march_7)
print(pcoa_plot_march_23)
print(pcoa_plot_april_12)
print(pcoa_plot_may_3)
print(pcoa_plot_june_9)

#Remove axis titles and figure titles from each plot and display them
pcoa_plot_march_7_clean <- pcoa_plot_march_7 + theme(axis.title = element_blank(), plot.title = element_blank())
pcoa_plot_march_23_clean <- pcoa_plot_march_23 + theme(axis.title = element_blank(), plot.title = element_blank())
pcoa_plot_april_12_clean <- pcoa_plot_april_12 + theme(axis.title = element_blank(), plot.title = element_blank())
pcoa_plot_may_3_clean <- pcoa_plot_may_3 + theme(axis.title = element_blank(), plot.title = element_blank())
pcoa_plot_june_9_clean <- pcoa_plot_june_9 + theme(axis.title = element_blank(), plot.title = element_blank())

#Print updated plots for each date
print(pcoa_plot_march_7_clean)
print(pcoa_plot_march_23_clean)
print(pcoa_plot_april_12_clean)
print(pcoa_plot_may_3_clean)
print(pcoa_plot_june_9_clean)

#Now make icons larger as well with conc
#Create chla_arrow_data for each date
chla_arrow_data_march_7_icon <- create_chla_arrow_data(combined_df_march_7[, c("Axis.1", "Axis.2")], combined_df_march_7$chla, axis_limits)
chla_arrow_data_march_237_icon <- create_chla_arrow_data(combined_df_march_23[, c("Axis.1", "Axis.2")], combined_df_march_23$chla, axis_limits)
chla_arrow_data_april_127_icon <- create_chla_arrow_data(combined_df_april_12[, c("Axis.1", "Axis.2")], combined_df_april_12$chla, axis_limits)
chla_arrow_data_may_37_icon <- create_chla_arrow_data(combined_df_may_3[, c("Axis.1", "Axis.2")], combined_df_may_3$chla, axis_limits)
chla_arrow_data_june_97_icon <- create_chla_arrow_data(combined_df_june_9[, c("Axis.1", "Axis.2")], combined_df_june_9$chla, axis_limits)

#Consistent shape and color assignments
shape_mapping <- c("Elastollan" = 16, "TPU181" = 17, "TPUFC2.1" = 18)
color_mapping <- c("Elastollan" = "darkorchid3", "TPU181" = "deeppink3", "TPUFC2.1" = "darkorange3")

#Function to create PCoA plot for a given combined data frame and arrow data
create_pcoa_plot <- function(combined_df, chla_arrow_data, title) {
  ggplot(combined_df, aes(x = Axis.1, y = Axis.2)) +
    geom_point(
      aes(
        color = plastic_type,
        shape = plastic_type,
        alpha = plastic_conc_numeric / 9,
        size = plastic_conc_numeric
      )
    ) +
    geom_text(aes(label = plastic_conc_numeric), size = 2, color = "black", vjust = 0.5, hjust = 0.5) +
    scale_shape_manual(values = shape_mapping) +
    scale_color_manual(values = color_mapping) +
    scale_alpha_continuous(range = c(0.3, 1)) +  # Ensure alpha is within a visible range
    scale_size_continuous(range = c(3, 10)) +    # Icons increase in size with concentration
    theme_minimal() +
    xlim(axis_limits$xlim) + ylim(axis_limits$ylim) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    geom_segment(
      data = chla_arrow_data,
      aes(x = Axis.1, y = Axis.2, xend = Axis.1_end, yend = Axis.2_end),
      arrow = arrow(length = unit(0.3, "cm")), color = "green4", size = 1
    ) +
    geom_text(
      data = chla_arrow_data,
      aes(x = Axis.1_end, y = Axis.2_end, label = "Chla"),
      color = "black", vjust = -0.5
    ) +
    labs(title = title)
}

#Create PCoA plots for each date
pcoa_plot_march_7_icon <- create_pcoa_plot(combined_df_march_7, chla_arrow_data_march_7, "March 7 PCoA Plot")
pcoa_plot_march_23_icon <- create_pcoa_plot(combined_df_march_23, chla_arrow_data_march_23, "March 23 PCoA Plot")
pcoa_plot_april_12_icon <- create_pcoa_plot(combined_df_april_12, chla_arrow_data_april_12, "April 12 PCoA Plot")
pcoa_plot_may_3_icon <- create_pcoa_plot(combined_df_may_3, chla_arrow_data_may_3, "May 3 PCoA Plot")
pcoa_plot_june_9_icon <- create_pcoa_plot(combined_df_june_9, chla_arrow_data_june_9, "June 9 PCoA Plot")

#Print updated plots for each date
print(pcoa_plot_march_7_icon)
print(pcoa_plot_march_23_icon)
print(pcoa_plot_april_12_icon)
print(pcoa_plot_may_3_icon)
print(pcoa_plot_june_9_icon)

#Remove axis titles and figure titles from each plot
pcoa_plot_march_7_icon <- pcoa_plot_march_7_icon + theme(axis.title = element_blank(), plot.title = element_blank())
pcoa_plot_march_23_icon <- pcoa_plot_march_23_icon + theme(axis.title = element_blank(), plot.title = element_blank())
pcoa_plot_april_12_icon <- pcoa_plot_april_12_icon + theme(axis.title = element_blank(), plot.title = element_blank())
pcoa_plot_may_3_icon <- pcoa_plot_may_3_icon + theme(axis.title = element_blank(), plot.title = element_blank())
pcoa_plot_june_9_icon <- pcoa_plot_june_9_icon + theme(axis.title = element_blank(), plot.title = element_blank())

#Print updated plots for each date
print(pcoa_plot_march_7_icon) #MANUSCRIPT [Fig. 3]
print(pcoa_plot_march_23_icon) #MANUSCRIPT [Fig. 3]
print(pcoa_plot_april_12_icon) #MANUSCRIPT [Fig. 3]
print(pcoa_plot_may_3_icon) #MANUSCRIPT [Fig. 3]
print(pcoa_plot_june_9_icon) #MANUSCRIPT [Fig. 3]


################################
###OK now start doing biplots###
################################

#Attach taxonomy table to phyloseq object
UniFrac_Phylo_Object <- merge_phyloseq(UniFrac_Phylo_Object, tax_table(as.matrix(taxa3.m)))

######## MARCH 7 ##########
#Subset phyloseq object for March 7  
march_7_phylo <- subset_samples(UniFrac_Phylo_Object, date == "March 7")

#Run PCoA ordination
set.seed(77)
ordination_march_7 <- ordinate(march_7_phylo, method = "PCoA", distance = "wunifrac")

#Extract ordination scores  
ordination_scores_march_7 <- as.data.frame(ordination_march_7$vectors[, 1:2])
ordination_scores_march_7$SampleID <- rownames(ordination_scores_march_7)

#Extract metadata and merge with ordination scores  
meta_march_7 <- as(sample_data(march_7_phylo), "data.frame")  
meta_march_7$SampleID <- rownames(meta_march_7)
ordination_metadata_march_7 <- left_join(meta_march_7, ordination_scores_march_7, by = "SampleID")

#Convert metadata to data frame  
ordination_metadata_march_7 <- as.data.frame(ordination_metadata_march_7)

#Define plastic groups  
ordination_metadata_march_7$plastic_group <- ifelse(
  ordination_metadata_march_7$plastic_type %in% c("TPU181", "TPUFC2.1"), 
  "Bioplastics", 
  "Elastollan"
)
ordination_metadata_march_7$plastic_group <- factor(ordination_metadata_march_7$plastic_group, 
                                                    levels = c("Elastollan", "Bioplastics"))

###############################################
################ ENVFIT #######################
###############################################

#Process OTU and taxonomy at Family level
otu <- otu_table(march_7_phylo)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
otu_table_march_7_family <- as.data.frame(as(otu, "matrix"))

tax_table_march_7_family <- as.data.frame(tax_table(march_7_phylo))
tax_table_march_7_family[] <- lapply(tax_table_march_7_family, trimws)

common_asvs <- intersect(rownames(tax_table_march_7_family), colnames(otu_table_march_7_family))
tax_table_march_7_family <- tax_table_march_7_family[common_asvs, , drop = FALSE]
otu_table_march_7_family <- otu_table_march_7_family[, common_asvs, drop = FALSE]

otu_table_march_7_family <- as.data.frame(t(otu_table_march_7_family))
otu_table_march_7_family$Family <- tax_table_march_7_family$Family
otu_table_march_7_family <- otu_table_march_7_family[!is.na(otu_table_march_7_family$Family) & 
                                                       otu_table_march_7_family$Family != "", ]
otu_table_march_7_family <- aggregate(. ~ Family, data = otu_table_march_7_family, FUN = sum)

#Set family as column names and transpose to sample × taxa for envfit
rownames(otu_table_march_7_family) <- otu_table_march_7_family$Family
otu_table_march_7_family$Family <- NULL
otu_table_march_7_family <- as.data.frame(t(otu_table_march_7_family))

#Process OTU and taxonomy at genus level
otu <- otu_table(march_7_phylo)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
otu_table_march_7_genus <- as.data.frame(as(otu, "matrix"))

tax_table_march_7_genus <- as.data.frame(tax_table(march_7_phylo))
tax_table_march_7_genus[] <- lapply(tax_table_march_7_genus, trimws)

common_asvs <- intersect(rownames(tax_table_march_7_genus), colnames(otu_table_march_7_genus))
tax_table_march_7_genus <- tax_table_march_7_genus[common_asvs, , drop = FALSE]
otu_table_march_7_genus <- otu_table_march_7_genus[, common_asvs, drop = FALSE]

otu_table_march_7_genus <- as.data.frame(t(otu_table_march_7_genus))
otu_table_march_7_genus$Genus <- tax_table_march_7_genus$Genus
otu_table_march_7_genus <- otu_table_march_7_genus[!is.na(otu_table_march_7_genus$Genus) & 
                                                     otu_table_march_7_genus$Genus != "", ]
otu_table_march_7_genus <- aggregate(. ~ Genus, data = otu_table_march_7_genus, FUN = sum)

rownames(otu_table_march_7_genus) <- otu_table_march_7_genus$Genus
otu_table_march_7_genus$Genus <- NULL
otu_table_march_7_genus <- as.data.frame(t(otu_table_march_7_genus))

#Calculate mean abundance
mean_abundance_families <- colMeans(otu_table_march_7_family)
mean_abundance_genera <- colMeans(otu_table_march_7_genus)

#Align samples for envfit
shared_samples_family <- intersect(rownames(ordination_scores_march_7), rownames(otu_table_march_7_family))
ordination_envfit_family <- ordination_scores_march_7[shared_samples_family, 1:2]
otu_envfit_family <- otu_table_march_7_family[shared_samples_family, ]

shared_samples_genus <- intersect(rownames(ordination_scores_march_7), rownames(otu_table_march_7_genus))
ordination_envfit_genus <- ordination_scores_march_7[shared_samples_genus, 1:2]
otu_envfit_genus <- otu_table_march_7_genus[shared_samples_genus, ]

#Remove non-variable families
otu_envfit_family <- otu_envfit_family[, apply(otu_envfit_family, 2, function(x) var(x, na.rm = TRUE) > 0)]

#Remove non-variable genera
otu_envfit_genus <- otu_envfit_genus[, apply(otu_envfit_genus, 2, function(x) var(x, na.rm = TRUE) > 0)]


#Run envfit and extract significant taxa (family)
set.seed(79)
envfit_march_7_family <- envfit(ordination_envfit_family, otu_envfit_family, permutations = 999)
sig_families_march_7 <- as.data.frame(envfit_march_7_family$vectors$arrows)
sig_families_march_7$p_value <- envfit_march_7_family$vectors$pvals
sig_families_march_7 <- sig_families_march_7[sig_families_march_7$p_value < 0.05, ]
sig_families_march_7$Family <- rownames(sig_families_march_7)
sig_families_march_7$mean_abundance <- mean_abundance_families[rownames(sig_families_march_7)]
sig_families_march_7 <- sig_families_march_7[order(-sig_families_march_7$mean_abundance), ]


#Run envfit and extract significant taxa (genus)
set.seed(67)
envfit_march_7_genus <- envfit(ordination_envfit_genus, otu_envfit_genus, permutations = 999)
sig_genera_march_7 <- as.data.frame(envfit_march_7_genus$vectors$arrows)
sig_genera_march_7$p_value <- envfit_march_7_genus$vectors$pvals
sig_genera_march_7 <- sig_genera_march_7[sig_genera_march_7$p_value < 0.05, ]
sig_genera_march_7$Genus <- rownames(sig_genera_march_7)
sig_genera_march_7$mean_abundance <- mean_abundance_genera[rownames(sig_genera_march_7)]
sig_genera_march_7 <- sig_genera_march_7[order(-sig_genera_march_7$mean_abundance), ]


#Select top 10 families and genera for biplot
top_10_families_march_7 <- head(sig_families_march_7, 10)
vectors_df_march_7_family <- data.frame(
  Axis.1 = top_10_families_march_7$Axis.1 * 0.2,  
  Axis.2 = top_10_families_march_7$Axis.2 * 0.2,
  Family = rownames(top_10_families_march_7)
)

top_10_genera_march_7 <- head(sig_genera_march_7, 10)
vectors_df_march_7_genus <- data.frame(
  Axis.1 = top_10_genera_march_7$Axis.1 * 0.2,  
  Axis.2 = top_10_genera_march_7$Axis.2 * 0.2,
  Genus = rownames(top_10_genera_march_7)
)

#Family biplot
biplot_march_7_family <- ggplot() +
  geom_point(data = ordination_metadata_march_7, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_march_7_family, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_march_7_family, aes(x = Axis.1, y = Axis.2, label = Family), 
                  color = "black", size = 3, max.overlaps = 15) +
  scale_color_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  labs(color = "Plastic type") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_blank(),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

#Genus biplot
biplot_march_7_genus <- ggplot() +
  geom_point(data = ordination_metadata_march_7, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_march_7_genus, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_march_7_genus, aes(x = Axis.1, y = Axis.2, label = Genus), 
                  color = "black", size = 3, max.overlaps = 15) +
  scale_color_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  labs(color = "Plastic type") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_blank(),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

#Print both plots
print(biplot_march_7_family)
print(biplot_march_7_genus)


######## MARCH 23 ##########
#Subset phyloseq object for March 23  
march_23_phylo <- subset_samples(UniFrac_Phylo_Object, date == "March 23")

#Run PCoA ordination
set.seed(44)
ordination_march_23 <- ordinate(march_23_phylo, method = "PCoA", distance = "wunifrac")

#Extract ordination scores  
ordination_scores_march_23 <- as.data.frame(ordination_march_23$vectors[, 1:2])
ordination_scores_march_23$SampleID <- rownames(ordination_scores_march_23)

#Extract metadata and merge with ordination scores  
meta_march_23 <- as(sample_data(march_23_phylo), "data.frame")  
meta_march_23$SampleID <- rownames(meta_march_23)
ordination_metadata_march_23 <- left_join(meta_march_23, ordination_scores_march_23, by = "SampleID")

#Convert metadata to data frame  
ordination_metadata_march_23 <- as.data.frame(ordination_metadata_march_23)

#Define plastic groups  
ordination_metadata_march_23$plastic_group <- ifelse(
  ordination_metadata_march_23$plastic_type %in% c("TPU181", "TPUFC2.1"), 
  "Bioplastics", 
  "Elastollan"
)
ordination_metadata_march_23$plastic_group <- factor(ordination_metadata_march_23$plastic_group, 
                                                     levels = c("Elastollan", "Bioplastics"))


#Process OTU and taxonomy at family level
otu <- otu_table(march_23_phylo)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
otu_table_march_23_family <- as.data.frame(as(otu, "matrix"))

tax_table_march_23_family <- as.data.frame(tax_table(march_23_phylo))
tax_table_march_23_family[] <- lapply(tax_table_march_23_family, trimws)

common_asvs <- intersect(rownames(tax_table_march_23_family), colnames(otu_table_march_23_family))
tax_table_march_23_family <- tax_table_march_23_family[common_asvs, , drop = FALSE]
otu_table_march_23_family <- otu_table_march_23_family[, common_asvs, drop = FALSE]

otu_table_march_23_family <- as.data.frame(t(otu_table_march_23_family))
otu_table_march_23_family$Family <- tax_table_march_23_family$Family
otu_table_march_23_family <- otu_table_march_23_family[!is.na(otu_table_march_23_family$Family) & 
                                                       otu_table_march_23_family$Family != "", ]
otu_table_march_23_family <- aggregate(. ~ Family, data = otu_table_march_23_family, FUN = sum)

#Set family as column names and transpose to sample × taxa for envfit
rownames(otu_table_march_23_family) <- otu_table_march_23_family$Family
otu_table_march_23_family$Family <- NULL
otu_table_march_23_family <- as.data.frame(t(otu_table_march_23_family))

#Process OTU and taxonomy at genus level
otu <- otu_table(march_23_phylo)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
otu_table_march_23_genus <- as.data.frame(as(otu, "matrix"))

tax_table_march_23_genus <- as.data.frame(tax_table(march_23_phylo))
tax_table_march_23_genus[] <- lapply(tax_table_march_23_genus, trimws)

common_asvs <- intersect(rownames(tax_table_march_23_genus), colnames(otu_table_march_23_genus))
tax_table_march_23_genus <- tax_table_march_23_genus[common_asvs, , drop = FALSE]
otu_table_march_23_genus <- otu_table_march_23_genus[, common_asvs, drop = FALSE]

otu_table_march_23_genus <- as.data.frame(t(otu_table_march_23_genus))
otu_table_march_23_genus$Genus <- tax_table_march_23_genus$Genus
otu_table_march_23_genus <- otu_table_march_23_genus[!is.na(otu_table_march_23_genus$Genus) & 
                                                     otu_table_march_23_genus$Genus != "", ]
otu_table_march_23_genus <- aggregate(. ~ Genus, data = otu_table_march_23_genus, FUN = sum)

rownames(otu_table_march_23_genus) <- otu_table_march_23_genus$Genus
otu_table_march_23_genus$Genus <- NULL
otu_table_march_23_genus <- as.data.frame(t(otu_table_march_23_genus))

#Calculate mean abundance
mean_abundance_families <- colMeans(otu_table_march_23_family)
mean_abundance_genera <- colMeans(otu_table_march_23_genus)

#Align samples for envfit
shared_samples_family <- intersect(rownames(ordination_scores_march_23), rownames(otu_table_march_23_family))
ordination_envfit_family <- ordination_scores_march_23[shared_samples_family, 1:2]
otu_envfit_family <- otu_table_march_23_family[shared_samples_family, ]

shared_samples_genus <- intersect(rownames(ordination_scores_march_23), rownames(otu_table_march_23_genus))
ordination_envfit_genus <- ordination_scores_march_23[shared_samples_genus, 1:2]
otu_envfit_genus <- otu_table_march_23_genus[shared_samples_genus, ]

#Remove non-variable families
otu_envfit_family <- otu_envfit_family[, apply(otu_envfit_family, 2, function(x) var(x, na.rm = TRUE) > 0)]

#Remove non-variable genera
otu_envfit_genus <- otu_envfit_genus[, apply(otu_envfit_genus, 2, function(x) var(x, na.rm = TRUE) > 0)]


#Run envfit and extract significant taxa (Family)
set.seed(55)
envfit_march_23_family <- envfit(ordination_envfit_family, otu_envfit_family, permutations = 999)
sig_families_march_23 <- as.data.frame(envfit_march_23_family$vectors$arrows)
sig_families_march_23$p_value <- envfit_march_23_family$vectors$pvals
sig_families_march_23 <- sig_families_march_23[sig_families_march_23$p_value < 0.05, ]
sig_families_march_23$Family <- rownames(sig_families_march_23)
sig_families_march_23$mean_abundance <- mean_abundance_families[rownames(sig_families_march_23)]
sig_families_march_23 <- sig_families_march_23[order(-sig_families_march_23$mean_abundance), ]


#Run envfit and extract significant taxa (Genus)
set.seed(66)
envfit_march_23_genus <- envfit(ordination_envfit_genus, otu_envfit_genus, permutations = 999)
sig_genera_march_23 <- as.data.frame(envfit_march_23_genus$vectors$arrows)
sig_genera_march_23$p_value <- envfit_march_23_genus$vectors$pvals
sig_genera_march_23 <- sig_genera_march_23[sig_genera_march_23$p_value < 0.05, ]
sig_genera_march_23$Genus <- rownames(sig_genera_march_23)
sig_genera_march_23$mean_abundance <- mean_abundance_genera[rownames(sig_genera_march_23)]
sig_genera_march_23 <- sig_genera_march_23[order(-sig_genera_march_23$mean_abundance), ]


#Select top 10 families and genera for biplot
top_10_families_march_23 <- head(sig_families_march_23, 10)
vectors_df_march_23_family <- data.frame(
  Axis.1 = top_10_families_march_23$Axis.1 * 0.2,  
  Axis.2 = top_10_families_march_23$Axis.2 * 0.2,
  Family = rownames(top_10_families_march_23)
)

top_10_genera_march_23 <- head(sig_genera_march_23, 10)
vectors_df_march_23_genus <- data.frame(
  Axis.1 = top_10_genera_march_23$Axis.1 * 0.2,  
  Axis.2 = top_10_genera_march_23$Axis.2 * 0.2,
  Genus = rownames(top_10_genera_march_23)
)

#Family biplot
biplot_march_23_family <- ggplot() +
  geom_point(data = ordination_metadata_march_23, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_march_23_family, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_march_23_family, aes(x = Axis.1, y = Axis.2, label = Family), 
                  color = "black", size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  labs(color = "Plastic type") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_blank(),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

#Genus biplot
biplot_march_23_genus <- ggplot() +
  geom_point(data = ordination_metadata_march_23, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_march_23_genus, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_march_23_genus, aes(x = Axis.1, y = Axis.2, label = Genus), 
                  color = "black", size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  labs(color = "Plastic type") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.title = element_blank(),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

#Print both plots
print(biplot_march_23_family)
print(biplot_march_23_genus)

######## APRIL 12 ##########
#Subset phyloseq object for April 12  
april_12_phylo <- subset_samples(UniFrac_Phylo_Object, date == "April 12")

#Run PCoA ordination
set.seed(99)
ordination_april_12 <- ordinate(april_12_phylo, method = "PCoA", distance = "wunifrac")

#Extract ordination scores
set.seed(88)
ordination_scores_april_12 <- as.data.frame(ordination_april_12$vectors[, 1:2])
ordination_scores_april_12$SampleID <- rownames(ordination_scores_april_12)

#Extract metadata and merge with ordination scores  
meta_april_12 <- as(sample_data(april_12_phylo), "data.frame")  
meta_april_12$SampleID <- rownames(meta_april_12)
ordination_metadata_april_12 <- left_join(meta_april_12, ordination_scores_april_12, by = "SampleID")

#Convert metadata to data frame  
ordination_metadata_april_12 <- as.data.frame(ordination_metadata_april_12)

#Define plastic groups  
ordination_metadata_april_12$plastic_group <- ifelse(
  ordination_metadata_april_12$plastic_type %in% c("TPU181", "TPUFC2.1"), 
  "Bioplastics", 
  "Elastollan"
)
ordination_metadata_april_12$plastic_group <- factor(ordination_metadata_april_12$plastic_group, 
                                                     levels = c("Elastollan", "Bioplastics"))

#Process OTU and taxonomy at family level
otu <- otu_table(april_12_phylo)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
otu_table_april_12_family <- as.data.frame(as(otu, "matrix"))

tax_table_april_12_family <- as.data.frame(tax_table(april_12_phylo))
tax_table_april_12_family[] <- lapply(tax_table_april_12_family, trimws)

common_asvs <- intersect(rownames(tax_table_april_12_family), colnames(otu_table_april_12_family))
tax_table_april_12_family <- tax_table_april_12_family[common_asvs, , drop = FALSE]
otu_table_april_12_family <- otu_table_april_12_family[, common_asvs, drop = FALSE]

otu_table_april_12_family <- as.data.frame(t(otu_table_april_12_family))
otu_table_april_12_family$Family <- tax_table_april_12_family$Family
otu_table_april_12_family <- otu_table_april_12_family[!is.na(otu_table_april_12_family$Family) & 
                                                         otu_table_april_12_family$Family != "", ]
otu_table_april_12_family <- aggregate(. ~ Family, data = otu_table_april_12_family, FUN = sum)

#Set family as column names and transpose to sample × taxa for envfit
rownames(otu_table_april_12_family) <- otu_table_april_12_family$Family
otu_table_april_12_family$Family <- NULL
otu_table_april_12_family <- as.data.frame(t(otu_table_april_12_family))

#Process OTU and taxonomy at genus level
otu <- otu_table(april_12_phylo)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
otu_table_april_12_genus <- as.data.frame(as(otu, "matrix"))

tax_table_april_12_genus <- as.data.frame(tax_table(april_12_phylo))
tax_table_april_12_genus[] <- lapply(tax_table_april_12_genus, trimws)

common_asvs <- intersect(rownames(tax_table_april_12_genus), colnames(otu_table_april_12_genus))
tax_table_april_12_genus <- tax_table_april_12_genus[common_asvs, , drop = FALSE]
otu_table_april_12_genus <- otu_table_april_12_genus[, common_asvs, drop = FALSE]

otu_table_april_12_genus <- as.data.frame(t(otu_table_april_12_genus))
otu_table_april_12_genus$Genus <- tax_table_april_12_genus$Genus
otu_table_april_12_genus <- otu_table_april_12_genus[!is.na(otu_table_april_12_genus$Genus) & 
                                                       otu_table_april_12_genus$Genus != "", ]
otu_table_april_12_genus <- aggregate(. ~ Genus, data = otu_table_april_12_genus, FUN = sum)

rownames(otu_table_april_12_genus) <- otu_table_april_12_genus$Genus
otu_table_april_12_genus$Genus <- NULL
otu_table_april_12_genus <- as.data.frame(t(otu_table_april_12_genus))

#Calculate mean abundance
mean_abundance_families <- colMeans(otu_table_april_12_family)
mean_abundance_genera <- colMeans(otu_table_april_12_genus)

#Align samples for envfit
shared_samples_family <- intersect(rownames(ordination_scores_april_12), rownames(otu_table_april_12_family))
ordination_envfit_family <- ordination_scores_april_12[shared_samples_family, 1:2]
otu_envfit_family <- otu_table_april_12_family[shared_samples_family, ]

shared_samples_genus <- intersect(rownames(ordination_scores_april_12), rownames(otu_table_april_12_genus))
ordination_envfit_genus <- ordination_scores_april_12[shared_samples_genus, 1:2]
otu_envfit_genus <- otu_table_april_12_genus[shared_samples_genus, ]

#Remove non-variable families
otu_envfit_family <- otu_envfit_family[, apply(otu_envfit_family, 2, function(x) var(x, na.rm = TRUE) > 0)]

#Remove non-variable genera
otu_envfit_genus <- otu_envfit_genus[, apply(otu_envfit_genus, 2, function(x) var(x, na.rm = TRUE) > 0)]


#Run envfit and extract significant taxa (Family)
set.seed(33)
envfit_april_12_family <- envfit(ordination_envfit_family, otu_envfit_family, permutations = 999)
sig_families_april_12 <- as.data.frame(envfit_april_12_family$vectors$arrows)
sig_families_april_12$p_value <- envfit_april_12_family$vectors$pvals
sig_families_april_12 <- sig_families_april_12[sig_families_april_12$p_value < 0.05, ]
sig_families_april_12$Family <- rownames(sig_families_april_12)
sig_families_april_12$mean_abundance <- mean_abundance_families[rownames(sig_families_april_12)]
sig_families_april_12 <- sig_families_april_12[order(-sig_families_april_12$mean_abundance), ]
print(sig_families_april_12)

#Run envfit and extract significant taxa (Genus)
set.seed(34)
envfit_april_12_genus <- envfit(ordination_envfit_genus, otu_envfit_genus, permutations = 999)
sig_genera_april_12 <- as.data.frame(envfit_april_12_genus$vectors$arrows)
sig_genera_april_12$p_value <- envfit_april_12_genus$vectors$pvals
sig_genera_april_12 <- sig_genera_april_12[sig_genera_april_12$p_value < 0.05, ]
sig_genera_april_12$Genus <- rownames(sig_genera_april_12)
sig_genera_april_12$mean_abundance <- mean_abundance_genera[rownames(sig_genera_april_12)]
sig_genera_april_12 <- sig_genera_april_12[order(-sig_genera_april_12$mean_abundance), ]


#Select top 10 families and genera for biplot
top_10_families_april_12 <- head(sig_families_april_12, 10)
vectors_df_april_12_family <- data.frame(
  Axis.1 = top_10_families_april_12$Axis.1 * 0.2,  
  Axis.2 = top_10_families_april_12$Axis.2 * 0.2,
  Family = rownames(top_10_families_april_12)
)

top_10_families_april_12$Family

top_10_genera_april_12 <- head(sig_genera_april_12, 10)
vectors_df_april_12_genus <- data.frame(
  Axis.1 = top_10_genera_april_12$Axis.1 * 0.2,  
  Axis.2 = top_10_genera_april_12$Axis.2 * 0.2,
  Genus = rownames(top_10_genera_april_12)
)

#Family biplot
biplot_april_12_family <- ggplot() +
  geom_point(data = ordination_metadata_april_12, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_april_12_family, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_april_12_family, aes(x = Axis.1, y = Axis.2, label = Family), 
                  color = "black", size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  labs(color = "Plastic type") +
  theme_minimal()

#Genus biplot
biplot_april_12_genus <- ggplot() +
  geom_point(data = ordination_metadata_april_12, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_april_12_genus, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_april_12_genus, aes(x = Axis.1, y = Axis.2, label = Genus), 
                  color = "black", size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  labs(color = "Plastic type") +
  theme_minimal()

#Print both plots
print(biplot_april_12_family)
print(biplot_april_12_genus)


######## MAY 3 ##########
#Subset phyloseq object for May 3  
may_3_phylo <- subset_samples(UniFrac_Phylo_Object, date == "May 3")

#Run PCoA ordination
set.seed(71)
ordination_may_3 <- ordinate(may_3_phylo, method = "PCoA", distance = "wunifrac")

#Extract ordination scores  
ordination_scores_may_3 <- as.data.frame(ordination_may_3$vectors[, 1:2])
ordination_scores_may_3$SampleID <- rownames(ordination_scores_may_3)

#Extract metadata and merge with ordination scores  
meta_may_3 <- as(sample_data(may_3_phylo), "data.frame")  
meta_may_3$SampleID <- rownames(meta_may_3)
ordination_metadata_may_3 <- left_join(meta_may_3, ordination_scores_may_3, by = "SampleID")

#Convert metadata to data frame  
ordination_metadata_may_3 <- as.data.frame(ordination_metadata_may_3)

#Define plastic groups  
ordination_metadata_may_3$plastic_group <- ifelse(
  ordination_metadata_may_3$plastic_type %in% c("TPU181", "TPUFC2.1"), 
  "Bioplastics", 
  "Elastollan"
)
ordination_metadata_may_3$plastic_group <- factor(ordination_metadata_may_3$plastic_group, 
                                                  levels = c("Elastollan", "Bioplastics"))

#Process OTU and taxonomy at family level
otu <- otu_table(may_3_phylo)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
otu_table_may_3_family <- as.data.frame(as(otu, "matrix"))

tax_table_may_3_family <- as.data.frame(tax_table(may_3_phylo))
tax_table_may_3_family[] <- lapply(tax_table_may_3_family, trimws)

common_asvs <- intersect(rownames(tax_table_may_3_family), colnames(otu_table_may_3_family))
tax_table_may_3_family <- tax_table_may_3_family[common_asvs, , drop = FALSE]
otu_table_may_3_family <- otu_table_may_3_family[, common_asvs, drop = FALSE]

otu_table_may_3_family <- as.data.frame(t(otu_table_may_3_family))
otu_table_may_3_family$Family <- tax_table_may_3_family$Family
otu_table_may_3_family <- otu_table_may_3_family[!is.na(otu_table_may_3_family$Family) & 
                                                         otu_table_may_3_family$Family != "", ]
otu_table_may_3_family <- aggregate(. ~ Family, data = otu_table_may_3_family, FUN = sum)

#Set family as column names and transpose to sample × taxa for envfit
rownames(otu_table_may_3_family) <- otu_table_may_3_family$Family
otu_table_may_3_family$Family <- NULL
otu_table_may_3_family <- as.data.frame(t(otu_table_may_3_family))

#Process OTU and taxonomy at genus level
otu <- otu_table(may_3_phylo)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
otu_table_may_3_genus <- as.data.frame(as(otu, "matrix"))

tax_table_may_3_genus <- as.data.frame(tax_table(may_3_phylo))
tax_table_may_3_genus[] <- lapply(tax_table_may_3_genus, trimws)

common_asvs <- intersect(rownames(tax_table_may_3_genus), colnames(otu_table_may_3_genus))
tax_table_may_3_genus <- tax_table_may_3_genus[common_asvs, , drop = FALSE]
otu_table_may_3_genus <- otu_table_may_3_genus[, common_asvs, drop = FALSE]

otu_table_may_3_genus <- as.data.frame(t(otu_table_may_3_genus))
otu_table_may_3_genus$Genus <- tax_table_may_3_genus$Genus
otu_table_may_3_genus <- otu_table_may_3_genus[!is.na(otu_table_may_3_genus$Genus) & 
                                                       otu_table_may_3_genus$Genus != "", ]
otu_table_may_3_genus <- aggregate(. ~ Genus, data = otu_table_may_3_genus, FUN = sum)

rownames(otu_table_may_3_genus) <- otu_table_may_3_genus$Genus
otu_table_may_3_genus$Genus <- NULL
otu_table_may_3_genus <- as.data.frame(t(otu_table_may_3_genus))

#Calculate mean abundance
mean_abundance_families <- colMeans(otu_table_may_3_family)
mean_abundance_genera <- colMeans(otu_table_may_3_genus)

#Align samples for envfit
shared_samples_family <- intersect(rownames(ordination_scores_may_3), rownames(otu_table_may_3_family))
ordination_envfit_family <- ordination_scores_may_3[shared_samples_family, 1:2]
otu_envfit_family <- otu_table_may_3_family[shared_samples_family, ]

shared_samples_genus <- intersect(rownames(ordination_scores_may_3), rownames(otu_table_may_3_genus))
ordination_envfit_genus <- ordination_scores_may_3[shared_samples_genus, 1:2]
otu_envfit_genus <- otu_table_may_3_genus[shared_samples_genus, ]

#Remove non-variable families
otu_envfit_family <- otu_envfit_family[, apply(otu_envfit_family, 2, function(x) var(x, na.rm = TRUE) > 0)]

#Remove non-variable genera
otu_envfit_genus <- otu_envfit_genus[, apply(otu_envfit_genus, 2, function(x) var(x, na.rm = TRUE) > 0)]


#Run envfit and extract significant taxa (Family)
set.seed(72)
envfit_may_3_family <- envfit(ordination_envfit_family, otu_envfit_family, permutations = 999)
sig_families_may_3 <- as.data.frame(envfit_may_3_family$vectors$arrows)
sig_families_may_3$p_value <- envfit_may_3_family$vectors$pvals
sig_families_may_3 <- sig_families_may_3[sig_families_may_3$p_value < 0.05, ]
sig_families_may_3$Family <- rownames(sig_families_may_3)
sig_families_may_3$mean_abundance <- mean_abundance_families[rownames(sig_families_may_3)]
sig_families_may_3 <- sig_families_may_3[order(-sig_families_may_3$mean_abundance), ]


#Run envfit and extract significant taxa (Genus)
set.seed(73)
envfit_may_3_genus <- envfit(ordination_envfit_genus, otu_envfit_genus, permutations = 999)
sig_genera_may_3 <- as.data.frame(envfit_may_3_genus$vectors$arrows)
sig_genera_may_3$p_value <- envfit_may_3_genus$vectors$pvals
sig_genera_may_3 <- sig_genera_may_3[sig_genera_may_3$p_value < 0.05, ]
sig_genera_may_3$Genus <- rownames(sig_genera_may_3)
sig_genera_may_3$mean_abundance <- mean_abundance_genera[rownames(sig_genera_may_3)]
sig_genera_may_3 <- sig_genera_may_3[order(-sig_genera_may_3$mean_abundance), ]


#Select top 10 families and genera for biplot
top_10_families_may_3 <- head(sig_families_may_3, 10)
vectors_df_may_3_family <- data.frame(
  Axis.1 = top_10_families_may_3$Axis.1 * 0.2,  
  Axis.2 = top_10_families_may_3$Axis.2 * 0.2,
  Family = rownames(top_10_families_may_3)
)

top_10_genera_may_3 <- head(sig_genera_may_3, 10)
vectors_df_may_3_genus <- data.frame(
  Axis.1 = top_10_genera_may_3$Axis.1 * 0.2,  
  Axis.2 = top_10_genera_may_3$Axis.2 * 0.2,
  Genus = rownames(top_10_genera_may_3)
)

#Family biplot
biplot_may_3_family <- ggplot() +
  geom_point(data = ordination_metadata_may_3, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_may_3_family, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_may_3_family, aes(x = Axis.1, y = Axis.2, label = Family), 
                  color = "black", size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  labs(color = "Plastic type") +
  theme_minimal()

#Genus biplot
biplot_may_3_genus <- ggplot() +
  geom_point(data = ordination_metadata_may_3, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_may_3_genus, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_may_3_genus, aes(x = Axis.1, y = Axis.2, label = Genus), 
                  color = "black", size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  labs(color = "Plastic type") +
  theme_minimal()

#Print both plots
print(biplot_may_3_family)
print(biplot_may_3_genus)

######## JUNE 9 ##########
#Phyloseq object for June 9  
june_9_phylo <- subset_samples(UniFrac_Phylo_Object, date == "June 9")

#Run PCoA ordination  
set.seed(78)
ordination_june_9 <- ordinate(june_9_phylo, method = "PCoA", distance = "wunifrac")

#Extract ordination scores  
ordination_scores_june_9 <- as.data.frame(ordination_june_9$vectors[, 1:2])
ordination_scores_june_9$SampleID <- rownames(ordination_scores_june_9)

#Extract metadata and merge with ordination scores  
meta_june_9 <- as(sample_data(june_9_phylo), "data.frame")  
meta_june_9$SampleID <- rownames(meta_june_9)
ordination_metadata_june_9 <- left_join(meta_june_9, ordination_scores_june_9, by = "SampleID")

#Convert metadata to data frame  
ordination_metadata_june_9 <- as.data.frame(ordination_metadata_june_9)

#Define plastic groups  
ordination_metadata_june_9$plastic_group <- ifelse(
  ordination_metadata_june_9$plastic_type %in% c("TPU181", "TPUFC2.1"), 
  "Bioplastics", 
  "Elastollan"
)
ordination_metadata_june_9$plastic_group <- factor(ordination_metadata_june_9$plastic_group, 
                                                   levels = c("Elastollan", "Bioplastics"))

#Process OTU and taxonomy at Family level
otu <- otu_table(june_9_phylo)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
otu_table_june_9_family <- as.data.frame(as(otu, "matrix"))

tax_table_june_9_family <- as.data.frame(tax_table(june_9_phylo))
tax_table_june_9_family[] <- lapply(tax_table_june_9_family, trimws)

common_asvs <- intersect(rownames(tax_table_june_9_family), colnames(otu_table_june_9_family))
tax_table_june_9_family <- tax_table_june_9_family[common_asvs, , drop = FALSE]
otu_table_june_9_family <- otu_table_june_9_family[, common_asvs, drop = FALSE]

otu_table_june_9_family <- as.data.frame(t(otu_table_june_9_family))
otu_table_june_9_family$Family <- tax_table_june_9_family$Family
otu_table_june_9_family <- otu_table_june_9_family[!is.na(otu_table_june_9_family$Family) & 
                                                         otu_table_june_9_family$Family != "", ]
otu_table_june_9_family <- aggregate(. ~ Family, data = otu_table_june_9_family, FUN = sum)

#Set family as column names and transpose to sample × taxa for envfit
rownames(otu_table_june_9_family) <- otu_table_june_9_family$Family
otu_table_june_9_family$Family <- NULL
otu_table_june_9_family <- as.data.frame(t(otu_table_june_9_family))

#Process OTU and taxonomy at genus level
otu <- otu_table(june_9_phylo)
if (taxa_are_rows(otu)) {
  otu <- t(otu)
}
otu_table_june_9_genus <- as.data.frame(as(otu, "matrix"))

tax_table_june_9_genus <- as.data.frame(tax_table(june_9_phylo))
tax_table_june_9_genus[] <- lapply(tax_table_june_9_genus, trimws)

common_asvs <- intersect(rownames(tax_table_june_9_genus), colnames(otu_table_june_9_genus))
tax_table_june_9_genus <- tax_table_june_9_genus[common_asvs, , drop = FALSE]
otu_table_june_9_genus <- otu_table_june_9_genus[, common_asvs, drop = FALSE]

otu_table_june_9_genus <- as.data.frame(t(otu_table_june_9_genus))
otu_table_june_9_genus$Genus <- tax_table_june_9_genus$Genus
otu_table_june_9_genus <- otu_table_june_9_genus[!is.na(otu_table_june_9_genus$Genus) & 
                                                       otu_table_june_9_genus$Genus != "", ]
otu_table_june_9_genus <- aggregate(. ~ Genus, data = otu_table_june_9_genus, FUN = sum)

rownames(otu_table_june_9_genus) <- otu_table_june_9_genus$Genus
otu_table_june_9_genus$Genus <- NULL
otu_table_june_9_genus <- as.data.frame(t(otu_table_june_9_genus))

#Calculate mean abundance
mean_abundance_families <- colMeans(otu_table_june_9_family)
mean_abundance_genera <- colMeans(otu_table_june_9_genus)

#Align samples for envfit
shared_samples_family <- intersect(rownames(ordination_scores_june_9), rownames(otu_table_june_9_family))
ordination_envfit_family <- ordination_scores_june_9[shared_samples_family, 1:2]
otu_envfit_family <- otu_table_june_9_family[shared_samples_family, ]

shared_samples_genus <- intersect(rownames(ordination_scores_june_9), rownames(otu_table_june_9_genus))
ordination_envfit_genus <- ordination_scores_june_9[shared_samples_genus, 1:2]
otu_envfit_genus <- otu_table_june_9_genus[shared_samples_genus, ]

#Remove non-variable families
otu_envfit_family <- otu_envfit_family[, apply(otu_envfit_family, 2, function(x) var(x, na.rm = TRUE) > 0)]

#Remove non-variable genera
otu_envfit_genus <- otu_envfit_genus[, apply(otu_envfit_genus, 2, function(x) var(x, na.rm = TRUE) > 0)]


#Run envfit and extract significant taxa (Family)
set.seed(79)
envfit_june_9_family <- envfit(ordination_envfit_family, otu_envfit_family, permutations = 999)
sig_families_june_9 <- as.data.frame(envfit_june_9_family$vectors$arrows)
sig_families_june_9$p_value <- envfit_june_9_family$vectors$pvals
sig_families_june_9 <- sig_families_june_9[sig_families_june_9$p_value < 0.05, ]
sig_families_june_9$Family <- rownames(sig_families_june_9)
sig_families_june_9$mean_abundance <- mean_abundance_families[rownames(sig_families_june_9)]
sig_families_june_9 <- sig_families_june_9[order(-sig_families_june_9$mean_abundance), ]


#Run envfit and extract significant taxa (Genus)
set.seed(80)
envfit_june_9_genus <- envfit(ordination_envfit_genus, otu_envfit_genus, permutations = 999)
sig_genera_june_9 <- as.data.frame(envfit_june_9_genus$vectors$arrows)
sig_genera_june_9$p_value <- envfit_june_9_genus$vectors$pvals
sig_genera_june_9 <- sig_genera_june_9[sig_genera_june_9$p_value < 0.05, ]
sig_genera_june_9$Genus <- rownames(sig_genera_june_9)
sig_genera_june_9$mean_abundance <- mean_abundance_genera[rownames(sig_genera_june_9)]
sig_genera_june_9 <- sig_genera_june_9[order(-sig_genera_june_9$mean_abundance), ]


#Select top 10 families and genera for biplot
top_10_families_june_9 <- head(sig_families_june_9, 10)
vectors_df_june_9_family <- data.frame(
  Axis.1 = top_10_families_june_9$Axis.1 * 0.2,  
  Axis.2 = top_10_families_june_9$Axis.2 * 0.2,
  Family = rownames(top_10_families_june_9)
)

top_10_genera_june_9 <- head(sig_genera_june_9, 10)
vectors_df_june_9_genus <- data.frame(
  Axis.1 = top_10_genera_june_9$Axis.1 * 0.2,  
  Axis.2 = top_10_genera_june_9$Axis.2 * 0.2,
  Genus = rownames(top_10_genera_june_9)
)

#Family biplot
biplot_june_9_family <- ggplot() +
  geom_point(data = ordination_metadata_june_9, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_june_9_family, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_june_9_family, aes(x = Axis.1, y = Axis.2, label = Family), 
                  color = "black", size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  labs(color = "Plastic type") +
  theme_minimal()

#Genus biplot
biplot_june_9_genus <- ggplot() +
  geom_point(data = ordination_metadata_june_9, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_june_9_genus, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_june_9_genus, aes(x = Axis.1, y = Axis.2, label = Genus), 
                  color = "black", size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  labs(color = "Plastic type") +
  theme_minimal()

#Print both plots
print(biplot_june_9_family)
print(biplot_june_9_genus)


#Print family level plots
print(biplot_march_7_family)
print(biplot_march_23_family)
print(biplot_april_12_family)
print(biplot_may_3_family)
print(biplot_june_9_family)

#Print genus level plots
print(biplot_march_7_genus)
print(biplot_march_23_genus)
print(biplot_april_12_genus)
print(biplot_may_3_genus)
print(biplot_june_9_genus)

#Function to add family level vectors to original PCoA plot
add_family_vectors <- function(pcoa_plot, vectors_df) {
  pcoa_plot + 
    geom_segment(
      data = vectors_df, 
      aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
      arrow = arrow(length = unit(0.2, "cm")), 
      color = "black"
    ) +
    geom_text_repel(  #So it doesn't override existing labels
      data = vectors_df, 
      aes(x = Axis.1, y = Axis.2, label = Family), 
      color = "black", size = 3,
      nudge_x = vectors_df$Axis.1 * 0.1,  
      nudge_y = vectors_df$Axis.2 * 0.1,  
      max.overlaps = 15
    )
}

#Function to add genus level vectors to an existing PCoA plot
add_genus_vectors <- function(pcoa_plot, vectors_df) {
  pcoa_plot + 
    geom_segment(
      data = vectors_df, 
      aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
      arrow = arrow(length = unit(0.2, "cm")), 
      color = "black"
    ) +
    geom_text_repel(
      data = vectors_df, 
      aes(x = Axis.1, y = Axis.2, label = Genus),  
      color = "black", size = 3, 
      nudge_x = vectors_df$Axis.1 * 0.1,  #Shift text away from arrows
      nudge_y = vectors_df$Axis.2 * 0.1,  
      max.overlaps = 15
    )
}

#Apply function to overlay vectors on each date specific PCoA plot
pcoa_plot_march_7_icon_family <- add_family_vectors(pcoa_plot_march_7_icon, vectors_df_march_7_family)
pcoa_plot_march_23_icon_family <- add_family_vectors(pcoa_plot_march_23_icon, vectors_df_march_23_family)
pcoa_plot_april_12_icon_family <- add_family_vectors(pcoa_plot_april_12_icon, vectors_df_april_12_family)
pcoa_plot_may_3_icon_family <- add_family_vectors(pcoa_plot_may_3_icon, vectors_df_may_3_family)
pcoa_plot_june_9_icon_family <- add_family_vectors(pcoa_plot_june_9_icon, vectors_df_june_9_family)

#Print updated family plots
print(pcoa_plot_march_7_icon_family)
print(pcoa_plot_march_23_icon_family)
print(pcoa_plot_april_12_icon_family) #MANUSCRIPT Extended Data Fig. 1]
print(pcoa_plot_may_3_icon_family) #MANUSCRIPT Extended Data Fig. 1]
print(pcoa_plot_june_9_icon_family) #MANUSCRIPT Extended Data Fig. 1]

#Print the list of each top 10
top_10_families_march_7$Family
top_10_families_march_23$Family
top_10_families_april_12$Family
top_10_families_may_3$Family
top_10_families_june_9$Family

#Apply function to overlay vectors on each date specific PCoA plot
pcoa_plot_march_7_icon_genus <- add_genus_vectors(pcoa_plot_march_7_icon, vectors_df_march_7_genus)
pcoa_plot_march_23_icon_genus <- add_genus_vectors(pcoa_plot_march_23_icon, vectors_df_march_23_genus)
pcoa_plot_april_12_icon_genus <- add_genus_vectors(pcoa_plot_april_12_icon, vectors_df_april_12_genus)
pcoa_plot_may_3_icon_genus <- add_genus_vectors(pcoa_plot_may_3_icon, vectors_df_may_3_genus)
pcoa_plot_june_9_icon_genus <- add_genus_vectors(pcoa_plot_june_9_icon, vectors_df_june_9_genus)

#Print updated genus plots
print(pcoa_plot_march_7_icon_genus)
print(pcoa_plot_march_23_icon_genus)
print(pcoa_plot_april_12_icon_genus)
print(pcoa_plot_may_3_icon_genus)
print(pcoa_plot_june_9_icon_genus)


######################################
###STOP##STOP##STOP###################
######################################

####DESEq2 test - FOR JUNE 9####
#load.packages("DESeq2")
#load.packages("indicspecies")
#load.packages("ggrepel")
#load.packages("tibble")
#library(DESeq2)
#library(indicspecies)
#library(ggrepel)
#library(tibble)

#Need to use raw data so need to recreate the subset data
#Extract and align sample metadata
sample_meta_df <- as.data.frame(BioPlastics_phylo@sam_data)
sample_meta_df <- sample_meta_df[sample_names(BioPlasticsASVtable1), ]
cat("Sample metadata aligned:", nrow(sample_meta_df), "samples\n")

#Set SampleID as rownames so it can match correctly
rownames(ordination_metadata_june_9) <- ordination_metadata_june_9$SampleID

#Reassign plastic_group now that rownames are aligned
sample_meta_df$plastic_group <- ordination_metadata_june_9[rownames(sample_meta_df), "plastic_group"]

#Rebuild the sample_data object
sample_meta_clean <- sample_data(sample_meta_df)

#Rebuild the phyloseq object
phylo_raw <- phyloseq(
  otu_table(BioPlasticsASVtable1, taxa_are_rows = FALSE),
  tax_table(BioPlastics_phylo@tax_table),
  sample_meta_clean,
  phy_tree(BioPlastics_phylo@phy_tree)
)

#Subset to June 9
june_9_phylo <- subset_samples(phylo_raw, date == "June_9")

#Check that plastic_group is correctly assigned
head(sample_data(june_9_phylo)$plastic_group)

#Filter out rare taxa (low read counts)
june_9_phylo_filtered <- prune_taxa(taxa_sums(june_9_phylo) > 10, june_9_phylo)

#Convert phyloseq object to DESeq2 format
dds <- phyloseq_to_deseq2(june_9_phylo_filtered, ~ plastic_group)

#Run DESeq2 analysis
dds <- DESeq(dds)

#Pull out results (Bioplastics vs. Elastollan)
res <- results(dds, contrast = c("plastic_group", "Bioplastics", "Elastollan"))

#Order results by significance (padj)
res <- res[order(res$padj, na.last = NA), ]

#Convert to a data frame and remove NAs
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "ASV") %>%
  filter(!is.na(padj))

#Merge taxonomic information
tax_table_june_9 <- as.data.frame(tax_table(june_9_phylo)) %>%
  rownames_to_column(var = "ASV")  #Make sure ASVs are stored as a column

#Merge DESeq2 results with taxonomy
res_df <- left_join(res_df, tax_table_june_9, by = "ASV")

#Filter significant taxa (padj)
sig_taxa <- res_df %>% filter(padj < 0.05)

#Print significant taxa
print(sig_taxa)

cat("Number of significantly different taxa:", nrow(sig_taxa), "\n")

sig_taxa <- as.data.frame(sig_taxa)  #Mure sure it's a dataframe
str(sig_taxa)
sig_taxa %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Bioplastics", "Elastollan")) %>%
  group_by(Direction) %>%
  summarise(Count = n())

#Print sig_taxa
print(sig_taxa$Family)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "DESeq2 Differential Abundance Analysis",
       x = "Log2 Fold Change (Bioplastics vs. Elastollan)",
       y = "-log10 Adjusted p-value") +
  theme_minimal()

top_taxa <- sig_taxa %>% arrange(padj) %>% head(15)  #Select top 15
print(top_taxa$Family)

ggplot(top_taxa, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"), labels = c("Elastollan", "Bioplastics")) +
  labs(title = "Top Differentially Abundant Taxa",
       x = "Family",
       y = "Log2 Fold Change (Bioplastics vs. Elastollan)",
       fill = "Higher in") +
  theme_minimal()



####DESEq2 test #### For May 3 ####
#Extract and align sample metadata
sample_meta_df <- as.data.frame(BioPlastics_phylo@sam_data)
sample_meta_df <- sample_meta_df[sample_names(BioPlasticsASVtable1), ]

#Add plastic_group from ordination metadata
#Set SampleID as rownames so it can match correctly
rownames(ordination_metadata_may_3) <- ordination_metadata_may_3$SampleID

#Reassign plastic_group now that rownames are aligned
sample_meta_df$plastic_group <- ordination_metadata_may_3[rownames(sample_meta_df), "plastic_group"]

#Rebuild the sample_data object
sample_meta_clean <- sample_data(sample_meta_df)

#Rebuild the phyloseq object
phylo_raw_may3 <- phyloseq(
  otu_table(BioPlasticsASVtable1, taxa_are_rows = FALSE),
  tax_table(BioPlastics_phylo@tax_table),
  sample_meta_clean,
  phy_tree(BioPlastics_phylo@phy_tree)
)

#Subset for May 3 samples
may_3_phylo <- subset_samples(phylo_raw_may3, date == "May_3")

#Check that plastic_group is correctly assigned
head(sample_data(may_3_phylo)$plastic_group)

#Filter out rare taxa (low read counts)
may_3_phylo_filtered <- prune_taxa(taxa_sums(may_3_phylo) > 10, may_3_phylo)

#Convert phyloseq object to DESeq2 format
dds <- phyloseq_to_deseq2(may_3_phylo_filtered, ~ plastic_group)

#Run DESeq2 analysis
dds <- DESeq(dds)

#Pull out results (Bioplastics vs. Elastollan)
res <- results(dds, contrast = c("plastic_group", "Bioplastics", "Elastollan"))

#Order results by significance (padj)
res <- res[order(res$padj, na.last = NA), ]

#Convert to a data frame and remove NAs
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "ASV") %>%
  filter(!is.na(padj))

#Merge taxonomic information
tax_table_may_3 <- as.data.frame(tax_table(may_3_phylo)) %>%
  rownames_to_column(var = "ASV")  #Make sure ASVs are stored as a column

#Merge DESeq2 results with taxonomy
res_df <- left_join(res_df, tax_table_may_3, by = "ASV")

#Filter significant taxa (padj)
sig_taxa <- res_df %>% filter(padj < 0.05)

#Print significant taxa
print(sig_taxa)

cat("Number of significantly different taxa:", nrow(sig_taxa), "\n")

sig_taxa <- as.data.frame(sig_taxa)  #Ensure it's a dataframe
str(sig_taxa)
sig_taxa %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Bioplastics", "Elastollan")) %>%
  group_by(Direction) %>%
  summarise(Count = n())

print(sig_taxa$Family)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "DESeq2 Differential Abundance Analysis",
       x = "Log2 Fold Change (Bioplastics vs. Elastollan)",
       y = "-log10 Adjusted p-value") +
  theme_minimal()

top_taxa <- sig_taxa %>% arrange(padj) %>% head(15)  # Select top 15
print(top_taxa$Family) 

ggplot(top_taxa, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"), labels = c("Elastollan", "Bioplastics")) +
  labs(title = "Top Differentially Abundant Taxa",
       x = "Family",
       y = "Log2 Fold Change (Bioplastics vs. Elastollan)",
       fill = "Higher in") +
  theme_minimal()


####DESEq2 test #### For April 12 ####
#Extract and align sample metadata
sample_meta_df <- as.data.frame(BioPlastics_phylo@sam_data)
sample_meta_df <- sample_meta_df[sample_names(BioPlasticsASVtable1), ]

#Add plastic_group from ordination metadata
#Set SampleID as rownames so we can match correctly
rownames(ordination_metadata_april_12) <- ordination_metadata_april_12$SampleID

#Reassign plastic_group now that rownames are aligned
sample_meta_df$plastic_group <- ordination_metadata_april_12[rownames(sample_meta_df), "plastic_group"]

#Rebuild the sample_data object
sample_meta_clean <- sample_data(sample_meta_df)

#Rebuild the phyloseq object
phylo_raw_april12 <- phyloseq(
  otu_table(BioPlasticsASVtable1, taxa_are_rows = FALSE),
  tax_table(BioPlastics_phylo@tax_table),
  sample_meta_clean,
  phy_tree(BioPlastics_phylo@phy_tree)
)

#Subset to April 12 samples
april_12_phylo <- subset_samples(phylo_raw_april12, date == "April_12")

#Check that plastic_group is correctly assigned
head(sample_data(april_12_phylo)$plastic_group)

#Filter out rare taxa (low read counts)
april_12_phylo_filtered <- prune_taxa(taxa_sums(april_12_phylo) > 10, april_12_phylo)

#Convert phyloseq object to DESeq2 format
dds <- phyloseq_to_deseq2(april_12_phylo_filtered, ~ plastic_group)

#Run DESeq2 analysis
dds <- DESeq(dds)

#Pull out results (Bioplastics vs. Elastollan)
res <- results(dds, contrast = c("plastic_group", "Bioplastics", "Elastollan"))

#Order results by significance (padj)
res <- res[order(res$padj, na.last = NA), ]

#Convert to a data frame and remove NAs
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "ASV") %>%
  filter(!is.na(padj))

#Merge taxonomic information
tax_table_april_12 <- as.data.frame(tax_table(april_12_phylo)) %>%
  rownames_to_column(var = "ASV")  # Ensure ASVs are stored as a column

#Merge DESeq2 results with taxonomy
res_df <- left_join(res_df, tax_table_april_12, by = "ASV")

#Filter significant taxa (padj)
sig_taxa <- res_df %>% filter(padj < 0.05)

#Print significant taxa
print(sig_taxa)

cat("Number of significantly different taxa:", nrow(sig_taxa), "\n")

sig_taxa <- as.data.frame(sig_taxa)  #Ensure it's a dataframe
str(sig_taxa)
sig_taxa %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Bioplastics", "Elastollan")) %>%
  group_by(Direction) %>%
  summarise(Count = n())

print(sig_taxa$Family)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "DESeq2 Differential Abundance Analysis",
       x = "Log2 Fold Change (Bioplastics vs. Elastollan)",
       y = "-log10 Adjusted p-value") +
  theme_minimal()

top_taxa <- sig_taxa %>% arrange(padj) %>% head(15)  # Select top 15
print(top_taxa$Family) 

ggplot(top_taxa, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"), labels = c("Elastollan", "Bioplastics")) +
  labs(title = "Top Differentially Abundant Taxa",
       x = "Family",
       y = "Log2 Fold Change (Bioplastics vs. Elastollan)",
       fill = "Higher in") +
  theme_minimal()


######################################
###STOP##STOP##STOP###################
######################################



####Richness over time#######

#Estimate richness (observed ASVs, Shannon)
richness_df <- estimate_richness(BioPlastics_phylo, measures = c("Observed", "Shannon"))

#Add sample metadata
richness_df <- cbind(richness_df, sample_data(BioPlastics_phylo))

#Adjust plastic_conc to factor for ordering
richness_df$plastic_conc <- factor(richness_df$plastic_conc, levels = sort(unique(richness_df$plastic_conc)))

#Create boxplot
richness_plot <- ggplot(richness_df, aes(x = date, y = Observed, fill = plastic_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = plastic_type), width = 0.2, size = 2, alpha = 0.7) +
  facet_grid(plastic_type ~ plastic_conc, scales = "free_x", space = "free") +
  scale_fill_manual(values = color_mapping) +
  scale_color_manual(values = color_mapping) +
  theme_minimal() +
  labs(
    x = "Date",
    y = "Observed Richness (ASVs)",
    fill = "Plastic Type",
    color = "Plastic Type"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title = element_text(face = "bold", size = 13),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(ncol = 3))

#Print plot
print(richness_plot)

dot_plot <- ggplot(richness_df, aes(x = plastic_conc, y = Observed, color = plastic_type)) +
  geom_point(size = 4, alpha = 0.8) +
  facet_wrap(~ date) +
  theme_minimal() +
  labs(x = "Plastic Concentration", y = "Observed ASV Richness") +
  scale_color_manual(values = color_mapping) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

#Print plot
print(dot_plot)


#Make sure date is in correct order
richness_df$date <- factor(richness_df$date, levels = c("March_7", "March_23", "April_12", "May_3", "June_9"))

#Order concentration for display (numeric or factor levels)
richness_df$plastic_conc <- factor(richness_df$plastic_conc, levels = sort(unique(richness_df$plastic_conc)))

#Plot
heatmap_plot <- ggplot(richness_df, aes(x = date, y = plastic_conc, fill = Observed)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Observed, 1)), color = "black", size = 3) +
  facet_wrap(~ plastic_type, nrow = 1) +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    x = "Date",
    y = "Plastic Concentration",
    fill = "Richness (ASVs)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    legend.position = "right"
  )

print(heatmap_plot)



######################################
###STOP##STOP##STOP###################
######################################

#Load in packages
#install.packages("indicspecies")
#BiocManager::install("DESeq2")
#library(indicspecies)
#library(DESeq2)


#### Phylum level analysis for 16S ####
#Subset the phyloseq object to focus on phylum level
BioP_phylum_16S <- tax_glom(BioPlastics_phylo, taxrank = "Phylum")

#So variables in sample_data are correctly formatted
sample_data_df_16S <- as.data.frame(sample_data(BioP_phylum_16S))  # Convert to data frame

#Convert variables to factors or numeric as required
sample_data_df_16S$plastic_type <- as.factor(sample_data_df_16S$plastic_type)
sample_data_df_16S$plastic_conc <- as.numeric(as.character(sample_data_df_16S$plastic_conc))  # Ensure numeric
sample_data_df_16S$date <- as.factor(sample_data_df_16S$date)

#Update the phyloseq object with the corrected sample data
sample_data(BioP_phylum_16S) <- sample_data(sample_data_df_16S)

#### Phylum Level ####
cat("\nRunning DESeq2 Analysis at Phylum Level (16S)...\n")
dds_phylum_16S <- phyloseq_to_deseq2(BioP_phylum_16S, ~ plastic_type * plastic_conc + date)
dds_phylum_16S <- DESeq(dds_phylum_16S)

#Get available terms in the model
cat("\nModel terms for Phylum Level (16S):\n")
print(resultsNames(dds_phylum_16S))

#Specify terms to extract
terms_to_extract_phylum_16S <- c("plastic_type_TPU181_vs_Elastollan",
                                 "plastic_type_TPUFC2.1_vs_Elastollan",
                                 "plastic_conc",
                                 "plastic_typeTPU181.plastic_conc",
                                 "plastic_typeTPUFC2.1.plastic_conc")

#Extract taxonomic names for phylum level
tax_names_phylum_16S <- as.data.frame(tax_table(BioP_phylum_16S))
tax_names_phylum_16S$ASV <- rownames(tax_names_phylum_16S)  # Add ASV as a column

#Function to extract and process significant results for a term
extract_significant_results_16S <- function(term, dds_object) {
  cat("\nProcessing term (16S):", term, "\n")
  results_phylum <- results(dds_object, name = term)
  
  #Filter significant results
  sig_results <- results_phylum[which(results_phylum$padj < 0.05), ]
  
  if (nrow(sig_results) > 0) {
    sig_results <- as.data.frame(sig_results)
    sig_results$effect <- term  #Add effect as a column
    return(sig_results)
  } else {
    return(NULL)  #Return NULL if no significant results
  }
}

#Apply function to all terms and combine results
sig_phylum_list_16S <- lapply(terms_to_extract_phylum_16S, extract_significant_results_16S, dds_object = dds_phylum_16S)
sig_phylum_list_16S <- Filter(Negate(is.null), sig_phylum_list_16S)  # Remove NULL elements

#Combine all significant results
sig_phylum_16S <- do.call(rbind, sig_phylum_list_16S)

#Merge with taxonomy
sig_phylum_16S <- merge(sig_phylum_16S, tax_names_phylum_16S, by.x = "row.names", by.y = "ASV")
colnames(sig_phylum_16S)[1] <- "ASV"  # Rename the first column for clarity

#Validate the merge result
cat("\nRows in sig_phylum_16S after merging:\n", nrow(sig_phylum_16S), "\n")
print(head(sig_phylum_16S))

#### Phylum Level Bar Plots ####

#Effect: plastic_type_TPU181_vs_Elastollan
if (!is.null(sig_phylum_16S) && nrow(sig_phylum_16S[sig_phylum_16S$effect == "plastic_type_TPU181_vs_Elastollan", ]) > 0) {
  plot_phylum_TPU181_16S <- ggplot(
    sig_phylum_16S[sig_phylum_16S$effect == "plastic_type_TPU181_vs_Elastollan", ],
    aes(x = reorder(get("Phylum"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Phylum",
      y = "Log2 Fold Change",
      title = "Differential Abundance - plastic_type_TPU181_vs_Elastollan (Phylum Level - 16S)"
    ) +
    theme_minimal()
  print(plot_phylum_TPU181_16S)
} else {
  cat("No data to plot for plastic_type_TPU181_vs_Elastollan (16S).\n")
}

#Effect: plastic_type_TPUFC2.1_vs_Elastollan
if (!is.null(sig_phylum_16S) && nrow(sig_phylum_16S[sig_phylum_16S$effect == "plastic_type_TPUFC2.1_vs_Elastollan", ]) > 0) {
  plot_phylum_TPUFC2_16S <- ggplot(
    sig_phylum_16S[sig_phylum_16S$effect == "plastic_type_TPUFC2.1_vs_Elastollan", ],
    aes(x = reorder(get("Phylum"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Phylum",
      y = "Log2 Fold Change",
      title = "Differential Abundance - plastic_type_TPUFC2.1_vs_Elastollan (Phylum Level - 16S)"
    ) +
    theme_minimal()
  print(plot_phylum_TPUFC2_16S)
} else {
  cat("No data to plot for plastic_type_TPUFC2.1_vs_Elastollan (16S).\n")
}

#Effect: plastic_conc
if (!is.null(sig_phylum_16S) && nrow(sig_phylum_16S[sig_phylum_16S$effect == "plastic_conc", ]) > 0) {
  plot_phylum_conc_16S <- ggplot(
    sig_phylum_16S[sig_phylum_16S$effect == "plastic_conc", ],
    aes(x = reorder(get("Phylum"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Phylum",
      y = "Log2 Fold Change",
      title = "Differential Abundance - plastic_conc (Phylum Level - 16S)"
    ) +
    theme_minimal()
  print(plot_phylum_conc_16S)
} else {
  cat("No data to plot for plastic_conc (16S).\n")
}

#Effect: plastic_typeTPU181.plastic_conc
if (!is.null(sig_phylum_16S) && nrow(sig_phylum_16S[sig_phylum_16S$effect == "plastic_typeTPU181.plastic_conc", ]) > 0) {
  plot_phylum_TPU181_conc_16S <- ggplot(
    sig_phylum_16S[sig_phylum_16S$effect == "plastic_typeTPU181.plastic_conc", ],
    aes(x = reorder(get("Phylum"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Phylum",
      y = "Log2 Fold Change",
      title = "Differential Abundance - plastic_typeTPU181.plastic_conc (Phylum Level - 16S)"
    ) +
    theme_minimal()
  print(plot_phylum_TPU181_conc_16S)
} else {
  cat("No data to plot for plastic_typeTPU181.plastic_conc (16S).\n")
}

#Effect: plastic_typeTPUFC2.1.plastic_conc
if (!is.null(sig_phylum_16S) && nrow(sig_phylum_16S[sig_phylum_16S$effect == "plastic_typeTPUFC2.1.plastic_conc", ]) > 0) {
  plot_phylum_TPUFC2_conc_16S <- ggplot(
    sig_phylum_16S[sig_phylum_16S$effect == "plastic_typeTPUFC2.1.plastic_conc", ],
    aes(x = reorder(get("Phylum"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Phylum",
      y = "Log2 Fold Change",
      title = "Differential Abundance - plastic_typeTPUFC2.1.plastic_conc (Phylum Level - 16S)"
    ) +
    theme_minimal()
  print(plot_phylum_TPUFC2_conc_16S)
} else {
  cat("No data to plot for plastic_typeTPUFC2.1.plastic_conc (16S).\n")
}



#### Family level analysis for 16S ####
#Subset the phyloseq object to focus on family level
BioP_family_16S <- tax_glom(BioPlastics_phylo, taxrank = "Family")

#So variables in sample_data are correctly formatted
sample_data_df_family_16S <- as.data.frame(sample_data(BioP_family_16S))  # Convert to data frame

#Convert variables to factors or numeric as required
sample_data_df_family_16S$plastic_type <- as.factor(sample_data_df_family_16S$plastic_type)
sample_data_df_family_16S$plastic_conc <- as.numeric(as.character(sample_data_df_family_16S$plastic_conc))  # Ensure numeric
sample_data_df_family_16S$date <- as.factor(sample_data_df_family_16S$date)

#Update the phyloseq object with the corrected sample data
sample_data(BioP_family_16S) <- sample_data(sample_data_df_family_16S)

#### Family Level ####
cat("\nRunning DESeq2 Analysis at Family Level (16S)...\n")
dds_family_16S <- phyloseq_to_deseq2(BioP_family_16S, ~ plastic_type * plastic_conc + date)
dds_family_16S <- DESeq(dds_family_16S)

#Get available terms in the model
cat("\nModel terms for Family Level (16S):\n")
print(resultsNames(dds_family_16S))

#Specify terms to extract
terms_to_extract_family_16S <- c("plastic_type_TPU181_vs_Elastollan",
                                 "plastic_type_TPUFC2.1_vs_Elastollan",
                                 "plastic_conc",
                                 "plastic_typeTPU181.plastic_conc",
                                 "plastic_typeTPUFC2.1.plastic_conc")

#Extract taxonomic names for family level
tax_names_family_16S <- as.data.frame(tax_table(BioP_family_16S))
tax_names_family_16S$ASV <- rownames(tax_names_family_16S)  # Add ASV as a column

#Function to extract and process significant results for a term
extract_significant_results_family_16S <- function(term, dds_object) {
  cat("\nProcessing term (16S - Family):", term, "\n")
  results_family <- results(dds_object, name = term)
  
  #Filter significant results
  sig_results <- results_family[which(results_family$padj < 0.05), ]
  
  if (nrow(sig_results) > 0) {
    sig_results <- as.data.frame(sig_results)
    sig_results$effect <- term  #Add effect as a column
    return(sig_results)
  } else {
    return(NULL)  #Return NULL if no significant results
  }
}

#Apply function to all terms and combine results
sig_family_list_16S <- lapply(terms_to_extract_family_16S, extract_significant_results_family_16S, dds_object = dds_family_16S)
sig_family_list_16S <- Filter(Negate(is.null), sig_family_list_16S)  # Remove NULL elements

#Combine all significant results
sig_family_16S <- do.call(rbind, sig_family_list_16S)

#Merge with taxonomy
sig_family_16S <- merge(sig_family_16S, tax_names_family_16S, by.x = "row.names", by.y = "ASV")
colnames(sig_family_16S)[1] <- "ASV"  # Rename the first column for clarity

#Validate the merge result
cat("\nRows in sig_family_16S after merging:\n", nrow(sig_family_16S), "\n")
print(head(sig_family_16S))

#### Family Level Bar Plots ####

#Effect: plastic_type_TPU181_vs_Elastollan
if (!is.null(sig_family_16S) && nrow(sig_family_16S[sig_family_16S$effect == "plastic_type_TPU181_vs_Elastollan", ]) > 0) {
  plot_family_TPU181_16S <- ggplot(
    sig_family_16S[sig_family_16S$effect == "plastic_type_TPU181_vs_Elastollan", ],
    aes(x = reorder(get("Family"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Family",
      y = "Log2 Fold Change",
      title = "Differential Abundance - plastic_type_TPU181_vs_Elastollan (Family Level - 16S)"
    ) +
    theme_minimal()
  print(plot_family_TPU181_16S)
} else {
  cat("No data to plot for plastic_type_TPU181_vs_Elastollan (16S - Family).\n")
}

#Effect: plastic_type_TPUFC2.1_vs_Elastollan
if (!is.null(sig_family_16S) && nrow(sig_family_16S[sig_family_16S$effect == "plastic_type_TPUFC2.1_vs_Elastollan", ]) > 0) {
  plot_family_TPUFC2_16S <- ggplot(
    sig_family_16S[sig_family_16S$effect == "plastic_type_TPUFC2.1_vs_Elastollan", ],
    aes(x = reorder(get("Family"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Family",
      y = "Log2 Fold Change",
      title = "Differential Abundance - plastic_type_TPUFC2.1_vs_Elastollan (Family Level - 16S)"
    ) +
    theme_minimal()
  print(plot_family_TPUFC2_16S)
} else {
  cat("No data to plot for plastic_type_TPUFC2.1_vs_Elastollan (16S - Family).\n")
}

#Effect: plastic_conc
if (!is.null(sig_family_16S) && nrow(sig_family_16S[sig_family_16S$effect == "plastic_conc", ]) > 0) {
  plot_family_conc_16S <- ggplot(
    sig_family_16S[sig_family_16S$effect == "plastic_conc", ],
    aes(x = reorder(get("Family"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Family",
      y = "Log2 Fold Change",
      title = "Differential Abundance - plastic_conc (Family Level - 16S)"
    ) +
    theme_minimal()
  print(plot_family_conc_16S)
} else {
  cat("No data to plot for plastic_conc (16S - Family).\n")
}

#Effect: plastic_typeTPU181.plastic_conc
if (!is.null(sig_family_16S) && nrow(sig_family_16S[sig_family_16S$effect == "plastic_typeTPU181.plastic_conc", ]) > 0) {
  plot_family_TPU181_conc_16S <- ggplot(
    sig_family_16S[sig_family_16S$effect == "plastic_typeTPU181.plastic_conc", ],
    aes(x = reorder(get("Family"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Family",
      y = "Log2 Fold Change",
      title = "Differential Abundance - plastic_typeTPU181.plastic_conc (Family Level - 16S)"
    ) +
    theme_minimal()
  print(plot_family_TPU181_conc_16S)
} else {
  cat("No data to plot for plastic_typeTPU181.plastic_conc (16S - Family).\n")
}

#Effect: plastic_typeTPUFC2.1.plastic_conc
if (!is.null(sig_family_16S) && nrow(sig_family_16S[sig_family_16S$effect == "plastic_typeTPUFC2.1.plastic_conc", ]) > 0) {
  plot_family_TPUFC2_conc_16S <- ggplot(
    sig_family_16S[sig_family_16S$effect == "plastic_typeTPUFC2.1.plastic_conc", ],
    aes(x = reorder(get("Family"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Family",
      y = "Log2 Fold Change",
      title = "Differential Abundance - plastic_typeTPUFC2.1.plastic_conc (Family Level - 16S)"
    ) +
    theme_minimal()
  print(plot_family_TPUFC2_conc_16S)
} else {
  cat("No data to plot for plastic_typeTPUFC2.1.plastic_conc (16S - Family).\n")
}

#######################################
#####Filter low count features#########

#Define a threshold for filtering#
count_threshold <- 10  #Minimum count across all samples

#Filter significant results to remove low count ASVs
sig_family_16S_filtered <- sig_family_16S[apply(counts(dds_family_16S, normalized = TRUE)[sig_family_16S$ASV, ], 1, max) > count_threshold, ]

cat("\nFiltered significant results:\n")
print(head(sig_family_16S_filtered))

#### Replot bar charts with filtered results ####

#Effect: plastic_type_TPU181_vs_Elastollan
if (!is.null(sig_family_16S_filtered) && nrow(sig_family_16S_filtered[sig_family_16S_filtered$effect == "plastic_type_TPU181_vs_Elastollan", ]) > 0) {
  plot_family_TPU181_16S_filtered <- ggplot(
    sig_family_16S_filtered[sig_family_16S_filtered$effect == "plastic_type_TPU181_vs_Elastollan", ],
    aes(x = reorder(get("Family"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Family",
      y = "Log2 Fold Change",
      title = "Differential Abundance (Filtered) - plastic_type_TPU181_vs_Elastollan (Family Level - 16S)"
    ) +
    theme_minimal()
  print(plot_family_TPU181_16S_filtered)
} else {
  cat("No data to plot for plastic_type_TPU181_vs_Elastollan (Filtered - Family).\n")
}

#Effect: plastic_type_TPUFC2.1_vs_Elastollan
if (!is.null(sig_family_16S_filtered) && nrow(sig_family_16S_filtered[sig_family_16S_filtered$effect == "plastic_type_TPUFC2.1_vs_Elastollan", ]) > 0) {
  plot_family_TPUFC2_16S_filtered <- ggplot(
    sig_family_16S_filtered[sig_family_16S_filtered$effect == "plastic_type_TPUFC2.1_vs_Elastollan", ],
    aes(x = reorder(get("Family"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Family",
      y = "Log2 Fold Change",
      title = "Differential Abundance (Filtered) - plastic_type_TPUFC2.1_vs_Elastollan (Family Level - 16S)"
    ) +
    theme_minimal()
  print(plot_family_TPUFC2_16S_filtered)
} else {
  cat("No data to plot for plastic_type_TPUFC2.1_vs_Elastollan (Filtered - Family).\n")
}

#Effect: plastic_conc
if (!is.null(sig_family_16S_filtered) && nrow(sig_family_16S_filtered[sig_family_16S_filtered$effect == "plastic_conc", ]) > 0) {
  plot_family_conc_16S_filtered <- ggplot(
    sig_family_16S_filtered[sig_family_16S_filtered$effect == "plastic_conc", ],
    aes(x = reorder(get("Family"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Family",
      y = "Log2 Fold Change",
      title = "Differential Abundance (Filtered) - plastic_conc (Family Level - 16S)"
    ) +
    theme_minimal()
  print(plot_family_conc_16S_filtered)
} else {
  cat("No data to plot for plastic_conc (Filtered - Family).\n")
}

#Effect: plastic_typeTPU181.plastic_conc
if (!is.null(sig_family_16S_filtered) && nrow(sig_family_16S_filtered[sig_family_16S_filtered$effect == "plastic_typeTPU181.plastic_conc", ]) > 0) {
  plot_family_TPU181_conc_16S_filtered <- ggplot(
    sig_family_16S_filtered[sig_family_16S_filtered$effect == "plastic_typeTPU181.plastic_conc", ],
    aes(x = reorder(get("Family"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Family",
      y = "Log2 Fold Change",
      title = "Differential Abundance (Filtered) - plastic_typeTPU181.plastic_conc (Family Level - 16S)"
    ) +
    theme_minimal()
  print(plot_family_TPU181_conc_16S_filtered)
} else {
  cat("No data to plot for plastic_typeTPU181.plastic_conc (Filtered - Family).\n")
}

#Effect: plastic_typeTPUFC2.1.plastic_conc
if (!is.null(sig_family_16S_filtered) && nrow(sig_family_16S_filtered[sig_family_16S_filtered$effect == "plastic_typeTPUFC2.1.plastic_conc", ]) > 0) {
  plot_family_TPUFC2_conc_16S_filtered <- ggplot(
    sig_family_16S_filtered[sig_family_16S_filtered$effect == "plastic_typeTPUFC2.1.plastic_conc", ],
    aes(x = reorder(get("Family"), log2FoldChange), y = log2FoldChange, fill = padj < 0.05)
  ) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
    labs(
      x = "Family",
      y = "Log2 Fold Change",
      title = "Differential Abundance (Filtered) - plastic_typeTPUFC2.1.plastic_conc (Family Level - 16S)"
    ) +
    theme_minimal()
  print(plot_family_TPUFC2_conc_16S_filtered)
} else {
  cat("No data to plot for plastic_typeTPUFC2.1.plastic_conc (Filtered - Family).\n")
}









