#This R script will import the exports from my QIIME2 pipeline script so that they can be combined into a single Phyloseq object for further processing, object extraction, and statistical analysis. 
#Lines beginning with (#) were not in the final script

#Clear work space if need
rm(list = ls())

#check R version
R.version.string

#installing absent packages
#install.packages("microbiome")
#BiocManager::install("microbiome") #might work better
#install.packages("BiocManager")
BiocManager::install("microbiome")
#install.packages("/Users/lab/vegetarian", repos = NULL, type = "source")
#BiocManager::install("phyloseq", force = TRUE)
install.packages("gdata")
install.packages("ecodist")
install.packages("ape")
install.packages("phytools")
install.packages("castor")
install.packages("doParallel")
install.packages("viridis")
install.packages("ggrepel")


#load packages; install first if library is vacant
library(phyloseq); packageVersion("phyloseq") #1.48.0
library(ggplot2); packageVersion("ggplot2") #3.5.1
library(microbiome)
library(vegetarian) 
library(gdata)
library(ecodist)
library(vegan)  
library(carData)
library(car) 
library(dplyr)
library(biomformat)
library(ape) #phylogenetic tools packages
library(maps)
library(phytools) #phylogenetic tools packages
library(Rcpp)
library(castor)
library(doParallel) #used for UniFrac
library(viridis)
library(ggrepel)
library(tidyverse)

set.seed(26) #permutational analyses to be reproducible

#18S data import
#set working directory and replace the path
setwd("C:/Users/DELL/Documents/R/projects/MP_18S")

#read in biom table
ASV_reads <- read_biom("18S_Phyto-ASVtable.biom") #uses biomformat package
ASV_table <- as.data.frame(as.matrix(biom_data(ASV_reads))) #shape shifting
otu_tab<-t(ASV_table) #Transpose cols and rows
dim(otu_tab)  #Review
colnames(otu_tab) #Review
rownames(otu_tab) #Review
df.ASV<-as.data.frame(otu_tab)
df.ASV[1:10,1:10]
rownames(df.ASV) #Review
colnames(df.ASV) #Review

#Reading in meta-data file
metadata<-read.csv("18S_Phyto_metadata.csv",header=TRUE)
dim(metadata) #150 samples x 11 columns

#Read in taxonomy table
taxa2<-read.csv("PR2_TaxonomyFinal.csv")
ph.headers<-c("Feature.ID" , "Domain" , "Super_Groups","Super_Kingdom","Kingdom","Phylum" ,"Class"  ,"Family" , "Genus" , "species") #These should match .csv file
colnames(taxa2)<-ph.headers
rownames(taxa2) 
rownames(taxa2)<-taxa2[,1]
rownames(taxa2) 
dim(taxa2) #5,594 x 10
taxa3<-taxa2[,2:10]
colnames(taxa3) 

#Read in phylogenetic tree
dated.18Stree<-read.tree(file="18S_Phyto_tree.nwk") #reading in tree, uses ape package
is.rooted(dated.18Stree) #TRUE
sample_names(dated.18Stree) #NULL
dated.18Stree$tip.label #For example "bf08d62b32cced86e829cba893bdf318" 

#Create phyloseq object
str(taxa3) #data.frame 291 x 7
taxa3.m<-as.matrix(taxa3) #VERY NECESSARY TO DO, DON'T SKIP. 
str(taxa3.m)
colnames(taxa3.m)
str(df.ASV) #data.frame 150 obs. of 6987 vars
str(metadata) #data.frame 150 obs. of 4 vars
colnames(metadata)

rownames(df.ASV)<-as.character(rownames(df.ASV))
as.character(rownames(df.ASV))==as.character(metadata[,1])  
colnames(df.ASV) #accession numbers
rownames(metadata)<-as.character(metadata[,1])
rownames(taxa3.m)<-as.character(rownames(taxa3.m)) 
samp.names<-as.character(metadata[,1])

#To set up sample names to match (originally marked NULL)
sample_names(df.ASV)<-samp.names
sample_names(metadata)<-samp.names
sample_names(taxa3.m)<-samp.names
sample_names(dated.18Stree)<-samp.names

############NEW################
# Check if rownames of df.ASV match sample names
identical(rownames(df.ASV), samp.names)  # Should be TRUE
# Check if rownames of metadata match sample names
identical(rownames(metadata), samp.names)  # Should be TRUE
# Check if column names of df.ASV match the ASV IDs in taxa3
identical(colnames(df.ASV), rownames(taxa3.m))  # Should be TRUE
# Check if tree tip labels match ASV IDs
identical(sort(dated.18Stree$tip.label), sort(colnames(df.ASV)))  # Should be TRUE


#This is the actual phyloseq object
BioP.phylo<-phyloseq(otu_table(df.ASV, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxa3.m), phy_tree(dated.18Stree))

# If using a phyloseq object already:
identical(sample_names(BioP.phylo), samp.names)  # Should be TRUE
#################NEW##############

#Replace "taxa names" to short hand so that they are easier to view in R
taxa_names(BioP.phylo) <- paste0("ASV", seq(ntaxa(BioP.phylo)))
BioP.phylo@otu_table[1:10,1:10]
dim(BioP.phylo@otu_table) #150 x 5594
rowSums(BioP.phylo@otu_table) 
colSums(BioP.phylo@otu_table)

#Examining taxonomic ranks to examine where chloroplasts and mitochondria are nested within
#Need to examine Unassigned / Unknown taxa and do blasts to see if they TRULY are unknown or are chloroplasts/host DNA in disguise. 

colnames(BioP.phylo@tax_table) #"Domain", "Super_Groups", "Super_Kingdom", "Kingdom", "Phylum", "Class", "Family", "Genus", "species"  
table(tax_table(BioP.phylo)[, "Domain"], exclude = NULL) #20 Bacteria, 4776 Eukaryota, 18 Eukaryota:mito, 1 Eukaryota:plas, 779 Unassigned ###remove all but the 4776
table(tax_table(BioP.phylo)[, "Super_Groups"], exclude = NULL)
table(tax_table(BioP.phylo)[, "Super_Kingdom"], exclude = NULL)
table(tax_table(BioP.phylo)[, "Kingdom"], exclude = NULL) #20 Bacteria, 779 Unassigned #Remove these two
table(tax_table(BioP.phylo)[, "Phylum"], exclude = NULL) #Examine 
table(tax_table(BioP.phylo)[, "Class"], exclude = NULL) #Examine 
#table(tax_table(BioP.phylo)[, "Order"], exclude = NULL) #No orders in 18S data to examine
table(tax_table(BioP.phylo)[, "Family"], exclude = NULL) #Examine
table(tax_table(BioP.phylo)[, "Genus"], exclude = NULL) 
table(tax_table(BioP.phylo)[, "species"], exclude = NULL)

p0<-subset_taxa(BioP.phylo,  !Domain %in% "Bacteria") #This function removes
table(tax_table(p0)[, "Domain"], exclude = NULL) #This line checks to see if it is in fact... removed. #Works! 

p0b<-subset_taxa(p0,  !Domain %in% "Eukaryota:mito")
table(tax_table(p0b)[, "Domain"], exclude = NULL) #works

p0c<-subset_taxa(p0b,  !Domain %in% "Eukaryota:plas")
table(tax_table(p0c)[, "Domain"], exclude = NULL) #works

p0d<-subset_taxa(p0c,  !Domain %in% "Unassigned")
table(tax_table(p0d)[, "Domain"], exclude = NULL) #works

#Now that we have just Eukaryota data that are not some mitochondrial sequences or plasmid (if you want to see what these are you can re run the script; avoid renaming the accession numbers as ASVs and use the accesion number to access the sequence in the sequence.qzv file once imported into view.qiime.org (control F accession number to find sequence and use blastn to verify what organism that taxa belongs to)); examine lower taxonomic classifications such as Phylum, order, and family. You're going to have to dig and make sure that there are not taxa present here that are not relevaent to your hypotheses. 
p1<-subset_taxa(p0d,  !Super_Groups %in% "Unassigned")
table(tax_table(p1)[, "Super_Groups"], exclude = NULL) 
p2<-subset_taxa(p1,  !Super_Kingdom %in% "Unassigned")
table(tax_table(p2)[, "Super_Kingdom"], exclude = NULL) 
p3a<-subset_taxa(p2,  !Kingdom %in% "Unassigned")
p3b<-subset_taxa(p3a,  !Kingdom %in% "Fungi") #And other variants containing 'Fungi'
table(tax_table(p3b)[, "Kingdom"], exclude = NULL) 

p4<-subset_taxa(p3b,  !Phylum %in% "Unassigned") 
table(tax_table(p4)[, "Phylum"], exclude = NULL) #works
p5<-subset_taxa(p4,  !Class %in% "Unassigned")
table(tax_table(p5)[, "Class"], exclude = NULL) 
p6<-subset_taxa(p5,  !Family %in% "Unassigned")
table(tax_table(p6)[, "Family"], exclude = NULL)


p6@sam_data
rowSums(p6@otu_table) #18S_JRD492, 18S_JRD493, and 18S_JRD495 are extremely low read depth in comparison to others.


#Remove samples with extremely low read depth
BioPlastics_phylo <- subset_samples(p6, sample_name != "18S_JRD492" & sample_name != "18S_JRD493" & sample_name != "18S_JRD495")

#Update ASV table after removing samples
BioPlasticsASVtable <- BioPlastics_phylo@otu_table

#Check dimensions after removing samples
dim(BioPlasticsASVtable) # 147 x 1491

#Sort the column sums to check the read abundance
sort(colSums(BioPlasticsASVtable), decreasing = FALSE) # Check for any low abundance ASVs

#Remove ASVs with zero abundance (singleton filtering)
BioPlasticsASVtable <- BioPlasticsASVtable[, colSums(BioPlasticsASVtable) > 0]

#Update the phyloseq object with the filtered ASV table
BioPlastics_phylo@otu_table <- BioPlasticsASVtable

#Final check
dim(BioPlasticsASVtable) #147 x 1489


sort(colSums(p6@otu_table),decreasing = FALSE) #Some singletons present, remove these too. 
dim(BioPlastics_phylo@otu_table) #147 x 1489 
dim(BioPlasticsASVtable) #147 x 1489 
min(colSums(BioPlasticsASVtable)) #Minimum is 1, we are OK
sort(colSums(BioPlastics_phylo@otu_table),decreasing = FALSE) 

#############################################
############Rarefaction step#################
#############################################

#Run multiple rarefaction - load library if vacant
if (!requireNamespace("metagMisc", quietly = TRUE)) {
  devtools::install_github("vmikk/metagMisc")
}
library(metagMisc)

#Find rarefaction depth
min_depth <- min(sample_sums(BioPlastics_phylo))  
print(min_depth) #18373

#Be sure format is correct
taxa_are_rows(BioPlastics_phylo)
if (!taxa_are_rows(BioPlastics_phylo)) {
  otu_table(BioPlastics_phylo) <- t(otu_table(BioPlastics_phylo))
  taxa_are_rows(BioPlastics_phylo) <- TRUE
}

#Perform multiple rarefactions and average
set.seed(28)
averaged_phylo <- phyloseq_mult_raref_avg(
  BioPlastics_phylo,
  SampSize = min_depth, #18373 (minimum depth object = min_depth)
  iter = 100,
  parallel = TRUE
)

#Create rarefied ASV-level phyloseq object (for UniFrac)
UniFrac_Phylo_Object <- phyloseq(otu_table(averaged_phylo, taxa_are_rows=FALSE),
                                 phy_tree(BioPlastics_phylo@phy_tree),
                                 tax_table(BioPlastics_phylo@tax_table),
                                 sample_data(BioPlastics_phylo@sam_data))

#Fix tree if needed
if (!ape::is.binary(phy_tree(UniFrac_Phylo_Object))) {
  phy_tree(UniFrac_Phylo_Object) <- ape::multi2di(phy_tree(UniFrac_Phylo_Object))
}

#Aggregate from rarefied ASV object (a family level object and a phylum level object)
BioPlastics_family <- tax_glom(UniFrac_Phylo_Object, taxrank = "Family")
BioPlastics_phylum <- tax_glom(UniFrac_Phylo_Object, taxrank = "Phylum")

#Weighted UniFrac
registerDoParallel(cores=4)
set.seed(71)
wdistUni <- UniFrac(UniFrac_Phylo_Object, weighted=TRUE, parallel=TRUE, fast=TRUE)

#Bray-Curtis from ASV
set.seed(72)
bray_phylo <- transform_sample_counts(UniFrac_Phylo_Object, function(x) x / sum(x))
bray_dist <- phyloseq::distance(bray_phylo, method = "bray")

#Bray-Curtis from Family
set.seed(73)
bray_family_phylo <- transform_sample_counts(BioPlastics_family, function(x) x / sum(x))
bray_dist_family <- phyloseq::distance(bray_family_phylo, method = "bray")

#Bray-Curtis from Phylum
set.seed(74)
bray_phylum_phylo <- transform_sample_counts(BioPlastics_phylum, function(x) x / sum(x))
bray_dist_phylum <- phyloseq::distance(bray_phylum_phylo, method = "bray")

#UniFrac for Family and Phylum
set.seed(75)
wdistUni_family <- UniFrac(BioPlastics_family, weighted=TRUE, parallel=TRUE, fast=TRUE)
set.seed(76)
wdistUni_phylum <- UniFrac(BioPlastics_phylum, weighted=TRUE, parallel=TRUE, fast=TRUE)


#Create a copy of BioPlastics_phylo for plotting
BioPlastics_plot <- BioPlastics_phylo


#Explore data visually
#Transforming phyloseq data for ggplot2 - WAIT for this one to process
physeq_df <- BioPlastics_plot %>%
  transform_sample_counts(function(x) x / sum(x)) %>% #Normalize the counts
  psmelt() %>%
  filter(!is.na(Phylum)) #So only rows with defined phylum are considered

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
  facet_grid(.~plastic_type, space="free", scales="free") #Adjust the facetting based on metadata

print(phyplot)


#Aggregating data by phylum, plastic_type, and plastic_conc
physeq_df2 <- physeq_df %>%
  group_by(date, Phylum, plastic_type, plastic_conc) %>%
  summarise(Abundance = sum(Abundance), .groups = 'drop')

#Creating the stacked bar plot with flipped facet display
phyplot2 <- ggplot(physeq_df2, aes(x=date, y=Abundance, fill=Phylum)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(plastic_type ~ plastic_conc, scales="free_x", space="free") +  #Flipped facet orientation
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Sample Date", y = "Relative Abundance", title = "Microbial Diversity by Plastic Type and Concentration [g/L] Over Time")

#Print the plot
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
  #Make it so plastic concentration is shown with 3 decimal places
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
  date = factor(date, levels = c(0, 16, 36, 56, 94))  #Treat date as discrete
  )

#Create stacked bar plot with updated labels and plasma color scheme
phyplot3 <- ggplot(physeq_df3, aes(x = date, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(plastic_type ~ plastic_conc, scales = "free_x", space = "free") +  #Flipped facet orientation
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for date
  scale_fill_viridis(option = "plasma", discrete = TRUE) +  #Use plasma color scheme
  #Update theme to bold titles, remove grids, and add box around each facet
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  
    axis.title.y = element_text(face = "bold", size = 13),  
    axis.text.x = element_text(size = 9, angle = 90, hjust = 0.5),  #Center align x-axis tick labels
    axis.text.y = element_text(size = 9),  #Increase size of y-axis tick labels
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),  #Bold and center the title
    legend.title = element_text(face = "bold", size = 10),  #Bold legend title
    legend.text = element_text(size = 7),  #Decreased size of legend text
    legend.key.size = unit(0.6, "lines"),  #Adjust spacing between legend items
    legend.key.height = unit(0.4, "lines"),  #Reduce height of legend keys
    panel.background = element_blank(),  #Remove grey background
    panel.grid.major = element_blank(),  #Remove major gridlines
    panel.grid.minor = element_blank(),  #Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  #Box around each plot with correct linewidth argument
  ) +
  guides(fill = guide_legend(ncol = 2)) +  #Adjust number of columns in legend
  labs(x = "Day Number", y = "Relative Abundance", title = "Eukaryotic diversity")

#Print the plot
print(phyplot3)

#Get unique groups in the 18S dataset
unique_phylum <- unique(physeq_df$Phylum)  
num_phylum <- length(unique_phylum)  #Count the number of unique groups

#Print the unique groups and their count
print(unique_phylum)
print(paste("Number of unique groups in the 18S dataset:", num_phylum))

###Continue updating these figure###

#Custom color palette with 42 colors
custom_col42 <- c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3",
                  "#114578","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0",
                  "#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0",
                  "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4",
                  "#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098",
                  "#F7F7C5","#784511","#A55E18","#D2781E","#E4913F","#EAAB6C","#F0C498",
                  "#781122","#A5182F","#D21E2C","#E43F5B","#EA6C81","#F098A7", "black")

#custom color palette with 27 colors
custom_col27 <- c("#781156", "#A51876", "#D21E96", "#E43FAD", "#EA6CC0", "#F098D3",  #Deep to light pink/purple hues
                  "#114578", "#185EA5", "#1E78D2", "#3F91E4", "#6CABEA", "#98C4F0",  #Purples and blues for balance
                  "#D2D21E", "#E4E43F", "#EAEA6C", "#F0F098",                      #Yellowish-orange to provide contrast
                  "#A55E18", "#D2781E", "#E4913F", "#EAAB6C", "#F0C498",            #Orange shades from dark to light
                  "#A5182F", "#D21E2C", "#E43F5B", "#EA6C81", "#F098A7", "black")            #Reds to pinks for additional variation

#Generate 294 distinct colors using colorRampPalette()
custom_col66 <- colorRampPalette(c("red", "blue", "green", "orange", "purple", "brown", "pink", "yellow", "cyan", "black"))(66)

#Stacked bar plot with custom colors for 27 phyla
phyplot4 <- ggplot(physeq_df3, aes(x = date, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(plastic_type ~ plastic_conc, scales = "free_x", space = "free") +  #Flipped facet orientation
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for date
  scale_fill_manual(values = custom_col66) +  #Use custom color palette with 27 colors
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
  guides(fill = guide_legend(ncol = 2)) +  #Adjust number of columns in legend
  labs(x = "Day Number", y = "Relative Abundance")  #Removed the title

#Print the plot
print(phyplot4)

#Move legend
phyplot4 <- ggplot(physeq_df3, aes(x = date, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +  #Increase bar width
  facet_grid(plastic_type ~ plastic_conc, scales = "free_x", space = "free") +  # Flipped facet orientation
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for date
  scale_fill_manual(values = custom_col66) +  #Use custom color palette with 27 colors
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

#Print the plot
print(phyplot4)


#Aggregating data by phylum, plastic_type, and plastic_conc and including family
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
  date = factor(date, levels = c(0, 16, 36, 56, 94))  #Treat date as discrete
  )

#Get unique groups in the 18S dataset
unique_family <- unique(physeq_df$Family)  #Replace physeq_df with the appropriate data frame name if needed
num_family <- length(unique_family)  #Count the number of unique groups

#Print the unique groups and their count
print(unique_family)
print(paste("Number of unique groups in the 18S dataset:", num_family))

#Generate 179 distinct colors using colorRampPalette()
custom_col179 <- colorRampPalette(c("red", "blue", "green", "orange", "purple", "brown", "pink", "yellow", "cyan", "black"))(294)

#Create bar plot faceted by Phylum, with families as fill
phyplot_family <- ggplot(physeq_df4, aes(x = date, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +  #Stack families for each date
  facet_wrap(~ Phylum, scales = "free_y") +  #Create a separate facet for each phylum
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for date
  scale_fill_manual(values = custom_col179) +  #Use custom color palette for families
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
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  # Box around each plot with correct linewidth argument
  ) +
  labs(x = "Day Number", y = "Relative Abundance", title = "Family-level Abundance for Each Phylum")  # Add plot title

#Print the plot
print(phyplot_family)


#Get unique values of plastic types and concentrations
plastic_types <- unique(physeq_df4$plastic_type)
plastic_concs <- unique(physeq_df4$plastic_conc)

#DONT RUN, TOO LONG TO GENERERATE
#Loop through each plastic type and concentration to create individual plots
#for (ptype in plastic_types) {
  #for (pconc in plastic_concs) {
    #Subset the data for the current plastic type and concentration
    #subset_data <- physeq_df4 %>%
      #filter(plastic_type == ptype, plastic_conc == pconc)
    
    #Create the bar plot for the current plastic type and concentration
    #phyplot_family <- ggplot(subset_data, aes(x = date, y = Abundance, fill = Family)) +
      #geom_bar(stat = "identity", position = "stack", width = 0.8) +  #Stack families for each date
      #facet_wrap(~ Phylum, scales = "free_y") +  #Create a separate facet for each phylum
      #scale_x_discrete(name = "Day Number") +  # Use discrete scale for date
      #scale_fill_manual(values = custom_col179) +  #Use custom color palette with 179 colors
      #Update theme to bold titles, remove grids, and add box around each facet
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

#Bar plot faceted by tank (plastic_type and plasic_conc) with families as stacked bars
phyplot_top_families <- ggplot(top_families_df, aes(x = date, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  #Stack families within each phylum for each date
  facet_grid(plastic_type ~ plastic_conc) +  #Create a separate facet for each tank (plastic type and concentration)
  scale_x_discrete(name = "Sampling Date") +  #Use discrete scale for dates
  scale_fill_manual(values = colorRampPalette(c("red", "blue", "green", "orange", "purple", "yellow", "cyan", "pink", "brown", "grey"))(length(unique(top_families_df$Family)))) +  # Use custom color palette for families
  #Update theme to bold titles, remove grids, and add box around each facet
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5),  # Tilt x-axis labels to reduce overlap
    axis.text.y = element_text(size = 7),  #Increase size of y-axis tick labels
    legend.position = "none",
    panel.background = element_blank(),  #Remove grey background
    panel.grid.major = element_blank(),  #Remove major gridlines
    panel.grid.minor = element_blank(),  #Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  # Box around each plot with correct linewidth argument
  ) +
  labs(
    x = "Sampling Date",
    y = "Relative Abundance",
    title = "Family-level Abundance of Top 3 Phyla for Each Tank"
  )  #Add plot title

#Print the plot
print(phyplot_top_families)

#Seperate figure for each of the 3 phyla
#Identify the top 3 phyla across the entire dataset
top_phyla_overall <- physeq_df4 %>%
  group_by(Phylum) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = 'drop') %>%
  arrange(desc(TotalAbundance)) %>%
  slice_max(order_by = TotalAbundance, n = 3) %>%
  pull(Phylum)

#Print the names of the top 3 phyla
print("Top 3 phyla in the 18S dataset:")
print(top_phyla_overall)

#Create separate data subsets for each of the top 3 phyla
physeq_df_phylum1 <- physeq_df4 %>% filter(Phylum == top_phyla_overall[1])
physeq_df_phylum2 <- physeq_df4 %>% filter(Phylum == top_phyla_overall[2])
physeq_df_phylum3 <- physeq_df4 %>% filter(Phylum == top_phyla_overall[3])

#Create separate plots for each of the top 3 phyla, with legends
#Plot for Phylum 1
phyplot_phylum1_Chlorophyceae <- ggplot(physeq_df_phylum1, aes(x = date, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  #Stack families for each date
  facet_grid(plastic_type ~ plastic_conc) +  #Create a separate facet for each tank
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for dates
  scale_fill_manual(values = colorRampPalette(c("red", "blue", "green", "orange", "purple", "yellow", "cyan", "pink", "brown", "grey"))(length(unique(physeq_df_phylum1$Family)))) +  # Use custom color palette for families
  #Update theme to bold titles, remove grids, add box around each facet, and add legend back in
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5),  # Center-align x-axis tick labels
    axis.text.y = element_text(size = 7),  #Increase size of y-axis tick labels
    legend.position = "right",  # Move the legend below the plot
    legend.title = element_text(face = "bold", size = 10),  #Bold legend title
    legend.text = element_text(size = 6),  #Customize legend text size
    legend.key.size = unit(0.6, "lines"),  #Adjust spacing between legend items
    legend.key.height = unit(0.4, "lines"),  #Reduce height of legend keys
    panel.background = element_blank(),  #Remove grey background
    panel.grid.major = element_blank(),  #Remove major gridlines
    panel.grid.minor = element_blank(),  #Remove minor gridlines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)  # Box around each plot with correct linewidth argument
  ) +
  guides(fill = guide_legend(ncol = 1)) +  #Adjust the number of columns in the legend for better fit
  labs(
    x = "Day Number",
    y = "Relative Abundance",
    title = paste("Family-level Abundance of Top Phylum:", top_phyla_overall[1])
  )  #Add plot title

#Plot for Phylum 2
phyplot_phylum2_Arthropoda <- ggplot(physeq_df_phylum2, aes(x = date, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  #Stack families for each date
  facet_grid(plastic_type ~ plastic_conc) +  #Create a separate facet for each tank
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for dates
  scale_fill_manual(values = colorRampPalette(c("red", "blue", "green", "orange", "purple", "yellow", "cyan", "pink", "brown", "grey"))(length(unique(physeq_df_phylum2$Family)))) +  # Use custom color palette for families
  #Update theme to bold titles, remove grids, add box around each facet, and add legend back in
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5),  # Center-align x-axis tick labels
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
    title = paste("Family-level Abundance of Top Phylum:", top_phyla_overall[2])
  )  #Add plot title

#Plot for Phylum 3
phyplot_phylum3_Chrysophyceae <- ggplot(physeq_df_phylum3, aes(x = date, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  #Stack families for each date
  facet_grid(plastic_type ~ plastic_conc) +  #Create a separate facet for each tank
  scale_x_discrete(name = "Day Number") +  #Use discrete scale for dates
  scale_fill_manual(values = colorRampPalette(c("red", "blue", "green", "orange", "purple", "yellow", "cyan", "pink", "brown", "grey"))(length(unique(physeq_df_phylum3$Family)))) +  #Use custom color palette for families
  #Update theme to bold titles, remove grids, add box around each facet, and add legend back in
  theme(
    axis.title.x = element_text(face = "bold", size = 13),  #Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 13),  #Bold y-axis title
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5),  # Center-align x-axis tick labels
    axis.text.y = element_text(size = 7),  #Increase size of y-axis tick labels
    legend.position = "right",  # Move the legend below the plot
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
print(phyplot_phylum1_Chlorophyceae)
print(phyplot_phylum2_Arthropoda)
print(phyplot_phylum3_Chrysophyceae)



###################################################
###########PERMANOVA analysis below################
###################################################

#Prepare metadata
#meta_df <- as(sample_data(BioPlastics_phylo), "data.frame")
meta_df <- as(sample_data(UniFrac_Phylo_Object), "data.frame") #On rarefied data 
meta_df$plastic_type <- as.factor(meta_df$plastic_type)
meta_df$plastic_conc <- as.numeric(as.character(meta_df$plastic_conc))
meta_df$date <- as.factor(meta_df$date)
unique(meta_df$date)
meta_df$date <- factor(meta_df$date, levels = c("March_7", "March_23", "April_12", "May_3", "June_9"))
#meta_df$date <- factor(meta_df$date, levels = c("March 7", "March 23", "April 12", "May 3", "June 9"))
#cat("Checking overlap between meta_df and bray_dist...\n")
#print(all(rownames(meta_df) %in% rownames(as.matrix(bray_dist))))

#Run global PERMANOVAs (Bray-Curtis and Weighted UniFrac at ASV, Phylum, and Family levels)
set.seed(91)
adonis2_bray_asv <- adonis2(bray_dist ~ plastic_type * plastic_conc * date + chla, data = meta_df, permutations = 1000, by = "terms")
set.seed(92)
adonis2_unifrac_asv <- adonis2(wdistUni ~ plastic_type * plastic_conc * date + chla, data = meta_df, permutations = 1000, by = "terms")
set.seed(93)
adonis2_bray_family <- adonis2(bray_dist_family ~ plastic_type * plastic_conc * date + chla, data = meta_df, permutations = 1000, by = "terms")
set.seed(94)
adonis2_unifrac_family <- adonis2(wdistUni_family ~ plastic_type * plastic_conc * date + chla, data = meta_df, permutations = 1000, by = "terms")
set.seed(95)
adonis2_bray_phylum <- adonis2(bray_dist_phylum ~ plastic_type * plastic_conc * date + chla, data = meta_df, permutations = 1000, by = "terms")
set.seed(96)
adonis2_unifrac_phylum <- adonis2(wdistUni_phylum ~ plastic_type * plastic_conc * date + chla, data = meta_df, permutations = 1000, by = "terms")

#Print PERMANOVA results
cat("ASV Level - Bray-Curtis\n"); print(adonis2_bray_asv)
cat("ASV Level - Weighted UniFrac\n"); print(adonis2_unifrac_asv)
cat("Family Level - Bray-Curtis\n"); print(adonis2_bray_family)
cat("Family Level - Weighted UniFrac\n"); print(adonis2_unifrac_family)
cat("Phylum Level - Bray-Curtis\n"); print(adonis2_bray_phylum)
cat("Phylum Level - Weighted UniFrac\n"); print(adonis2_unifrac_phylum)

#Summarize and compare models
models <- list(
  ASV_Bray = adonis2_bray_asv,
  ASV_UniFrac = adonis2_unifrac_asv,
  Family_Bray = adonis2_bray_family,
  Family_UniFrac = adonis2_unifrac_family,
  Phylum_Bray = adonis2_bray_phylum,
  Phylum_UniFrac = adonis2_unifrac_phylum
)

#Pull out R2 and sig interactions for each model into a table
summarize_adonis2 <- function(model, name) {
  aov_df <- as.data.frame(model)
  aov_df$term <- rownames(aov_df)
  total_r2 <- sum(aov_df$R2[!(aov_df$term %in% c("Residual", "Total"))])
  n <- sum(aov_df$Df)
  p <- sum(!(aov_df$term %in% c("Residual", "Total")))
  adj_r2 <- 1 - ((1 - total_r2) * (n - 1) / (n - p - 1))
  significant_terms <- aov_df %>%
    filter(!(term %in% c("Residual", "Total")) & `Pr(>F)` <= 0.05) %>%
    pull(term) %>%
    paste(collapse = ", ")
  tibble(
    Model = name,
    Total_R2 = round(total_r2, 3),
    Adjusted_R2 = round(adj_r2, 3),
    Significant_Terms = ifelse(nchar(significant_terms) > 0, significant_terms, "None")
  )
}

#Create and print summary table
summary_table <- bind_rows(lapply(names(models), function(n) summarize_adonis2(models[[n]], n)))
print(summary_table)
View(summary_table)
#write.csv(summary_table, file = "18S_Distance_Matrix_comparisons_table_global_PERMANOVAs.csv", row.names = FALSE) #Write .csv into directory if desired


#Plot R2 summary
summary_long <- summary_table %>%
  pivot_longer(cols = c(Total_R2, Adjusted_R2), names_to = "Metric", values_to = "R2")

summary_long$Model <- factor(summary_long$Model, levels = summary_table$Model)

ggplot(summary_long, aes(x = Model, y = R2, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("#4B9CD3", "#FFA500")) +
  labs(title = "Total and Adjusted R2 Across PERMANOVA Models",
       x = "Model", y = "R2 Value", fill = "Metric") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#PERMDISP (homogeneity of dispersion) checks
dispersion_checks <- function(dist_obj, grouping_var, label) {
  disp <- betadisper(dist_obj, grouping_var)
  cat(paste0("\nPERMDISP for ", label, ":\n"))
  print(anova(disp))
  set.seed(26)  #Makes it reproducible
  print(permutest(disp))
  print(TukeyHSD(disp))
}

#Check for each model created
dispersion_checks(bray_dist, meta_df$plastic_type, "ASV - Bray-Curtis")
dispersion_checks(wdistUni, meta_df$plastic_type, "ASV - UniFrac")
dispersion_checks(bray_dist_family, meta_df$plastic_type, "Family - Bray-Curtis")
dispersion_checks(wdistUni_family, meta_df$plastic_type, "Family - UniFrac")
dispersion_checks(bray_dist_phylum, meta_df$plastic_type, "Phylum - Bray-Curtis")
dispersion_checks(wdistUni_phylum, meta_df$plastic_type, "Phylum - UniFrac")

####PERMANOVA for each date separately in 18S dataset####

#List of dates
dates <- c("March_7", "March_23", "April_12", "May_3", "June_9")
#dates <- c("March 7", "March 23", "April 12", "May 3", "June 9")

#Convert distance matrices to as.matrix for subsetting
bray_asv_matrix     <- as.matrix(bray_dist)
unifrac_asv_matrix  <- as.matrix(wdistUni)
bray_fam_matrix     <- as.matrix(bray_dist_family)
unifrac_fam_matrix  <- as.matrix(wdistUni_family)
bray_phy_matrix     <- as.matrix(bray_dist_phylum)
unifrac_phy_matrix  <- as.matrix(wdistUni_phylum)

#Initiate result storage
results_bray_asv        <- list()
results_unifrac_asv     <- list()
results_bray_family     <- list()
results_unifrac_family  <- list()
results_bray_phylum     <- list()
results_unifrac_phylum  <- list()

#Loop through each date
for (d in dates) {
  subset_df <- meta_df[meta_df$date == d, ]
  sample_ids <- rownames(subset_df)
  
  #Subset distance matrices
  subset_bray_asv     <- as.dist(bray_asv_matrix[sample_ids, sample_ids])
  subset_unifrac_asv  <- as.dist(unifrac_asv_matrix[sample_ids, sample_ids])
  subset_bray_fam     <- as.dist(bray_fam_matrix[sample_ids, sample_ids])
  subset_unifrac_fam  <- as.dist(unifrac_fam_matrix[sample_ids, sample_ids])
  subset_bray_phy     <- as.dist(bray_phy_matrix[sample_ids, sample_ids])
  subset_unifrac_phy  <- as.dist(unifrac_phy_matrix[sample_ids, sample_ids])
  
  #ASV Level
  set.seed(26 + which(dates == d))
  res_bray_asv <- adonis2(subset_bray_asv ~ plastic_type * plastic_conc + chla, data = subset_df, permutations = 1000, by = "terms")
  set.seed(27 + which(dates == d))
  res_unifrac_asv <- adonis2(subset_unifrac_asv ~ plastic_type * plastic_conc + chla, data = subset_df, permutations = 1000, by = "terms")
  
  #Family Level
  set.seed(28 + which(dates == d))
  res_bray_fam <- adonis2(subset_bray_fam ~ plastic_type * plastic_conc + chla, data = subset_df, permutations = 1000, by = "terms")
  set.seed(29 + which(dates == d))
  res_unifrac_fam <- adonis2(subset_unifrac_fam ~ plastic_type * plastic_conc + chla, data = subset_df, permutations = 1000, by = "terms")
  
  #Phylum Level
  set.seed(30 + which(dates == d))
  res_bray_phy <- adonis2(subset_bray_phy ~ plastic_type * plastic_conc + chla, data = subset_df, permutations = 1000, by = "terms")
  set.seed(31 + which(dates == d))
  res_unifrac_phy <- adonis2(subset_unifrac_phy ~ plastic_type * plastic_conc + chla, data = subset_df, permutations = 1000, by = "terms")
  
  #Store results
  results_bray_asv[[d]]       <- res_bray_asv
  results_unifrac_asv[[d]]    <- res_unifrac_asv
  results_bray_family[[d]]    <- res_bray_fam
  results_unifrac_family[[d]] <- res_unifrac_fam
  results_bray_phylum[[d]]    <- res_bray_phy
  results_unifrac_phylum[[d]] <- res_unifrac_phy
  
  #Print results for each date
  cat(paste0("\n--- ", d, " ---\n"))
  cat("ASV - Bray-Curtis:\n"); print(res_bray_asv)
  cat("ASV - UniFrac:\n"); print(res_unifrac_asv)
  cat("Family - Bray-Curtis:\n"); print(res_bray_fam)
  cat("Family - UniFrac:\n"); print(res_unifrac_fam)
  cat("Phylum - Bray-Curtis:\n"); print(res_bray_phy)
  cat("Phylum - UniFrac:\n"); print(res_unifrac_phy)
}


#Function to print results by category
print_results_by_type <- function(result_list, label) {
  cat(paste0("\n### ", label, " ###\n"))
  for (d in names(result_list)) {
    cat(paste0("\n--- ", d, " ---\n"))
    print(result_list[[d]])
  }
}

#Print all results
print_results_by_type(results_bray_asv,       "ASV - Bray-Curtis")
print_results_by_type(results_unifrac_asv,    "ASV - UniFrac")
print_results_by_type(results_bray_family,    "Family - Bray-Curtis")
print_results_by_type(results_unifrac_family, "Family - UniFrac")
print_results_by_type(results_bray_phylum,    "Phylum - Bray-Curtis")
print_results_by_type(results_unifrac_phylum, "Phylum - UniFrac")



################################
####18S PCoA Ordination Code####
################################

#Confirm metadata and phyloseq object alignment
all(rownames(meta_df) %in% sample_names(UniFrac_Phylo_Object))
all(rownames(meta_df) %in% sample_names(BioPlastics_family))
all(rownames(meta_df) %in% sample_names(BioPlastics_phylum))

#Add metadata to phyloseq objects
sample_data(UniFrac_Phylo_Object) <- sample_data(meta_df)
sample_data(BioPlastics_family)   <- sample_data(meta_df)
sample_data(BioPlastics_phylum)   <- sample_data(meta_df)

#Subset by date
asv_physeq_18S    <- subset_samples(UniFrac_Phylo_Object, date %in% c("March_7", "March_23", "April_12", "May_3", "June_9")) #ASV level
family_physeq_18S <- subset_samples(BioPlastics_family,   date %in% c("March_7", "March_23", "April_12", "May_3", "June_9")) #Family level
phylum_physeq_18S <- subset_samples(BioPlastics_phylum,   date %in% c("March_7", "March_23", "April_12", "May_3", "June_9")) #Phylum level

#Remove taxa with zero counts
asv_physeq_18S    <- prune_taxa(taxa_sums(asv_physeq_18S) > 0, asv_physeq_18S)
family_physeq_18S <- prune_taxa(taxa_sums(family_physeq_18S) > 0, family_physeq_18S)
phylum_physeq_18S <- prune_taxa(taxa_sums(phylum_physeq_18S) > 0, phylum_physeq_18S)

#Ordinate - PCoA
#Bray-Curtis
set.seed(26)
pcoa_bray_asv_18S    <- ordinate(physeq = asv_physeq_18S,    method = "PCoA", distance = "bray")
set.seed(27)
pcoa_bray_family_18S <- ordinate(physeq = family_physeq_18S, method = "PCoA", distance = "bray")
set.seed(28)
pcoa_bray_phylum_18S <- ordinate(physeq = phylum_physeq_18S, method = "PCoA", distance = "bray")

#Weighted Unifrac
set.seed(29)
pcoa_unifrac_asv_18S    <- ordinate(physeq = asv_physeq_18S,    method = "PCoA", distance = "wunifrac")
set.seed(30)
pcoa_unifrac_family_18S <- ordinate(physeq = family_physeq_18S, method = "PCoA", distance = "wunifrac")
set.seed(31)
pcoa_unifrac_phylum_18S <- ordinate(physeq = phylum_physeq_18S, method = "PCoA", distance = "wunifrac")

#Format date labels for display
meta_df$date <- gsub("_", " ", meta_df$date)
meta_df$date <- factor(meta_df$date, levels = c("March 7", "March 23", "April 12", "May 3", "June 9"))

#Define plot aesthetics
custom_shapes <- c("Elastollan" = 16, "TPU181" = 17, "TPUFC2.1" = 18)
custom_colors <- c("March 7" = "goldenrod2", "March 23" = "darkorchid3", 
                   "April 12" = "deeppink3", "May 3" = "darkorange3", "June 9" = "lightblue3")

#Plotting function
plot_pcoa_18S <- function(physeq, ordination, title) {
  sample_data(physeq) <- sample_data(meta_df)
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

#Generate plots for each distance metric
plot_bray_asv_18S    <- plot_pcoa_18S(asv_physeq_18S,    pcoa_bray_asv_18S,    "Bray-Curtis PCoA - ASV Level (18S)")
plot_bray_family_18S <- plot_pcoa_18S(family_physeq_18S, pcoa_bray_family_18S, "Bray-Curtis PCoA - Family Level (18S)")
plot_bray_phylum_18S <- plot_pcoa_18S(phylum_physeq_18S, pcoa_bray_phylum_18S, "Bray-Curtis PCoA - Phylum Level (18S)")

plot_unifrac_asv_18S    <- plot_pcoa_18S(asv_physeq_18S,    pcoa_unifrac_asv_18S,    "Weighted UniFrac PCoA - ASV Level (18S)")
plot_unifrac_family_18S <- plot_pcoa_18S(family_physeq_18S, pcoa_unifrac_family_18S, "Weighted UniFrac PCoA - Family Level (18S)")
plot_unifrac_phylum_18S <- plot_pcoa_18S(phylum_physeq_18S, pcoa_unifrac_phylum_18S, "Weighted UniFrac PCoA - Phylum Level (18S)")

#Display plots
print(plot_bray_asv_18S)
print(plot_bray_family_18S)
print(plot_bray_phylum_18S)

print(plot_unifrac_asv_18S)
print(plot_unifrac_family_18S)
print(plot_unifrac_phylum_18S)


#############################################################################
###Ok moving forward with Weighted UniFrac - ASV level - for further analysis
#############################################################################

###OK try to do a three panel by plastic type###
#Make sure plastic_type is a character in metadata
meta_df$plastic_type <- as.character(meta_df$plastic_type)
sample_data(UniFrac_Phylo_Object) <- sample_data(meta_df)

#Custom shape mapping for plastic_type
custom_shapes <- c("Elastollan" = 16, "TPU181" = 17, "TPUFC2.1" = 18)

#Custom color mapping for date
custom_colors <- c("March 7" = "goldenrod2", "March 23" = "darkorchid3", 
                   "April 12" = "deeppink3", "May 3" = "darkorange3", "June 9" = "lightblue3")

#Subset the phyloseq object by each plastic type
elastollan_phylo <- subset_samples(UniFrac_Phylo_Object, plastic_type == "Elastollan")
elastollan_phylo <- prune_taxa(taxa_sums(elastollan_phylo) > 0, elastollan_phylo)

tpu181_phylo <- subset_samples(UniFrac_Phylo_Object, plastic_type == "TPU181")
tpu181_phylo <- prune_taxa(taxa_sums(tpu181_phylo) > 0, tpu181_phylo)

tpufc2_phylo <- subset_samples(UniFrac_Phylo_Object, plastic_type == "TPUFC2.1")
tpufc2_phylo <- prune_taxa(taxa_sums(tpufc2_phylo) > 0, tpufc2_phylo)

#Function to create PCoA plot for a given subsetted phyloseq object and plastic type
create_pcoa_plot <- function(phylo_obj, plastic_type_label) {
  set.seed(26)
  #Perform PCoA using Weighted UniFrac distance
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

#Update the plots
#Define manual concentration labels directly
#Assuming concentrations are in the original phyloseq metadata as "plastic_conc"
manual_concentration_labels <- data.frame(
  plastic_conc = c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385),
  plastic_conc_label = as.character(0:9)  # Labels from 0 to 9
)

#Function to update PCoA plot by adding annotations without modifying meta_df
update_pcoa_plot_manual_labels <- function(pcoa_plot, phylo_obj, ordination) {
  # Extract metadata for the given phyloseq object
  local_meta_df <- as(sample_data(phylo_obj), "data.frame")
  local_meta_df$SampleID <- rownames(local_meta_df)  # Add SampleID to use for merging
  
  #Extract ordination coordinates as a data frame
  ordination_df <- as.data.frame(ordination$vectors[, 1:2])
  ordination_df$SampleID <- rownames(ordination_df)  # Add SampleID to use for merging
  colnames(ordination_df)[1:2] <- c("Axis.1", "Axis.2")  # Rename columns for easy plotting
  
  #Merge ordination data with metadata
  combined_df <- left_join(local_meta_df, ordination_df, by = "SampleID")
  
  #Merge with manual labels for plastic concentration
  combined_df <- left_join(combined_df, manual_concentration_labels, by = "plastic_conc")
  
  #Update the existing plot with annotations (no lines)
  updated_plot <- pcoa_plot +
    geom_point(data = combined_df, aes(x = Axis.1, y = Axis.2), size = 5, alpha = 0.8) + # Reduce size of points slightly
    geom_text(data = combined_df, aes(x = Axis.1, y = Axis.2, label = plastic_conc_label), size = 4, color = "black", vjust = 0.5, hjust = 0.5) + # Annotations with numbers 0-9
    theme(
      plot.title = element_blank(),   # Remove title
      legend.position = "none"        # Remove legends
    )
  
  return(updated_plot)
}

#Update each PCoA plot with manual labels
pcoa_plot_elastollan_updated <- update_pcoa_plot_manual_labels(pcoa_plot_elastollan, elastollan_phylo, ordinate(elastollan_phylo, method = "PCoA", distance = "wunifrac"))
pcoa_plot_tpu181_updated <- update_pcoa_plot_manual_labels(pcoa_plot_tpu181, tpu181_phylo, ordinate(tpu181_phylo, method = "PCoA", distance = "wunifrac"))
pcoa_plot_tpufc2_updated <- update_pcoa_plot_manual_labels(pcoa_plot_tpufc2, tpufc2_phylo, ordinate(tpufc2_phylo, method = "PCoA", distance = "wunifrac"))

#display the updated plots
print(pcoa_plot_elastollan_updated)
print(pcoa_plot_tpu181_updated)
print(pcoa_plot_tpufc2_updated)

#Ensure `custom_colors` and `manual_concentration_labels` are correctly defined
custom_colors <- c(
  "March 7" = "goldenrod2",   # Day 0
  "March 23" = "darkorchid3",  # Day 16
  "April 12" = "deeppink3",  # Day 36
  "May 3" = "darkorange3",     # Day 56
  "June 9" = "lightblue3"     # Day 94
)

manual_concentration_labels <- data.frame(
  plastic_conc = c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385),
  plastic_conc_label = as.character(0:9)  # Labels from 0 to 9
)

#Function to update PCoA plots with transparency and size based on plastic concentration
update_pcoa_plot_transparency <- function(pcoa_plot, phylo_obj, ordination) {
  #Extract metadata for the given phyloseq object
  local_meta_df <- as(sample_data(phylo_obj), "data.frame")
  local_meta_df$SampleID <- rownames(local_meta_df)  # Add SampleID to use for merging
  
  #Extract ordination coordinates as a data frame
  ordination_df <- as.data.frame(ordination$vectors[, 1:2])
  ordination_df$SampleID <- rownames(ordination_df)  # Add SampleID to use for merging
  colnames(ordination_df)[1:2] <- c("Axis.1", "Axis.2")  # Rename columns for easy plotting
  
  #Merge ordination data with metadata
  combined_df <- left_join(local_meta_df, ordination_df, by = "SampleID")
  
  #Merge with manual labels for plastic concentration
  combined_df <- left_join(combined_df, manual_concentration_labels, by = "plastic_conc")
  
  #Ensure plastic_conc_label is numeric
  combined_df$plastic_conc_label <- as.numeric(combined_df$plastic_conc_label)
  
  #Update the existing plot with transparency and size scaled by concentration
  updated_plot <- ggplot(combined_df, aes(x = Axis.1, y = Axis.2)) +
    geom_point(
      aes(
        size = plastic_conc_label,
        alpha = plastic_conc_label / max(plastic_conc_label, na.rm = TRUE),
        color = date,   # Retain coloring by date
        shape = plastic_type  # Retain shapes by plastic type
      ),
      show.legend = FALSE
    ) +
    geom_text(
      aes(label = plastic_conc_label),
      size = 2,  # Text size for annotations
      color = "black",
      vjust = 0.5,
      hjust = 0.5
    ) +
    scale_size_continuous(
      range = c(3, 10),  # Adjust point sizes based on concentration
      guide = "none"     # Suppress size from legend
    ) +
    scale_alpha_continuous(
      range = c(0.3, 1), # Make transparency differences more dramatic
      guide = "none"     # Suppress alpha from legend
    ) +
    scale_color_manual(
      values = custom_colors,
      breaks = names(custom_colors),
      labels = c("0", "16", "36", "56", "94")  # Replace labels with day numbers
    ) +
    scale_shape_manual(
      values = c("Elastollan" = 16, "TPU181" = 17, "TPUFC2.1" = 18)
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",  # Remove all legends
      plot.title = element_blank(),  # Optionally remove title
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1)
    )
  
  return(updated_plot)
}

#Rebuild each PCoA plot for 18S without any legends, adding transparency
set.seed(42)
pcoa_plot_elastollan_transparency <- update_pcoa_plot_transparency(
  pcoa_plot_elastollan,
  elastollan_phylo,
  ordinate(elastollan_phylo, method = "PCoA", distance = "wunifrac")
)

set.seed(43)
pcoa_plot_tpu181_transparency <- update_pcoa_plot_transparency(
  pcoa_plot_tpu181,
  tpu181_phylo,
  ordinate(tpu181_phylo, method = "PCoA", distance = "wunifrac")
)

set.seed(44)
pcoa_plot_tpufc2_transparency <- update_pcoa_plot_transparency(
  pcoa_plot_tpufc2,
  tpufc2_phylo,
  ordinate(tpufc2_phylo, method = "PCoA", distance = "wunifrac")
)

#Display the updated plots for 18S
print(pcoa_plot_elastollan_transparency)
print(pcoa_plot_tpu181_transparency)
print(pcoa_plot_tpufc2_transparency)

#####Put them on the same scale###
# Manually set x and y axis limits
x_limits <- c(-0.5, 0.25)
y_limits <- c(-0.5, 0.6)

# Rebuild each plot with fixed axis limits
pcoa_plot_elastollan_transparency <- pcoa_plot_elastollan_transparency +
  xlim(x_limits) +
  ylim(y_limits)

pcoa_plot_tpu181_transparency <- pcoa_plot_tpu181_transparency +
  xlim(x_limits) +
  ylim(y_limits)

pcoa_plot_tpufc2_transparency <- pcoa_plot_tpufc2_transparency +
  xlim(x_limits) +
  ylim(y_limits)

# Display the updated plots
print(pcoa_plot_elastollan_transparency)
print(pcoa_plot_tpu181_transparency)
print(pcoa_plot_tpufc2_transparency)




####################################################
#Ok now lets do each date separately with chla vector
####################################################

#Define manual concentration labels directly
unique(meta_df$date)

manual_concentration_labels <- data.frame(
  plastic_conc = c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385),
  plastic_conc_label = as.character(0:9)  # Labels from 0 to 9
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
  set.seed(42)
  # Perform PCoA using weighted UniFrac distance
  date_pcoa <- ordinate(physeq = phylo_obj, method = "PCoA", distance = "wunifrac")
  
  #Extract metadata for the given phyloseq object
  local_meta_df <- as(sample_data(phylo_obj), "data.frame")
  local_meta_df$SampleID <- rownames(local_meta_df)  # Add SampleID to use for merging
  
  #Extract ordination coordinates as a data frame from the PCoA object
  ordination_df <- as.data.frame(date_pcoa$vectors[, 1:2])
  ordination_df$SampleID <- rownames(ordination_df)  # Add SampleID to use for merging
  colnames(ordination_df)[1:2] <- c("Axis.1", "Axis.2")  # Rename columns for easy plotting
  
  #Merge ordination data with metadata for the current date
  combined_df <- left_join(local_meta_df, ordination_df, by = "SampleID")
  
  #Merge with manual labels for plastic concentration for the current date
  combined_df <- left_join(combined_df, manual_concentration_labels, by = "plastic_conc")
  
  #Fit the environmental variable (`chla`) to the ordination using `envfit`
  ord_scores <- ordination_df[, c("Axis.1", "Axis.2")]  # Extract the ordination scores for envfit
  envfit_result <- envfit(ord_scores, local_meta_df$chla, permutations = 1000)
  
  #Create PCoA plot for the current date
  pcoa_plot <- ggplot(combined_df, aes(x = Axis.1, y = Axis.2, color = plastic_type)) +
    geom_point(size = 5, alpha = 0.8) +  # Points with defined size
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

#Display the plots for each date
print(pcoa_plot_march_7)
print(pcoa_plot_march_23)
print(pcoa_plot_april_12)
print(pcoa_plot_may_3)
print(pcoa_plot_june_9)

#update these plots
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


###Insert simplified code
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
expand_factor <- 0.1  # Expands by 10% to ensure enough space for arrows
x_range <- x_max - x_min
y_range <- y_max - y_min

axis_limits <- list(
  xlim = c(x_min - expand_factor * x_range, x_max + expand_factor * x_range),
  ylim = c(y_min - expand_factor * y_range, y_max + expand_factor * y_range)
)

#Function to create combined data frame for a given phyloseq object
create_combined_df <- function(phylo_obj) {
  set.seed(26)
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

#Create `chla_arrow_data` for each date
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
    geom_text(aes(label = plastic_conc_numeric), size = 4, color = "black", vjust = 0.5, hjust = 0.5) +
    scale_shape_manual(values = shape_mapping) +
    scale_color_manual(values = color_mapping) +
    scale_alpha_continuous(range = c(0.3, 1)) +  # Ensure alpha is within a visible range
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

#Manually define labels
manual_concentration_labels <- data.frame(
  plastic_conc = c(0.000, 0.004, 0.008, 0.013, 0.023, 0.041, 0.072, 0.126, 0.220, 0.385),
  plastic_conc_numeric = 0:9  # Numeric values to be used for plotting size and alpha
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

#Create PCoA plots for each date
pcoa_plot_march_7 <- create_pcoa_plot(combined_df_march_7, chla_arrow_data_march_7, "March 7 PCoA Plot")
pcoa_plot_march_23 <- create_pcoa_plot(combined_df_march_23, chla_arrow_data_march_23, "March 23 PCoA Plot")
pcoa_plot_april_12 <- create_pcoa_plot(combined_df_april_12, chla_arrow_data_april_12, "April 12 PCoA Plot")
pcoa_plot_may_3 <- create_pcoa_plot(combined_df_may_3, chla_arrow_data_may_3, "May 3 PCoA Plot")
pcoa_plot_june_9 <- create_pcoa_plot(combined_df_june_9, chla_arrow_data_june_9, "June 9 PCoA Plot")

#Display the updated plots for each date
print(pcoa_plot_march_7)
print(pcoa_plot_march_23)
print(pcoa_plot_april_12)
print(pcoa_plot_may_3)
print(pcoa_plot_june_9)

#Remove axis figure titles from each plot and display them
pcoa_plot_march_7 <- pcoa_plot_march_7 + theme(plot.title = element_blank())
pcoa_plot_march_23 <- pcoa_plot_march_23 + theme(plot.title = element_blank())
pcoa_plot_april_12 <- pcoa_plot_april_12 + theme(plot.title = element_blank())
pcoa_plot_may_3 <- pcoa_plot_may_3 + theme(plot.title = element_blank())
pcoa_plot_june_9 <- pcoa_plot_june_9 + theme(plot.title = element_blank())

#Display the updated plots for each date
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

#Display the updated plots for each date
print(pcoa_plot_march_7_clean)
print(pcoa_plot_march_23_clean)
print(pcoa_plot_april_12_clean)
print(pcoa_plot_may_3_clean)
print(pcoa_plot_june_9_clean)

#Now make icons larger as well with palstic_conc
#Create `chla_arrow_data` for each date
chla_arrow_data_march_7_icon <- create_chla_arrow_data(combined_df_march_7[, c("Axis.1", "Axis.2")], combined_df_march_7$chla, axis_limits)
chla_arrow_data_march_237_icon <- create_chla_arrow_data(combined_df_march_23[, c("Axis.1", "Axis.2")], combined_df_march_23$chla, axis_limits)
chla_arrow_data_april_127_icon <- create_chla_arrow_data(combined_df_april_12[, c("Axis.1", "Axis.2")], combined_df_april_12$chla, axis_limits)
chla_arrow_data_may_37_icon <- create_chla_arrow_data(combined_df_may_3[, c("Axis.1", "Axis.2")], combined_df_may_3$chla, axis_limits)
chla_arrow_data_june_97_icon <- create_chla_arrow_data(combined_df_june_9[, c("Axis.1", "Axis.2")], combined_df_june_9$chla, axis_limits)

#Keep same shape and color assignments
shape_mapping <- c("Elastollan" = 16, "TPU181" = 17, "TPUFC2.1" = 18)
color_mapping <- c("Elastollan" = "darkorchid3", "TPU181" = "deeppink3", "TPUFC2.1" = "darkorange3")

#Function to create PCoA plot for a given combined_df and arrow data
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

#Display the updated plots for each date
print(pcoa_plot_march_7_icon)
print(pcoa_plot_march_23_icon)
print(pcoa_plot_april_12_icon)
print(pcoa_plot_may_3_icon)
print(pcoa_plot_june_9_icon)

#Remove axis titles/figure titles from each plot
pcoa_plot_march_7_icon <- pcoa_plot_march_7_icon + theme(axis.title = element_blank(), plot.title = element_blank())
pcoa_plot_march_23_icon <- pcoa_plot_march_23_icon + theme(axis.title = element_blank(), plot.title = element_blank())
pcoa_plot_april_12_icon <- pcoa_plot_april_12_icon + theme(axis.title = element_blank(), plot.title = element_blank())
pcoa_plot_may_3_icon <- pcoa_plot_may_3_icon + theme(axis.title = element_blank(), plot.title = element_blank())
pcoa_plot_june_9_icon <- pcoa_plot_june_9_icon + theme(axis.title = element_blank(), plot.title = element_blank())

#Display updated plots for each date
print(pcoa_plot_march_7_icon)
print(pcoa_plot_march_23_icon)
print(pcoa_plot_april_12_icon)
print(pcoa_plot_may_3_icon)
print(pcoa_plot_june_9_icon)






#Attach taxonomy table to Phyloseq object
UniFrac_Phylo_Object <- merge_phyloseq(UniFrac_Phylo_Object, tax_table(as.matrix(taxa3.m)))

#Bi Plots for 10 most abundant sig taxa
################ MARCH 7 ################
#Check levels of date
sample_data(UniFrac_Phylo_Object)$date <- as.character(sample_data(UniFrac_Phylo_Object)$date)
unique(sample_data(UniFrac_Phylo_Object)$date)
#Subset Phyloseq Object for March 7
march_7_phylo <- subset_samples(UniFrac_Phylo_Object, date == "March 7")

#Run PCoA ordination
set.seed(2)
ordination_march_7 <- ordinate(march_7_phylo, method = "PCoA", distance = "wunifrac")


#Extract ordination scores
ordination_scores_march_7 <- as.data.frame(ordination_march_7$vectors[, 1:2])
ordination_scores_march_7$SampleID <- rownames(ordination_scores_march_7)

#Extract metadata & merge with ordination scores
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


### ----- Process OTU + Taxonomy: FAMILY level -----
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
otu_table_march_7_family <- otu_table_march_7_family[
  !is.na(otu_table_march_7_family$Family) & otu_table_march_7_family$Family != "", ]
otu_table_march_7_family <- aggregate(. ~ Family, data = otu_table_march_7_family, FUN = sum)

rownames(otu_table_march_7_family) <- otu_table_march_7_family$Family
otu_table_march_7_family$Family <- NULL
otu_table_march_7_family <- as.data.frame(t(otu_table_march_7_family))

### ----- Process OTU + Taxonomy: GENUS level -----
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
otu_table_march_7_genus <- otu_table_march_7_genus[
  !is.na(otu_table_march_7_genus$Genus) & otu_table_march_7_genus$Genus != "", ]
otu_table_march_7_genus <- aggregate(. ~ Genus, data = otu_table_march_7_genus, FUN = sum)

rownames(otu_table_march_7_genus) <- otu_table_march_7_genus$Genus
otu_table_march_7_genus$Genus <- NULL
otu_table_march_7_genus <- as.data.frame(t(otu_table_march_7_genus))

### ----- ENVFIT & SIGNIFICANT TAXA -----

# Compute mean abundance
mean_abundance_families <- colMeans(otu_table_march_7_family)
mean_abundance_genera <- colMeans(otu_table_march_7_genus)

# Align samples
shared_samples_family <- intersect(rownames(ordination_scores_march_7), rownames(otu_table_march_7_family))
ordination_envfit_family <- ordination_scores_march_7[shared_samples_family, 1:2]
otu_envfit_family <- otu_table_march_7_family[shared_samples_family, ]

shared_samples_genus <- intersect(rownames(ordination_scores_march_7), rownames(otu_table_march_7_genus))
ordination_envfit_genus <- ordination_scores_march_7[shared_samples_genus, 1:2]
otu_envfit_genus <- otu_table_march_7_genus[shared_samples_genus, ]

# Filter out non-variable taxa
otu_envfit_family <- otu_envfit_family[, apply(otu_envfit_family, 2, function(x) var(x, na.rm = TRUE) > 0)]
otu_envfit_genus  <- otu_envfit_genus[,  apply(otu_envfit_genus, 2,  function(x) var(x, na.rm = TRUE) > 0)]

# Run envfit: FAMILY
set.seed(26)
envfit_march_7_family <- envfit(ordination_envfit_family, otu_envfit_family, permutations = 999)
sig_families_march_7 <- as.data.frame(envfit_march_7_family$vectors$arrows)
sig_families_march_7$p_value <- envfit_march_7_family$vectors$pvals
sig_families_march_7 <- sig_families_march_7[sig_families_march_7$p_value < 0.05, ]
sig_families_march_7$Family <- rownames(sig_families_march_7)
sig_families_march_7$mean_abundance <- mean_abundance_families[rownames(sig_families_march_7)]
sig_families_march_7 <- sig_families_march_7[order(-sig_families_march_7$mean_abundance), ]
print(sig_families_march_7)

# Run envfit: GENUS
set.seed(27)
envfit_march_7_genus <- envfit(ordination_envfit_genus, otu_envfit_genus, permutations = 999)
sig_genera_march_7 <- as.data.frame(envfit_march_7_genus$vectors$arrows)
sig_genera_march_7$p_value <- envfit_march_7_genus$vectors$pvals
sig_genera_march_7 <- sig_genera_march_7[sig_genera_march_7$p_value < 0.05, ]
sig_genera_march_7$Genus <- rownames(sig_genera_march_7)
sig_genera_march_7$mean_abundance <- mean_abundance_genera[rownames(sig_genera_march_7)]
sig_genera_march_7 <- sig_genera_march_7[order(-sig_genera_march_7$mean_abundance), ]
print(sig_genera_march_7)

#Select Top 10 families and genera for biplot
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

#Print plots
print(biplot_march_7_family)
print(biplot_march_7_genus)

################ MARCH 23 ################
#Subset phyloseq object for March 7
march_23_phylo <- subset_samples(UniFrac_Phylo_Object, date == "March 23")

#Run PCoA ordination
set.seed(3)
ordination_march_23 <- ordinate(march_23_phylo, method = "PCoA", distance = "wunifrac")


#Extract ordination scores
ordination_scores_march_23 <- as.data.frame(ordination_march_23$vectors[, 1:2])
ordination_scores_march_23$SampleID <- rownames(ordination_scores_march_23)

#Extract metadata & merge with ordination scores
meta_march_23 <- as(sample_data(march_23_phylo), "data.frame")
meta_march_23$SampleID <- rownames(meta_march_23)
ordination_metadata_march_23 <- left_join(meta_march_23, ordination_scores_march_23, by = "SampleID")

#Convert metadata to data frame  
ordination_metadata_march_23 <- as.data.frame(ordination_metadata_march_23)

#define plastic groups  
ordination_metadata_march_23$plastic_group <- ifelse(
  ordination_metadata_march_23$plastic_type %in% c("TPU181", "TPUFC2.1"), 
  "Bioplastics", 
  "Elastollan"
)
ordination_metadata_march_23$plastic_group <- factor(ordination_metadata_march_23$plastic_group, 
                                                    levels = c("Elastollan", "Bioplastics"))


### ----- Process OTU + Taxonomy: FAMILY level -----
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
otu_table_march_23_family <- otu_table_march_23_family[
  !is.na(otu_table_march_23_family$Family) & otu_table_march_23_family$Family != "", ]
otu_table_march_23_family <- aggregate(. ~ Family, data = otu_table_march_23_family, FUN = sum)

rownames(otu_table_march_23_family) <- otu_table_march_23_family$Family
otu_table_march_23_family$Family <- NULL
otu_table_march_23_family <- as.data.frame(t(otu_table_march_23_family))

### ----- Process OTU + Taxonomy: GENUS level -----
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
otu_table_march_23_genus <- otu_table_march_23_genus[
  !is.na(otu_table_march_23_genus$Genus) & otu_table_march_23_genus$Genus != "", ]
otu_table_march_23_genus <- aggregate(. ~ Genus, data = otu_table_march_23_genus, FUN = sum)

rownames(otu_table_march_23_genus) <- otu_table_march_23_genus$Genus
otu_table_march_23_genus$Genus <- NULL
otu_table_march_23_genus <- as.data.frame(t(otu_table_march_23_genus))

### ----- ENVFIT & SIGNIFICANT TAXA -----

# Compute mean abundance
mean_abundance_families <- colMeans(otu_table_march_23_family)
mean_abundance_genera <- colMeans(otu_table_march_23_genus)

# Align samples
shared_samples_family <- intersect(rownames(ordination_scores_march_23), rownames(otu_table_march_23_family))
ordination_envfit_family <- ordination_scores_march_23[shared_samples_family, 1:2]
otu_envfit_family <- otu_table_march_23_family[shared_samples_family, ]

shared_samples_genus <- intersect(rownames(ordination_scores_march_23), rownames(otu_table_march_23_genus))
ordination_envfit_genus <- ordination_scores_march_23[shared_samples_genus, 1:2]
otu_envfit_genus <- otu_table_march_23_genus[shared_samples_genus, ]

# Filter out non-variable taxa
otu_envfit_family <- otu_envfit_family[, apply(otu_envfit_family, 2, function(x) var(x, na.rm = TRUE) > 0)]
otu_envfit_genus  <- otu_envfit_genus[,  apply(otu_envfit_genus, 2,  function(x) var(x, na.rm = TRUE) > 0)]

# Run envfit: FAMILY
set.seed(30)
envfit_march_23_family <- envfit(ordination_envfit_family, otu_envfit_family, permutations = 999)
sig_families_march_23 <- as.data.frame(envfit_march_23_family$vectors$arrows)
sig_families_march_23$p_value <- envfit_march_23_family$vectors$pvals
sig_families_march_23 <- sig_families_march_23[sig_families_march_23$p_value < 0.05, ]
sig_families_march_23$Family <- rownames(sig_families_march_23)
sig_families_march_23$mean_abundance <- mean_abundance_families[rownames(sig_families_march_23)]
sig_families_march_23 <- sig_families_march_23[order(-sig_families_march_23$mean_abundance), ]
print(sig_families_march_23)

# Run envfit: GENUS
set.seed(31)
envfit_march_23_genus <- envfit(ordination_envfit_genus, otu_envfit_genus, permutations = 999)
sig_genera_march_23 <- as.data.frame(envfit_march_23_genus$vectors$arrows)
sig_genera_march_23$p_value <- envfit_march_23_genus$vectors$pvals
sig_genera_march_23 <- sig_genera_march_23[sig_genera_march_23$p_value < 0.05, ]
sig_genera_march_23$Genus <- rownames(sig_genera_march_23)
sig_genera_march_23$mean_abundance <- mean_abundance_genera[rownames(sig_genera_march_23)]
sig_genera_march_23 <- sig_genera_march_23[order(-sig_genera_march_23$mean_abundance), ]
print(sig_genera_march_23)

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
biplot_march_23_genus <- ggplot() +
  geom_point(data = ordination_metadata_march_23, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_march_23_genus, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_march_23_genus, aes(x = Axis.1, y = Axis.2, label = Genus), 
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
print(biplot_march_23_family)
print(biplot_march_23_genus)


################ APRIL 12 ################
#Subset phyloseq object for March 7
april_12_phylo <- subset_samples(UniFrac_Phylo_Object, date == "April 12")

#Run PCoA ordination
set.seed(44)
ordination_april_12 <- ordinate(april_12_phylo, method = "PCoA", distance = "wunifrac")


#Extract ordination scores
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


### ----- Process OTU + Taxonomy: FAMILY level -----
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
otu_table_april_12_family <- otu_table_april_12_family[
  !is.na(otu_table_april_12_family$Family) & otu_table_april_12_family$Family != "", ]
otu_table_april_12_family <- aggregate(. ~ Family, data = otu_table_april_12_family, FUN = sum)

rownames(otu_table_april_12_family) <- otu_table_april_12_family$Family
otu_table_april_12_family$Family <- NULL
otu_table_april_12_family <- as.data.frame(t(otu_table_april_12_family))

### ----- Process OTU + Taxonomy: GENUS level -----
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
otu_table_april_12_genus <- otu_table_april_12_genus[
  !is.na(otu_table_april_12_genus$Genus) & otu_table_april_12_genus$Genus != "", ]
otu_table_april_12_genus <- aggregate(. ~ Genus, data = otu_table_april_12_genus, FUN = sum)

rownames(otu_table_april_12_genus) <- otu_table_april_12_genus$Genus
otu_table_april_12_genus$Genus <- NULL
otu_table_april_12_genus <- as.data.frame(t(otu_table_april_12_genus))

### ----- ENVFIT & SIGNIFICANT TAXA -----

# Compute mean abundance
mean_abundance_families <- colMeans(otu_table_april_12_family)
mean_abundance_genera <- colMeans(otu_table_april_12_genus)

# Align samples
shared_samples_family <- intersect(rownames(ordination_scores_april_12), rownames(otu_table_april_12_family))
ordination_envfit_family <- ordination_scores_april_12[shared_samples_family, 1:2]
otu_envfit_family <- otu_table_april_12_family[shared_samples_family, ]

shared_samples_genus <- intersect(rownames(ordination_scores_april_12), rownames(otu_table_april_12_genus))
ordination_envfit_genus <- ordination_scores_april_12[shared_samples_genus, 1:2]
otu_envfit_genus <- otu_table_april_12_genus[shared_samples_genus, ]

# Filter out non-variable taxa
otu_envfit_family <- otu_envfit_family[, apply(otu_envfit_family, 2, function(x) var(x, na.rm = TRUE) > 0)]
otu_envfit_genus  <- otu_envfit_genus[,  apply(otu_envfit_genus, 2,  function(x) var(x, na.rm = TRUE) > 0)]

# Run envfit: FAMILY
set.seed(45)
envfit_april_12_family <- envfit(ordination_envfit_family, otu_envfit_family, permutations = 999)
sig_families_april_12 <- as.data.frame(envfit_april_12_family$vectors$arrows)
sig_families_april_12$p_value <- envfit_april_12_family$vectors$pvals
sig_families_april_12 <- sig_families_april_12[sig_families_april_12$p_value < 0.05, ]
sig_families_april_12$Family <- rownames(sig_families_april_12)
sig_families_april_12$mean_abundance <- mean_abundance_families[rownames(sig_families_april_12)]
sig_families_april_12 <- sig_families_april_12[order(-sig_families_april_12$mean_abundance), ]
print(sig_families_april_12)

# Run envfit: GENUS
set.seed(46)
envfit_april_12_genus <- envfit(ordination_envfit_genus, otu_envfit_genus, permutations = 999)
sig_genera_april_12 <- as.data.frame(envfit_april_12_genus$vectors$arrows)
sig_genera_april_12$p_value <- envfit_april_12_genus$vectors$pvals
sig_genera_april_12 <- sig_genera_april_12[sig_genera_april_12$p_value < 0.05, ]
sig_genera_april_12$Genus <- rownames(sig_genera_april_12)
sig_genera_april_12$mean_abundance <- mean_abundance_genera[rownames(sig_genera_april_12)]
sig_genera_april_12 <- sig_genera_april_12[order(-sig_genera_april_12$mean_abundance), ]
print(sig_genera_april_12)

#Select top 10 families and genera for biplot
top_10_families_april_12 <- head(sig_families_april_12, 10)
vectors_df_april_12_family <- data.frame(
  Axis.1 = top_10_families_april_12$Axis.1 * 0.2,  
  Axis.2 = top_10_families_april_12$Axis.2 * 0.2,
  Family = rownames(top_10_families_april_12)
)

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
biplot_april_12_genus <- ggplot() +
  geom_point(data = ordination_metadata_april_12, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_april_12_genus, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_april_12_genus, aes(x = Axis.1, y = Axis.2, label = Genus), 
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
print(biplot_april_12_family)
print(biplot_april_12_genus)


################ MAY 3 ################
#Subset phyloseq object for March 7
may_3_phylo <- subset_samples(UniFrac_Phylo_Object, date == "May 3")

#Run PCoA ordination
set.seed(47)
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


### ----- Process OTU + Taxonomy: FAMILY level -----
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
otu_table_may_3_family <- otu_table_may_3_family[
  !is.na(otu_table_may_3_family$Family) & otu_table_may_3_family$Family != "", ]
otu_table_may_3_family <- aggregate(. ~ Family, data = otu_table_may_3_family, FUN = sum)

rownames(otu_table_may_3_family) <- otu_table_may_3_family$Family
otu_table_may_3_family$Family <- NULL
otu_table_may_3_family <- as.data.frame(t(otu_table_may_3_family))

### ----- Process OTU + Taxonomy: GENUS level -----
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
otu_table_may_3_genus <- otu_table_may_3_genus[
  !is.na(otu_table_may_3_genus$Genus) & otu_table_may_3_genus$Genus != "", ]
otu_table_may_3_genus <- aggregate(. ~ Genus, data = otu_table_may_3_genus, FUN = sum)

rownames(otu_table_may_3_genus) <- otu_table_may_3_genus$Genus
otu_table_may_3_genus$Genus <- NULL
otu_table_may_3_genus <- as.data.frame(t(otu_table_may_3_genus))

### ----- ENVFIT & SIGNIFICANT TAXA -----

# Compute mean abundance
mean_abundance_families <- colMeans(otu_table_may_3_family)
mean_abundance_genera <- colMeans(otu_table_may_3_genus)

# Align samples
shared_samples_family <- intersect(rownames(ordination_scores_may_3), rownames(otu_table_may_3_family))
ordination_envfit_family <- ordination_scores_may_3[shared_samples_family, 1:2]
otu_envfit_family <- otu_table_may_3_family[shared_samples_family, ]

shared_samples_genus <- intersect(rownames(ordination_scores_may_3), rownames(otu_table_may_3_genus))
ordination_envfit_genus <- ordination_scores_may_3[shared_samples_genus, 1:2]
otu_envfit_genus <- otu_table_may_3_genus[shared_samples_genus, ]

# Filter out non-variable taxa
otu_envfit_family <- otu_envfit_family[, apply(otu_envfit_family, 2, function(x) var(x, na.rm = TRUE) > 0)]
otu_envfit_genus  <- otu_envfit_genus[,  apply(otu_envfit_genus, 2,  function(x) var(x, na.rm = TRUE) > 0)]

# Run envfit: FAMILY
set.seed(48)
envfit_may_3_family <- envfit(ordination_envfit_family, otu_envfit_family, permutations = 999)
sig_families_may_3 <- as.data.frame(envfit_may_3_family$vectors$arrows)
sig_families_may_3$p_value <- envfit_may_3_family$vectors$pvals
sig_families_may_3 <- sig_families_may_3[sig_families_may_3$p_value < 0.05, ]
sig_families_may_3$Family <- rownames(sig_families_may_3)
sig_families_may_3$mean_abundance <- mean_abundance_families[rownames(sig_families_may_3)]
sig_families_may_3 <- sig_families_may_3[order(-sig_families_may_3$mean_abundance), ]
print(sig_families_may_3)

# Run envfit: GENUS
set.seed(49)
envfit_may_3_genus <- envfit(ordination_envfit_genus, otu_envfit_genus, permutations = 999)
sig_genera_may_3 <- as.data.frame(envfit_may_3_genus$vectors$arrows)
sig_genera_may_3$p_value <- envfit_may_3_genus$vectors$pvals
sig_genera_may_3 <- sig_genera_may_3[sig_genera_may_3$p_value < 0.05, ]
sig_genera_may_3$Genus <- rownames(sig_genera_may_3)
sig_genera_may_3$mean_abundance <- mean_abundance_genera[rownames(sig_genera_may_3)]
sig_genera_may_3 <- sig_genera_may_3[order(-sig_genera_may_3$mean_abundance), ]
print(sig_genera_may_3)

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
biplot_may_3_genus <- ggplot() +
  geom_point(data = ordination_metadata_may_3, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_may_3_genus, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_may_3_genus, aes(x = Axis.1, y = Axis.2, label = Genus), 
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
print(biplot_may_3_family)
print(biplot_may_3_genus)

################JUNE 9################
#Subset phyloseq object for June 9
june_9_phylo <- subset_samples(UniFrac_Phylo_Object, date == "June 9")

#Run PCoA ordination
set.seed(50)
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


### ----- Process OTU + Taxonomy: FAMILY level -----
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
otu_table_june_9_family <- otu_table_june_9_family[
  !is.na(otu_table_june_9_family$Family) & otu_table_june_9_family$Family != "", ]
otu_table_june_9_family <- aggregate(. ~ Family, data = otu_table_june_9_family, FUN = sum)

rownames(otu_table_june_9_family) <- otu_table_june_9_family$Family
otu_table_june_9_family$Family <- NULL
otu_table_june_9_family <- as.data.frame(t(otu_table_june_9_family))

### ----- Process OTU + Taxonomy: GENUS level -----
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
otu_table_june_9_genus <- otu_table_june_9_genus[
  !is.na(otu_table_june_9_genus$Genus) & otu_table_june_9_genus$Genus != "", ]
otu_table_june_9_genus <- aggregate(. ~ Genus, data = otu_table_june_9_genus, FUN = sum)

rownames(otu_table_june_9_genus) <- otu_table_june_9_genus$Genus
otu_table_june_9_genus$Genus <- NULL
otu_table_june_9_genus <- as.data.frame(t(otu_table_june_9_genus))

### ----- ENVFIT & SIGNIFICANT TAXA -----

# Compute mean abundance
mean_abundance_families <- colMeans(otu_table_june_9_family)
mean_abundance_genera <- colMeans(otu_table_june_9_genus)

# Align samples
shared_samples_family <- intersect(rownames(ordination_scores_june_9), rownames(otu_table_june_9_family))
ordination_envfit_family <- ordination_scores_june_9[shared_samples_family, 1:2]
otu_envfit_family <- otu_table_june_9_family[shared_samples_family, ]

shared_samples_genus <- intersect(rownames(ordination_scores_june_9), rownames(otu_table_june_9_genus))
ordination_envfit_genus <- ordination_scores_june_9[shared_samples_genus, 1:2]
otu_envfit_genus <- otu_table_june_9_genus[shared_samples_genus, ]

# Filter out non-variable taxa
otu_envfit_family <- otu_envfit_family[, apply(otu_envfit_family, 2, function(x) var(x, na.rm = TRUE) > 0)]
otu_envfit_genus  <- otu_envfit_genus[,  apply(otu_envfit_genus, 2,  function(x) var(x, na.rm = TRUE) > 0)]

# Run envfit: FAMILY
set.seed(22)
envfit_june_9_family <- envfit(ordination_envfit_family, otu_envfit_family, permutations = 999)
sig_families_june_9 <- as.data.frame(envfit_june_9_family$vectors$arrows)
sig_families_june_9$p_value <- envfit_june_9_family$vectors$pvals
sig_families_june_9 <- sig_families_june_9[sig_families_june_9$p_value < 0.05, ]
sig_families_june_9$Family <- rownames(sig_families_june_9)
sig_families_june_9$mean_abundance <- mean_abundance_families[rownames(sig_families_june_9)]
sig_families_june_9 <- sig_families_june_9[order(-sig_families_june_9$mean_abundance), ]
print(sig_families_june_9)

# Run envfit: GENUS
set.seed(23)
envfit_june_9_genus <- envfit(ordination_envfit_genus, otu_envfit_genus, permutations = 999)
sig_genera_june_9 <- as.data.frame(envfit_june_9_genus$vectors$arrows)
sig_genera_june_9$p_value <- envfit_june_9_genus$vectors$pvals
sig_genera_june_9 <- sig_genera_june_9[sig_genera_june_9$p_value < 0.05, ]
sig_genera_june_9$Genus <- rownames(sig_genera_june_9)
sig_genera_june_9$mean_abundance <- mean_abundance_genera[rownames(sig_genera_june_9)]
sig_genera_june_9 <- sig_genera_june_9[order(-sig_genera_june_9$mean_abundance), ]
print(sig_genera_june_9)

#Select Top 10 families and genera for biplot
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

list(top_10_families_june_9)

#Family biplot
biplot_june_9_family <- ggplot() +
  geom_point(data = ordination_metadata_june_9, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_june_9_family, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_june_9_family, aes(x = Axis.1, y = Axis.2, label = Family), 
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
biplot_june_9_genus <- ggplot() +
  geom_point(data = ordination_metadata_june_9, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_june_9_genus, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_june_9_genus, aes(x = Axis.1, y = Axis.2, label = Genus), 
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
print(biplot_june_9_family)
print(biplot_june_9_genus)


##################################################################
###############For subset of plots to JUST look at ALGAE##########
##################################################################

###Define algal families (May 3 + June 9) ###CONFIRM IF USED#########
algal_families <- c("Chlamydomonadales_X", "Cyanidiales_X", "Chaetopeltidaceae", 
                    "Chrysochromulinaceae", "Microthamniales_X", "Euglenaceae", 
                    "Klebsormidiophyceae_XX", "Xanthophyceae_XX", "Prasiolales_X", 
                    "Hemiselmidaceae")


#Filter significant algae from ordination vectors May 3
algal_sig_families_may3 <- sig_families_may_3 %>% filter(Family %in% algal_families)

#Filter significant algae from ordination vectors June 9
algal_sig_families_june9 <- sig_families_june_9 %>% filter(Family %in% algal_families)

#Create vector data for algae in ordination biplot May 3
vectors_df_may_3_family_algae <- data.frame(
  Axis.1 = algal_sig_families_may3$Axis.1 * 0.2,  # Scale for visibility
  Axis.2 = algal_sig_families_may3$Axis.2 * 0.2,
  Family = algal_sig_families_may3$Family
)

#Create vector data for algae in ordination biplot June 9
vectors_df_june_9_family_algae <- data.frame(
  Axis.1 = algal_sig_families_june9$Axis.1 * 0.2,  
  Axis.2 = algal_sig_families_june9$Axis.2 * 0.2,
  Family = algal_sig_families_june9$Family
)

biplot_may_3_algae <- ggplot() +
  geom_point(data = ordination_metadata_may_3, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_may_3_family_algae, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_may_3_family_algae, aes(x = Axis.1, y = Axis.2, label = Family), 
                  color = "black", size = 3, max.overlaps = 15) +
  scale_color_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  labs(color = "Plastic type", title = "Algal Families in Ordination (May 3)") +
  theme_minimal()

biplot_june_9_algae <- ggplot() +
  geom_point(data = ordination_metadata_june_9, aes(x = Axis.1, y = Axis.2, color = plastic_group), size = 3) +
  geom_segment(data = vectors_df_june_9_family_algae, aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = vectors_df_june_9_family_algae, aes(x = Axis.1, y = Axis.2, label = Family), 
                  color = "black", size = 3, max.overlaps = 15) +
  scale_color_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  labs(color = "Plastic type", title = "Algal Families in Ordination (June 9)") +
  theme_minimal()

#Print new ordination plots for algae only
print(biplot_may_3_algae)
print(biplot_june_9_algae)

#####################################################
#####################################################

#Family Level
print(biplot_march_7_family)
print(biplot_march_23_family)
print(biplot_april_12_family)
print(biplot_may_3_family)
print(biplot_june_9_family)
print(biplot_may_3_algae)
print(biplot_june_9_algae)

#Genus Level
print(biplot_march_7_genus)
print(biplot_march_23_genus)
print(biplot_april_12_genus)
print(biplot_may_3_genus)
print(biplot_june_9_genus)


###Try a more effective function to keep labels out of the way
# Function to add family-level vectors and labels to an existing PCoA plot
add_family_vectors <- function(pcoa_plot, vectors_df) {
  pcoa_plot + 
    geom_segment(
      data = vectors_df, 
      aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
      arrow = arrow(length = unit(0.2, "cm")), 
      color = "black"
    ) +
    geom_text_repel(
      data = vectors_df, 
      aes(x = Axis.1 * 1.1, y = Axis.2 * 1.1, label = Family), 
      color = "black", size = 3,
      nudge_x = vectors_df$Axis.1 * 0.1,
      nudge_y = vectors_df$Axis.2 * 0.1,
      box.padding = 0.5,         # Add padding around the labels
      point.padding = 0.5,       # Add spacing between labels and plot points
      force = 2,                 # Stronger repulsion between overlapping labels
      force_pull = 0.1,          # Limit how far labels can be pulled from arrow tips
      max.overlaps = Inf,        # Don't drop any overlapping labels
      min.segment.length = 0.1   # Ensures the segment is drawn
    )
}


#Function to add genus level vectors to existing PCoA plot
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
      nudge_x = vectors_df$Axis.1 * 0.1,  # Shift text slightly away from arrows
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
pcoa_plot_may_3_icon_family_algae <- add_family_vectors(pcoa_plot_may_3_icon, vectors_df_may_3_family_algae)
pcoa_plot_june_9_icon_family_algae <- add_family_vectors(pcoa_plot_june_9_icon, vectors_df_june_9_family_algae)


#Display the updated plots, family
print(pcoa_plot_march_7_icon_family)
print(pcoa_plot_march_23_icon_family)
print(pcoa_plot_april_12_icon_family)
print(pcoa_plot_may_3_icon_family)
print(pcoa_plot_june_9_icon_family)
print(pcoa_plot_may_3_icon_family_algae)
print(pcoa_plot_june_9_icon_family_algae)

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

#Display the updated plots, genus
print(pcoa_plot_march_7_icon_genus)
print(pcoa_plot_march_23_icon_genus)
print(pcoa_plot_april_12_icon_genus)
print(pcoa_plot_may_3_icon_genus)
print(pcoa_plot_june_9_icon_genus)




############################
####DESEq2 test ############ 
############################
library('DESeq2')
library('indicspecies')
library('ggrepel')
library('tibble')

unique(sample_data(BioPlastics_phylo)$date)

#Rerun these subsets with the non rarefied data for DESeq2 
march_7_phylo <- subset_samples(BioPlastics_phylo, date == "March_7")
march_23_phylo <- subset_samples(BioPlastics_phylo, date == "March_23")
april_12_phylo <- subset_samples(BioPlastics_phylo, date == "April_12")
may_3_phylo <- subset_samples(BioPlastics_phylo, date == "May_3")
june_9_phylo <- subset_samples(BioPlastics_phylo, date == "June_9")

###########
#FOR JUNE 9
###########

#plastic_group exists in sample metadata
sample_data(june_9_phylo)$plastic_group <- ordination_metadata_june_9$plastic_group

#Filter out rare taxa (low read counts)
june_9_phylo_filtered <- prune_taxa(taxa_sums(june_9_phylo) > 10, june_9_phylo)

#Convert phyloseq object to DESeq2 format
dds <- phyloseq_to_deseq2(june_9_phylo_filtered, ~ plastic_group)

#Run DESeq2 analysis
dds <- DESeq(dds)

#Extract results (Bioplastics vs Elastollan)
res <- results(dds, contrast = c("plastic_group", "Bioplastics", "Elastollan"))


#Order results by significance (padj)
res <- res[order(res$padj, na.last = NA), ]

#Convert to a data.frame and remove NAs
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "ASV") %>%
  filter(!is.na(padj))

#Merge taxonomic information
tax_table_june_9 <- as.data.frame(tax_table(june_9_phylo)) %>%
  rownames_to_column(var = "ASV")  # Ensure ASVs are stored as a column

#Merge DESeq2 results with taxonomy
res_df <- left_join(res_df, tax_table_june_9, by = "ASV")

#Filter significant taxa (padj < 0.05)
sig_taxa <- res_df %>% filter(padj < 0.05)

#Print sig taxa
print(sig_taxa)

cat("Number of significantly different taxa:", nrow(sig_taxa), "\n")

sig_taxa <- as.data.frame(sig_taxa)  # Ensure it's a dataframe
str(sig_taxa)
sig_taxa %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Bioplastics", "Elastollan")) %>%
  group_by(Direction) %>%
  summarise(Count = n())

print(sig_taxa$Family)
sig_taxa %>% select(Family, log2FoldChange) %>% print()

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "DESeq2 Differential Abundance Analysis",
       x = "Log2 Fold Change (Bioplastics vs. Elastollan)",
       y = "-log10 Adjusted p-value") +
  theme_minimal()

top_taxa <- sig_taxa %>% arrange(padj) %>% head(15)  # Select top 15

june_9_log_fold_plot <- ggplot(top_taxa, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"), labels = c("Elastollan", "Bioplastics")) +
  labs(title = "Top Differentially Abundant Taxa June 9",
       x = "Family",
       y = "Log2 Fold Change (Bioplastics vs. Elastollan)",
       fill = "Higher in") +
  theme_minimal()

print(june_9_log_fold_plot)

####Update formatting####
june_9_log_fold_plot <- ggplot(top_taxa, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"), labels = c("Elastollan", "Bioplastics")) +
  labs(
    x = "Family",
    y = "Log2 Fold Change (Bioplastics vs. Elastollan)",
    fill = "Higher in"
  ) +
  theme_classic() +  # Use classic theme for better control
  theme(
    panel.grid = element_blank(),           # Remove grid lines
    plot.title = element_blank(),           # Remove the main title
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Add a box
  )

print(june_9_log_fold_plot)

#########################
####DESEq2 test FOR May 3
#########################

#Ensure plastic_group exists in sample metadata
sample_data(may_3_phylo)$plastic_group <- ordination_metadata_may_3$plastic_group

#Filter out rare taxa (low read counts)
may_3_phylo_filtered <- prune_taxa(taxa_sums(may_3_phylo) > 10, may_3_phylo)

#Convert phyloseq object to DESeq2 format
dds <- phyloseq_to_deseq2(may_3_phylo_filtered, ~ plastic_group)

#Run DESeq2 analysis
dds <- DESeq(dds)

#Extract results (Bioplastics vs. Elastollan)
res <- results(dds, contrast = c("plastic_group", "Bioplastics", "Elastollan"))

#Order results by significance (adjusted p-value)
res <- res[order(res$padj, na.last = NA), ]

#Convert to a data frame and remove NAs
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "ASV") %>%
  filter(!is.na(padj))

#Merge taxonomic information
tax_table_may_3 <- as.data.frame(tax_table(may_3_phylo)) %>%
  rownames_to_column(var = "ASV")  # Ensure ASVs are stored as a column

#Merge DESeq2 results with taxonomy
res_df <- left_join(res_df, tax_table_may_3, by = "ASV")

#Filter significant taxa (padj < 0.05)
sig_taxa <- res_df %>% filter(padj < 0.05)

#Print significant taxa
print(sig_taxa)

cat("Number of significantly different taxa:", nrow(sig_taxa), "\n")

sig_taxa <- as.data.frame(sig_taxa)  # Ensure it's a dataframe
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

may_3_log_fold_plot <- ggplot(top_taxa, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"), labels = c("Elastollan", "Bioplastics")) +
  labs(title = "Top Differentially Abundant Taxa May 3",
       x = "Family",
       y = "Log2 Fold Change (Bioplastics vs. Elastollan)",
       fill = "Higher in") +
  theme_minimal()

print(may_3_log_fold_plot)

###Update format
may_3_log_fold_plot <- ggplot(top_taxa, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"), labels = c("Elastollan", "Bioplastics")) +
  labs(
    x = "Family",
    y = "Log2 Fold Change (Bioplastics vs. Elastollan)",
    fill = "Higher in"
  ) +
  theme_classic() +  # Use classic theme for better control
  theme(
    panel.grid = element_blank(),           # Remove grid lines
    plot.title = element_blank(),           # Remove the main title
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Add a box
  )

print(may_3_log_fold_plot)





##########################################################################
###################### TOP 6 PLASTIC CONCENTRATIONS ######################
#########################   May 3 & June 9   #############################
##########################################################################

###Identify Top 6 plastic concentrations
#Select the top 6 highest plastic concentrations across both datasets
top_6_concentrations <- ordination_metadata_may_3 %>%
  arrange(desc(plastic_conc)) %>%
  distinct(plastic_conc) %>%
  head(6) %>%
  pull(plastic_conc)

#Subset phyloseq objects for May 3 and June 9
may_3_phylo_top6 <- subset_samples(may_3_phylo, plastic_conc %in% top_6_concentrations)
june_9_phylo_top6 <- subset_samples(june_9_phylo, plastic_conc %in% top_6_concentrations)

###Run DESeq2 for May 3 (top 6 plastic conc)**
#Chec if metadata matches the subset before assigning plastic_group
ordination_metadata_may_3_top6 <- ordination_metadata_may_3 %>% 
  filter(SampleID %in% sample_names(may_3_phylo_top6))  

#Assign plastic_group to the subset
sample_data(may_3_phylo_top6)$plastic_group <- ordination_metadata_may_3_top6$plastic_group  

#Filter out rare taxa
may_3_phylo_filtered_top6 <- prune_taxa(taxa_sums(may_3_phylo_top6) > 10, may_3_phylo_top6)

#Convert to DESeq2 format
dds_may3_top6 <- phyloseq_to_deseq2(may_3_phylo_filtered_top6, ~ plastic_group)

#Run DESeq2
dds_may3_top6 <- DESeq(dds_may3_top6)

#Extract results
res_may3_top6 <- results(dds_may3_top6, contrast = c("plastic_group", "Bioplastics", "Elastollan"))

#Filter significant taxa
sig_taxa_may3_top6 <- as.data.frame(res_may3_top6) %>%
  rownames_to_column(var = "ASV") %>%
  filter(!is.na(padj) & padj < 0.05)

#Merge taxonomy information
tax_table_may_3 <- as.data.frame(tax_table(may_3_phylo)) %>%
  rownames_to_column(var = "ASV")
sig_taxa_may3_top6 <- left_join(sig_taxa_may3_top6, tax_table_may_3, by = "ASV")

#Print summary
print(sig_taxa_may3_top6)

#Volcano plot
may_3_top6_volcano <- ggplot(sig_taxa_may3_top6, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "DESeq2 Analysis (Top 6 Plastic Concentrations - May 3)",
       x = "Log2 Fold Change (Bioplastics vs. Elastollan)",
       y = "-log10 Adjusted p-value") +
  theme_minimal()

print(may_3_top6_volcano)



###Run DESeq2 for June 9 (Top 6 Plastic Concentrations)###
#Metadata matches the subset before assigning plastic_group
ordination_metadata_june_9_top6 <- ordination_metadata_june_9 %>% 
  filter(SampleID %in% sample_names(june_9_phylo_top6))  

#Assign plastic_group to the subset
sample_data(june_9_phylo_top6)$plastic_group <- ordination_metadata_june_9_top6$plastic_group  

#Filter out rare taxa
june_9_phylo_filtered_top6 <- prune_taxa(taxa_sums(june_9_phylo_top6) > 10, june_9_phylo_top6)

#Convert to DESeq2 format
dds_june9_top6 <- phyloseq_to_deseq2(june_9_phylo_filtered_top6, ~ plastic_group)

#Run DESeq2
dds_june9_top6 <- DESeq(dds_june9_top6)

#Extract results
res_june9_top6 <- results(dds_june9_top6, contrast = c("plastic_group", "Bioplastics", "Elastollan"))

#Filter significant taxa
sig_taxa_june9_top6 <- as.data.frame(res_june9_top6) %>%
  rownames_to_column(var = "ASV") %>%
  filter(!is.na(padj) & padj < 0.05)

#Merge taxonomy information
tax_table_june_9 <- as.data.frame(tax_table(june_9_phylo)) %>%
  rownames_to_column(var = "ASV")
sig_taxa_june9_top6 <- left_join(sig_taxa_june9_top6, tax_table_june_9, by = "ASV")

#Print summary
print(sig_taxa_june9_top6)

#Volcano plot
june_9_top6_volcano <- ggplot(sig_taxa_june9_top6, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "DESeq2 Analysis (Top 6 Plastic Concentrations - June 9)",
       x = "Log2 Fold Change (Bioplastics vs. Elastollan)",
       y = "-log10 Adjusted p-value") +
  theme_minimal()

print(june_9_top6_volcano)

###Bar Plot for May 3 (Top 6 Plastic Concentrations)###
may_3_top6_bar_plot <- sig_taxa_may3_top6 %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Bioplastics", "Elastollan")) %>%
  ggplot(aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = Direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Elastollan" = "red", "Bioplastics" = "blue")) +
  labs(title = "Top Differentially Abundant Taxa (Top 6 Plastic Concentrations - May 3)",
       x = "Taxa (Family Level)",
       y = "Log2 Fold Change",
       fill = "Higher in") +
  theme_minimal()

print(may_3_top6_bar_plot)

#Update formatting
may_3_top6_bar_plot <- sig_taxa_may3_top6 %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Bioplastics", "Elastollan")) %>%
  ggplot(aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = Direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Elastollan" = "red", "Bioplastics" = "blue")) +
  labs(
    x = "Taxa (Family Level)",
    y = "Log2 Fold Change"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(may_3_top6_bar_plot)

#June 9 top 6 conc
june_9_top6_bar_plot <- sig_taxa_june9_top6 %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Bioplastics", "Elastollan")) %>%
  ggplot(aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = Direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Elastollan" = "red", "Bioplastics" = "blue")) +
  labs(title = "Top Differentially Abundant Taxa (Top 6 Plastic Concentrations - June 9)",
       x = "Taxa (Family Level)",
       y = "Log2 Fold Change",
       fill = "Higher in") +
  theme_minimal()

print(june_9_top6_bar_plot)

#Update format
june_9_top6_bar_plot <- sig_taxa_june9_top6 %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Bioplastics", "Elastollan")) %>%
  ggplot(aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = Direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Elastollan" = "red", "Bioplastics" = "blue")) +
  labs(
    x = "Taxa (Family Level)",
    y = "Log2 Fold Change"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(june_9_top6_bar_plot)

###Just list for copying###
#Extract and print taxa sig higher in Elastollan for May 3
elastollan_taxa_may3 <- sig_taxa_may3_top6 %>%
  filter(log2FoldChange < 0) %>%
  dplyr::select(Family)

print("Taxa significantly higher in Elastollan (May 3):")
print(elastollan_taxa_may3)

#Extract and print taxa sig higher in Elastollan for June 9
elastollan_taxa_june9 <- sig_taxa_june9_top6 %>%
  filter(log2FoldChange < 0) %>%
  dplyr::select(Family)

print("Taxa significantly higher in Elastollan (June 9):")
print(elastollan_taxa_june9)





###STOP#####STOP####STOP####STOP###
###################################
#### Not used exploration code ####
###################################
###STOP#####STOP####STOP####STOP###



#############################################################
#####WHICH TAXA are driving chla in elastollan tanks on june_9########

#Subset for petroleum plastic tanks (Elastollan)
petroleum_tanks <- subset_samples(june_9_phylo, plastic_group == "Elastollan")

#Extract OTU table and metadata
otu_table_petroleum <- as.data.frame(as(otu_table(petroleum_tanks), "matrix"))
meta_petroleum <- as(sample_data(petroleum_tanks), "data.frame")  # Convert to data frame

#Fix SampleID issue
meta_petroleum <- meta_petroleum %>% rownames_to_column(var = "SampleID")

#Merge chla Data with all OTU abundance, family level
chla_otu_data <- meta_petroleum %>%
  select(SampleID, chla) %>%
  inner_join(otu_table_june_9_family, by = "SampleID")

#Check if merge Worked
print(dim(chla_otu_data))  # Should return a valid number of rows and columns
print(head(chla_otu_data))  # Preview merged dataset

#Remove families with zero variance before correlation analysis
chla_otu_data_filtered <- chla_otu_data %>%
  select(-SampleID) %>%
  select(where(~ var(.x, na.rm = TRUE) > 0))  # Keep only families with variance

#Compute Spearman correlations b/w families and chla
chla_correlations <- chla_otu_data_filtered %>%
  summarise(across(everything(), ~ cor(.x, chla, method = "spearman", use = "pairwise.complete.obs"), 
                   .names = "cor_{.col}"))

#Check results
print(chla_correlations)

#Convert correlations to long format for visualization
chla_correlations_long <- pivot_longer(
  chla_correlations, 
  cols = everything(), 
  names_to = "Family", 
  values_to = "Correlation"
)

#Clean up family names (remove "cor_" prefix)
chla_correlations_long$Family <- gsub("cor_", "", chla_correlations_long$Family)

#Select top correlated families
top_chla_drivers <- chla_correlations_long %>%
  arrange(desc(abs(Correlation))) %>%  # Sort by absolute correlation
  filter(abs(Correlation) > 0.5)  # Threshold (adjust if needed)

#Visualize the top chla correlated families
ggplot(top_chla_drivers, aes(x = reorder(Family, Correlation), y = Correlation)) +
  geom_col(fill = ifelse(top_chla_drivers$Correlation > 0, "darkgreen", "darkred")) +
  coord_flip() +
  labs(title = "Families Correlated with Chl-a in Petroleum Tanks",
       x = "Family",
       y = "Spearman Correlation with Chl-a") +
  theme_minimal()

#Print and review
print(top_chla_drivers, n = 36)




#################################
## TRY algae in potrilium tanks##

#FAMILY Filter metadata for Elastollan (petroleum-based plastic)
petroleum_tanks <- ordination_metadata_june_9 %>%
  filter(plastic_group == "Elastollan")

#Check column names before calling distinct()
colnames(otu_table_june_9_family)

otu_table_june_9_family <- otu_table_june_9_family %>%
  rownames_to_column(var = "SampleID") %>%
  distinct(SampleID, .keep_all = TRUE)

#Ensure SampleID is unique
otu_table_june_9_family <- otu_table_june_9_family %>%
  distinct(SampleID, .keep_all = TRUE)

if (!"SampleID" %in% colnames(otu_table_june_9_family)) {
  stop("Error: 'SampleID' column is missing from `otu_table_june_9_family`.")
}

if (!"SampleID" %in% colnames(petroleum_tanks)) {
  stop("Error: 'SampleID' column is missing from `petroleum_tanks`.")
}

#Merge chla with family-level OTU abundance
chla_otu_data <- petroleum_tanks %>%
  select(SampleID, chla) %>%
  inner_join(otu_table_june_9_family, by = "SampleID")

print(chla_otu_data)

#Remove genera with zero variance
chla_otu_data_filtered <- chla_otu_data %>%
  select(-SampleID) %>%
  select(where(~ var(.x, na.rm = TRUE) > 0))

print(chla_otu_data_filtered)

#Recompute correlation
chla_correlations <- chla_otu_data_filtered %>%
  summarise(across(everything(), ~ cor(.x, chla, use = "pairwise.complete.obs"), .names = "cor_{.col}"))

print(chla_correlations)

#Convert to long format for easy visualization
chla_correlations_long <- pivot_longer(chla_correlations, cols = everything(), names_to = "Family", values_to = "Correlation")
print(chla_correlations_long)

#Filter for strong positive correlations
top_algal_drivers <- chla_correlations_long %>%
  arrange(desc(Correlation)) %>%
  filter(Correlation > 0.5)  # Adjust threshold as needed

print(top_algal_drivers)

#Extract significant genera from envfit()
sig_family <- sig_families_june_9$Family

#Find overlap with chla-correlated genera
key_algal_drivers <- top_algal_drivers %>%
  filter(Family %in% sig_family)

print(key_algal_drivers)

ggplot(top_algal_drivers, aes(x = reorder(Family, Correlation), y = Correlation)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(title = "Top Algal Genera Driving Chlorophyll-a in Petroleum Tanks",
       x = "Family",
       y = "Correlation with Chla") +
  theme_minimal()

print(top_algal_drivers)

print(nrow(key_algal_drivers))
print(key_algal_drivers)

biplot_june_9_family_algae <- biplot_june_9_family +
  geom_text(data = key_algal_drivers, aes(label = Family), color = "red", size = 4, fontface = "bold")

print(biplot_june_9_family_algae)

####################################################
#Family filter metadata for Elastollan (petroleum-based plastic)
petroleum_tanks <- ordination_metadata_june_9 %>%
  filter(plastic_group == "Elastollan")

#Ensure SampleID is unique
otu_table_june_9_genus <- otu_table_june_9_genus %>%
  distinct(SampleID, .keep_all = TRUE)

#Ensure SampleID exists in OTU table
if (!"SampleID" %in% colnames(otu_table_june_9_family)) {
  otu_table_june_9_family <- otu_table_june_9_family %>%
    rownames_to_column(var = "SampleID")
}

#Ensure SampleID formats match
petroleum_tanks$SampleID <- trimws(tolower(petroleum_tanks$SampleID))
otu_table_june_9_family$SampleID <- trimws(tolower(otu_table_june_9_family$SampleID))

#Check if SampleID values exist in both datasets
table(petroleum_tanks$SampleID %in% otu_table_june_9_family$SampleID)

#Merge chla with family-level OTU abundance
chla_otu_data_family <- petroleum_tanks %>%
  select(SampleID, chla) %>%
  inner_join(otu_table_june_9_family, by = "SampleID")


#Remove genera with zero variance
chla_otu_data_filtered <- chla_otu_data %>%
  select(-SampleID) %>%
  select(where(~ var(.x, na.rm = TRUE) > 0))

#Recompute correlation
chla_correlations <- chla_otu_data_filtered %>%
  summarise(across(everything(), ~ cor(.x, chla, use = "pairwise.complete.obs"), .names = "cor_{.col}"))


#Convert to long format for easy visualization
chla_correlations_long <- pivot_longer(chla_correlations, cols = everything(), names_to = "Family", values_to = "Correlation")

#Filter for strong positive correlations
top_algal_drivers <- chla_correlations_long %>%
  arrange(desc(Correlation)) %>%
  filter(Correlation > 0.5)  # Adjust threshold as needed

print(top_algal_drivers)

#Extract significant genera from envfit()
sig_genera <- sig_families_june_9$Family

#Find overlap with chla correlated genera
key_algal_drivers <- top_algal_drivers %>%
  filter(Family %in% sig_genera)

print(key_algal_drivers)

ggplot(top_algal_drivers, aes(x = reorder(Family, Correlation), y = Correlation)) +
  geom_col(fill = "seagreen3") +
  coord_flip() +
  labs(
    title = "Top Algal Families Driving Chlorophyll-a in Petroleum Tanks",
    subtitle = "Based on correlation with chla levels",
    x = "Family",
    y = "Correlation with Chlorophyll-a"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(face = "bold", color = "darkblue"))

print(top_algal_drivers)

biplot_june_9_family_algae <- biplot_june_9_family +
  geom_text(data = key_algal_drivers, aes(label = Family), color = "red", size = 4, fontface = "bold")

print(biplot_june_9_family_algae)


#######################################
#####Diversity Metrics Explore
# 1. Estimate Shannon diversity
alpha_df <- estimate_richness(BioPlastics_phylo, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)

# 2. Extract metadata and merge
meta_df <- as(sample_data(BioPlastics_phylo), "data.frame")
meta_df$SampleID <- rownames(meta_df)

# 3. Combine diversity metrics and metadata
meta_alpha <- left_join(alpha_df, meta_df, by = "SampleID")

# 4. Clean up date labels for display
meta_alpha$date <- gsub("_", " ", meta_alpha$date)
meta_alpha$date <- factor(meta_alpha$date, levels = c("March 7", "March 23", "April 12", "May 3", "June 9"))

ggplot(meta_alpha, aes(x = date, y = Shannon, fill = plastic_conc)) +
  geom_col(position = position_dodge2(preserve = "single"), color = "black") +
  facet_wrap(~plastic_type, nrow = 1) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Shannon Diversity per Tank", y = "Shannon Index", x = "Date", fill = "Plastic Conc (g/L)") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggplot(meta_alpha, aes(x = plastic_conc, y = Shannon, color = plastic_type)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~date) +
  labs(title = "Shannon Diversity vs. Plastic Concentration", y = "Shannon Index", x = "Plastic Concentration (g/L)") +
  theme_minimal(base_size = 13)






















































#####RICHNESS WITH PLASTIC TYPE
#Ensure the date column is in character format
sample_data(BioPlastics_phylo)$date <- as.character(sample_data(BioPlastics_phylo)$date)

process_18S_data <- function(phylo_obj, date) {
  #Subset phyloseq object
  phylo_subset <- subset_samples(phylo_obj, date == date)
  
  #Run PCoA ordination
  ordination <- ordinate(phylo_subset, method = "PCoA", distance = "wunifrac")
  
  #Extract ordination scores
  ordination_scores <- as.data.frame(ordination$vectors[, 1:2])
  ordination_scores$SampleID <- rownames(ordination_scores)
  
  #Extract metadata and merge with ordination scores
  meta <- as.data.frame(sample_data(phylo_subset))
  meta$SampleID <- rownames(meta)
  ordination_metadata <- left_join(meta, ordination_scores, by = "SampleID")
  
  #Process OTU and taxonomy tables, family level
  otu_table <- as.data.frame(as(otu_table(phylo_subset), "matrix"))
  tax_table <- as.data.frame(tax_table(phylo_subset))
  
  #Trim spaces
  tax_table[] <- lapply(tax_table, trimws)
  
  #Keep only ASVs present in both tables
  common_asvs <- intersect(rownames(tax_table), colnames(otu_table))
  tax_table <- tax_table[common_asvs, , drop = FALSE]
  otu_table <- otu_table[, common_asvs, drop = FALSE]
  
  #Transpose and aggregate by family
  otu_table <- as.data.frame(t(otu_table))
  otu_table$Family <- tax_table$Family
  otu_table <- otu_table[!is.na(otu_table$Family) & otu_table$Family != "", ]
  otu_table <- aggregate(. ~ Family, data = otu_table, FUN = sum)
  
  #Set family as row names and transpose back
  rownames(otu_table) <- otu_table$Family
  otu_table$Family <- NULL
  otu_table <- as.data.frame(t(otu_table))
  
  #Calculate richness (number of unique families)
  meta$richness <- rowSums(otu_table > 0)
  
  #Assign plastic type groups
  meta$plastic_group <- ifelse(meta$plastic_type %in% c("TPU181", "TPUFC2.1"), "Bioplastics", "Elastollan")
  meta$plastic_group <- factor(meta$plastic_group, levels = c("Elastollan", "Bioplastics"))
  
  #Return processed metadata with richness
  return(meta)
}

dates <- c("March_7", "March_23", "April_12", "May_3", "June_9")

#Run the function for each date and combine results
richness_results <- lapply(dates, function(d) process_18S_data(BioPlastics_phylo, d))
richness_df <- do.call(rbind, richness_results)

#View results
print(head(richness_df))

ggplot(richness_df, aes(x = plastic_group, y = richness, fill = plastic_group)) +
  geom_boxplot() +
  facet_wrap(~ date) +
  labs(y = "Algal Richness (Unique Families)", x = "Plastic Type") +
  scale_fill_manual(values = c("Bioplastics" = "green3", "Elastollan" = "darkorchid3")) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))







