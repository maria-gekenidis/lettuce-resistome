TEST
## CREDITS and many thanks go to Mahendra Mariadassou of whose scripts I took a
## lot of inspiration and commands. I also acknowledge the Genetic Diversity
## Centre (Zurich) and especially Jean-Claude Walser for their support in
## microbiome analysis.


####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  Get started (script NR.1)                 ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
### Project : p660
### Data Typ: MiSeq PE300 16S
### Run(s)  : 4 (run180504 / run180516 / run190117 / run190227)
### User(s) : Maria Stergiou (maria.stergiou@agroscope.admin.ch)
### === === === === === === === === === === === === === === === === ===
#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

# Define paths to use
idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)

## R and Bioconductor libraries 
library(ape)
library(DESeq2)
library(dplyr)
library(gapminder)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(labdsv)
library(microbiome)
library(phyloseq)
library(phyloseq.extended)
library(plotly)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(scales)
library(tidyverse)
library(vegan)


## Verify working directory and list all files
setwd("working-directory-of-your-choice")
getwd()
list.files(getwd())


#### Data import ####
otufile.ZOTU_c99      <- file.path(idir, "p660_run180504_run180516_run190117_run190227_16S_ZOTU_c99_Count_Sintax_Silva.txt")  # this file is the combined OTU and TAX
mapfile               <- file.path(idir, "SWL_sample_data.txt")                                                               # metadata
treefile.ZOTU_c99     <- file.path(idir, "p660_run180504_run180516_run190117_run190227_16S_ZOTU_c99_MSA.tre")
refseqfile.ZOTU_c99   <- file.path(idir, "p660_run180504_run180516_run190117_run190227_16S_ZOTU_c99.fa")                      # optional


## Import data into phyloseq
d.ZOTU_c99 <- import_qiime(otufilename = otufile.ZOTU_c99, mapfilename = mapfile, treefilename = treefile.ZOTU_c99)
d.ZOTU_c99
d.ZOTU99 <- d.ZOTU_c99


#### Filtering data according to your needs ####
d.ZOTU99 <- subset_samples(d.ZOTU99, ID != "L5DIIb")
d.ZOTU99 <- prune_taxa(taxa_sums(d.ZOTU99) > 0, d.ZOTU99)
d.ZOTU99

# Filter out chloroplasts and mitochondria
d.ZOTU99 <- subset_taxa(d.ZOTU99,  (Class != "Chloroplast") | is.na(Class))
d.ZOTU99 <- subset_taxa(d.ZOTU99,  (Family != "Mitochondria") | is.na(Family))
d.ZOTU99


#### Data adjustments by Maria Gekenidis-Stergiou ####

## At this point it is useful for further processing to write "unknown" into all
## empty (<NA>) cells of the taxonomy table and then return it to your phyloseq
## object.
tax.clean <- data.frame(tax_table(d.ZOTU99))

levels(tax.clean$Kingdom) <- c(levels(tax.clean$Kingdom), "unknown")    # add level "unknown" to Kingdom-levels
tax.clean$Kingdom[is.na(tax.clean$Kingdom)] <- "unknown"                # insert "unknown" into all NA cells

levels(tax.clean$Phylum) <- c(levels(tax.clean$Phylum), "unknown")
tax.clean$Phylum[is.na(tax.clean$Phylum)] <- "unknown"

levels(tax.clean$Class) <- c(levels(tax.clean$Class), "unknown")
tax.clean$Class[is.na(tax.clean$Class)] <- "unknown"

levels(tax.clean$Order) <- c(levels(tax.clean$Order), "unknown")
tax.clean$Order[is.na(tax.clean$Order)] <- "unknown"

levels(tax.clean$Family) <- c(levels(tax.clean$Family), "unknown")
tax.clean$Family[is.na(tax.clean$Family)] <- "unknown"

levels(tax.clean$Genus) <- c(levels(tax.clean$Genus), "unknown")
tax.clean$Genus[is.na(tax.clean$Genus)] <- "unknown"

levels(tax.clean$Species) <- c(levels(tax.clean$Species), "unknown")
tax.clean$Species[is.na(tax.clean$Species)] <- "unknown"

tax_table(d.ZOTU99) <- as.matrix(tax.clean)


# Filter out unassigned Bacteria and Archaea
d.ZOTU99 <- subset_taxa(d.ZOTU99,  (Phylum != "unknown"))
d.ZOTU99


## Change "Time" in the metadata table from integer to factor
class(d.ZOTU99@sam_data[["Time"]])
d.ZOTU99@sam_data[["Time"]] <- factor(d.ZOTU99@sam_data[["Time"]])
class(d.ZOTU99@sam_data[["Time"]])


###
#### Raw data exploration ####
###

## Count Range
range(sample_sums(d.ZOTU99))

## Miscellaneous information about my phyloseq object
summarize_phyloseq(d.ZOTU99)
ntaxa(d.ZOTU99)
nsamples(d.ZOTU99)
taxa_names(d.ZOTU99)
sample_names(d.ZOTU99)
taxa_sums(d.ZOTU99)
sample_sums(d.ZOTU99)
sample_sums(d.ZOTU99)[c("L0ZI", "L0ZII", "S0ZI", "S0ZII")]
rank_names(d.ZOTU99)
sample_variables(d.ZOTU99)
get_taxa(d.ZOTU99, "S0ZI")                    # count table
get_sample(d.ZOTU99, "ZOTU50")             # count table
get_variable(d.ZOTU99, c("Time", "System"))   # meta info


## Sequencing Depth Plot
plot(sample_sums(d.ZOTU99), xaxt = "n", xlab = "",
     ylab = "Fragment-Length Counts",
     pch = 19,
     cex = c(sample_sums(d.ZOTU99)/mean(sample_sums(d.ZOTU99))),
     col = rgb(0.5,0.5,0.5,alpha = 0.3),
     main = "Raw Counts per Sample / Sequencing Depth",
     ylim = range(0,sample_sums(d.ZOTU99)*1.01))
mycolors <- rep("gray",length(sample_sums(d.ZOTU99)))
mycolors[get_variable(d.ZOTU99, "System") == "L"] <- "green"
mycolors[get_variable(d.ZOTU99, "System") == "S"] <- "red"
mycolors[get_variable(d.ZOTU99, "System") == "W"] <- "blue"
mycolors[get_variable(d.ZOTU99, "System") == "M"] <- "orange"
points(sample_sums(d.ZOTU99), pch = 3, col = mycolors, cex = 1.5)
v <- sample_names(d.ZOTU99)
axis(side = 1, at = seq(1,length(sample_sums(d.ZOTU99))), labels = v, tck = -.02, las = 2, col.axis = 1)
abline(h = mean(sample_sums(d.ZOTU99)), col = "lightgreen")
abline(h = median(sample_sums(d.ZOTU99)), col = "lightblue")
legend(120, 620000, legend = c("L","S","W","M"),
       col = c("green", "red", "blue","orange"), lty = 1, cex = 0.5, lwd = 2)


## Alpha-Diversity-Plots
plot_richness(d.ZOTU99, x = "System", measures = c("Observed","Chao1", "Shannon"), color = "System") +
  geom_point(size = 5, alpha = 0.7)


###
#### Subsetting for later analyses ####
###

## Simplifying / subsetting your phyloseq object
d.ZOTU99_S <- subset_samples(d.ZOTU99, System == "S")
d.ZOTU99_L <- subset_samples(d.ZOTU99, System == "L")
d.ZOTU99_W <- subset_samples(d.ZOTU99, System == "W")
d.ZOTU99_M <- subset_samples(d.ZOTU99, System == "M")

## After subsetting, there might be OTUs without counts.
# Number of OTUs without counts:
sum(taxa_sums(d.ZOTU99_S) == 0)
sum(taxa_sums(d.ZOTU99_L) == 0)
sum(taxa_sums(d.ZOTU99_W) == 0)
sum(taxa_sums(d.ZOTU99_M) == 0)

# Or percentage of OTUs without counts:
100/(sum(taxa_sums(d.ZOTU99_S) == 0) + sum(taxa_sums(d.ZOTU99_S) > 0))*sum(taxa_sums(d.ZOTU99_S) == 0)            # 58%
100/(sum(taxa_sums(d.ZOTU99_L) == 0) + sum(taxa_sums(d.ZOTU99_L) > 0))*sum(taxa_sums(d.ZOTU99_L) == 0)            # 48%
100/(sum(taxa_sums(d.ZOTU99_W) == 0) + sum(taxa_sums(d.ZOTU99_W) > 0))*sum(taxa_sums(d.ZOTU99_W) == 0)            # 34%
100/(sum(taxa_sums(d.ZOTU99_M) == 0) + sum(taxa_sums(d.ZOTU99_M) > 0))*sum(taxa_sums(d.ZOTU99_M) == 0)            # 82%

# Remove these OTUs from each subset:
d.ZOTU99_S_no0 <- prune_taxa(taxa_sums(d.ZOTU99_S) > 0, d.ZOTU99_S)
d.ZOTU99_L_no0 <- prune_taxa(taxa_sums(d.ZOTU99_L) > 0, d.ZOTU99_L)
d.ZOTU99_W_no0 <- prune_taxa(taxa_sums(d.ZOTU99_W) > 0, d.ZOTU99_W)
d.ZOTU99_M_no0 <- prune_taxa(taxa_sums(d.ZOTU99_M) > 0, d.ZOTU99_M)


## Alpha-diversity for single systems
plot_richness(d.ZOTU99_L_no0, x = "Time",
              measures = c("Observed", "Shannon", "Simpson"),
              color = "Treatment")
plot_richness(d.ZOTU99_S_no0, x = "Time",
              measures = c("Observed", "Shannon", "Simpson"),
              color = "Treatment")
plot_richness(d.ZOTU99_W_no0, x = "Time",
              measures = c("Observed", "Shannon", "Simpson"),
              color = "Treatment")


## Export filtered and subsetted phyloseq objects for downstream scripts
RData.filename <- "GDC_ZOTU99-Silva_1.RData"
output <- file.path(odir, RData.filename)

save(d.ZOTU99,
     d.ZOTU99_S_no0,
     d.ZOTU99_L_no0,
     d.ZOTU99_W_no0,
     d.ZOTU99_M_no0,
     file = output)



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  RAREFY/TRANSFORM COUNTS (script NR.2)     ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

# Define paths to use
idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
RData.filename <- "GDC_ZOTU99-Silva_1.RData"
input <- file.path(odir, RData.filename)

load(input)

rm(RData.filename, input)


#### Rarefaction of phyloseq objects
### Full physeq object ####
## Quick approach to determine optimal rarefaction level is to look at sample_sums
max(sample_sums(d.ZOTU99))
min(sample_sums(d.ZOTU99))
sort(sample_sums(d.ZOTU99))    # to see which are the lowest depths and pick a cutoff


## Rarefy at a certain level.
# Choose a depth from the above output which is reasonable (not loosing too many
# samples nor depth):
rare_cutoff <- 22470
d.ZOTU99.rare_custom <- rarefy_even_depth(d.ZOTU99, sample.size = rare_cutoff, rngseed = 20200228)
sample_sums(d.ZOTU99.rare_custom)[1:5]


## Without defining a cutoff for rarefaction, all samples will be rarefied to
## the lowest depth sample. In the present data, two L6 samples have much lower
## depth than the rest > remove t6.
d.ZOTU99_not6 <- subset_samples(d.ZOTU99, Time != "6")
d.ZOTU99_not6 <- prune_taxa(taxa_sums(d.ZOTU99_not6) > 0, d.ZOTU99_not6)
min(sample_sums(d.ZOTU99_not6))

## The above value being OK for rarefaction, I rarefy without defining a cutoff:
d.ZOTU99_not6.rare <- rarefy_even_depth(d.ZOTU99_not6, rngseed = 20200228)
sample_sums(d.ZOTU99_not6.rare)[1:5]


### Lettuce physeq object ####
## Quick approach to determine optimal rarefaction level is to look at sample_sums
max(sample_sums(d.ZOTU99_L_no0))
min(sample_sums(d.ZOTU99_L_no0))
sort(sample_sums(d.ZOTU99_L_no0))   # to see which are the lowest depths and pick a cutoff


## Rarefy at a certain level.
# choose a depth from the above output which is reasonable (not loosing too many
# samples nor depth):
rare_cutoff <- 26497
d.ZOTU99_L_no0.rare_custom <- rarefy_even_depth(d.ZOTU99_L_no0, sample.size = rare_cutoff, rngseed = 20200228)
sample_sums(d.ZOTU99_L_no0.rare_custom)[1:5]


## Without defining a cutoff for rarefaction, all samples will be rarefied to
## the lowest depth sample. In the present data, two L6 samples have much lower
## depth than the rest > remove t6.
d.ZOTU99_L_no0_not6 <- subset_samples(d.ZOTU99_L_no0, Time != "6")
d.ZOTU99_L_no0_not6 <- prune_taxa(taxa_sums(d.ZOTU99_L_no0_not6) > 0, d.ZOTU99_L_no0_not6)
min(sample_sums(d.ZOTU99_L_no0_not6))

## The above value being OK for rarefaction, I rarefy without defining a cutoff:
d.ZOTU99_L_no0_not6.rare <- rarefy_even_depth(d.ZOTU99_L_no0_not6, rngseed = 20200228)
sample_sums(d.ZOTU99_L_no0_not6.rare)[1:5]


### Soil physeq object ####
## Quick approach to determine optimal rarefaction level is to look at sample_sums
max(sample_sums(d.ZOTU99_S_no0))
min(sample_sums(d.ZOTU99_S_no0))
sort(sample_sums(d.ZOTU99_S_no0))   # to see which are the lowest depths and pick a cutoff


## Rarefy at a certain level.
# choose a depth from the above output which is reasonable (not loosing too many
# samples nor depth):
rare_cutoff <- 65963
d.ZOTU99_S_no0.rare_custom <- rarefy_even_depth(d.ZOTU99_S_no0, sample.size = rare_cutoff, rngseed = 20200228)
sample_sums(d.ZOTU99_S_no0.rare_custom)[1:5]


## Without defining a cutoff for rarefaction, all samples will be rarefied to
## the lowest depth sample. In the present data, two L6 samples have much lower
## depth than the rest > remove t6.
d.ZOTU99_S_no0_not6 <- subset_samples(d.ZOTU99_S_no0, Time != "6")
d.ZOTU99_S_no0_not6 <- prune_taxa(taxa_sums(d.ZOTU99_S_no0_not6) > 0, d.ZOTU99_S_no0_not6)
min(sample_sums(d.ZOTU99_S_no0_not6))

## The above value being OK for rarefaction, I rarefy without defining a cutoff:
d.ZOTU99_S_no0_not6.rare <- rarefy_even_depth(d.ZOTU99_S_no0_not6, rngseed = 20200228)
sample_sums(d.ZOTU99_S_no0_not6.rare)[1:5]


### Water physeq object ####
## Quick approach to determine optimal rarefaction level is to look at sample_sums
max(sample_sums(d.ZOTU99_W_no0))
min(sample_sums(d.ZOTU99_W_no0))
sort(sample_sums(d.ZOTU99_W_no0))   # to see which are the lowest depths and pick a cutoff


## Rarefy at a certain level.
# choose a depth from the above output which is reasonable (not loosing too many
# samples nor depth):
rare_cutoff <- 22470
d.ZOTU99_W_no0.rare_custom <- rarefy_even_depth(d.ZOTU99_W_no0, sample.size = rare_cutoff, rngseed = 20200228)
sample_sums(d.ZOTU99_W_no0.rare_custom)[1:5]


## Without defining a cutoff for rarefaction, all samples will be rarefied to
## the lowest depth sample. In the present data, two L6 samples have much lower
## depth than the rest > remove t6.
d.ZOTU99_W_no0_not6 <- subset_samples(d.ZOTU99_W_no0, Time != "6")
d.ZOTU99_W_no0_not6 <- prune_taxa(taxa_sums(d.ZOTU99_W_no0_not6) > 0, d.ZOTU99_W_no0_not6)
min(sample_sums(d.ZOTU99_W_no0_not6))

## The above value being OK for rarefaction, I rarefy without defining a cutoff:
d.ZOTU99_W_no0_not6.rare <- rarefy_even_depth(d.ZOTU99_W_no0_not6, rngseed = 20200228)
sample_sums(d.ZOTU99_W_no0_not6.rare)[1:5]


### Manure physeq object ####
## Quick approach to determine optimal rarefaction level is to look at sample_sums
max(sample_sums(d.ZOTU99_M_no0))
min(sample_sums(d.ZOTU99_M_no0))
sort(sample_sums(d.ZOTU99_M_no0))   # to see which are the lowest depths and pick a cutoff


## Rarefy at a certain level.
# choose a depth from the above output which is reasonable (not loosing too many
# samples nor depth):
rare_cutoff <- 89996
d.ZOTU99_M_no0.rare_custom <- rarefy_even_depth(d.ZOTU99_M_no0, sample.size = rare_cutoff, rngseed = 20200228)
sample_sums(d.ZOTU99_M_no0.rare_custom)[1:3]



#### Export rarefied physeq objects for downstream scripts ####

RData.filename <- "GDC_ZOTU99-Silva_2.RData"
output <- file.path(odir, RData.filename)

rm(d.ZOTU99,
   d.ZOTU99_L_no0,
   d.ZOTU99_M_no0,
   d.ZOTU99_S_no0,
   d.ZOTU99_W_no0,
   idir, odir, rare_cutoff, RData.filename)

save.image(file = output)



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  COMPOSITION PLOTS (script NR.3)          ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

# Define paths to use
idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
RData.filename <- "GDC_ZOTU99-Silva_1.RData"
input <- file.path(odir, RData.filename)

load(input)

rm(RData.filename, input)

## Create an additional physeq object omitting lettuce
d.ZOTU99_SWM <- subset_samples(d.ZOTU99, System != "L")
sum(taxa_sums(d.ZOTU99_SWM) == 0)
d.ZOTU99_SWM_no0 <- prune_taxa(taxa_sums(d.ZOTU99_SWM) > 0, d.ZOTU99_SWM)

## Create an additional physeq object omitting water to keep only S and M
d.ZOTU99_SM <- subset_samples(d.ZOTU99_SWM_no0, System != "W")
sum(taxa_sums(d.ZOTU99_SM) == 0)
d.ZOTU99_SM_no0 <- prune_taxa(taxa_sums(d.ZOTU99_SM) > 0, d.ZOTU99_SM)


#### Plot composition rather than raw counts (custom function) ####

## Note that x-axis labels can be removed by setting axis.text.x =
## element_blank() vs. axis.text.x = element_text(angle = 90, vjust = 0.5, hjust
## = 1)

## Note aesthaetics when displaying several systems: replace facet_grid by
## facet_wrap command (as used with water or manure) to get same width for all
## systems

## Note that some samples visibly don't reach 100% if Archeae are not displayed.

## So if you don't wish to focus on certain subgroup such as Kingdom = Bacteria
## use the following:
pL <- phyloseq.extended::plot_composition(d.ZOTU99_L_no0, NULL, NULL, "Phylum", numberOfTaxa = 10, fill = "Phylum") 
pL <- pL + facet_grid(~Time, scales = "free_x", space = "free_x") +
  ggtitle("") +
  guides(fill = guide_legend(nrow = 11)) +
  labs(x = "Sample ID", y = "Relative Abundance") +
  theme(plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 0.5),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5,"cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold")) +
  scale_y_continuous(expand = c(0, 0))
plot(pL)

pSM <- phyloseq.extended::plot_composition(d.ZOTU99_SM_no0, NULL, NULL, "Phylum", numberOfTaxa = 10, fill = "Phylum") 
pSM <- pSM + facet_grid(~System, scales = "free_x", space = "free_x") +
  ggtitle("") +
  guides(fill = guide_legend(nrow = 11)) +
  labs(x = "Sample ID", y = "Relative Abundance") +
  theme(plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 0.5),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5,"cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold")) +
  scale_y_continuous(expand = c(0, 0))
plot(pSM)

pW <- phyloseq.extended::plot_composition(d.ZOTU99_W_no0, NULL, NULL, "Phylum", numberOfTaxa = 10, fill = "Phylum") 
pW <- pW +
  facet_wrap(~Time, scales = "free_x", nrow = 1) +
  ggtitle("") +
  guides(fill = guide_legend(nrow = 11)) +
  labs(x = "Sample ID", y = "Relative Abundance") +
  theme(plot.title = element_text(color = "black", face = "bold", size = 18, hjust = 0.5),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5,"cm"),
        legend.margin = margin(0, 21, 0, 5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        strip.text.x = element_text(size = 10, face = "bold")) +
  scale_y_continuous(expand = c(0, 0))
plot(pW)


ggpubr::ggarrange(pL, pSM, pW, labels = "AUTO", common.legend = FALSE,
                  legend = "right", ncol = 1, nrow = 3)


# Up to this point, each plot contains the top10 taxa, ordered by abundance. If
# we want to combine several plots in one figure, however, taxa occurring in
# several plots might have different colours in each plot. The following remedy
# to the problem was proposed by ***Mahendra Mariadassou
# <mahendra.mariadassou@inrae.fr>.***

## combine all phyla and put "Other" at the end (several steps needed when > 2 plots)
L_SM_phyla <- union(pL$data$Phylum, pSM$data$Phylum) |>  
  setdiff("Other") |>  
  sort() |>  
  union("Other")
L_SM_phyla

all_phyla <- union(L_SM_phyla, pW$data$Phylum) |>  
  setdiff("Other") |>  
  sort() |>  
  union("Other")
all_phyla


#And the corresponding palette, created using the brewer pal (to mimic what's done
#in the package, but you can tweak the colors if you want)
my_palette <- c(
  RColorBrewer::brewer.pal(length(all_phyla) - 3, "Paired"), ## colors from pre-existing palette
  "#FFEC8B", "#B15928", "black"                              ## custom colors for last phyla & "Other"
)
names(my_palette) <- all_phyla
my_palette


#We now enforce the correct phylum level and order in all plots
pL$data$Phylum <- factor(pL$data$Phylum, levels = all_phyla)
pSM$data$Phylum <- factor(pSM$data$Phylum, levels = all_phyla)
pW$data$Phylum <- factor(pW$data$Phylum, levels = all_phyla)


#We now change the fill scale of each plot:
pL_new <- pL + scale_fill_manual(values = my_palette, drop = FALSE) + guides(fill = guide_legend(ncol = 1))
pSM_new <- pSM + scale_fill_manual(values = my_palette, drop = FALSE) + guides(fill = guide_legend(ncol = 1))
pW_new <- pW + scale_fill_manual(values = my_palette, drop = FALSE) + guides(fill = guide_legend(ncol = 1))


#Before combining them in the final plot:
ggpubr::ggarrange(pL_new, pSM_new, pW_new, labels = "AUTO", common.legend = TRUE,
          legend = "right", ncol = 1, nrow = 3) %>%
  ggpubr::ggexport(filename = "relative_abund_all_MS.png", width = 4500, height = 6000, res = 600)



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  ALPHA DIVERSITIES (script NR.4)           ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

# Define paths to use
idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
RData.filename <- "GDC_ZOTU99-Silva_1.RData"
input <- file.path(odir, RData.filename)

load(input)


RData.filename <- "GDC_ZOTU99-Silva_2.RData"
input <- file.path(odir, RData.filename)

load(input)


rm(RData.filename, input, output)


#### Plot richness of the different sample types using 4 diversities ####
## Richness represented by boxplots. NOTE: by deleting geom_boxplot(...) the
## boxes disappear and only points remain, which is more meaningful in my case
## of 3 points per box!
# Lettuce ####
rL <- plot_richness(d.ZOTU99_L_no0, 
                    color = "Treatment", 
                    measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                    x = "Time") +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) + 
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, lettuce") +
  xlab("Time (weeks)")
plot(rL)

rL.rare <- plot_richness(d.ZOTU99_L_no0.rare_custom, 
                         color = "Treatment", 
                         measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                         x = "Time") +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) + 
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, lettuce rarefied") +
  xlab("Time (weeks)")
plot(rL.rare)

rL_not6 <- plot_richness(d.ZOTU99_L_no0_not6, 
                         color = "Treatment", 
                         measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                         x = "Time") +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) + 
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, lettuce") +
  xlab("Time (weeks)")
plot(rL_not6)

rL_not6.rare <- plot_richness(d.ZOTU99_L_no0_not6.rare, 
                              color = "Treatment", 
                              measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                              x = "Time") +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) + 
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, lettuce rarefied") +
  xlab("Time (weeks)")
plot(rL_not6.rare)


# Soil ####
rS <- plot_richness(d.ZOTU99_S_no0, 
                    color = "Treatment", 
                    measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                    x = "Time") +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) + 
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, soil") +
  xlab("Time (weeks)")
plot(rS)

rS.rare <- plot_richness(d.ZOTU99_S_no0.rare_custom, 
                         color = "Treatment", 
                         measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                         x = "Time") +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) + 
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, soil rarefied") +
  xlab("Time (weeks)")
plot(rS.rare)

rS_not6 <- plot_richness(d.ZOTU99_S_no0_not6, 
                         color = "Treatment", 
                         measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                         x = "Time") +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) + 
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, soil") +
  xlab("Time (weeks)")
plot(rS_not6)

rS_not6.rare <- plot_richness(d.ZOTU99_S_no0_not6.rare, 
                              color = "Treatment", 
                              measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                              x = "Time") +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) + 
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, soil rarefied") +
  xlab("Time (weeks)")
plot(rS_not6.rare)


# Water: with t0 ####
rW <- plot_richness(d.ZOTU99_W_no0, 
                    color = "Treatment", 
                    measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                    x = "Time") +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) + 
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, water") +
  xlab("Time (weeks)")
plot(rW)

rW.rare <- plot_richness(d.ZOTU99_W_no0.rare_custom, 
                         color = "Treatment", 
                         measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                         x = "Time") +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) + 
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, water rarefied") +
  xlab("Time (weeks)")
plot(rW.rare)

rW_not6 <- plot_richness(d.ZOTU99_W_no0_not6, 
                         color = "Treatment", 
                         measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                         x = "Time") +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) + 
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, water") +
  xlab("Time (weeks)")
plot(rW_not6)

rW_not6.rare <- plot_richness(d.ZOTU99_W_no0_not6.rare, 
                              color = "Treatment", 
                              measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
                              x = "Time") +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) + 
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, water rarefied") +
  xlab("Time (weeks)") +
  xlab("Time (weeks)")
plot(rW_not6.rare)


# Water: without t0 (plants weren't irrigated with it) ####

d.ZOTU99_W_not6_not0 <-
  subset_samples(d.ZOTU99_W_no0_not6, Time != "0")
d.ZOTU99_W_not6_not0 <-
  prune_taxa(taxa_sums(d.ZOTU99_W_not6_not0) > 0, d.ZOTU99_W_not6_not0)
d.ZOTU99_W_not6_not0.rare <-
  rarefy_even_depth(d.ZOTU99_W_not6_not0, rngseed = 20200228)

rW_not6_not0 <- plot_richness(
  d.ZOTU99_W_not6_not0,
  color = "Treatment",
  measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
  x = "Time"
) +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) +
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, water") +
  xlab("Time (weeks)")
plot(rW_not6_not0)


rW_not6_not0.rare <- plot_richness(
  d.ZOTU99_W_not6_not0.rare,
  color = "Treatment",
  measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
  x = "Time"
) +
  theme_bw() + geom_boxplot(aes(fill = Treatment), position =  position_dodge2(preserve = "single")) +
  geom_point() + theme(axis.text.x = element_text(vjust = 0.5)) +
  ggtitle("GDC: ZOTU99_Silva, water rarefied") +
  xlab("Time (weeks)")
plot(rW_not6_not0.rare)



#### ANOVA: estimate of Time and Treatment (and depth) on richness ####
## For ANOVA in my case it makes sense to omit t0 and look only at t1 - t5. I
## observed that including t0 makes Time and Treatment look collinear as
## indicated by a large kappa(aov.model), (probably because 0 and Z are the same
## thing and very different from the other times and treatments...) Was done for
## water above, therefore here for lettuce and soil. Remember always rarefy
## AFTER filtering

d.ZOTU99_L_not6_not0 <- subset_samples(d.ZOTU99_L_no0_not6, Time != "0")
d.ZOTU99_L_not6_not0 <- prune_taxa(taxa_sums(d.ZOTU99_L_not6_not0) > 0, d.ZOTU99_L_not6_not0)
d.ZOTU99_L_not6_not0.rare <- rarefy_even_depth(d.ZOTU99_L_not6_not0, rngseed = 20200228)

d.ZOTU99_S_not6_not0 <- subset_samples(d.ZOTU99_S_no0_not6, Time != "0")
d.ZOTU99_S_not6_not0 <- prune_taxa(taxa_sums(d.ZOTU99_S_not6_not0) > 0, d.ZOTU99_S_not6_not0)
d.ZOTU99_S_not6_not0.rare <- rarefy_even_depth(d.ZOTU99_S_not6_not0, rngseed = 20200228)



## If required update sample_data (metadata) ####
## Do this for all L, S, W
## If you forget to make time a factor with levels, ANOVA results will change
## drastically...
filename <- "NEWEST_SWL_sample_data.txt"
input <- file.path(idir, filename)

new.metadata <- read.table(input, header = T)

rownames(new.metadata) <- new.metadata$ID
physeq <- d.ZOTU99_L_not6_not0.rare
sample_data(physeq) <- new.metadata

class(physeq@sam_data[["Time"]])
physeq@sam_data[["Time"]] <- factor(physeq@sam_data[["Time"]])
class(physeq@sam_data[["Time"]])

sample_variables(physeq)

d.ZOTU99_L_not6_not0.rare <- physeq


## ------- don't do the uneven depth part, just to see depth's significance ------- ##
## Prepare data (uneven depth, all time points) ####
alpha.diversity <- estimate_richness(d.ZOTU99_L_no0, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
data.ZOTU99_L <- cbind(sample_data(d.ZOTU99_L_no0), alpha.diversity)
data.ZOTU99_L$Depth <- sample_sums(d.ZOTU99_L_no0)

## Perform ANOVA on observed richness, on data with _uneven_ depth
d.ZOTU99_L_no0.richness.anova <- aov(Observed ~ Depth + Time + Treatment, data.ZOTU99_L)
summary(d.ZOTU99_L_no0.richness.anova)


## Since depth is highly significant, and if you don't want to correct for Depth,
## compute alpha diversities on rarefied samples! 

## Prepare data (even depth, selected time points) ####
alpha.diversity <- estimate_richness(d.ZOTU99_L_not6_not0.rare, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
data.ZOTU99_L_not6_not0.rare <- cbind(sample_data(d.ZOTU99_L_not6_not0.rare), alpha.diversity)
data.ZOTU99_L_not6_not0.rare$Depth <- sample_sums(d.ZOTU99_L_not6_not0.rare)

alpha.diversity <- estimate_richness(d.ZOTU99_S_not6_not0.rare, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
data.ZOTU99_S_not6_not0.rare <- cbind(sample_data(d.ZOTU99_S_not6_not0.rare), alpha.diversity)
data.ZOTU99_S_not6_not0.rare$Depth <- sample_sums(d.ZOTU99_S_not6_not0.rare)

alpha.diversity <- estimate_richness(d.ZOTU99_W_not6_not0.rare, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
data.ZOTU99_W_not6_not0.rare <- cbind(sample_data(d.ZOTU99_W_not6_not0.rare), alpha.diversity)
data.ZOTU99_W_not6_not0.rare$Depth <- sample_sums(d.ZOTU99_W_not6_not0.rare)



## ANOVA on observed richness, data with _even_ depth ####
## TASK-1: change between Observed / Chao1 / Shannon / InvSimpson
## TASK-2: having time as integer or factor changes outcome excessively...! Do
##         both! Time in weeks vs. days however makes no difference!
## NOTE-1: changing the sequence of the variables completely changes
##         significance (e.g. Humid_7d in S from *** when placed last it becomes ns when
##         placed first...)
## NOTE-2: as opposed to kruskal.test, here you can combine several factors
d.ZOTU99_L_not6_not0.rare.richness.anova <- aov(InvSimpson ~ Time + Treatment, data.ZOTU99_L_not6_not0.rare)
summary(d.ZOTU99_L_not6_not0.rare.richness.anova)

d.ZOTU99_S_not6_not0.rare.richness.anova <- aov(InvSimpson ~ Time + Treatment, data.ZOTU99_S_not6_not0.rare)
summary(d.ZOTU99_S_not6_not0.rare.richness.anova)

d.ZOTU99_W_not6_not0.rare.richness.anova <- aov(InvSimpson ~ Treatment + Time, data.ZOTU99_W_not6_not0.rare)
summary(d.ZOTU99_W_not6_not0.rare.richness.anova)


#### Posthoc testing (as suggested by Jean-claude): ####
## In the following comparisons, switch between Observed / Chao1 / Shannon /
## InvSimpson. Combine sample data and diversity measure.
df <- data.ZOTU99_L_not6_not0.rare

## Kruskal-Wallis test
kruskal.test(InvSimpson ~ Time, df) # 


## Nemenyi Test for Multiple Comparisons (post hoc analysis)
# “Tukey” or “Chisq”
# p.adjust.method can be chosen (default: none)
PMCMRplus::kwAllPairsNemenyiTest(x = df$InvSimpson, g = df$Time, dist = "Tukey")

## Dunn Test for Multiple Comparisons  (post hoc analysis)
# bh : Benjamini and Hochberg
# p.adjust.method bonferroni (default: holm)
PMCMRplus::kwAllPairsDunnTest(InvSimpson ~ Time, data = df, method = "bh")


save(d.ZOTU99_L_not6_not0, d.ZOTU99_L_not6_not0.rare,
     d.ZOTU99_S_not6_not0, d.ZOTU99_S_not6_not0.rare,
     d.ZOTU99_W_not6_not0, d.ZOTU99_W_not6_not0.rare,
     rL_not6, rS_not6, rW_not6_not0,
     file = paste(odir, "GDC_ZOTU99-Silva_4.RData", sep = "/"))




####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  BETA DIVERSITIES (script NR.5)            ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

# Define paths to use
idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
RData.filename <- "GDC_ZOTU99-Silva_2.RData"
input <- file.path(odir, RData.filename)

load(input)

rm(RData.filename, input, output)


#### Compute beta diversity distances ####
## distances are computed with distance()

## Available distances are
## Bray-Curtis  : "bray"
## Jaccard      : "cc"
## Unifrac (UF) : "unifrac"
## Weighted UF  : "wunifrac"

dist.bc  <- phyloseq::distance(d.ZOTU99_L_no0.rare_custom, method = "bray")
dist.jac <- phyloseq::distance(d.ZOTU99_L_no0.rare_custom, method = "cc")
dist.uf  <- phyloseq::distance(d.ZOTU99_L_no0.rare_custom, method = "unifrac")
dist.wuf <- phyloseq::distance(d.ZOTU99_L_no0.rare_custom, method = "wunifrac")


## Different orders could reveal different structures
## Task: Test different variables of interest such as Time or Treatment
SampleOrder <- levels(reorder(sample_names(d.ZOTU99_L_no0.rare_custom),
                              as.numeric(get_variable(d.ZOTU99_L_no0.rare_custom,
                              "Time"))))


#### Visualize distance matrices with heatmaps ####
plot_dist_as_heatmap(dist.bc, order = SampleOrder, show.names = TRUE) + ggtitle("Bray Curtis; ordered by time")
plot_dist_as_heatmap(dist.jac, order = SampleOrder, show.names = TRUE) + ggtitle("Jaccard; ordered by time")
plot_dist_as_heatmap(dist.uf, order = SampleOrder, show.names = TRUE) + ggtitle("Unifrac; ordered by time")
plot_dist_as_heatmap(dist.wuf, order = SampleOrder, show.names = TRUE) + ggtitle("Weighted Unifrac; ordered by time")


#### Compute ordination objects and plot ####
#### NOTE: You can switch between ordination methods: PCoA (or MDS) / NMDS

### Lettuce only ####
p.bc <- plot_ordination(d.ZOTU99_L_no0.rare_custom,
                        ordinate(d.ZOTU99_L_no0.rare_custom, "PCoA", "bray"),
                        color = "Time") +
  theme_classic() + theme(axis.text=element_text(size=13), axis.title=element_text(size=15, face = "bold"),
                          legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=13),
                          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
                          legend.margin=margin(0,30,0,0)) +
  stat_ellipse(aes(group = Time), level = 0.95) + geom_point(size = 3)
plot(p.bc)


p.jac <- plot_ordination(d.ZOTU99_L_no0.rare_custom,
                         ordinate(d.ZOTU99_L_no0.rare_custom, "PCoA", "cc"),
                         color = "Time") +
  theme_classic() + theme(axis.text=element_text(size=13), axis.title=element_text(size=15, face = "bold"),
                          legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=13),
                          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
                          legend.margin=margin(0,30,0,0)) +
  stat_ellipse(aes(group = Time), level = 0.95) + geom_point(size = 3)
plot(p.jac)


p.uf <- plot_ordination(d.ZOTU99_L_no0.rare_custom,
                        ordinate(d.ZOTU99_L_no0.rare_custom, "PCoA", "unifrac"),
                        color = "Time") +
  theme_classic() + theme(axis.text=element_text(size=13), axis.title=element_text(size=15, face = "bold"),
                          legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=13),
                          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                          legend.margin=margin(0,30,0,0)) +
  stat_ellipse(aes(group = Time), level = 0.95) + geom_point(size = 3)
plot(p.uf)


p.wuf <- plot_ordination(d.ZOTU99_L_no0.rare_custom,
                         ordinate(d.ZOTU99_L_no0.rare_custom, "PCoA", "wunifrac"),
                         color = "Time") +
  theme_classic() + theme(axis.text=element_text(size=13), axis.title=element_text(size=15, face = "bold"),
                          legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=13),
                          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                          legend.margin=margin(0,30,0,0)) +
  stat_ellipse(aes(group = Time), level = 0.95) + geom_point(size = 3)
plot(p.wuf)


## Combine all four plots into one for publication
ggpubr::ggarrange(p.jac, p.bc, p.uf, p.wuf,
                  labels = c("A", "B", "C", "D"),
                  font.label = list(size = 20),
                  ncol = 2, nrow = 2)



### All systems (manure, soil, lettuce, water) ####
custom.col <- c("#FF69B4", "#C4961A", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

p.bc <- plot_ordination(d.ZOTU99.rare_custom,
                        ordinate(d.ZOTU99.rare_custom, "PCoA", "bray"),
                        color = "Time", shape = "System") +
  theme_classic() + theme(axis.text=element_text(size=13), axis.title=element_text(size=15, face = "bold"),
                          legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=13),
                          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
                          legend.margin=margin(0,30,0,0)) + geom_point(size = 3, stroke = 1.5) + scale_colour_manual(values = custom.col)

plot(p.bc)


p.jac <- plot_ordination(d.ZOTU99.rare_custom,
                        ordinate(d.ZOTU99.rare_custom, "PCoA", "cc"),
                        color = "Time", shape = "System") +
  theme_classic() + theme(axis.text=element_text(size=13), axis.title=element_text(size=15, face = "bold"),
                          legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=13),
                          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
                          legend.margin=margin(0,30,0,0)) + geom_point(size = 3, stroke = 1.5) + scale_colour_manual(values = custom.col)
plot(p.jac)


p.uf <- plot_ordination(d.ZOTU99.rare_custom,
                         ordinate(d.ZOTU99.rare_custom, "PCoA", "unifrac"),
                         color = "Time", shape = "System") +
  theme_classic() + theme(axis.text=element_text(size=13), axis.title=element_text(size=15, face = "bold"),
                          legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=13),
                          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
                          legend.margin=margin(0,30,0,0)) + geom_point(size = 3, stroke = 1.5) + scale_colour_manual(values = custom.col)
plot(p.uf)


p.wuf <- plot_ordination(d.ZOTU99.rare_custom,
                        ordinate(d.ZOTU99.rare_custom, "PCoA", "wunifrac"),
                        color = "Time", shape = "System") +
  theme_classic() + theme(axis.text=element_text(size=13), axis.title=element_text(size=15, face = "bold"),
                          legend.title=element_text(size=14, face = "bold"), legend.text=element_text(size=13),
                          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
                          legend.margin=margin(0,30,0,0)) + geom_point(size = 3, stroke = 1.5) + scale_colour_manual(values = custom.col)
plot(p.wuf)


## Combine all four plots into one for publication
ggpubr::ggarrange(p.jac, p.bc, p.uf, p.wuf,
                  labels = c("A", "B", "C", "D"),
                  font.label = list(size = 20),
                  ncol = 2, nrow = 2)



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  CLUSTERING (script NR.6)                  ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

# Define paths to use
idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
RData.filename <- "GDC_ZOTU99-Silva_2.RData"
input <- file.path(odir, RData.filename)

load(input)

rm(RData.filename, input, output)


## If required update sample_data ####
## Task: Do this for all L, S, W
filename <- "NEWEST_SWL_sample_data.txt"
input <- file.path(idir, filename)

new.metadata <- read.table(input, header = T)

rownames(new.metadata) <- new.metadata$ID
physeq <- d.ZOTU99_L_no0_not6.rare
sample_data(physeq) <- new.metadata

class(physeq@sam_data[["Time"]])
physeq@sam_data[["Time"]] <- factor(physeq@sam_data[["Time"]])
class(physeq@sam_data[["Time"]])

sample_variables(physeq)

d.ZOTU99_L_no0_not6.rare <- physeq


## If required add new columns to sample_data ####

## Manually add new columns to your sample_data
sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_Tunn_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_Tunn_3h > 15,
         "high",
         "low")

sample_data(d.ZOTU99_L_no0_not6.rare)$Humid_3h_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Humid_3h > 75,
         "humid",
         "normal")

sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_2m_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_2m_7d > 19,
         "warm",
         "normal")
sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_2m_binary[sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_2m_7d < 16] <-
  "cool"

sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_5cm_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_5cm_7d > 19,
         "warm",
         "normal")

sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_soil_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_soil_7d > 19,
         "warm",
         "normal")

sample_data(d.ZOTU99_L_no0_not6.rare)$Humid_7d_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Humid_7d > 75,
         "humid",
         "normal")

sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_intens_binary <-
  ifelse(
    sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_intens_7d > 10,
    "strong",
    "intermediate")
sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_intens_binary[sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_intens_7d < 2.5] <-
  "weak"

sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_7d > 5,
         "strong",
         "intermediate")
sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_binary[sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_7d < 1.5] <-
  "weak"

sample_data(d.ZOTU99_L_no0_not6.rare)$Radiation_binary <-
  ifelse(
    sample_data(d.ZOTU99_L_no0_not6.rare)$Radiation_7d > 6000,
    "strong",
    "intermediate")
sample_data(d.ZOTU99_L_no0_not6.rare)$Radiation_binary[sample_data(d.ZOTU99_L_no0_not6.rare)$Radiation_7d < 5000] <-
  "weak"



####  CLUSTERING ####
## Clustering is done with following general syntax
## Task: Switch between ward.D2 / single / complete

par(mfrow = c(2, 2))


clust.jac  <- plot_clust(d.ZOTU99_L_no0_not6.rare,
                         dist = "cc",
                         method = "ward.D2",
                         color = "Time",
                         title = "clustering tree: ward.D2 / Jaccard") 

clust.bc  <- plot_clust(d.ZOTU99_L_no0_not6.rare,
                        dist = "bray",
                        method = "ward.D2",
                        color = "Time",
                        title = "clustering tree: ward.D2 / Bray-Curtis")

clust.uf  <- plot_clust(d.ZOTU99_L_no0_not6.rare,
                        dist = "unifrac",
                        method = "ward.D2",
                        color = "Time",
                        title = "clustering tree: ward.D2 / Unifrac") 

clust.wuf <- plot_clust(d.ZOTU99_L_no0_not6.rare,
                        dist = "wunifrac",
                        method = "ward.D2",
                        color = "Time",
                        title = "clustering tree: ward.D2 / wUnifrac") 

par(mfrow = c(1, 1))



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  HEATMAPS (script NR.7)                    ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

# Define paths to use
idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
RData.filename <- "GDC_ZOTU99-Silva_2.RData"
input <- file.path(odir, RData.filename)

load(input)

rm(RData.filename, input, output)


## If required update sample_data ####
## Task: Do this for all L, S, W
filename <- "NEWEST_SWL_sample_data.txt"
input <- file.path(idir, filename)

new.metadata <- read.table(input, header = T)

rownames(new.metadata) <- new.metadata$ID
physeq <- d.ZOTU99_L_no0_not6.rare
sample_data(physeq) <- new.metadata

class(physeq@sam_data[["Time"]])
physeq@sam_data[["Time"]] <- factor(physeq@sam_data[["Time"]])
class(physeq@sam_data[["Time"]])

sample_variables(physeq)

d.ZOTU99_L_no0_not6.rare <- physeq


## If required add new columns to sample_data ####

## Manually add a a new column to your sample_data
sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_Tunn_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_Tunn_3h > 15,
         "high",
         "low")

sample_data(d.ZOTU99_L_no0_not6.rare)$Humid_3h_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Humid_3h > 75,
         "humid",
         "normal")

sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_2m_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_2m_7d > 19,
         "warm",
         "normal")
sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_2m_binary[sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_2m_7d < 16] <-
  "cool"

sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_5cm_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_5cm_7d > 19,
         "warm",
         "normal")

sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_soil_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Temp_soil_7d > 19,
         "warm",
         "normal")

sample_data(d.ZOTU99_L_no0_not6.rare)$Humid_7d_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Humid_7d > 75,
         "humid",
         "normal")

sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_intens_binary <-
  ifelse(
    sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_intens_7d > 10,
    "strong",
    "intermediate")
sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_intens_binary[sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_intens_7d < 2.5] <-
  "weak"

sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_binary <-
  ifelse(sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_7d > 5,
         "strong",
         "intermediate")
sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_binary[sample_data(d.ZOTU99_L_no0_not6.rare)$Precip_7d < 1.5] <-
  "weak"

sample_data(d.ZOTU99_L_no0_not6.rare)$Radiation_binary <-
  ifelse(
    sample_data(d.ZOTU99_L_no0_not6.rare)$Radiation_7d > 6000,
    "strong",
    "intermediate")
sample_data(d.ZOTU99_L_no0_not6.rare)$Radiation_binary[sample_data(d.ZOTU99_L_no0_not6.rare)$Radiation_7d < 5000] <-
  "weak"

#### Heatmap generation ####
## Base heatmap with adapted color scale and clustered by Time; but too many taxa
p <-
  plot_heatmap(
    d.ZOTU99_L_no0_not6.rare,
    low = "yellow",
    high = "red",
    na.value = "white"
  ) +
  facet_grid( ~ Time, scales = "free_x")
plot(p)


## Heatmap for top X taxa only
## Task: switch between physeq and physeq.rare, it makes a difference! ***Mahendra is using non-rarefied for top_200***
top_200      <-
  names(sort(taxa_sums(d.ZOTU99_L_no0_not6.rare), decreasing = TRUE)[1:200])
small.physeq <- prune_taxa(top_200, d.ZOTU99_L_no0_not6.rare)
p            <-
  plot_heatmap(
    small.physeq,
    low = "yellow",
    high = "red",
    na.value = "white"
  ) +
  facet_grid( ~ Time, scales = "free_x")
plot(p)


## You can choose a custom taxa order in the heatmap
## Group OTUs in a clustering tree using prevalence and specificity patterns
taxa_description <-
  estimate_prevalence(d.ZOTU99_L_no0_not6.rare,
                      group = "Time",
                      format = "wide")
taxa_distance    <- dist(taxa_description)
taxa_clust       <- hclust(taxa_distance, method = "ward.D2")

otu.order        <- taxa_clust$labels[taxa_clust$order]

## Order taxa according to the ordering in the tree
p                <-
  plot_heatmap(
    d.ZOTU99_L_no0_not6.rare,
    taxa.order = otu.order,
    low = "yellow",
    high = "red",
    na.value = "white"
  ) +
  facet_grid( ~ Time, scales = "free_x")
plot(p)



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  MULTIVARIATE ANOVA (script NR.8)          ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

# Define paths to use
idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
RData.filename <- "GDC_ZOTU99-Silva_2.RData"
input <- file.path(odir, RData.filename)

load(input)

rm(RData.filename, input, output)


#### Multivariate ANOVA with adonis ####
metadata <- as(sample_data(d.ZOTU99_L_no0.rare_custom), "data.frame")


## First examine if within-groups dispersion is homogeneous (betadisper should be NS).
##      What is this doing?
##        - by using betadisper (and subsequently permutesting the model), we
##        test for homogeneous dispersions within groups, which is needed to 
##        consider ADONIS results as valid (homogeneous dispersion is pre-
##        requisite for adonis)! That is, if betadisper is *significant*
##        you should NOT do ADONIS!
##        - if betadisper is NOT significant, then you can do ADONIS to see
##        whether community structure differs *significantly*

# Task-1: adjust accordingly to use either Time or Treatment
groups  <- metadata[,"Treatment"]

ord.bc  <- phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "bray")
ord.cc  <- phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "cc")
ord.uf  <- phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "unifrac")
ord.wuf <- phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "wunifrac")

mod.bc  <- betadisper(ord.bc, groups)
mod.cc  <- betadisper(ord.cc, groups)
mod.uf  <- betadisper(ord.uf, groups)
mod.wuf <- betadisper(ord.wuf, groups)

permutest(mod.bc, perm = 9999, pairwise = T)
permutest(mod.cc, perm = 9999, pairwise = T)
permutest(mod.uf, perm = 9999, pairwise = T)
permutest(mod.wuf, perm = 9999, pairwise = T)

# If dispersion is different between groups (i.e. above test is sig.), the
# following plots can give you an impression:
# For exact comparison of group pairs, look at TukeyHSD which shows you which
# pairs have sig. different dispersion.
plot(mod.bc)
plot(mod.cc)
plot(mod.uf)
plot(mod.wuf)

boxplot(mod.bc)
boxplot(mod.cc)
boxplot(mod.uf)
boxplot(mod.wuf)

# TukeyHSD as post-hoc test, but if betadisper is significant (meaning the
# groups have HETEROGENEOUS dispersions), the ADONIS results are not
# trustworthy, since ADONIS' assumption of homogeneous dispersion is NOT met.
mod.bc.HSD  <- TukeyHSD(mod.bc)
mod.cc.HSD  <- TukeyHSD(mod.cc)
mod.uf.HSD  <- TukeyHSD(mod.uf)
mod.wuf.HSD <- TukeyHSD(mod.wuf)

mod.bc.HSD
plot(mod.bc.HSD)

mod.cc.HSD
plot(mod.cc.HSD)

mod.uf.HSD
plot(mod.uf.HSD)

mod.wuf.HSD
plot(mod.wuf.HSD)


## Second, you can now do adonis FOR DISTANCE MEASURES WITH non-significant
## BETADISPER, to see if groups differ:
## Task-2: Test different variables, one by one
adonis2(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "bray")     ~ Treatment,
  data = metadata,
  perm = 9999
)
adonis2(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "cc")       ~ Treatment,
  data = metadata,
  perm = 9999
)
adonis2(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "unifrac")  ~ Treatment,
  data = metadata,
  perm = 9999
)
adonis2(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "wunifrac") ~ Treatment,
  data = metadata,
  perm = 9999
)

## Task-3: Test different variables in combination if interested
adonis2(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "bray")     ~ Time + Treatment,
  data = metadata,
  perm = 9999
)
adonis2(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "cc")       ~ Time + Treatment,
  data = metadata,
  perm = 9999
)
adonis2(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "unifrac")  ~ Time + Treatment,
  data = metadata,
  perm = 9999
)
adonis2(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "wunifrac") ~ Time + Treatment,
  data = metadata,
  perm = 9999
)


#### Multivariate ANOVA with PAIRWISE adonis ####

## to get an insight into which pairs of groups significantly differ, one can
## use pairwise.adonis:
# import it:
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni',reduce=NULL,perm=999)
{
  
  co <- combn(unique(as.character(factors)),2)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  
  for(elem in 1:ncol(co)){
    if(inherits(x, 'dist')){
      x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                      factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    }
    
    else  (
      if (sim.function == 'daisy'){
        x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
      } 
      else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    )
    
    ad <- adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])],
                 permutations = perm);
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c(Df,ad$aov.tab[1,1])
    SumsOfSqs <- c(SumsOfSqs, ad$aov.tab[1,2])
    F.Model <- c(F.Model,ad$aov.tab[1,4]);
    R2 <- c(R2,ad$aov.tab[1,5]);
    p.value <- c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
  
  if(!is.null(reduce)){
    pairw.res <- subset (pairw.res, grepl(reduce,pairs))
    pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
    
    sig = c(rep('',length(pairw.res$p.adjusted)))
    sig[pairw.res$p.adjusted <= 0.1] <-'.'
    sig[pairw.res$p.adjusted <= 0.05] <-'*'
    sig[pairw.res$p.adjusted <= 0.01] <-'**'
    sig[pairw.res$p.adjusted <= 0.001] <-'***'
    pairw.res <- data.frame(pairw.res[,1:7],sig)
  }
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
} 


# define the vector 'factors' which will tell the function e.g. which treatment
# each data entry belongs to
factors_treat <- sample_data(d.ZOTU99_L_no0.rare_custom)$Treatment
factors_time  <- sample_data(d.ZOTU99_L_no0.rare_custom)$Time


# run it (beware of defaults like permutation level 999 and p-value correction
# 'bonferroni')
pairwise.adonis(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "bray"),
  factors_treat,
  perm = 9999,
  p.adjust.m = 'BH'
)
pairwise.adonis(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "cc"),
  factors_treat,
  perm = 9999,
  p.adjust.m = 'BH'
)
pairwise.adonis(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "unifrac"),
  factors_treat,
  perm = 9999,
  p.adjust.m = 'BH'
)
pairwise.adonis(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "wunifrac"),
  factors_treat,
  perm = 9999,
  p.adjust.m = 'BH'
)

pairwise.adonis(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "bray"),
  factors_time,
  perm = 9999,
  p.adjust.m = 'BH'
)
pairwise.adonis(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "cc"),
  factors_time,
  perm = 9999,
  p.adjust.m = 'BH'
)
pairwise.adonis(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "unifrac"),
  factors_time,
  perm = 9999,
  p.adjust.m = 'BH'
)
pairwise.adonis(
  phyloseq::distance(d.ZOTU99_L_no0.rare_custom, "wunifrac"),
  factors_time,
  perm = 9999,
  p.adjust.m = 'BH'
)



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  DIFFERENTIAL ABUNDANCE (script NR.9)      ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

# Define paths to use
idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
RData.filename <- "GDC_ZOTU99-Silva_1.RData"
input <- file.path(odir, RData.filename)

load(input)

RData.filename <- "GDC_ZOTU99-Silva_4.RData"
input <- file.path(odir, RData.filename)

load(input)

rm(RData.filename, input)


#### Prepare datasets ####
# To compare only time points w1 to w5, select these two time points and rarefy
# data to even depth
d.ZOTU99_L_t15 <- subset_samples(d.ZOTU99_L_no0, Time %in% c(1, 5))
d.ZOTU99_L_t15 <- prune_taxa(taxa_sums(d.ZOTU99_L_t15) > 0, d.ZOTU99_L_t15)
d.ZOTU99_L_t15.rare <- rarefy_even_depth(d.ZOTU99_L_t15, rngseed = 20200311)

## Define physeq object for investigating effect of manure and river water
## application
physeq <- d.ZOTU99_L_not6_not0.rare
physeq@sam_data[["Time"]] <- factor(physeq@sam_data[["Time"]])
class(physeq@sam_data[["Time"]])
sample_variables(physeq)

## Add new columns to sample_data, to enable DESeq2 pairwise comparisons
sample_data(physeq)$Manure[sample_data(physeq)$Treatment %in% c("B", "C")] <-
  TRUE
sample_data(physeq)$Manure[sample_data(physeq)$Treatment %in% c("A", "D")] <-
  FALSE

sample_data(physeq)$River[sample_data(physeq)$Treatment %in% c("A", "B")] <-
  TRUE
sample_data(physeq)$River[sample_data(physeq)$Treatment %in% c("C", "D")] <-
  FALSE

sample_data(d.ZOTU99_L_t15.rare)$Plantlets[sample_data(d.ZOTU99_L_t15.rare)$Time == 1] <-
  TRUE
sample_data(d.ZOTU99_L_t15.rare)$Plantlets[sample_data(d.ZOTU99_L_t15.rare)$Time == 5] <-
  FALSE


#### DESeq2 ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# First, agglomerate taxa at the taxonomic level of your choice.                    #
#  Note that once agglomerated, data need to be re-imported before agglomerating    #
#  them at another level, especially when going from higher to lower ranks (e.g.    #
#  order > family).                                                                 #
#                                                                                   #
# Second, convert physeq- to deseq2-object and run DESeq2                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## Task: Choose a level at which you want to agglomerate [must be one of
## rank_names(physeq)].
## Note: in [,1:x] x must be adapted accordingly, also in later command cbind().

# For effect of river water or manure
physeq                          <- phyloseq::tax_glom(physeq, taxrank = "Family")
tax_table(physeq)               <- tax_table(physeq)[,1:5]

# For effect of time
d.ZOTU99_L_t15.rare             <-
  phyloseq::tax_glom(d.ZOTU99_L_t15.rare, taxrank = "Family")
tax_table(d.ZOTU99_L_t15.rare)  <-
  tax_table(d.ZOTU99_L_t15.rare)[, 1:5]

## Convert physeq to deseq2 object
cds.man <- phyloseq_to_deseq2(physeq, ~ Manure)
cds.riv <- phyloseq_to_deseq2(physeq, ~ River)
cds.t   <- phyloseq_to_deseq2(d.ZOTU99_L_t15.rare, ~ Plantlets)

dds.man <- DESeq2::DESeq(cds.man)
dds.riv <- DESeq2::DESeq(cds.riv)
dds.t   <- DESeq2::DESeq(cds.t)


## Since we only have a *binary factors*, we can use the following syntax to
## format the log2-fold change from the fitted model.
results.man <- DESeq2::results(dds.man, name = "ManureTRUE", tidy = TRUE)
results.riv <- DESeq2::results(dds.riv, name = "RiverTRUE", tidy = TRUE)
results.t   <- DESeq2::results(dds.t, name = "PlantletsTRUE", tidy = TRUE)

## Rename the column called "row" to "OTU"
da.otus.man <- results.man %>% rename(OTU = row)
da.otus.riv <- results.riv %>% rename(OTU = row)
da.otus.t   <- results.t   %>% rename(OTU = row)


## Select only OTUs with significant padj (default for correction is BH = FDR;
## other correction methods must be specified). Add taxonomic information in
## addition to ZOTU-numbers.

## Task: choose padj to your taste

taxtable   <- data.frame(tax_table(physeq))
taxtable15 <- data.frame(tax_table(d.ZOTU99_L_t15.rare))


da.otus.man <- subset(da.otus.man, padj < 0.05)
dim(da.otus.man)
da.otus.man <- cbind(da.otus.man, taxtable[rownames(taxtable) %in% da.otus.man$OTU,])
da.otus.man <- arrange(da.otus.man, log2FoldChange)
head(da.otus.man, 5)

da.otus.riv <- subset(da.otus.riv, padj < 0.01)
dim(da.otus.riv)
da.otus.riv <- cbind(da.otus.riv, taxtable[rownames(taxtable) %in% da.otus.riv$OTU,])
da.otus.riv <- arrange(da.otus.riv, log2FoldChange)
head(da.otus.riv, 5)

da.otus.t   <- subset(da.otus.t, padj < 0.01)
dim(da.otus.t)
da.otus.t   <- cbind(da.otus.t, taxtable15[rownames(taxtable15) %in% da.otus.t$OTU,])
da.otus.t   <- arrange(da.otus.t, log2FoldChange)
head(da.otus.t, 5)




#### Enriched taxa (DA: differentially abundant) ####
## We're going to add a "da.class" rank to the taxonomy table. That rank tells
## us whether:
## - a taxon is not DA: "none"
## - a taxon is DA and enriched e.g. after manure application: "With Manure"
## - a taxon is DA and enriched e.g. without manure application: "No Manure"
## We first create the "da.class" vector. By default, all taxa are set to "None"
## and we're going to change the status of ZOTUs appearing in the da.otus data.frame

da.class        <- rep("None", ntaxa(physeq))
names(da.class) <- taxa_names(physeq)

manure.otus              <-
  subset(da.otus.man, log2FoldChange > 0)$OTU ## set of DA taxa enriched with manure?
no.manure.otus           <-
  subset(da.otus.man, log2FoldChange < 0)$OTU ## set of DA taxa enriched without manure?
da.class[manure.otus]    <- "With Manure"
da.class[no.manure.otus] <- "No Manure"


da.class.t15        <- rep("None", ntaxa(d.ZOTU99_L_t15.rare))
names(da.class.t15) <- taxa_names(d.ZOTU99_L_t15.rare)

plantlets.otus               <-
  subset(da.otus.t, log2FoldChange > 0)$OTU ## set of DA taxa enriched on plantlets
salad.otus                   <-
  subset(da.otus.t, log2FoldChange < 0)$OTU ## set of DA taxa enriched on salad at harvest
da.class.t15[plantlets.otus] <- "Plantlets"
da.class.t15[salad.otus]     <- "Salad at harvest"


tax_table(physeq) <- 
  cbind(tax_table(physeq)[, 1:5], da.class)
head(tax_table(physeq)) ## Updated!

tax_table(d.ZOTU99_L_t15.rare) <-
  cbind(tax_table(d.ZOTU99_L_t15.rare)[, 1:5], da.class.t15)
head(tax_table(d.ZOTU99_L_t15.rare)) ## Updated!


## With this additional "rank", we can go back to composition plots and assess
## the composition within a class of DA taxa (here enriched after manure or river)

## Define count_to_prop function
count_to_prop <- function(x) { return( x / sum(x) )}


## Plot top N differentially abundant taxa (here top 10 are plotted)
## Task-1: Adjust the taxonomic rank to what you chose above for agglomeration
## Task-2: Adapt the number of taxa to be plotted in plot_composition()

# Effect of manure application:
manure_fraction <- physeq %>% 
  transform_sample_counts(fun = count_to_prop) %>% ## compositional data 
  subset_taxa(da.class == "No Manure")              ## of DA taxa only 
pM1 <-
  plot_composition(manure_fraction,
                   NULL,
                   NULL,
                   "Family",
                   fill = "Family",
                   5,
                   taxaOrder = "name")
pM1 <-
  pM1 + facet_wrap( ~ Manure, scales = "free_x", nrow = 1) + ylim(0, 1)
plot(pM1)

manure_fraction <- physeq %>% 
  transform_sample_counts(fun = count_to_prop) %>% ## compositional data 
  subset_taxa(da.class == "With Manure")              ## of DA taxa only 
pM2 <-
  plot_composition(manure_fraction,
                   NULL,
                   NULL,
                   "Family",
                   fill = "Family",
                   5,
                   taxaOrder = "name")
pM2 <-
  pM2 + facet_wrap( ~ Manure, scales = "free_x", nrow = 1) + ylim(0, 1)
plot(pM2)


# Effect of time (young [week1] vs. mature [week5] lettuce)
plantlet_fraction <- d.ZOTU99_L_t15.rare %>% 
  transform_sample_counts(fun = count_to_prop) %>% ## compositional data 
  subset_taxa(da.class.t15 == "Plantlets")  ## of DA taxa only 
pL1 <-
  plot_composition(plantlet_fraction, NULL, NULL, "Family", fill = "Family", 5)
pL1 <-
  pL1 + facet_wrap( ~ Time, scales = "free_x", nrow = 1) + ylim(0, 1) +
  scale_fill_brewer(palette = "Paired") +
  theme(
    axis.text.y = element_text(
      colour = "black",
      size = 12,
      face = "bold"
    ),
    axis.text.x = element_text(
      colour = "black",
      face = "bold",
      size = 12,
      hjust = 0.9,
      vjust = 0.5
    ),
    legend.text = element_text(
      size = 12,
      face = "bold.italic",
      colour = "black"
    ),
    legend.position = "right",
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      face = "bold",
      size = 16,
      colour = "black",
      margin = margin(r = 10)
    ),
    legend.title = element_text(
      size = 12,
      colour = "black",
      face = "bold"
    ),
    plot.title = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5),
    legend.key = element_blank()
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
plot(pL1)

plantlet_fraction <- d.ZOTU99_L_t15.rare %>% 
  transform_sample_counts(fun = count_to_prop) %>% ## compositional data 
  subset_taxa(da.class.t15 == "Salad at harvest")  ## of DA taxa only 
pL5 <-
  plot_composition(plantlet_fraction, NULL, NULL, "Family", fill = "Family", 5)
pL5 <-
  pL5 + facet_wrap( ~ Time, scales = "free_x", nrow = 1) + ylim(0, 1) +
  scale_fill_brewer(palette = "Paired") +
  theme(
    axis.text.y = element_text(
      colour = "black",
      size = 12,
      face = "bold"
    ),
    axis.text.x = element_text(
      colour = "black",
      face = "bold",
      size = 12,
      hjust = 0.9,
      vjust = 0.5
    ),
    legend.margin = margin(r = 25),
    legend.text = element_text(
      size = 12,
      face = "bold.italic",
      colour = "black"
    ),
    legend.position = "right",
    axis.title.x = element_text(
      face = "bold",
      size = 16,
      colour = "black"
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 16,
      colour = "black",
      margin = margin(r = 10)
    ),
    legend.title = element_text(
      size = 12,
      colour = "black",
      face = "bold"
    ),
    plot.title = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5),
    legend.key = element_blank()
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
plot(pL5)


ggpubr::ggarrange(pL1, pL5,
                  labels = "AUTO", common.legend = FALSE,
                  font.label = list(size = 20, color = "black", face = "bold", family = NULL),
                  legend = "right",
                  ncol = 1, nrow = 2) %>%
  ggpubr::ggexport(filename = "DA_microbiome_L1-vs-L5.png", width = 5000, height = 5000, res = 600)



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  Rel. Abundance numbers (script NR.10)     ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

# Define paths to use
idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
RData.filename <- "GDC_ZOTU99-Silva_1.RData"
input <- file.path(odir, RData.filename)

load(input)

rm(RData.filename, input)


#### Make numbers relative (chose L, M, S, _or_ W data) ####
# Task-1: change the phyloseq object between *_L_*, *_M_*, *_S_*, or *_W_*
ps <- tax_glom(d.ZOTU99_L_no0, "Phylum")
ps0 <- transform_sample_counts(ps, function(x)100 * x / sum(x))

# Task-2: Then use either option (a) or option (b):
# (a) for S,L,W: merge samples if wanted to get means, e.g. mean of all L1, L2,
# L3 ... samples
ps1 <- merge_samples(ps0, "Sys_Time")
ps2 <- transform_sample_counts(ps1, function(x)100 * x / sum(x))
otu_table(ps2) <- t(otu_table(ps2)) # transpose the matrix otus !!!
OTUg <- otu_table(ps2)
TAXg <- tax_table(ps2)[,"Phylum"]
GTable <- merge(TAXg, OTUg, by = 0, all = TRUE)
GTable$Row.names = NULL
GTable$Mean = rowMeans(GTable[,-c(1)], na.rm = TRUE)
GTable_L <- GTable[order(GTable$Mean, decreasing = TRUE),]
head(GTable_L, 10)

# (b) for MANURE use these commands, since there is no mean to be calculated...
ps1 <- merge_samples(ps0, "Sys_Time")
ps2 <- transform_sample_counts(ps1, function(x)100 * x / sum(x))
otu_table(ps2) <- t(otu_table(ps2)) # transpose the matrix otus !!!
OTUg <- otu_table(ps2)
TAXg <- tax_table(ps2)[,"Phylum"]
GTable <- merge(TAXg, OTUg, by = 0, all = TRUE)
GTable$Row.names = NULL
GTable_M <- GTable[order(GTable$M0, decreasing = TRUE),]
head(GTable_M, 10)

# export data frames into csv files
write.csv(GTable_L,"rel.abund_phyla_L.csv", row.names = FALSE)
write.csv(GTable_M,"rel.abund_phyla_M.csv", row.names = FALSE)
write.csv(GTable_S,"rel.abund_phyla_S.csv", row.names = FALSE)
write.csv(GTable_W,"rel.abund_phyla_W.csv", row.names = FALSE)



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  Venn Diagram                              ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
## Based on results from deepARG ####
##  identity    70%
##  e-value     1e-10
##  coverage    60%
##  probability 0.7
venn.plot_1 <- draw.quad.venn(
  area1 = 256,
  area2 = 95,
  area3 = 102,
  area4 = 141,
  n12 = 41,
  n13 = 94,
  n14 = 126,
  n23 = 39,
  n24 = 35,
  n34 = 69,
  n123 = 34,
  n124 = 33,
  n134 = 68,
  n234 = 31,
  n1234 = 30,
  category = c("lettuce", "manure", "soil", "water"),
  fill = c("lime green", "orange", "saddle brown", "blue"),
  lty = "dashed",
  cex = 7,
  cat.cex = 7,
  cat.col = c("lime green", "orange", "saddle brown", "blue"),
  margin = 0.02
)

# Writing to file
png(
  filename = "QuadVenn_all_deepARG_70-0.7.png",
  width = 1600,
  height = 1700,
  units = "px"
)
grid.draw(venn.plot_1);
dev.off()



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  PCoA clustering à la AGRsOAP v2.0         ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
Data.filename_16S <- "stagetwo_all.normalize_16s.subtype.tab.txt"
Data.filename_cell <- "stagetwo_all.normalize_cellnumber.subtype.tab.txt"
Data.filename_ppm <- "stagetwo_all.ppm.subtype.txt"

input_16S <- file.path(idir, Data.filename_16S)
input_cell <- file.path(idir, Data.filename_cell)
input_ppm <- file.path(idir, Data.filename_ppm)

rm(Data.filename_16S, Data.filename_cell, Data.filename_ppm)

#### Generate PCoA plots ####
## The following vector must be adapted accordingly!
cols <-
  c(rep("Lettuce", 9), "Manure", rep("Soil", 9), rep("Water", 2)) #add example data titles if included


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# NOTE-1: if using *.mergesubtype.* tables as mothertables, which contain            #
# __apart from my samples__ also some example data from public DBs, the              #
# information of the first row (aminoglycoside__aac(2')-I) will get lost (there      #
# seems to be some bug which I can't fix because the numbers for this ARG are        #
# missing for the example data).                                                     #
#                                                                                    #                         
# NOTE-2: PARAMETERS to tweak:                                                       #                           
#         (1) type of table to use as mothertable: c(16s, cell). There are also      #
#             the "gene" tables which, how- ever are very complex, containing ALL    #
#             information incl. type and subtype. Therefore don't use these tables   #
#             for PCoA plots.                                                        #
#         (2) method to use for PCoA plot: c(bray, jaccard). There are more but      #
#             I'm not sure if meaningful. ARGsOAP team was using bray.               #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


# 16S normalized ####
mothertable  <-
  read.table(
    file = input_16S,
    sep = "\t",
    header = T,
    row.names = 1,
    quote = "",
    stringsAsFactors = FALSE
  )
mothertable  <-
  mothertable[which(rowSums(mothertable) != 0), ] # remove zero lines
mothertable  <- t(mothertable)
samples.orig <-
  rownames(mothertable)                           # recording original sample names
mothertable  <-
  mothertable[which(rowSums(mothertable) != 0), ] # would remove samples with only zeroes
dim(mothertable)


# update cols in case certain samples have been filtered out.
cols <- cols[match(rownames(mothertable), samples.orig)]

png("stage2_16S_subtype_bray.png", res = 300, width = 2200, height = 2400)

vd       <- vegdist(mothertable, method = "bray")
vd.pco   <- pco(vd, k = 10)
pcoadata <- data.frame(vd.pco$points[, 1], vd.pco$points[, 2], cols)
Environment <- pcoadata$cols
pc1n     <- vd.pco$eig[1] / sum(vd.pco$eig)
pc2n     <- vd.pco$eig[2] / sum(vd.pco$eig)
xl       <- paste("PCoA1 (", (pc1n * 10000) %/% 100, "%)", sep = "")
yl       <- paste("PCoA2 (", (pc2n * 10000) %/% 100, "%)", sep = "")

p <-
  ggplot(pcoadata,
         aes(x = vd.pco.points...1., y = vd.pco.points...2., color = Environment)) +
  geom_point(size = 3, alpha = .8) +
  scale_color_manual(values = c("#00B81F", "#F8766D", "#C59900", "#592EFF")) +
  geom_text(
    aes(label = samples.orig),
    hjust = -0.6,
    vjust = 0,
    key_glyph = "point"
  )

p + theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.background = element_rect(fill = alpha('white', 0.2))
) + labs(fill = "Samples subtype") + xlab(xl) + ylab(yl) + xlim(-0.6, 0.6)
dev.off()



# Cell number normalized ####
mothertable <-
  read.table(
    file = input_cell,
    sep = "\t",
    header = T,
    row.names = 1,
    quote = "",
    stringsAsFactors = FALSE
  )
mothertable <-
  mothertable[which(rowSums(mothertable) != 0), ] # remove zero lines
mothertable <- t(mothertable)
samples.orig <- rownames(mothertable)
mothertable <- mothertable[which(rowSums(mothertable) != 0), ]
dim(mothertable)

# update cols in case certain samples have been filtered out.
cols <- cols[match(rownames(mothertable), samples.orig)]

png(
  "stage2_cell_subtype_bray.png",
  res = 300,
  width = 2200,
  height = 2400
)
vd <- vegdist(mothertable, method = "bray")
vd.pco <- pco(vd, k = 10)
pcoadata <- data.frame(vd.pco$points[, 1], vd.pco$points[, 2], cols)
Environment <- pcoadata$cols
pc1n <- vd.pco$eig[1] / sum(vd.pco$eig)
pc2n <- vd.pco$eig[2] / sum(vd.pco$eig)
xl <- paste("PCoA1 (", (pc1n * 10000) %/% 100, "%)", sep = "")
yl <- paste("PCoA2 (", (pc2n * 10000) %/% 100, "%)", sep = "")
p <-
  ggplot(pcoadata,
         aes(x = vd.pco.points...1., y = vd.pco.points...2., color = Environment)) +
  geom_point(size = 3, alpha = .8) +
  scale_color_manual(values = c("#00B81F", "#F8766D", "#C59900", "#592EFF")) +
  geom_text(
    aes(label = samples.orig),
    hjust = -0.6,
    vjust = 0,
    key_glyph = "point"
  )

p + theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.background = element_rect(fill = alpha('white', 0.2))
) + labs(fill = "Samples subtype") + xlab(xl) + ylab(yl) + xlim(-0.6, 0.8)
dev.off()


# ppm normalized ####
mothertable <-
  read.table(
    file = input_ppm,
    sep = "\t",
    header = T,
    row.names = 1,
    quote = "",
    stringsAsFactors = FALSE
  )
mothertable <-
  mothertable[which(rowSums(mothertable) != 0), ] # remove zero lines
mothertable <- t(mothertable)
samples.orig <- rownames(mothertable)
mothertable <- mothertable[which(rowSums(mothertable) != 0), ]
dim(mothertable)

# update cols in case certain samples have been filtered out.
cols <- cols[match(rownames(mothertable), samples.orig)]

png(
  "stage2_ppm_subtype_bray.png",
  res = 300,
  width = 2200,
  height = 2400
)
vd <- vegdist(mothertable, method = "bray")
vd.pco <- pco(vd, k = 10)
pcoadata <- data.frame(vd.pco$points[, 1], vd.pco$points[, 2], cols)
Environment <- pcoadata$cols
pc1n <- vd.pco$eig[1] / sum(vd.pco$eig)
pc2n <- vd.pco$eig[2] / sum(vd.pco$eig)
xl <- paste("PCoA1 (", (pc1n * 10000) %/% 100, "%)", sep = "")
yl <- paste("PCoA2 (", (pc2n * 10000) %/% 100, "%)", sep = "")
p <-
  ggplot(pcoadata,
         aes(x = vd.pco.points...1., y = vd.pco.points...2., color = Environment)) +
  geom_point(size = 3, alpha = .8) +
  scale_color_manual(values = c("#00B81F", "#F8766D", "#C59900", "#592EFF")) +
  geom_text(
    aes(label = samples.orig),
    hjust = -0.6,
    vjust = 0,
    key_glyph = "point"
  )

p + theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.background = element_rect(fill = alpha('white', 0.2))
) + labs(fill = "Samples subtype") + xlab(xl) + ylab(yl) + xlim(-0.6, 0.8)
dev.off()


## To combine several plots, save them in different variables (not like above
## all called p) #### Then use the following:

p.ex <-
  ggpubr::ggarrange(
    p.bray,
    p.jac,
    labels = "AUTO",
    common.legend = TRUE,
    legend = "right",
    font.label = list(size = 15, face = "bold"),
    ncol = 2,
    nrow = 1
  )

ggpubr::ggexport(
  p.ex,
  filename = "stage2_NORMALIZATION-METHOD_subtype_DIST-MEASURE.png",
  width = 6000,
  height = 2800,
  res = 600
)


####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  Non-metric multidimensional scaling       ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
#### Notes ####
### https://jkzorz.github.io/2019/06/06/NMDS.html
# (1) Interesting fact about NMDS: is iterative, repeating itself until
#     it finds the best solution. This means that each time you produce
#     an NMDS plot from scratch it may look slightly different!
# (2) Being non-metric, NMDS is rank-based: instead of using the actual
#     values to calculate distances, it uses ranks. So for example, instead
#     of saying that sample A is 5 points away from sample B, and 10 from
#     sample C, you would instead say that: sample A is the “1st” most close
#     sample to B, and sample C is the “2nd” most close.
# (3) NMDS uses a distance matrix as an input. Of course, you can choose
#     different distances such as Bray-Curtis (presence-absence and abundance)
#     or Jaccard (presence-absence only).


#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
Data.filename_RFAM        <- "RFAM_CPM-t.csv"
Data.filename_BacMet      <- "BacMet_CPM-t.csv"
Data.filename_BacMet.RFAM <- "BacMet-RFAM_CPM-t.csv"

input_RFAM        <- file.path(idir, Data.filename_RFAM)
input_BacMet      <- file.path(idir, Data.filename_BacMet)
input_BacMet.RFAM <- file.path(idir, Data.filename_BacMet.RFAM)

rm(Data.filename_RFAM, Data.filename_BacMet, Data.filename_BacMet.RFAM)


#### Load data ####
# complete set containing either/or/and RFAM and BacMet annotations:
d.pc_RFAM         <- read.csv(input_RFAM, sep = ";", header = T)
d.pc_BacMet       <- read.csv(input_BacMet, sep = ";", header = T)
d.pc_BacMet.RFAM  <- read.csv(input_BacMet.RFAM, sep = ";", header = T)


#### Process data ####
# Make a matrix containing only ARG abundance data (no sample names)
d.com_RFAM <- d.pc_RFAM[,8:ncol(d.pc_RFAM)]
m.com_RFAM <- as.matrix(d.com_RFAM)

d.com_BacMet <- d.pc_BacMet[,8:ncol(d.pc_BacMet)]
m.com_BacMet <- as.matrix(d.com_BacMet)

d.com_BacMet.RFAM <- d.pc_BacMet.RFAM[,8:ncol(d.pc_BacMet.RFAM)]
m.com_BacMet.RFAM <- as.matrix(d.com_BacMet.RFAM)


# Generate raw NMDS plot ####
# NOTE that k = 3 means you want 3 dimensions to be returned; default is 2; more dimensions
# help to reduce stress (see below)
# distance measures available include
#     "manhattan", "euclidean", "canberra",
#     "clark", "bray", "kulczynski", "jaccard",
#     "gower", "altGower", "morisita", "horn",
#     "mountford", "raup", "binomial", "chao",
#     "cao" or "mahalanobis"

set.seed(20200729)
nmds.RFAM <- metaMDS(m.com_RFAM, k = 3, distance = "bray")
nmds.RFAM

nmds.BacMet <- metaMDS(m.com_BacMet, k = 3, distance = "bray")
nmds.BacMet

nmds.BacMet.RFAM <- metaMDS(m.com_BacMet.RFAM, k = 3, distance = "bray")
nmds.BacMet.RFAM

# Look at stress value; ideally <0.2; if 0 that means you have a strong outlier
# which you might want to remove.


# Plot the base NMDS plot if you want, though not very informative
plot(nmds.RFAM)
plot(nmds.BacMet)
plot(nmds.BacMet.RFAM)

# Not very informative, therefore export results to plot in ggplot
data.scores.RFAM                  <- as.data.frame(scores(nmds.RFAM$points))
data.scores.RFAM$Sample           <- d.pc_RFAM$SID
data.scores.RFAM$System           <- d.pc_RFAM$system
data.scores.RFAM$Time             <- d.pc_RFAM$time
data.scores.RFAM$Treat            <- d.pc_RFAM$treat
data.scores.RFAM$SysTime          <- d.pc_RFAM$sys.time
data.scores.RFAM$SysTreat         <- d.pc_RFAM$sys.treat
data.scores.RFAM$SysTimeTreat     <- d.pc_RFAM$sys.time.treat
head(data.scores.RFAM)

data.scores.BacMet                  <- as.data.frame(scores(nmds.BacMet$points))
data.scores.BacMet$Sample           <- d.pc_BacMet$SID
data.scores.BacMet$System           <- d.pc_BacMet$system
data.scores.BacMet$Time             <- d.pc_BacMet$time
data.scores.BacMet$Treat            <- d.pc_BacMet$treat
data.scores.BacMet$SysTime          <- d.pc_BacMet$sys.time
data.scores.BacMet$SysTreat         <- d.pc_BacMet$sys.treat
data.scores.BacMet$SysTimeTreat     <- d.pc_BacMet$sys.time.treat
head(data.scores.BacMet)

data.scores.BacMet.RFAM               <- as.data.frame(scores(nmds.BacMet.RFAM$points))
data.scores.BacMet.RFAM$Sample        <- d.pc_BacMet.RFAM$SID
data.scores.BacMet.RFAM$System        <- d.pc_BacMet.RFAM$system
data.scores.BacMet.RFAM$Time          <- d.pc_BacMet.RFAM$time
data.scores.BacMet.RFAM$Treat         <- d.pc_BacMet.RFAM$treat
data.scores.BacMet.RFAM$SysTime       <- d.pc_BacMet.RFAM$sys.time
data.scores.BacMet.RFAM$SysTreat      <- d.pc_BacMet.RFAM$sys.treat
data.scores.BacMet.RFAM$SysTimeTreat  <- d.pc_BacMet.RFAM$sys.time.treat
head(data.scores.BacMet.RFAM)


# Plot data with ggplot ####
# Adjust plot titles according to what input data and what k-value was used!
# To adjust
#   - position of data labels inside the plot, use hjust = 2, vjust = -0.5 inside geom_text()
#   - text of data labels inside the plot, modify aes(label = Treat) accordingly
# For colour options available in Brewer, consult http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually
# Nice colour options for gradients include Blues, Greens, Oranges; also nice: Pastel1

k3.nmds_RFAM <- ggplot(data.scores.RFAM, aes(x = MDS1, y = MDS2)) +
  geom_point(size = 7, aes(shape = System, colour = Time)) +
  scale_color_brewer(palette = "Oranges") +
  scale_shape_manual(values = c(18, 17, 15, 19)) +
  geom_text(aes(label = Treat)) +
  guides(shape = guide_legend(override.aes = list(size = 4)),
         colour = guide_legend(override.aes = list(size = 4))) +
  theme(
    axis.text.y = element_text(
      colour = "black",
      size = 12,
      face = "bold"
    ),
    axis.text.x = element_text(
      colour = "black",
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    legend.position = "right",
    axis.title.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(
      face = "bold",
      size = 14,
      colour = "black"
    ),
    legend.title = element_text(
      size = 12,
      colour = "black",
      face = "bold"
    ),
    plot.title = element_text(
      size = 15,
      colour = "black",
      face = "bold"
    ),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 1),
    legend.key = element_blank()
  ) +
  labs(
    x = "NMDS1",
    colour = "Time",
    y = "NMDS2",
    shape = "System"
  )


k3.nmds_BacMet <- ggplot(data.scores.BacMet, aes(x = MDS1, y = MDS2)) +
  geom_point(size = 7, aes(shape = System, colour = Time)) +
  scale_color_brewer(palette = "Oranges") +
  scale_shape_manual(values = c(18, 17, 15, 19)) +
  geom_text(aes(label = Treat)) +
  guides(shape = guide_legend(override.aes = list(size = 4)),
         colour = guide_legend(override.aes = list(size = 4))) +
  theme(
    axis.text.y = element_text(
      colour = "black",
      size = 12,
      face = "bold"
    ),
    axis.text.x = element_text(
      colour = "black",
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    legend.position = "right",
    axis.title.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(
      face = "bold",
      size = 14,
      colour = "black"
    ),
    legend.title = element_text(
      size = 12,
      colour = "black",
      face = "bold"
    ),
    plot.title = element_text(
      size = 15,
      colour = "black",
      face = "bold"
    ),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 1),
    legend.key = element_blank()
  ) +
  labs(
    x = "NMDS1",
    colour = "Time",
    y = "NMDS2",
    shape = "System"
  )

k3.nmds_BacMet.RFAM <- ggplot(data.scores.BacMet.RFAM, aes(x = MDS1, y = MDS2)) +
  geom_point(size = 7, aes(shape = System, colour = Time)) +
  scale_color_brewer(palette = "Oranges") +
  scale_shape_manual(values = c(18, 17, 15, 19)) +
  geom_text(aes(label = Treat)) +
  guides(shape = guide_legend(override.aes = list(size = 4)),
         colour = guide_legend(override.aes = list(size = 4))) +
  theme(
    axis.text.y = element_text(
      colour = "black",
      size = 12,
      face = "bold"
    ),
    axis.text.x = element_text(
      colour = "black",
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    legend.position = "right",
    axis.title.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(
      face = "bold",
      size = 14,
      colour = "black"
    ),
    legend.title = element_text(
      size = 12,
      colour = "black",
      face = "bold"
    ),
    plot.title = element_text(
      size = 15,
      colour = "black",
      face = "bold"
    ),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 1),
    legend.key = element_blank()
  ) +
  labs(
    x = "NMDS1",
    colour = "Time",
    y = "NMDS2",
    shape = "System"
  )


k3.nmds_RFAM
k3.nmds_BacMet
k3.nmds_BacMet.RFAM


# For plots with significance ellipses and stress value indicated
k3.nmds.ell.90_RFAM_bray          <-
  k3.nmds_RFAM + stat_ellipse(aes(group = SysTime), level = 0.9) + labs(subtitle = "stress = 0.045")
k3.nmds.ell.90_RFAM_bray

k3.nmds.ell.90_BacMet_bray        <-
  k3.nmds_BacMet + stat_ellipse(aes(group = SysTime), level = 0.9) + labs(subtitle = "stress = 0.044")
k3.nmds.ell.90_BacMet_bray

k3.nmds.ell.90_RFAM.BacMet_bray   <-
  k3.nmds_BacMet.RFAM + stat_ellipse(aes(group = SysTime), level = 0.9) + labs(subtitle = "stress = 0.039")
k3.nmds.ell.90_RFAM.BacMet_bray


# arrange several plots together including A,B,C ... labeling

ggarrange(k3.nmds.ell.90_RFAM_bray,
          k3.nmds.ell.90_BacMet_bray,
          k3.nmds.ell.90_RFAM.BacMet_bray,
          labels = "AUTO", common.legend = TRUE,
          font.label = list(size = 20, color = "black", face = "bold", family = NULL),
          legend = "right",
          ncol = 1, nrow = 3) %>%
  ggexport(filename = "BacMet-RFAM_CPM_Bray.png", width = 4800, height = 8000, res = 600)



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  Calypso & Procrustes Export (script NR.12)####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

# Define paths to use
idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
RData.filename <- "GDC_ZOTU99-Silva_1.RData"
input <- file.path(odir, RData.filename)

load(input)

rm(RData.filename, input)


#### Select samples ####
## Samples used in shotgun approach
d.ZOTU99 <- subset_samples(d.ZOTU99, Time != "2")
d.ZOTU99 <- subset_samples(d.ZOTU99, Time != "3")
d.ZOTU99 <- subset_samples(d.ZOTU99, Time != "4")
d.ZOTU99 <- subset_samples(d.ZOTU99, Time != "6")
d.ZOTU99 <- subset_samples(d.ZOTU99, Sys_Time != "W0")

d.ZOTU99 <- prune_taxa(taxa_sums(d.ZOTU99) > 0, d.ZOTU99)
d.ZOTU99


#### Aggregate technical replicates ####
## ...to get rid of technical replicates I,II,III and thereby match the shotgun
## data, where no tech. reps exist

## With the merge_samples function, the abundance values of merged samples are
## *summed*, so make sure to do any pre-processing to account for differences in
## sequencing effort before merging, unless you want seq. depth effect:

## Quick approach to determine optimal rarefaction level is to look at sample_sums
max(sample_sums(d.ZOTU99))
min(sample_sums(d.ZOTU99))
sort(sample_sums(d.ZOTU99))  # to see which are the lowest depths and pick a cutoff

## To rarefy at a certain level:
rare_cutoff <- 32647         # choose number from above output
d.ZOTU99.rare <-
  rarefy_even_depth(d.ZOTU99, sample.size = rare_cutoff, rngseed = 20211112)
sample_sums(d.ZOTU99.rare)[1:5]

d.ZOTU99.rare.m <- merge_samples(d.ZOTU99.rare, "Sys_Time_Treat")
print(sample_data(d.ZOTU99.rare.m)[, "Sys_Time_Treat"])


## note that
## (1) sample_data has got completely messed up by sample merging -> replace it:

new.SD <-
  read.csv("path-to-your-input-files/SWL_SHOT_sample_data.txt",
           header = TRUE,
           sep = "\t")
row.names(new.SD) <- new.SD$ID
sample_data(d.ZOTU99.rare.m) <- new.SD

## (2) in the otu-table now suddenly taxa (=ZOTUs) are columns instead of rows ->
## correct this:
otu_table(d.ZOTU99.rare.m) <- t(otu_table(d.ZOTU99.rare.m))


## check what the merge has done based on top10 abundance ZOTUs (summing up OTU
## abundances?)

OTUnames10     <- names(sort(taxa_sums(d.ZOTU99.rare), TRUE)[1:10])
ZOTU10         <- prune_taxa(OTUnames10,  d.ZOTU99.rare)

mZOTU10        <- prune_taxa(OTUnames10, d.ZOTU99.rare.m)
manure_samples <-
  sample_names(subset(sample_data(d.ZOTU99.rare), System == "M"))
print(manure_samples)

otu_table(ZOTU10, taxa_are_rows = T)[, manure_samples]
rowSums(otu_table(ZOTU10, taxa_are_rows = T)[, manure_samples])
otu_table(mZOTU10)[, "M0Z"]


#### Agglomerate taxa ####
## To choose a rank for agglomerating use one of the rank_names(physeq)
## [1] "Kingdom"
## [2] "Phylum" 
## [3] "Class"
## [4] "Order"   
## [5] "Family"
## [6] "Genus"
## [7] "Species"

d.ZOTU99.rare.m.glom7 <-
  tax_glom(d.ZOTU99.rare.m,
           taxrank = rank_names(d.ZOTU99.rare.m)[7],
           NArm = TRUE)
d.ZOTU99.rare.m.glom6 <-
  tax_glom(d.ZOTU99.rare.m,
           taxrank = rank_names(d.ZOTU99.rare.m)[6],
           NArm = TRUE)


#### Export for ResistoXplorer (http://www.resistoxplorer.no/faces/home.xhtml)
## for this, refer to script NR.13, option (2)



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  Export OTU plus TAX (script NR.13)        ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
## Export OTU plus TAX for conversion to .biom using biom convert ####
## and subsequent upload to Calypso

## Replace accordingly:
# name of filtered phyloseq object; in my case:
#  - d.ZOTU99        for Calypso
#  - d.ZOTU99.rare.m for ResistoXplorer (get it by running script NR.12)
# name of output file

otu.tax.table <-
  data.frame(
    "#OTU ID" = rownames(phyloseq::otu_table(d.ZOTU99.rare.m)@.Data),
    phyloseq::otu_table(d.ZOTU99.rare.m)@.Data,
    phyloseq::tax_table(d.ZOTU99.rare.m)@.Data,
    check.names = FALSE
  )
rownames(otu.tax.table) <- NULL
head(otu.tax.table, 10)


## Format of taxonomy information depends on downstream platform ####
## EITHER (1) semantics needed for Calypso:
otu.tax.table$Kingdom <- paste0(' k__', otu.tax.table$Kingdom)
otu.tax.table$Phylum <- paste0(' p__', otu.tax.table$Phylum)
otu.tax.table$Class <- paste0(' c__', otu.tax.table$Class)
otu.tax.table$Order <- paste0(' o__', otu.tax.table$Order)
otu.tax.table$Family <- paste0(' f__', otu.tax.table$Family)
otu.tax.table$Genus <- paste0(' g__', otu.tax.table$Genus)
otu.tax.table$Species <- paste0(' s__', otu.tax.table$Species)
head(otu.tax.table, 10)

otu.tax.table <- tidyr::unite(otu.tax.table, taxonomy,
                              Kingdom:Species, sep = ";")

otu.tax.table$taxonomy <- gsub("unknown", "", otu.tax.table$taxonomy)
head(otu.tax.table, 10)


## OR (2) semantics needed for ResistoXplorer:
otu.tax.table <- tidyr::unite(otu.tax.table, taxonomy,
                              Kingdom:Species, sep = "|")

otu.tax.table$taxonomy <- gsub("unknown", "", otu.tax.table$taxonomy)
head(otu.tax.table, 10)


## Export data ####
write.table(otu.tax.table,
            row.names = FALSE,
            sep = "\t",
            file = "path-to-your-output-files/noMiCh_noUA_t015_GDC_ZOTU99.rare.m_Silva.txt")



####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####
###  Procrustes rotations                      ####
####-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-####

#### SETUP ####

## clean/reset environment 
rm(list = ls()) 

## Init
linked_path_input <- "path-to-your-input-files"
linked_path_output <- "path-to-your-output-files"

idir <- file.path(linked_path_input)
odir <- file.path(linked_path_output)

rm(linked_path_input, linked_path_output)


#### Data import ####
Data.filename_Phylum       <- "R_Phylum-Feature_RFAM-BacMet_PCoA_Bray_L.csv"
Data.filename_Class        <- "R_Class-Feature_RFAM-BacMet_PCoA_Bray_L.csv"
Data.filename_Order        <- "R_Order-Feature_RFAM-BacMet_PCoA_Bray_L.csv"
Data.filename_Family       <- "R_Family-Feature_RFAM-BacMet_PCoA_Bray_L.csv"
Data.filename_Genus        <- "R_Genus-Feature_RFAM-BacMet_PCoA_Bray_L.csv"
Data.filename_Feature      <- "R_Feature-Feature_RFAM-BacMet_PCoA_Bray_L.csv"

input_Phylum         <- file.path(idir, Data.filename_Phylum)
input_Class          <- file.path(idir, Data.filename_Class)
input_Order          <- file.path(idir, Data.filename_Order)
input_Family         <- file.path(idir, Data.filename_Family)
input_Genus          <- file.path(idir, Data.filename_Genus)
input_Feature        <- file.path(idir, Data.filename_Feature)

rm(Data.filename_Phylum,
   Data.filename_Class,
   Data.filename_Order,
   Data.filename_Family,
   Data.filename_Genus,
   Data.filename_Feature)

p1 <- read.csv(input_Phylum, sep = ";", header = TRUE)
p2 <- read.csv(input_Class, sep = ";", header = TRUE)
p3 <- read.csv(input_Order, sep = ";", header = TRUE)
p4 <- read.csv(input_Family, sep = ";", header = TRUE)
p5 <- read.csv(input_Genus, sep = ";", header = TRUE)
p6 <- read.csv(input_Feature, sep = ";", header = TRUE)


#### Graphs ####
## Create graphs
theme_set(theme_bw(16))

plot1 <- p1 %>%
  ggplot(aes(dim1, dim2)) +
  geom_point(aes(shape = Data, color = Time), cex = 5) +
  scale_shape_manual(values = c(0, 1)) +
  geom_line(aes(group = paired)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  xlab("PCoA1") +
  ylab("PCoA2") +
  ylim(NA, 0.39) +
  geom_text(aes(label = treatment, vjust = -1.2)) +
  theme(panel.grid = element_blank(), axis.line = element_line(), panel.border = element_blank())

plot2 <- p2 %>%
  ggplot(aes(dim1, dim2)) +
  geom_point(aes(shape = Data, color = Time), cex = 5) +
  scale_shape_manual(values = c(0, 1)) +
  geom_line(aes(group = paired)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  xlab("PCoA1") +
  ylab("PCoA2") +
  ylim(NA, 0.39) +
  geom_text(aes(label = treatment, vjust = -1.2)) +
  theme(panel.grid = element_blank(), axis.line = element_line(), panel.border = element_blank())

plot3 <- p3 %>%
  ggplot(aes(dim1, dim2)) +
  geom_point(aes(shape = Data, color = Time), cex = 5) +
  scale_shape_manual(values = c(0, 1)) +
  geom_line(aes(group = paired)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  xlab("PCoA1") +
  ylab("PCoA2") +
  ylim(NA, 0.39) +
  geom_text(aes(label = treatment, vjust = -1.2)) +
  theme(panel.grid = element_blank(), axis.line = element_line(), panel.border = element_blank())

plot4 <- p4 %>%
  ggplot(aes(dim1, dim2)) +
  geom_point(aes(shape = Data, color = Time), cex = 5) +
  scale_shape_manual(values = c(0, 1)) +
  geom_line(aes(group = paired)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  xlab("PCoA1") +
  ylab("PCoA2") +
  ylim(NA, 0.39) +
  geom_text(aes(label = treatment, vjust = -1.2)) +
  theme(panel.grid = element_blank(), axis.line = element_line(), panel.border = element_blank())

plot5 <- p5 %>%
  ggplot(aes(dim1, dim2)) +
  geom_point(aes(shape = Data, color = Time), cex = 5) +
  scale_shape_manual(values = c(0, 1)) +
  geom_line(aes(group = paired)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  xlab("PCoA1") +
  ylab("PCoA2") +
  ylim(NA, 0.39) +
  geom_text(aes(label = treatment, vjust = -1.2)) +
  theme(panel.grid = element_blank(), axis.line = element_line(), panel.border = element_blank())

plot6 <- p6 %>%
  ggplot(aes(dim1, dim2)) +
  geom_point(aes(shape = Data, color = Time), cex = 5) +
  scale_shape_manual(values = c(0, 1)) +
  geom_line(aes(group = paired)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  xlab("PCoA1") +
  ylab("PCoA2") +
  ylim(NA, 0.39) +
  geom_text(aes(label = treatment, vjust = -1.2)) +
  theme(panel.grid = element_blank(), axis.line = element_line(), panel.border = element_blank())


## Combine and export
ggpubr::ggarrange(plot1, plot2, plot3, plot4, plot5, plot6,
                  labels = "AUTO", common.legend = TRUE, legend = "right",
                  font.label = list(size = 20, face = "bold"),
                  ncol = 2, nrow = 3) %>%
  ggpubr::ggexport(filename = "Procrustes_RFAM-BacMet_CPM_PCoA-Bray_all.png",
                   width = 7000, height = 9000, res = 600)
