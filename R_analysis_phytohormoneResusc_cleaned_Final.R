#######################################################################################
####Bacterial community Analysis for Drought Experiment using Bean and Switchgrass#####
#######################################################################################

#######################################################################################
####Navigate to HPCC remote address "/mnt/research/ShadeLab/WorkingSpace/Bandopadhyay_WorkingSpace/Bean_SW_drought_rhizobiome/FinalCode_Clean_ForGit"
#######################################################################################

#######################################################################################
#####Contaminant, mitochondria, chloroplast removal####################################
#######################################################################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(grid)
library(reshape2)
library(grid)
#install.packages("vctrs")
#update.packages("tidyverse")
#update.packages("scales")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

library(decontam)

otu = read.csv("otu_table_or_99.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
tax = read.csv("taxonomy-rep-seqs-or-99-edited.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("sample-metadata-decontam.csv", sep=",", row.names=1)
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 17610 taxa and 717 samples ]
sample_data() Sample Data:       [ 717 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 17610 taxa by 7 taxonomic ranks ]

sample_names(phyloseq_merged)

##check column names of taxonomy table
colnames(tax_table(phyloseq_merged))

##remove mito, chloroplast and archaea, eukarya
phyloseq_merged_clean <- phyloseq_merged %>%
  subset_taxa(
    Domain == "d__Bacteria" &
      Family   != "f__Chloroplast" &
      Family  != "f__Mitochondria"
  )
phyloseq_merged_clean

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 17331 taxa and 717 samples ]
sample_data() Sample Data:       [ 717 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 17331 taxa by 7 taxonomic ranks ]

head(sample_data(phyloseq_merged_clean))

##inspect library sizes
df <- as.data.frame(sample_data(phyloseq_merged_clean)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phyloseq_merged_clean)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
librarysize<- ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point() ##Supplemental Figure S2
ggsave(filename = "librarySize.tiff", plot = librarysize,
       width = 17,
       height = 15, units = c("cm"),
       dpi = 300)

##identify contaminants by prevalence. in this method, the distribution 
#of the frequency of each sequence feature as a function of the 
#prevalence is used to identify contaminants. In this method, 
#the prevalence (presence/absence across samples) of each sequence feature 
#in true positive samples is compared to the prevalence in negative controls 
#to identify contaminants.

sample_data(phyloseq_merged_clean)$is.neg <- sample_data(phyloseq_merged_clean)$Sample_or_Control == "Control Sample"
contamdf.prev0.5 <- isContaminant(phyloseq_merged_clean, method="prevalence", neg="is.neg", threshold=0.5)
Warning messages:
  1: In .is_contaminant(seqtab, conc = conc, neg = neg, method = method,  :
                          Removed 3 samples with zero total counts (or frequency).
  2: In .is_contaminant(seqtab, conc = conc, neg = neg, method = method,  :
                                                Removed 3 samples with zero total counts (or frequency).

table(contamdf.prev0.5$contaminant)
##FALSE  TRUE 
##17283    48    

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(phyloseq_merged_clean, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev0.5$contaminant)
prevalence<- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") ##Supplemental Figure S3

ggsave(filename = "prevalence.tiff", plot = prevalence,
       width = 17,
       height = 15, units = c("cm"),
       dpi = 300)

write.csv(df.pa, "contaminant-table-dna-cdna.csv")

write.csv(contamdf.prev0.5, "contaminant-prev-0.5-dna-cdna.csv")

##removing contaminants from phyloseq object
df.pa <- read.csv("contaminant-prev-0.5-dna-cdna.csv")
View(df.pa)
subset.df <- subset(df.pa, contaminant== "FALSE")
View(subset.df)
keep.taxa <- as.vector(subset.df$X)
keep.taxa

phyloseq_merged_clean_decontam <- prune_taxa(keep.taxa, phyloseq_merged_clean)
phyloseq_merged_clean_decontam

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 17283 taxa and 717 samples ]
sample_data() Sample Data:       [ 717 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 17283 taxa by 7 taxonomic ranks ]

##export otu table out of phyloseq object 

OTU1 = as(otu_table(phyloseq_merged_clean_decontam), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

write.csv(OTUdf, "OTU_clean.csv")

#remove neg controls from otu table-- might add more rows which are rowsums=0, read back in to rarefy reads to 15K

otu = read.csv("OTU_clean.csv", sep=",", row.names=1)#with no contaminants, mito, chloroplast, neg control 
tax = read.csv("taxonomy-rep-seqs-or-99-edited.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("sample-metadata-decontam.csv", sep=",", row.names=1)

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 17283 taxa and 666 samples ]
sample_data() Sample Data:       [ 666 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 17283 taxa by 7 taxonomic ranks ]


#############################################
#####rarefaction curves-- all samples#######
############################################

otu1 = read.csv("OTU_clean.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
tax1 = read.csv("taxonomy-rep-seqs-or-99-edited.csv", sep=",", row.names=1)
tax1 = as.matrix(tax1)
#otu1 = as.matrix(otu1)
metadata1 = read.csv("sample-metadata-decontam.csv", sep=",", row.names=1)
OTU1 = otu_table(otu1, taxa_are_rows = TRUE)
TAX1 = tax_table(tax1)
meta1 = sample_data(metadata1)

phyloseq_merged1 = phyloseq(OTU1, TAX1, meta1)
phyloseq_merged1


phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 17283 taxa and 666 samples ]
sample_data() Sample Data:       [ 666 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 17283 taxa by 7 taxonomic ranks ]

#col<- cyan4
library(vegan)
library(scales)
library(grid)
library(tidyverse)
rarecurve(t(as.data.frame(otu_table(phyloseq_merged1))), label = FALSE, col="cyan4") ## Supplemental Figure S7

###############################
sample_df <- data.frame(sample_data(phyloseq_merged))
View(sample_df)
write.csv(sample_df, "samples_pre_rarefied.csv")

library(vegan)
#rarecurve(t(otu_table(phyloseq_merged)), step=50)
phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged, sample.size = 12000, rngseed = TRUE, trimOTUs=FALSE)

`set.seed(TRUE)` was used to initialize repeatable random subsampling.
Please record this for your records so others can reproduce.
Try `set.seed(TRUE); .Random.seed` for the full vector
...
77 samples removedbecause they contained fewer reads than `sample.size`.
Up to first five removed samples are: 
  
  DNABE0X1DNABE0X2DNABE0X3DNABE0Y1DNABE0Y2

##7 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: DNA_129DNA_146DNA_163DNA_164DNA_245

sample_df <- data.frame(sample_data(phyloseq_rarefy))

View(sample_df)
write.csv(sample_df, "samples_rarefied.csv")

##use anti-join to select non-present samples in rarefied object. change sample id header to SampleID
#This join is like eg: for df1-df2, it selects all rows from df1 that are not present in df2.
df1<- read.csv("samples_pre_rarefied.csv")
df2 <- read.csv("samples_rarefied.csv")
df= df1 %>% anti_join(df2,by="SampleID")

View(df)

write.csv(df, "samples-lost-rarefaction.csv")

OTU2 = as(otu_table(phyloseq_rarefy), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, "OTU_clean_noneg_rarefied12k.csv")

#######################################################################
###### Divide into DNA and cDNA dataset for rna/dna ratio calculations 
#######################################################################
library(tidyverse)
require(dplyr) 
require(tibble)
dna_cdna<- read.csv("OTU_clean_noneg_rarefied12k_removedsamples.csv", row.names=1) ##new rarefied and cleaned OTU table, with 15k reads per sample, and exact match of DNA RNA samples
dna_cdna = dna_cdna[rowSums(dna_cdna[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(dna_cdna)
dna_cdna ##otus reduced to 13,138
dna_cdna1 <-as_tibble(dna_cdna, rownames="OTUID") 
dna_cdna1
has_rownames(dna_cdna1)
rownames(dna_cdna1)
dna<- dna_cdna1 %>% select(OTUID, starts_with("DNA"))
View(dna) #324 samples
cdna <- dna_cdna1 %>% select(OTUID, starts_with("RNA"))
View(cdna) #324 samples
write.csv(dna, "DNA.csv")
write.csv(cdna, "RNA.csv") ##removed the first column after writing to csv

##already removed all negative samples from DNA and CDNA dataset
##check that all rowsums are >1, remove any rowsums=0 

#subset DNA dataframe to check only those OTUs that appear in at least one sample
DNA<- read.csv("DNA.csv", row.names=1)
View(DNA)## 13138 otus
dnanozero = DNA[rowSums(DNA[])>0, ,drop=FALSE]
dnanozero
View(dnanozero)##6343 otus
write.csv(dnanozero, "dna_noneg_nozero.csv")

dnanozero["Total"] <- rowSums(dnanozero)
dnanozero
dnanozero1 = dnanozero %>% filter(dnanozero$Total== 0)

View(dnanozero1) ##nothing here, all zero totals have been successfully removed


##subset cDNA dataframe to check only those OTUs that appear in at least one sample
CDNA<- read.csv("RNA.csv", row.names=1)
View(CDNA)## 13138 otus
cdnanozero = CDNA[rowSums(CDNA[])>0, ,drop=FALSE]
cdnanozero
View(cdnanozero)## 8950 otus
write.csv(cdnanozero, "rna_noneg_nozero.csv")


##merge to see overlap in OTUs between dna and cdna datasets

A = read.csv("dna_noneg_nozero.csv")

B = read.csv("rna_noneg_nozero.csv")

C<- inner_join(A, B, by = "X") ##2155 taxa across 648 columns, join keeps both sets of data from A and B
y <- plyr::match_df(A, B, on = "X") ##2155 taxa across 324 columns, match_df only keeps DNA columns (values from left table A) 
View(y)
View(C)

##if negative samples removed, then 2155 taxa match between dna and cdna datasets.



#######################################################################################
###### Use DNA and cDNA dataset for rna/dna ratio calculations ########################
#######################################################################################

library(tidyverse)

##making the headers match between dna and rna datasets

dna<- read.csv("DNAcopy.csv", header=TRUE, row.names = 1)
cdna<- read.csv("RNAcopy.csv", header=TRUE, row.names = 1)
View(dna)
View(cdna)
#dna.long<-pivot_longer(dna, !OTUID, names_to = "sample_id", values_to = "Abundance")
#View(dna.long)

dna1<- dna %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="dna", -ASV)#gather columns into key value pairs
View(dna1)
cdna1<- cdna %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="cdna", -ASV)
View(cdna1)
merge.df = full_join(dna1, cdna1, by=c("ASV", "SampleID"))
View(merge.df)


dna_wide <- merge.df %>% select(ASV, dna, SampleID) %>%pivot_wider(names_from="SampleID", values_from="dna")  
View(dna_wide) 

cdna_wide <- merge.df %>% select(ASV, cdna, SampleID) %>%pivot_wider(names_from="SampleID", values_from="cdna")  
View(cdna_wide) 

write.csv(dna_wide, "dna_wide.csv")
write.csv(cdna_wide, "cdna_wide.csv")

dna<- read.csv("DNA_mat_noheader_final.csv", header=FALSE) ##remove OTU names (rownames) and header row with sample names to create matrix for calculation of ratios
cdna <- read.csv("RNA_mat_noheader_final.csv", header=FALSE) ##remove OTU names (rownames) and header row with sample names to create matrix for calculation of ratios
View(dna)
View(cdna)


##WITHOUT phantom taxa
cdna.mat <- as.matrix(cdna)
dna.mat<- as.matrix(dna)
merge.mat<- cdna.mat/dna.mat
View(merge.mat)

merge.mat[!is.finite(merge.mat)] <- 0
View(merge.mat)


write.csv(merge.mat, "ratio_nophantoms.csv")


###Method 1 for ratio calculation
cdna.mat <- as.matrix(cdna)

View(cdna.mat)
dna.mat<- as.matrix(dna)

new.mat= ifelse (dna.mat == 0 & cdna.mat>0, 100, cdna.mat/dna.mat) ##this works
View(new.mat)

write.csv(new.mat, "cdna-dna-ratio-method1.csv")


##selecting ratios greater than equal to 1

new.mat[new.mat < 1] <- 0 
View(new.mat)

new.mat[!is.finite(new.mat)] <- 0
View(new.mat)

write.csv(new.mat, "cdna-dna-allfinite-greaterthanEqone_method1.csv")

#####Method 2 for calculating ratios

dna.mat<- as.matrix(dna)
dna.mat[dna.mat == 0] <- 1

View(dna.mat)

cdna.mat <- as.matrix(cdna)

View(cdna.mat)
merge.mat<- cdna.mat/dna.mat


View(merge.mat)
write.csv(dna.mat, "dna_mat_method2_final.csv")##need this for replacing with dna abun
write.csv(merge.mat, "cdna-dna-ratio-method2_final.csv")

##selecting ratios greater than equal to 1

merge.mat[merge.mat < 1] <- 0 
View(merge.mat)


merge.mat[!is.finite(merge.mat)] <- 0
View(merge.mat)

write.csv(merge.mat, "cdna-dna-allfinite-greaterthanEqone_method2_final.csv") 


##replacing ratios greater than equal to 1 with abund values in **dna table**, using method 2

cdna_ratio<- read.csv("cdna-dna-allfinite-greaterthanEqone_method2_final.csv") ##remove first column in csv file
View(cdna_ratio)
dna_abun <- read.csv("dna_mat_method2_final.csv") ##use this to filter to dna abun #remove first column in csv file
View(dna_abun)
cdna_ratio_mat<- as.matrix(cdna_ratio)
View(cdna_ratio_mat)
dna_abun_mat <- as.matrix(dna_abun)
View(dna_abun_mat)

mask <- cdna_ratio_mat > 0
mask

cdna_ratio_mat[mask] <- dna_abun_mat[mask]
cdna_ratio_mat
View(cdna_ratio_mat)


write.csv(cdna_ratio_mat, "finalactivecdna_method2_changing_to_dnaabun_final.csv")


##remove rowsums equals zero
cdna_active<- read.csv("finalactivecdna_method2_changing_to_dnaabun_final.csv", row.names=1) ##put otu feature IDs into first column

activefinalnozero = cdna_active[rowSums(cdna_active[])>0, ,drop=FALSE]

View(activefinalnozero) ##8887



write.csv(activefinalnozero, "activefinalnozero_method2_changing_to_dnaabun_final.csv")


####################################################################################################
####################### Microbiome analysis of active community data using phyloseq ################
####################################################################################################

otu = read.csv("activefinalnozero_method2_changing_to_dnaabun_final.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 

any(rowSums(otu[])<1)
#FALSE
any(rowSums(otu[])<=1)
#TRUE, singletons present, KEEPing singletons


otu1 = otu[rowSums(otu[])>1, ,drop=FALSE]
otu1 ##no singletons
View(otu1) ##9962 OTUs in this dataset with no singletons (singletons defined here as OTUs present only in one sample in one copy, rowsums==1)
##We will keep the singletons

tax = read.csv("taxonomy-rep-seqs-or-99-edited.csv", sep=",", row.names=1)
tax = as.matrix(tax)

##changing to edaphic metadata for constrained ordinations
metadata = read.csv("sample-metadata-edaphic-constrained.csv", sep=",", row.names=1)

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)
data1<- data.frame(sample_data(metadata))
View(data1)

phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged


##after correcting tables
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 324 samples ]
sample_data() Sample Data:       [ 324 samples by 14 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]
> 


data<- data.frame(sample_data(phyloseq_merged))
View(data)

any(taxa_sums(phyloseq_merged) <= 1) #true

any(taxa_sums(phyloseq_merged) < 1) #false


if (!requireNamespace("BiocManager", quietly = TRUE)) # Install BiocManager if not already installed. 
  install.packages("BiocManager")

BiocManager::install("microbiome")

# From GitHub
devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
microbiome::summarize_phyloseq(phyloseq_merged)

[[1]]
[1] "1] Min. number of reads = 66"

[[2]]
[1] "2] Max. number of reads = 9153"

[[3]]
[1] "3] Total number of reads = 523836"

[[4]]
[1] "4] Average number of reads = 1616.77777777778"

[[5]]
[1] "5] Median number of reads = 1456"

[[6]]
[1] "7] Sparsity = 0.983030038346532"

[[7]]
[1] "6] Any OTU sum to 1 or less? YES"

[[8]]
[1] "8] Number of singletons = 6538"

[[9]]
[1] "9] Percent of OTUs that are singletons \n        (i.e. exactly one read detected across all samples)73.074773667151"

[[10]]
[1] "10] Number of sample variables are: 7"

[[11]]
[1] "Plant"               "Treatment"           "Timepoint"           "Treatment_Timepoint" "Mesocosm"           
[6] "Nucleic_Acid"        "Sample_or_Control"  



# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(phyloseq_merged))
sum(sample_sum_df)
##523836  (with singletons)
# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

smin <- min(sample_sums(phyloseq_merged))
smin # 86
smean <- mean(sample_sums(phyloseq_merged))
smean # 3510.91
smax <- max(sample_sums(phyloseq_merged))
smax # 8318



#### Rarefaction curves

library(scales)
library(vegan)
library(dplyr)
rarecurve(t(otu_table(phyloseq_merged)), step=20, label = FALSE, col="cyan4") ##Supplemental Figure S10

###################################################
## Bacterial active community microbiome analysis
###################################################


########## Stacked barplots ##############

metadata<- read.csv("sample-metadata-edaphic-constrained.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones, this file includes all samples, bean and switchgrass
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
phyloseq_merged_final

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 324 samples ]
sample_data() Sample Data:       [ 324 samples by 14 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]


phyloseq_class <- phyloseq_merged_final %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."

View(phyloseq_class)



devtools::install_github("doehm/evoPalette")
library(evoPalette)
launch_evo_palette()

palette_new41 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#fcfc00","#ffff9e","#ffb7ef","#fa7efc","#ae09ea","#521899",
                  "#a0fffc","#1e0047", "#CBD588", "#599861", "#508578","#FFD700","#FFA500",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD",
                  "#D14285")


install.packages("ggh4x")


# install.packages("devtools")
devtools::install_github("teunbrand/ggh4x")


View(phyloseq_class)

##this is for barplot
phyloseq_class$Treatment <- recode(phyloseq_class$Treatment, 
                            'Field-Soil (PreHormone)' = 'Field-Soil',
                            'Predry (PreHormone)' = 'Predry',
                            'Postdry (PreHormone)' = 'Postdry',
                            'Water-acclimated (PreHormone)' ='Water-acclimated')


phyloseq_class$Treatment_new = factor(phyloseq_class$Treatment, levels=c("Field-Soil","Predry","Postdry", "Water-acclimated", "Abscisic-acid","Salicylic-acid","Methanol-control", "Water-control"))
phyloseq_class$Timepoint_new =factor(phyloseq_class$Timepoint, levels = c("PreHormone", "1-Day", "7-Day", "14-Day"))
View(phyloseq_class)


library(ggh4x)
remotes::install_gitlab('r-packages/rock')

library(rock)

install.packages('preregr')
library(preregr)

phyloseq_class_new<- phyloseq_class %>% filter(Treatment!="Abscisic-acid")
test<- phyloseq_class_new %>% filter(Treatment=="Abscisic-acid")
View(test)
#View(phyloseq_class_new)


phyloseq_class_new<- phyloseq_class %>% filter(Treatment!="Salicylic-acid")

devtools::install_github("teunbrand/ggh4x")
plot_sa<- ggplot(phyloseq_class_new, aes(x=Sample, y=Abundance, fill=Class)) + 
  #facet_grid2(Plant~Treatment_new+Timepoint_new,scales="free", independent="x")+
  theme(strip.text.x = element_text(size = 30, face=c("bold")))+ 
  facet_wrap2(Plant~Treatment_new+Timepoint_new,scales="free_x",nrow=2)+
  theme(strip.text.y = element_text(size = 35))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(palette_new41))+
  #ggtitle("Microbial Composition of treatments as compared to baseline and controls for \nA. Abscisic acid")+
  ggtitle("B. Salicylic acid")+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=35),
                                     axis.title.y=element_text(size=35),
                                     plot.title=element_text(size=45, face="bold")) + 
  theme(legend.title = element_text(size=35))+ theme(legend.text = element_text(size=30))
  
plot_sa 

plot_aba<- ggplot(phyloseq_class_new, aes(x=Sample, y=Abundance, fill=Class)) + 
  #facet_grid2(Plant~Treatment_new+Timepoint_new,scales="free", independent="x")+
  theme(strip.text.x = element_text(size = 30, face=c("bold")))+ 
  facet_wrap2(Plant~Treatment_new+Timepoint_new,scales="free_x",nrow=2)+
  theme(strip.text.y = element_text(size = 35))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(palette_new41))+
  ggtitle("Microbial Composition of treatments as compared to baseline and controls for \nA. Abscisic acid")+
  #ggtitle("B. Salicylic acid")+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=35),
                                     axis.title.y=element_text(size=35),
                                     plot.title=element_text(size=45, face="bold")) + 
  theme(legend.title = element_text(size=35))+ theme(legend.text = element_text(size=30))
plot_aba 

figure = ggpubr::ggarrange(plot_aba, plot_sa, 
                           #labels = c("A", "B"),
                           ncol = 1, nrow = 2,
                           common.legend = T,
                           legend = c("right"))


ggsave(filename = "barplot_class_allsamples_hormone_resus_SA_02.21.25.tiff", plot = plot_sa, ##CHANGED FILE NAMES AS NEEDED
       width = 160,
       height = 50, units = c("cm"),
       dpi = 300, limitsize = FALSE)
ggsave(filename = "barplot_class_allsamples_hormone_resus_ABA_02.21.25.tiff", plot = plot_aba, ##CHANGED FILE NAMES AS NEEDED
       width = 160,
       height = 50, units = c("cm"),
       dpi = 300, limitsize = FALSE)

ggsave(filename = "barplot_class_allsamples_hormone_resus_ABA_SA_02.21.25.tiff", figure, ##CHANGED FILE NAMES AS NEEDED
       width = 160,
       height = 100, units = c("cm"),
       dpi = 300, limitsize = FALSE)



################################################
######### Principal Coordinate analysis #########
#################################################
###pcoa plots by SA and ABA for bean and sw
##############

metadata<- read.csv("sample-metadata-edaphic-constrained.csv")
metadata<- metadata %>% filter(Plant=="Bean", Treatment!="Salicylic-acid")##change the plant and tretament for bean and SW aba and sa plots

keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
phyloseq_merged_final


##bean aba
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 126 samples ]
sample_data() Sample Data:       [ 126 samples by 14 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]


## bean sa
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 127 samples ]
sample_data() Sample Data:       [ 127 samples by 14 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]


##sw aba
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 126 samples ]
sample_data() Sample Data:       [ 126 samples by 14 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]


##sw sa

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 126 samples ]
sample_data() Sample Data:       [ 126 samples by 14 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]


####constrained ordination without mesocosm effect
##bean SA remained the same ordination with and without mesocosm effect included
##SW SA changed so saved both versions

sample_df <- data.frame(sample_data(phyloseq_merged_final)) 
View(sample_df)
cap_ord <- ordinate(
  physeq = phyloseq_merged_final, 
  method = "CAP",
  distance = "bray",
  formula = ~ Treatment + Timepoint)


pcoa_bean<-plot_ordination(
  physeq = phyloseq_merged_final,
  ordination = cap_ord,
  color = "Treatment",
  shape = "Timepoint",
  title = "Constrained ordination - Bean"
) + 
  geom_point(aes(color = Treatment, shape=Timepoint, size=Timepoint)) +
  scale_shape_manual(values=c(16,16,16,16), breaks=c('PreHormone', '1-Day', '7-Day', '14-Day'))+
  scale_color_manual(name="Treatment", limits = c("Field-Soil (PreHormone)", "Predry (PreHormone)", "Postdry (PreHormone)", "Water-acclimated (PreHormone)", "Abscisic-acid", "Methanol-control", "Water-control"), values = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                                                                                                                                                                                                               "#0072B2", "#D55E00", "#CC79A7")) +
  scale_fill_manual(name="Treatment", limits = c("Field-Soil (PreHormone)", "Predry (PreHormone)", "Postdry (PreHormone)", "Water-acclimated (PreHormone)", "Abscisic-acid", "Methanol-control", "Water-control"), values = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                                                                                                                                                                                                              "#0072B2", "#D55E00", "#CC79A7")) +
  scale_size_manual(values=c(2.5,4.2,5.5,6.5), breaks=c('PreHormone', '1-Day', '7-Day', '14-Day'))+
  theme(plot.title = element_text(hjust = 0.5))

pcoa_bean

ggsave(filename="pcoa-sw-sa-withOUT-mesocosmeffect-constrained.TIFF", plot=pcoa_bean, width=6.8, height=6, units="in", dpi=720)
ggsave(filename="pcoa-sw-aba-withOUT-mesocosmeffect-constrained.TIFF", plot=pcoa_bean, width=6.8, height=6, units="in", dpi=720)
ggsave(filename="pcoa-bean-sa-withOUT-mesocosmeffect-constrained.TIFF", plot=pcoa_bean, width=6.8, height=6, units="in", dpi=720)
ggsave(filename="pcoa-bean-aba-withOUT-mesocosmeffect-constrained.TIFF", plot=pcoa_bean, width=6.8, height=6, units="in", dpi=720)

#######################################
##### Statistical anslyses ############
#######################################

##Community differences by crop type permanova

metadata<- read.csv("sample-metadata-decontam.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
phyloseq_merged_final

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 324 samples ]
sample_data() Sample Data:       [ 324 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]

set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_final, method = "bray") 
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged_final)) 
View(sample_df)
set.seed(1)
library(vegan)
adonis2(phyloseq_bray ~ Plant, data = sample_df)

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ Plant, data = sample_df)
Df SumOfSqs      R2      F Pr(>F)    
Plant      1    9.393 0.08274 29.045  0.001 ***
  Residual 322  104.137 0.91726                  
Total    323  113.530 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



##permanova for sw 
##by treatment
metadata<- read.csv("sample-metadata-switchgrass.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones
metadata<- metadata%>% filter(Mesocosm!="SW0")
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
phyloseq_merged_final

##switchgrass
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 162 samples ]
sample_data() Sample Data:       [ 162 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]

##switchgrass minus SW0

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 154 samples ]
sample_data() Sample Data:       [ 154 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]

set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_final, method = "bray") 
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged_final)) 
View(sample_df)
set.seed(1)
library(vegan)
adonis2(phyloseq_bray ~ Treatment, data = sample_df)


##WIHTOUT SW 0
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ Treatment, data = sample_df)
Df SumOfSqs      R2      F Pr(>F)    
Treatment   4    8.076 0.16249 7.2269  0.001 ***
  Residual  149   41.626 0.83751                  
Total     153   49.702 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##WITH SW0

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ Treatment, data = sample_df)
Df SumOfSqs      R2      F Pr(>F)    
Treatment   7   11.995 0.21932 6.1805  0.001 ***
  Residual  154   42.698 0.78068                  
Total     161   54.694 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


library(pairwiseAdonis)
set.seed(1)
pairwise.adonis2(phyloseq_bray ~ Treatment, data = sample_df)

         


##by time (Timepoint_new)
set.seed(1)
adonis2(phyloseq_bray ~ Treatment_Timepoint, data = sample_df)

##WITHOUT SW0
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ Treatment_Timepoint, data = sample_df)
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  12   18.884 0.37994 7.1998  0.001 ***
  Residual            141   30.818 0.62006                  
Total               153   49.702 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> 

##WITH SW0
  
  Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ Treatment_Timepoint, data = sample_df)
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  15   22.803 0.41693 6.9599  0.001 ***
  Residual            146   31.890 0.58307                  
Total               161   54.694 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


set.seed(1)
pairwise.adonis2(phyloseq_bray ~ Treatment_Timepoint, data = sample_df)
##water acclimated and 1-D SA/ABA are not significantly differing in active community in switchgrass
##switchgrass does not change its soil community when exposed to stress signals, stabilizing?

##WITH SW0


$parent_call
[1] "phyloseq_bray ~ Treatment_Timepoint , strata = Null , permutations 999"

$`01 - Field Soil_vs_02 - Predry`
Df SumOfSqs      R2      F Pr(>F)
Treatment_Timepoint  1  0.43867 0.34983 1.6141    0.2
Residual             3  0.81531 0.65017              
Total                4  1.25398 1.00000              

$`01 - Field Soil_vs_05 - Methanol control (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3056 0.29672 5.0628  0.006 **
  Residual            12   3.0946 0.70328                 
Total               13   4.4002 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_06 - ABA (1D)`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1   1.7678 0.5955 19.139  0.002 **
  Residual            13   1.2008 0.4045                 
Total               14   2.9686 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2476 0.25036 4.3416  0.002 **
  Residual            13   3.7357 0.74964                 
Total               14   4.9834 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_10 - ABA (7D)`
Df SumOfSqs      R2     F Pr(>F)   
Treatment_Timepoint  1   1.7984 0.62687 21.84  0.004 **
  Residual            13   1.0705 0.37313                
Total               14   2.8689 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3056 0.27646 4.9672  0.005 **
  Residual            13   3.4170 0.72354                 
Total               14   4.7226 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3880 0.31586 6.0021  0.003 **
  Residual            13   3.0063 0.68414                 
Total               14   4.3943 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2506 0.24983 4.3294  0.004 **
  Residual            13   3.7553 0.75017                 
Total               14   5.0059 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3646 0.29446 5.4257  0.003 **
  Residual            13   3.2695 0.70554                 
Total               14   4.6341 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_15 - SA (14D)`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1   1.3576 0.2946 5.4293  0.002 **
  Residual            13   3.2507 0.7054                 
Total               14   4.6083 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.4255 0.36483 6.8927  0.004 **
  Residual            12   2.4818 0.63517                 
Total               13   3.9073 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_04 - Water acclimated`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1   1.2494 0.2526 4.3936  0.004 **
  Residual            13   3.6970 0.7474                 
Total               14   4.9464 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2802 0.26576 4.7055  0.006 **
  Residual            13   3.5370 0.73424                 
Total               14   4.8172 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3010 0.27478 4.9256  0.008 **
  Residual            13   3.4337 0.72522                 
Total               14   4.7346 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)
Treatment_Timepoint  1  0.97363 0.52461 4.4142    0.1
Residual             4  0.88227 0.47539              
Total                5  1.85590 1.00000              

$`02 - Predry_vs_05 - Methanol control (1D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.1419 0.30044 4.7243  0.014 *
  Residual            11   2.6589 0.69956                
Total               12   3.8008 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_06 - ABA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1  1.46658 0.65717 23.003  0.009 **
  Residual            12  0.76509 0.34283                 
Total               13  2.23167 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1    1.098 0.24966 3.9926  0.018 *
  Residual            12    3.300 0.75034                
Total               13    4.398 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.4857 0.70064 28.085  0.014 *
  Residual            12   0.6348 0.29936                
Total               13   2.1205 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.1430 0.27714 4.6007   0.01 **
  Residual            12   2.9813 0.72286                 
Total               13   4.1243 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_12 - Water control (7D)`
Df SumOfSqs      R2     F Pr(>F)  
Treatment_Timepoint  1   1.2024 0.31869 5.613  0.015 *
  Residual            12   2.5706 0.68131               
Total               13   3.7730 1.00000               
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.0973 0.24844 3.9668  0.011 *
  Residual            12   3.3196 0.75156                
Total               13   4.4169 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.1715 0.29249 4.9609  0.015 *
  Residual            12   2.8339 0.70751                
Total               13   4.0054 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_15 - SA (14D)`
Df SumOfSqs      R2     F Pr(>F)  
Treatment_Timepoint  1   1.1708 0.29374 4.991  0.014 *
  Residual            12   2.8150 0.70626               
Total               13   3.9858 1.00000               
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.2354 0.37648 6.6418  0.014 *
  Residual            11   2.0461 0.62352                
Total               12   3.2815 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.1048 0.25305 4.0652  0.014 *
  Residual            12   3.2613 0.74695                
Total               13   4.3661 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.1257 0.26632 4.3558  0.006 **
  Residual            12   3.1013 0.73368                 
Total               13   4.2270 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_9 - Methanol control (7D)`
Df SumOfSqs     R2    F Pr(>F)   
Treatment_Timepoint  1   1.1417 0.2758 4.57  0.009 **
  Residual            12   2.9980 0.7242               
Total               13   4.1397 1.0000               
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`02 - Predry_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)
Treatment_Timepoint  1  0.96288 0.68315 6.4683    0.1
Residual             3  0.44658 0.31685              
Total                4  1.40947 1.00000              

$`05 - Methanol control (1D)_vs_06 - ABA (1D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.5106 0.33164 10.42  0.001 ***
  Residual            21   3.0444 0.66836                 
Total               22   4.5550 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)
Treatment_Timepoint  1   0.4458 0.07398 1.6778  0.192
Residual            21   5.5793 0.92602              
Total               22   6.0251 1.00000              

$`05 - Methanol control (1D)_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1699 0.42681 15.637  0.001 ***
  Residual            21   2.9141 0.57319                  
Total               22   5.0840 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.7491 0.12465 2.9905   0.03 *
  Residual            21   5.2606 0.87535                
Total               22   6.0097 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2     F Pr(>F)   
Treatment_Timepoint  1   1.1430 0.19072 4.949  0.009 **
  Residual            21   4.8499 0.80928                
Total               22   5.9929 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.4541 0.20617 5.454  0.001 ***
  Residual            21   5.5989 0.79383                 
Total               22   7.0530 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.2782 0.30823 9.3569  0.001 ***
  Residual            21   5.1131 0.69177                  
Total               22   7.3914 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.8633 0.26781 7.6809  0.001 ***
  Residual            21   5.0943 0.73219                  
Total               22   6.9575 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.4388 0.24961 6.6526  0.001 ***
  Residual            20   4.3254 0.75039                  
Total               21   5.7641 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.8862 0.13789 3.3588  0.032 *
  Residual            21   5.5406 0.86211                
Total               22   6.4267 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.5651 0.09505 2.2057  0.083 .
Residual            21   5.3806 0.90495                
Total               22   5.9457 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2080 0.18627 4.8071  0.008 **
  Residual            21   5.2772 0.81373                 
Total               22   6.4853 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3717 0.33477 6.0388  0.007 **
  Residual            12   2.7259 0.66523                 
Total               13   4.0976 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.6647 0.31114 9.937  0.001 ***
  Residual            22   3.6855 0.68886                 
Total               23   5.3502 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.9817 0.66013 42.731  0.001 ***
  Residual            22   1.0203 0.33987                  
Total               23   3.0020 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6392 0.32745 10.711  0.001 ***
  Residual            22   3.3667 0.67255                  
Total               23   5.0060 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.3104 0.43869 17.194  0.001 ***
  Residual            22   2.9561 0.56131                  
Total               23   5.2665 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.7227 0.42359 16.167  0.001 ***
  Residual            22   3.7051 0.57641                  
Total               23   6.4278 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.6958 0.45575 18.422  0.001 ***
  Residual            22   3.2193 0.54425                  
Total               23   5.9151 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   3.1349 0.49483 21.549  0.001 ***
  Residual            22   3.2005 0.50517                  
Total               23   6.3354 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_16 - Water control (14D)`
Df SumOfSqs    R2      F Pr(>F)    
Treatment_Timepoint  1   2.4807 0.505 21.424  0.001 ***
  Residual            21   2.4316 0.495                  
Total               22   4.9122 1.000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7399 0.32301 10.497  0.001 ***
  Residual            22   3.6468 0.67699                  
Total               23   5.3867 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   0.9568 0.21531 6.0367  0.001 ***
  Residual            22   3.4868 0.78469                  
Total               23   4.4435 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.5931 0.43388 16.861  0.001 ***
  Residual            22   3.3834 0.56612                  
Total               23   5.9766 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_03 - Postdry`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1  2.05502 0.7118 32.108  0.005 **
  Residual            13  0.83205 0.2882                 
Total               14  2.88707 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.5115 0.41398 15.541  0.001 ***
  Residual            22   3.5552 0.58602                  
Total               23   6.0667 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.9449 0.13801 3.5222  0.031 *
  Residual            22   5.9017 0.86199                
Total               23   6.8466 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.1235 0.16985 4.5012  0.025 *
  Residual            22   5.4911 0.83015                
Total               23   6.6145 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2     F Pr(>F)   
Treatment_Timepoint  1   1.4307 0.18651 5.044  0.003 **
  Residual            22   6.2400 0.81349                
Total               23   7.6707 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.3421 0.28928 8.9545  0.001 ***
  Residual            22   5.7543 0.71072                  
Total               23   8.0964 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_15 - SA (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.9129 0.2501 7.3374  0.001 ***
  Residual            22   5.7354 0.7499                  
Total               23   7.6483 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2     F Pr(>F)   
Treatment_Timepoint  1   1.6023 0.24393 6.775  0.003 **
  Residual            21   4.9665 0.75607                
Total               22   6.5688 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)
Treatment_Timepoint  1   0.5878 0.08683 2.0918   0.13
Residual            22   6.1817 0.91317              
Total               23   6.7695 1.00000              

$`08 - Water control (1D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.5524 0.08403 2.0182   0.09 .
Residual            22   6.0217 0.91597                
Total               23   6.5741 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2714 0.17684 4.7262   0.01 **
  Residual            22   5.9184 0.82316                 
Total               23   7.1898 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2020 0.26308 4.6411  0.004 **
  Residual            13   3.3670 0.73692                 
Total               14   4.5691 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6464 0.33718 11.192  0.001 ***
  Residual            22   3.2365 0.66282                  
Total               23   4.8829 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.6110 0.48024 20.328  0.001 ***
  Residual            22   2.8258 0.51976                  
Total               23   5.4368 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.8579 0.34198 11.434  0.001 ***
  Residual            22   3.5748 0.65802                  
Total               23   5.4326 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_14 - ABA (14D)`
Df SumOfSqs     R2      F Pr(>F)  
Treatment_Timepoint  1   0.6749 0.1793 4.8063  0.015 *
  Residual            22   3.0890 0.8207                
Total               23   3.7639 1.0000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_15 - SA (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.5136 0.3302 10.846  0.001 ***
  Residual            22   3.0702 0.6698                  
Total               23   4.5837 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.6478 0.53501 24.162  0.001 ***
  Residual            21   2.3013 0.46499                  
Total               22   4.9491 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_04 - Water acclimated`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   2.7077 0.43503 16.94  0.001 ***
  Residual            22   3.5165 0.56497                 
Total               23   6.2242 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1889 0.39473 14.347  0.001 ***
  Residual            22   3.3565 0.60527                  
Total               23   5.5454 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.3656 0.29567 9.2354  0.001 ***
  Residual            22   3.2531 0.70433                  
Total               23   4.6188 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1  2.08497 0.74818 38.624  0.005 **
  Residual            13  0.70176 0.25182                 
Total               14  2.78673 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.3342 0.20505 5.6748  0.001 ***
  Residual            22   5.1723 0.79495                  
Total               23   6.5065 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3688 0.18776 5.0855  0.003 **
  Residual            22   5.9213 0.81224                 
Total               23   7.2900 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7248 0.24088 6.9809  0.001 ***
  Residual            22   5.4355 0.75912                  
Total               23   7.1603 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6064 0.22873 6.5244  0.001 ***
  Residual            22   5.4167 0.77127                  
Total               23   7.0230 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_16 - Water control (14D)`
Df SumOfSqs     R2     F Pr(>F)    
Treatment_Timepoint  1   1.6655 0.2638 7.525  0.001 ***
  Residual            21   4.6478 0.7362                 
Total               22   6.3132 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2157 0.17174 4.5617   0.01 **
  Residual            22   5.8629 0.82826                 
Total               23   7.0786 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.9489 0.14266 3.6607  0.015 *
  Residual            22   5.7029 0.85734                
Total               23   6.6519 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3232 0.19114 5.1988  0.006 **
  Residual            22   5.5996 0.80886                 
Total               23   6.9229 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.3592 0.30839 5.7966  0.001 ***
  Residual            13   3.0482 0.69161                  
Total               14   4.4074 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.2472 0.18455 4.979  0.001 ***
  Residual            22   5.5106 0.81545                 
Total               23   6.7578 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.2106 0.30552 9.6784  0.001 ***
  Residual            22   5.0249 0.69448                  
Total               23   7.2354 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.8349 0.26822 8.0637  0.001 ***
  Residual            22   5.0060 0.73178                  
Total               23   6.8409 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)
Treatment_Timepoint  1   0.3551 0.07734 1.7602  0.183
Residual            21   4.2371 0.92266              
Total               22   4.5923 1.00000              

$`12 - Water control (7D)_vs_04 - Water acclimated`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1   1.3224 0.1952 5.3361  0.008 **
  Residual            22   5.4523 0.8048                 
Total               23   6.7747 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3154 0.19908 5.4683  0.008 **
  Residual            22   5.2923 0.80092                 
Total               23   6.6077 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2132 0.18949 5.1436  0.008 **
  Residual            22   5.1890 0.81051                 
Total               23   6.4021 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.5010 0.36269 7.3981  0.001 ***
  Residual            13   2.6376 0.63731                  
Total               14   4.1386 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.9499 0.14128 3.6195  0.017 *
  Residual            22   5.7738 0.85872                
Total               23   6.7238 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.5642 0.08928 2.1566  0.051 .
Residual            22   5.7550 0.91072                
Total               23   6.3191 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6466 0.24826 6.9351  0.001 ***
  Residual            21   4.9861 0.75174                  
Total               22   6.6327 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6076 0.20587 5.7033  0.001 ***
  Residual            22   6.2013 0.79413                  
Total               23   7.8089 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6547 0.21501 6.0259  0.001 ***
  Residual            22   6.0413 0.78499                  
Total               23   7.6960 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)
Treatment_Timepoint  1   0.4432 0.06945 1.6419  0.195
Residual            22   5.9379 0.93055              
Total               23   6.3811 1.00000              

$`13 - Methanol control (14D)_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.4322 0.29721 5.4978  0.002 **
  Residual            13   3.3866 0.70279                 
Total               14   4.8188 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.9871 0.15777 4.1212  0.011 *
  Residual            22   5.2692 0.84223                
Total               23   6.2563 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.2969 0.33792 10.718  0.001 ***
  Residual            21   4.5003 0.66208                  
Total               22   6.7973 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.0935 0.26809 8.0584  0.001 ***
  Residual            22   5.7155 0.73191                  
Total               23   7.8091 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1607 0.28002 8.5565  0.001 ***
  Residual            22   5.5555 0.71998                  
Total               23   7.7162 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.1972 0.18005 4.8308  0.005 **
  Residual            22   5.4522 0.81995                 
Total               23   6.6494 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.6258 0.35916 7.2858  0.005 **
  Residual            13   2.9008 0.64084                 
Total               14   4.5266 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1376 0.32295 10.017  0.001 ***
  Residual            21   4.4815 0.67705                  
Total               22   6.6191 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.0576 0.26535 7.9463  0.001 ***
  Residual            22   5.6967 0.73465                  
Total               23   7.7543 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1394 0.27871 8.5009  0.001 ***
  Residual            22   5.5367 0.72129                  
Total               23   7.6760 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_9 - Methanol control (7D)`
Df SumOfSqs     R2      F Pr(>F)  
Treatment_Timepoint  1   0.7241 0.1176 2.9321   0.03 *
  Residual            22   5.4333 0.8824                
Total               23   6.1575 1.0000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.5599 0.35119 7.0366  0.003 **
  Residual            13   2.8820 0.64881                 
Total               14   4.4419 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`16 - Water control (14D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6414 0.24986 6.9948  0.001 ***
  Residual            21   4.9278 0.75014                  
Total               22   6.5691 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`16 - Water control (14D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.6250 0.25419 7.1574  0.002 **
  Residual            21   4.7678 0.74581                 
Total               22   6.3928 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`16 - Water control (14D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.5634 0.25104 7.0389  0.001 ***
  Residual            21   4.6644 0.74896                  
Total               22   6.2279 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`16 - Water control (14D)_vs_03 - Postdry`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1   1.6078 0.4321 9.1305  0.005 **
  Residual            12   2.1131 0.5679                 
Total               13   3.7208 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.6269 0.09485 2.3053   0.09 .
Residual            22   5.9829 0.90515                
Total               23   6.6099 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3805 0.19015 5.1655  0.006 **
  Residual            22   5.8796 0.80985                 
Total               23   7.2601 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.1052 0.24929 4.3169  0.003 **
  Residual            13   3.3282 0.75071                 
Total               14   4.4335 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3970 0.19631 5.3736  0.005 **
  Residual            22   5.7196 0.80369                 
Total               23   7.1167 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2209 0.27816 5.0094  0.004 **
  Residual            13   3.1682 0.72184                 
Total               14   4.3891 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`9 - Methanol control (7D)_vs_03 - Postdry`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3691 0.30878 5.8072  0.004 **
  Residual            13   3.0649 0.69122                 
Total               14   4.4341 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

attr(,"class")
[1] "pwadstrata" "list" 



##by mesocosm

set.seed(1)
adonis2(phyloseq_bray ~ Mesocosm, data = sample_df)


##WITHOUT SW0
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ Mesocosm, data = sample_df)
Df SumOfSqs      R2      F Pr(>F)   
Mesocosm   3    2.668 0.05367 2.8359  0.002 **
  Residual 150   47.035 0.94633                 
Total    153   49.702 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

set.seed(1)
pairwise.adonis2(phyloseq_bray ~ Mesocosm, data = sample_df)


##by mesocoms without SW0

$parent_call
[1] "phyloseq_bray ~ Mesocosm , strata = Null , permutations 999"

$SW1_vs_SW2
Df SumOfSqs      R2      F Pr(>F)  
Mesocosm  1   0.8527 0.03335 2.5875  0.033 *
  Residual 75  24.7155 0.96665                
Total    76  25.5681 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$SW1_vs_SW3
Df SumOfSqs      R2      F Pr(>F)    
Mesocosm  1   1.7107 0.06803 5.4018  0.001 ***
  Residual 74  23.4353 0.93197                  
Total    75  25.1460 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$SW1_vs_SW4
Df SumOfSqs      R2      F Pr(>F)    
Mesocosm  1   1.8269 0.07113 5.7431  0.001 ***
  Residual 75  23.8571 0.92887                  
Total    76  25.6840 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$SW2_vs_SW3
Df SumOfSqs      R2      F Pr(>F)
Mesocosm  1   0.3568 0.01516 1.1547   0.31
Residual 75  23.1774 0.98484              
Total    76  23.5342 1.00000              

$SW2_vs_SW4
Df SumOfSqs      R2      F Pr(>F)
Mesocosm  1   0.3439 0.01436 1.1075  0.344
Residual 76  23.5992 0.98564              
Total    77  23.9431 1.00000              

$SW3_vs_SW4
Df SumOfSqs      R2      F Pr(>F)
Mesocosm  1   0.2623 0.01161 0.8813  0.543
Residual 75  22.3191 0.98839              
Total    76  22.5813 1.00000              

attr(,"class")
[1] "pwadstrata" "list"

##mesocosms are different from each other, block on mesocosms?



##bean crop 


metadata<- read.csv("sample-metadata-bean.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones
metadata<- metadata %>% filter(Mesocosm!="BE0")


keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
phyloseq_merged_final

#bean
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 162 samples ]
sample_data() Sample Data:       [ 162 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]

##bean minus BE0

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 153 samples ]
sample_data() Sample Data:       [ 153 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]

set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_final, method = "bray") 
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged_final)) 
View(sample_df)
set.seed(1)
library(vegan)
set.seed(1)
adonis2(phyloseq_bray ~ Treatment, data = sample_df)

##WITH BE0
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ Treatment, data = sample_df)
Df SumOfSqs      R2      F Pr(>F)    
Treatment   7   12.436 0.25153 7.3932  0.001 ***
  Residual  154   37.007 0.74847                  
Total     161   49.443 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


library(pairwiseAdonis)
set.seed(1)
pairwise.adonis2(phyloseq_bray ~ Treatment, data = sample_df)




set.seed(1)
adonis2(phyloseq_bray ~ Treatment_Timepoint, data = sample_df)

##WITH BE0
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ Treatment_Timepoint, data = sample_df)
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  15   19.538 0.39517 6.3593  0.001 ***
  Residual            146   29.905 0.60483                  
Total               161   49.443 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



set.seed(1)
pairwise.adonis2(phyloseq_bray ~ Treatment_Timepoint, data = sample_df)
$parent_call
[1] "phyloseq_bray ~ Treatment_Timepoint , strata = Null , permutations 999"

$`03 - Postdry_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7336 0.33537 7.5689  0.001 ***
  Residual            15   3.4357 0.66463                  
Total               16   5.1693 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_05 - Methanol control (1D)`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1   1.8786 0.3857 9.4182  0.002 **
  Residual            15   2.9920 0.6143                 
Total               16   4.8706 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_06 - ABA (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.0158 0.42322 11.006  0.001 ***
  Residual            15   2.7472 0.57678                  
Total               16   4.7631 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_07 - SA (1D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.7755 0.3388 7.6861  0.001 ***
  Residual            15   3.4651 0.6612                  
Total               16   5.2406 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_09 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.8619 0.37856 9.1375  0.001 ***
  Residual            15   3.0565 0.62144                  
Total               16   4.9185 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3240 0.21685 4.1535  0.002 **
  Residual            15   4.7813 0.78315                 
Total               16   6.1053 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.5884 0.28256 5.9076  0.001 ***
  Residual            15   4.0331 0.71744                  
Total               16   5.6214 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.9698 0.41143 10.486  0.001 ***
  Residual            15   2.8178 0.58857                  
Total               16   4.7876 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3834 0.26309 4.6412  0.002 **
  Residual            13   3.8749 0.73691                 
Total               14   5.2583 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7165 0.32648 6.7862  0.001 ***
  Residual            14   3.5411 0.67352                  
Total               15   5.2576 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7024 0.31245 6.8167  0.001 ***
  Residual            15   3.7461 0.68755                  
Total               16   5.4485 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.2155 0.50711 15.433  0.001 ***
  Residual            15   2.1534 0.49289                  
Total               16   4.3690 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.5794 0.28406 5.9514  0.001 ***
  Residual            15   3.9807 0.71594                  
Total               16   5.5601 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_01 - Field Soil`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1  0.76033 0.33911 2.5656  0.055 .
Residual             5  1.48179 0.66089                
Total                6  2.24211 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`03 - Postdry_vs_02 - Predry`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1  0.74555 0.33653 2.5361  0.034 *
  Residual             5  1.46986 0.66347                
Total                6  2.21542 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_05 - Methanol control (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.6756 0.14756 3.8082  0.002 **
  Residual            22   3.9032 0.85244                 
Total               23   4.5788 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_06 - ABA (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   0.6572 0.15228 3.9518  0.001 ***
  Residual            22   3.6584 0.84772                  
Total               23   4.3155 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.7886 0.15268 3.9643  0.005 **
  Residual            22   4.3762 0.84732                 
Total               23   5.1648 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_09 - Methanol control (7D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.1158 0.2195 6.1871  0.001 ***
  Residual            22   3.9677 0.7805                  
Total               23   5.0835 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.9960 0.14891 3.8491  0.018 *
  Residual            22   5.6924 0.85109                
Total               23   6.6884 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.3056 0.20891 5.8097  0.001 ***
  Residual            22   4.9442 0.79109                  
Total               23   6.2498 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_12 - Water control (7D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.4004 0.27302 8.262  0.001 ***
  Residual            22   3.7289 0.72698                 
Total               23   5.1293 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.4865 0.23699 6.2119  0.001 ***
  Residual            20   4.7861 0.76301                  
Total               21   6.2726 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1694 0.20801 5.5155  0.001 ***
  Residual            21   4.4523 0.79199                  
Total               22   5.6217 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.5048 0.24421 7.1085  0.001 ***
  Residual            22   4.6572 0.75579                  
Total               23   6.1620 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7613 0.36497 12.644  0.001 ***
  Residual            22   3.0646 0.63503                  
Total               23   4.8258 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.6148 0.11164 2.7648  0.041 *
  Residual            22   4.8919 0.88836                
Total               23   5.5067 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_01 - Field Soil`
Df SumOfSqs      R2     F Pr(>F)  
Treatment_Timepoint  1   1.1602 0.32652 5.818  0.019 *
  Residual            12   2.3929 0.67348               
Total               13   3.5531 1.00000               
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_02 - Predry`
Df SumOfSqs     R2      F Pr(>F)  
Treatment_Timepoint  1   1.1606 0.3277 5.8492  0.017 *
  Residual            12   2.3810 0.6723                
Total               13   3.5416 1.0000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_06 - ABA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.4698 0.12751 3.2153  0.002 **
  Residual            22   3.2147 0.87249                 
Total               23   3.6845 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.5442 0.12156 3.0444   0.01 **
  Residual            22   3.9325 0.87844                 
Total               23   4.4767 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_09 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.4045 0.10297 2.5254  0.006 **
  Residual            22   3.5240 0.89703                 
Total               23   3.9285 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.0996 0.17321 4.6091  0.002 **
  Residual            22   5.2488 0.82679                 
Total               23   6.3484 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.2285 0.21444 6.0054  0.001 ***
  Residual            22   4.5005 0.78556                  
Total               23   5.7290 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6104 0.32895 10.784  0.001 ***
  Residual            22   3.2853 0.67105                  
Total               23   4.8957 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1538 0.20993 5.3142  0.001 ***
  Residual            20   4.3424 0.79007                  
Total               21   5.4962 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1774 0.22704 6.1683  0.001 ***
  Residual            21   4.0086 0.77296                  
Total               22   5.1860 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1230 0.21043 5.8633  0.001 ***
  Residual            22   4.2135 0.78957                  
Total               23   5.3365 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1010 0.44495 17.636  0.001 ***
  Residual            22   2.6209 0.55505                  
Total               23   4.7219 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.8078 0.15369 3.9952  0.002 **
  Residual            22   4.4482 0.84631                 
Total               23   5.2560 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_01 - Field Soil`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.2556 0.39178 7.7297  0.015 *
  Residual            12   1.9492 0.60822                
Total               13   3.2048 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_02 - Predry`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1   1.2627 0.3946 7.8216  0.009 **
  Residual            12   1.9373 0.6054                 
Total               13   3.2001 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.5217 0.12393 3.1121  0.016 *
  Residual            22   3.6877 0.87607                
Total               23   4.2094 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_09 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.0041 0.23442 6.7364  0.001 ***
  Residual            22   3.2792 0.76558                  
Total               23   4.2833 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_10 - ABA (7D)`
Df SumOfSqs      R2     F Pr(>F)  
Treatment_Timepoint  1   0.9187 0.15512 4.039  0.012 *
  Residual            22   5.0040 0.84488               
Total               23   5.9227 1.00000               
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_11 - SA (7D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   0.9203 0.1778 4.7575  0.001 ***
  Residual            22   4.2557 0.8222                  
Total               23   5.1760 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.6736 0.35502 12.11  0.001 ***
  Residual            22   3.0405 0.64498                 
Total               23   4.7141 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7564 0.30003 8.5728  0.001 ***
  Residual            20   4.0976 0.69997                  
Total               21   5.8540 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_14 - ABA (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.0272 0.2144 5.7312  0.001 ***
  Residual            21   3.7638 0.7856                  
Total               22   4.7910 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_15 - SA (14D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.6696 0.29611 9.255  0.001 ***
  Residual            22   3.9687 0.70389                 
Total               23   5.6383 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_16 - Water control (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   2.0479 0.4629 18.961  0.001 ***
  Residual            22   2.3761 0.5371                  
Total               23   4.4240 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   0.8596 0.16978 4.4988  0.001 ***
  Residual            22   4.2034 0.83022                  
Total               23   5.0630 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_01 - Field Soil`
Df SumOfSqs      R2     F Pr(>F)   
Treatment_Timepoint  1   1.2849 0.42982 9.046   0.01 **
  Residual            12   1.7045 0.57018                
Total               13   2.9893 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_02 - Predry`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2896 0.43245 9.1434  0.005 **
  Residual            12   1.6925 0.56755                 
Total               13   2.9822 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_09 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   0.7874 0.16458 4.3341  0.001 ***
  Residual            22   3.9970 0.83542                  
Total               23   4.7845 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_10 - ABA (7D)`
Df SumOfSqs     R2      F Pr(>F)  
Treatment_Timepoint  1   0.8968 0.1355 3.4483  0.026 *
  Residual            22   5.7218 0.8645                
Total               23   6.6186 1.0000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.6083 0.10897 2.6906  0.054 .
Residual            22   4.9736 0.89103                
Total               23   5.5818 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.8696 0.33221 10.944  0.001 ***
  Residual            22   3.7583 0.66779                  
Total               23   5.6279 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.2818 0.21023 5.3239  0.001 ***
  Residual            20   4.8154 0.78977                  
Total               21   6.0972 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.2199 0.21396 5.7164  0.001 ***
  Residual            21   4.4817 0.78604                  
Total               22   5.7016 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1413 0.19584 5.3576  0.001 ***
  Residual            22   4.6866 0.80416                  
Total               23   5.8279 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.4502 0.44195 17.423  0.001 ***
  Residual            22   3.0939 0.55805                  
Total               23   5.5442 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_08 - Water control (1D)`
Df SumOfSqs     R2      F Pr(>F)  
Treatment_Timepoint  1   0.7668 0.1348 3.4277  0.019 *
  Residual            22   4.9212 0.8652                
Total               23   5.6880 1.0000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_01 - Field Soil`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.1637 0.32452 5.7651  0.025 *
  Residual            12   2.4223 0.67548                
Total               13   3.5860 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_02 - Predry`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.1655 0.32594 5.8025  0.022 *
  Residual            12   2.4104 0.67406                
Total               13   3.5759 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.2154 0.18616 5.0324  0.001 ***
  Residual            22   5.3133 0.81384                  
Total               23   6.5286 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_11 - SA (7D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.5339 0.2515 7.3921  0.001 ***
  Residual            22   4.5650 0.7485                  
Total               23   6.0989 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6233 0.32642 10.661  0.001 ***
  Residual            22   3.3498 0.67358                  
Total               23   4.9731 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.7080 0.13842 3.2132  0.003 **
  Residual            20   4.4069 0.86158                 
Total               21   5.1149 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1460 0.21957 5.9083  0.001 ***
  Residual            21   4.0731 0.78043                  
Total               22   5.2191 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1003 0.20459 5.6586  0.001 ***
  Residual            22   4.2780 0.79541                  
Total               23   5.3784 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.2718 0.45828 18.612  0.001 ***
  Residual            22   2.6854 0.54172                  
Total               23   4.9572 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.9310 0.17103 4.5389  0.002 **
  Residual            22   4.5127 0.82897                 
Total               23   5.4437 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_01 - Field Soil`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.2417 0.38143 7.3995  0.012 *
  Residual            12   2.0137 0.61857                
Total               13   3.2555 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_02 - Predry`
Df SumOfSqs     R2      F Pr(>F)  
Treatment_Timepoint  1   1.2484 0.3841 7.4838  0.011 *
  Residual            12   2.0018 0.6159                
Total               13   3.2503 1.0000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.7562 0.10732 2.6449  0.032 *
  Residual            22   6.2898 0.89268                
Total               23   7.0460 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6692 0.24751 7.2364  0.001 ***
  Residual            22   5.0745 0.75249                  
Total               23   6.7437 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.0833 0.15014 3.5334  0.006 **
  Residual            20   6.1317 0.84986                 
Total               21   7.2149 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)
Treatment_Timepoint  1   0.4851 0.07721 1.7571  0.119
Residual            21   5.7979 0.92279              
Total               22   6.2830 1.00000              

$`10 - ABA (7D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3101 0.17915 4.8016  0.008 **
  Residual            22   6.0028 0.82085                 
Total               23   7.3129 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1787 0.33066 10.868  0.001 ***
  Residual            22   4.4102 0.66934                  
Total               23   6.5889 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.9995 0.13811 3.5254  0.013 *
  Residual            22   6.2375 0.86189                
Total               23   7.2370 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_01 - Field Soil`
Df SumOfSqs     R2   F Pr(>F)  
Treatment_Timepoint  1   0.9658 0.2053 3.1  0.012 *
  Residual            12   3.7385 0.7947             
Total               13   4.7043 1.0000             
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_02 - Predry`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.9596 0.20477 3.0899  0.015 *
  Residual            12   3.7266 0.79523                
Total               13   4.6862 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_12 - Water control (7D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   2.0228 0.3186 10.286  0.001 ***
  Residual            22   4.3263 0.6814                  
Total               23   6.3491 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.4611 0.21347 5.4281  0.001 ***
  Residual            20   5.3834 0.78653                  
Total               21   6.8445 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_14 - ABA (14D)`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1   1.0838 0.1767 4.5072  0.002 **
  Residual            21   5.0496 0.8233                 
Total               22   6.1334 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2559 0.19291 5.2584  0.002 **
  Residual            22   5.2545 0.80709                 
Total               23   6.5105 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.5340 0.40898 15.224  0.001 ***
  Residual            22   3.6619 0.59102                  
Total               23   6.1959 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_08 - Water control (1D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.4707 0.2113 5.8942  0.001 ***
  Residual            22   5.4892 0.7887                  
Total               23   6.9599 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_01 - Field Soil`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.0827 0.26582 4.3448  0.022 *
  Residual            12   2.9903 0.73418                
Total               13   4.0729 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_02 - Predry`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.0863 0.26726 4.3768  0.025 *
  Residual            12   2.9783 0.73274                
Total               13   4.0647 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.4819 0.26227 7.1104  0.001 ***
  Residual            20   4.1681 0.73773                  
Total               21   5.6500 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7109 0.30853 9.3699  0.001 ***
  Residual            21   3.8344 0.69147                  
Total               22   5.5452 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.9147 0.32158 10.428  0.001 ***
  Residual            22   4.0393 0.67842                  
Total               23   5.9540 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_16 - Water control (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1  0.43142 0.1499 3.8792  0.001 ***
  Residual            22  2.44666 0.8501                  
Total               23  2.87808 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   0.9853 0.18735 5.0719  0.001 ***
  Residual            22   4.2740 0.81265                  
Total               23   5.2593 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_01 - Field Soil`
Df SumOfSqs    R2      F Pr(>F)  
Treatment_Timepoint  1   1.2592 0.415 8.5127  0.012 *
  Residual            12   1.7750 0.585                
Total               13   3.0342 1.000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_02 - Predry`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.2658 0.41792 8.6156  0.016 *
  Residual            12   1.7631 0.58208                
Total               13   3.0290 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.2788 0.20725 4.9673  0.001 ***
  Residual            19   4.8915 0.79275                  
Total               20   6.1703 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.1511 0.18425 4.5174  0.002 **
  Residual            20   5.0964 0.81575                 
Total               21   6.2475 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_16 - Water control (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.9666 0.3595 11.226  0.001 ***
  Residual            20   3.5038 0.6405                  
Total               21   5.4704 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2029 0.18409 4.5126  0.003 **
  Residual            20   5.3311 0.81591                 
Total               21   6.5339 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_01 - Field Soil`
Df SumOfSqs      R2     F Pr(>F)  
Treatment_Timepoint  1   1.0077 0.26243 3.558  0.037 *
  Residual            10   2.8321 0.73757               
Total               11   3.8398 1.00000               
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_02 - Predry`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.0142 0.26451 3.5963  0.029 *
  Residual            10   2.8202 0.73549                
Total               11   3.8344 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_15 - SA (14D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.5773 0.24879 6.955  0.001 ***
  Residual            21   4.7626 0.75121                 
Total               22   6.3400 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_16 - Water control (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   2.1685 0.4062 14.366  0.001 ***
  Residual            21   3.1700 0.5938                  
Total               22   5.3385 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.4104 0.22011 5.9269  0.001 ***
  Residual            21   4.9973 0.77989                  
Total               22   6.4077 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_01 - Field Soil`
Df SumOfSqs     R2      F Pr(>F)  
Treatment_Timepoint  1   1.1319 0.3118 4.9836  0.017 *
  Residual            11   2.4984 0.6882                
Total               12   3.6303 1.0000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_02 - Predry`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.1337 0.31317 5.0156  0.016 *
  Residual            11   2.4864 0.68683                
Total               12   3.6202 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.4703 0.42262 16.103  0.001 ***
  Residual            22   3.3749 0.57738                  
Total               23   5.8453 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2942 0.19922 5.4731  0.002 **
  Residual            22   5.2022 0.80078                 
Total               23   6.4964 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_01 - Field Soil`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.1196 0.29288 4.9701  0.009 **
  Residual            12   2.7033 0.70712                 
Total               13   3.8229 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_02 - Predry`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.1305 0.29581 5.0408  0.012 *
  Residual            12   2.6914 0.70419                
Total               13   3.8219 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`16 - Water control (14D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.5598 0.30174 9.5068  0.001 ***
  Residual            22   3.6096 0.69826                  
Total               23   5.1694 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`16 - Water control (14D)_vs_01 - Field Soil`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3663 0.55161 14.762  0.007 **
  Residual            12   1.1106 0.44839                 
Total               13   2.4769 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`16 - Water control (14D)_vs_02 - Predry`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.3723 0.55536 14.988  0.013 *
  Residual            12   1.0987 0.44464                
Total               13   2.4710 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_01 - Field Soil`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.0681 0.26662 4.3627  0.011 *
  Residual            12   2.9380 0.73338                
Total               13   4.0061 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_02 - Predry`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.0718 0.26809 4.3955  0.009 **
  Residual            12   2.9260 0.73191                 
Total               13   3.9978 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`01 - Field Soil_vs_02 - Predry`
Df SumOfSqs      R2     F Pr(>F)
Treatment_Timepoint  1  0.25989 0.37831 1.217 0.3333
Residual             2  0.42708 0.62169             
Total                3  0.68697 1.00000             

attr(,"class")
[1] "pwadstrata" "list" 
     
##BY MESOCOSM
set.seed(1)
adonis2(phyloseq_bray ~ Mesocosm, data = sample_df)

##without BE0

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ Mesocosm, data = sample_df)
Df SumOfSqs      R2      F Pr(>F)
Mesocosm   3    1.183 0.02712 1.3846  0.112
Residual 149   42.418 0.97288              
Total    152   43.601 1.00000   



set.seed(1)
pairwise.adonis2(phyloseq_bray ~ Mesocosm, data = sample_df)
     




#######bean  and switchgrass samples combined mantel correlations###########
  
  
# Calculate bray curtis distance matrix

phyloseq_merged ##324 samples and 8947 OTUS

set.seed(1)
phyloseq_bray <- phyloseq::distance(phyloseq_merged, method = "bray")
phyloseq_bray

library("ape")
set.seed(1)
random_tree = rtree(ntaxa(phyloseq_merged), rooted=TRUE, tip.label=taxa_names(phyloseq_merged))
plot(random_tree)
phyloseq_tree = merge_phyloseq(phyloseq_merged, random_tree)
phyloseq_tree

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 324 samples ]
sample_data() Sample Data:       [ 324 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 8887 tips and 8886 internal nodes ]

##########REDO WITH SEQUENCE DATA
###using rep-seqs-dn-99 for generating tree in qiime2 instead of random tree
tree_bean_sw<- read.tree("tree_bean_sw.nwk")

phyloseq_tree = merge_phyloseq(phyloseq_merged, tree_bean_sw)
phyloseq_tree

#calculate weighted unifrac distance matrix
phyloseq_wunifrac <- phyloseq::distance(phyloseq_tree, method = "wunifrac" )

###########################################################
#calculate weighted unifrac distance matrix
phyloseq_wunifrac <- phyloseq::distance(phyloseq_tree, method = "wunifrac" )

##mantel tests to test correlation between wunifrac and BC (https://jkzorz.github.io/2019/07/08/mantel-test.html)
library(scales)
library(reshape2)
library(grid)
set.seed(1)
mantel_cor = mantel(phyloseq_bray, phyloseq_wunifrac, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_cor


Mantel statistic based on Spearmans rank correlation rho 

Call:
mantel(xdis = phyloseq_bray, ydis = phyloseq_wunifrac, method = "spearman",      permutations = 9999, na.rm = TRUE) 

Mantel statistic r: 0.6455 
      Significance: 1e-04 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0340 0.0439 0.0527 0.0622 
Permutation: free
Number of permutations: 9999


mantel_cor = mantel(phyloseq_bray, phyloseq_wunifrac, method = "pearson", permutations = 9999, na.rm = TRUE)
mantel_cor

Mantel statistic based on Pearsons product-moment correlation 

Call:
mantel(xdis = phyloseq_bray, ydis = phyloseq_wunifrac, method = "pearson",      permutations = 9999, na.rm = TRUE) 

Mantel statistic r: 0.7223 
      Significance: 1e-04 

Upper quantiles of permutations (null model):
   90%    95%  97.5%    99% 
0.0277 0.0360 0.0441 0.0526 
Permutation: free
Number of permutations: 9999


##permanova with wunifrac, comparing bean and switchgrass

sample_df <- data.frame(sample_data(phyloseq_merged)) 
View(sample_df)
set.seed(1)
adonis2(phyloseq_wunifrac ~ Plant, data = sample_df)

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_wunifrac ~ Plant, data = sample_df)
Df SumOfSqs      R2      F Pr(>F)    
Plant      1    4.907 0.08804 31.084  0.001 ***
  Residual 322   50.831 0.91196                  
Total    323   55.738 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##plot correlation

plot<- plot(phyloseq_bray,phyloseq_wunifrac,pch=16,cex=0.5,col="black",bty="l")



####################################################
########## Interaction Permanovas ##################
####################################################

##Switchgrass interaction (using permanova)


#metadata<- read.csv("sample-metadata-sw-no-sw0-no-watacclim.csv")
metadata<- read.csv("sample-metadata-sw-no-sw0.csv")
#metadata <- metadata %>% filter (Treatment!="Salicylic-acid")
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_sw <- prune_samples(keep.samples, phyloseq_merged) ##changed to _merged from _rarefy
phyloseq_sw

##without SW0 mesocosm containing predry, field and postdry soils
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 154 samples ]
sample_data() Sample Data:       [ 154 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]


set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_sw, method = "bray")
phyloseq_bray


# Adonis test
sample_df <- data.frame(sample_data(phyloseq_sw))
View(sample_df)
set.seed(1)

perm <- how(nperm = 999)
setBlocks(perm) <- with(sample_df, Mesocosm)

# Run adonis2 (without SW0)
set.seed(1)
adonis2(formula = phyloseq_bray ~ Treatment*Timepoint,
                    permutations = perm, data = sample_df)


Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Blocks:  with(sample_df, Mesocosm) 
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ Treatment * Timepoint, data = sample_df, permutations = perm)
Df SumOfSqs      R2       F Pr(>F)    
Treatment             4    8.076 0.16249  9.2372  0.001 ***
  Timepoint             2    5.589 0.11245 12.7855  0.001 ***
  Treatment:Timepoint   6    5.219 0.10501  3.9797  0.001 ***
  Residual            141   30.818 0.62006                   
Total               153   49.702 1.00000                   
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




# Run pairwiseadonis2 (without SW0)
set.seed(1)
pairwise.adonis2(phyloseq_bray ~ Treatment_Timepoint,
         data = sample_df, strata='Mesocosm')

$parent_call
[1] "phyloseq_bray ~ Treatment_Timepoint , strata = Mesocosm , permutations 999"

$`05 - Methanol control (1D)_vs_06 - ABA (1D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.5106 0.33164 10.42  0.001 ***
  Residual            21   3.0444 0.66836                 
Total               22   4.5550 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)
Treatment_Timepoint  1   0.4458 0.07398 1.6778  0.169
Residual            21   5.5793 0.92602              
Total               22   6.0251 1.00000              

$`05 - Methanol control (1D)_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1699 0.42681 15.637  0.001 ***
  Residual            21   2.9141 0.57319                  
Total               22   5.0840 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.7491 0.12465 2.9905  0.029 *
  Residual            21   5.2606 0.87535                
Total               22   6.0097 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2     F Pr(>F)  
Treatment_Timepoint  1   1.1430 0.19072 4.949  0.012 *
  Residual            21   4.8499 0.80928               
Total               22   5.9929 1.00000               
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.4541 0.20617 5.454  0.001 ***
  Residual            21   5.5989 0.79383                 
Total               22   7.0530 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.2782 0.30823 9.3569  0.001 ***
  Residual            21   5.1131 0.69177                  
Total               22   7.3914 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.8633 0.26781 7.6809  0.001 ***
  Residual            21   5.0943 0.73219                  
Total               22   6.9575 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.4388 0.24961 6.6526  0.001 ***
  Residual            20   4.3254 0.75039                  
Total               21   5.7641 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.8862 0.13789 3.3588  0.003 **
  Residual            21   5.5406 0.86211                 
Total               22   6.4267 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.5651 0.09505 2.2057  0.003 **
  Residual            21   5.3806 0.90495                 
Total               22   5.9457 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.2080 0.18627 4.8071  0.001 ***
  Residual            21   5.2772 0.81373                  
Total               22   6.4853 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.6647 0.31114 9.937  0.001 ***
  Residual            22   3.6855 0.68886                 
Total               23   5.3502 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.9817 0.66013 42.731  0.001 ***
  Residual            22   1.0203 0.33987                  
Total               23   3.0020 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6392 0.32745 10.711  0.001 ***
  Residual            22   3.3667 0.67255                  
Total               23   5.0060 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.3104 0.43869 17.194  0.001 ***
  Residual            22   2.9561 0.56131                  
Total               23   5.2665 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.7227 0.42359 16.167  0.001 ***
  Residual            22   3.7051 0.57641                  
Total               23   6.4278 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.6958 0.45575 18.422  0.001 ***
  Residual            22   3.2193 0.54425                  
Total               23   5.9151 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   3.1349 0.49483 21.549  0.001 ***
  Residual            22   3.2005 0.50517                  
Total               23   6.3354 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_16 - Water control (14D)`
Df SumOfSqs    R2      F Pr(>F)    
Treatment_Timepoint  1   2.4807 0.505 21.424  0.001 ***
  Residual            21   2.4316 0.495                  
Total               22   4.9122 1.000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7399 0.32301 10.497  0.001 ***
  Residual            22   3.6468 0.67699                  
Total               23   5.3867 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   0.9568 0.21531 6.0367  0.001 ***
  Residual            22   3.4868 0.78469                  
Total               23   4.4435 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.5931 0.43388 16.861  0.001 ***
  Residual            22   3.3834 0.56612                  
Total               23   5.9766 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.5115 0.41398 15.541  0.001 ***
  Residual            22   3.5552 0.58602                  
Total               23   6.0667 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.9449 0.13801 3.5222  0.019 *
  Residual            22   5.9017 0.86199                
Total               23   6.8466 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   1.1235 0.16985 4.5012  0.015 *
  Residual            22   5.4911 0.83015                
Total               23   6.6145 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2     F Pr(>F)   
Treatment_Timepoint  1   1.4307 0.18651 5.044  0.003 **
  Residual            22   6.2400 0.81349                
Total               23   7.6707 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.3421 0.28928 8.9545  0.001 ***
  Residual            22   5.7543 0.71072                  
Total               23   8.0964 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_15 - SA (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.9129 0.2501 7.3374  0.001 ***
  Residual            22   5.7354 0.7499                  
Total               23   7.6483 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2     F Pr(>F)   
Treatment_Timepoint  1   1.6023 0.24393 6.775  0.002 **
  Residual            21   4.9665 0.75607                
Total               22   6.5688 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.5878 0.08683 2.0918  0.054 .
Residual            22   6.1817 0.91317                
Total               23   6.7695 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.5524 0.08403 2.0182  0.054 .
Residual            22   6.0217 0.91597                
Total               23   6.5741 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`08 - Water control (1D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2714 0.17684 4.7262  0.002 **
  Residual            22   5.9184 0.82316                 
Total               23   7.1898 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6464 0.33718 11.192  0.001 ***
  Residual            22   3.2365 0.66282                  
Total               23   4.8829 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.6110 0.48024 20.328  0.001 ***
  Residual            22   2.8258 0.51976                  
Total               23   5.4368 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.8579 0.34198 11.434  0.001 ***
  Residual            22   3.5748 0.65802                  
Total               23   5.4326 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_14 - ABA (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   0.6749 0.1793 4.8063  0.001 ***
  Residual            22   3.0890 0.8207                  
Total               23   3.7639 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_15 - SA (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.5136 0.3302 10.846  0.001 ***
  Residual            22   3.0702 0.6698                  
Total               23   4.5837 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.6478 0.53501 24.162  0.001 ***
  Residual            21   2.3013 0.46499                  
Total               22   4.9491 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_04 - Water acclimated`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   2.7077 0.43503 16.94  0.001 ***
  Residual            22   3.5165 0.56497                 
Total               23   6.2242 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1889 0.39473 14.347  0.001 ***
  Residual            22   3.3565 0.60527                  
Total               23   5.5454 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.3656 0.29567 9.2354  0.001 ***
  Residual            22   3.2531 0.70433                  
Total               23   4.6188 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3342 0.20505 5.6748  0.002 **
  Residual            22   5.1723 0.79495                 
Total               23   6.5065 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.3688 0.18776 5.0855  0.001 ***
  Residual            22   5.9213 0.81224                  
Total               23   7.2900 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7248 0.24088 6.9809  0.001 ***
  Residual            22   5.4355 0.75912                  
Total               23   7.1603 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6064 0.22873 6.5244  0.001 ***
  Residual            22   5.4167 0.77127                  
Total               23   7.0230 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_16 - Water control (14D)`
Df SumOfSqs     R2     F Pr(>F)    
Treatment_Timepoint  1   1.6655 0.2638 7.525  0.001 ***
  Residual            21   4.6478 0.7362                 
Total               22   6.3132 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2157 0.17174 4.5617  0.003 **
  Residual            22   5.8629 0.82826                 
Total               23   7.0786 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.9489 0.14266 3.6607  0.009 **
  Residual            22   5.7029 0.85734                 
Total               23   6.6519 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.3232 0.19114 5.1988  0.001 ***
  Residual            22   5.5996 0.80886                  
Total               23   6.9229 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.2472 0.18455 4.979  0.001 ***
  Residual            22   5.5106 0.81545                 
Total               23   6.7578 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.2106 0.30552 9.6784  0.001 ***
  Residual            22   5.0249 0.69448                  
Total               23   7.2354 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.8349 0.26822 8.0637  0.001 ***
  Residual            22   5.0060 0.73178                  
Total               23   6.8409 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)
Treatment_Timepoint  1   0.3551 0.07734 1.7602  0.214
Residual            21   4.2371 0.92266              
Total               22   4.5923 1.00000              

$`12 - Water control (7D)_vs_04 - Water acclimated`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1   1.3224 0.1952 5.3361  0.003 **
  Residual            22   5.4523 0.8048                 
Total               23   6.7747 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3154 0.19908 5.4683  0.002 **
  Residual            22   5.2923 0.80092                 
Total               23   6.6077 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2132 0.18949 5.1436   0.01 **
  Residual            22   5.1890 0.81051                 
Total               23   6.4021 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.9499 0.14128 3.6195  0.008 **
  Residual            22   5.7738 0.85872                 
Total               23   6.7238 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.5642 0.08928 2.1566  0.069 .
Residual            22   5.7550 0.91072                
Total               23   6.3191 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6466 0.24826 6.9351  0.001 ***
  Residual            21   4.9861 0.75174                  
Total               22   6.6327 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6076 0.20587 5.7033  0.001 ***
  Residual            22   6.2013 0.79413                  
Total               23   7.8089 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.6547 0.21501 6.0259  0.002 **
  Residual            22   6.0413 0.78499                 
Total               23   7.6960 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.4432 0.06945 1.6419  0.095 .
Residual            22   5.9379 0.93055                
Total               23   6.3811 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   0.9871 0.15777 4.1212  0.001 ***
  Residual            22   5.2692 0.84223                  
Total               23   6.2563 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.2969 0.33792 10.718  0.001 ***
  Residual            21   4.5003 0.66208                  
Total               22   6.7973 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.0935 0.26809 8.0584  0.001 ***
  Residual            22   5.7155 0.73191                  
Total               23   7.8091 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1607 0.28002 8.5565  0.001 ***
  Residual            22   5.5555 0.71998                  
Total               23   7.7162 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1972 0.18005 4.8308  0.001 ***
  Residual            22   5.4522 0.81995                  
Total               23   6.6494 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1376 0.32295 10.017  0.001 ***
  Residual            21   4.4815 0.67705                  
Total               22   6.6191 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.0576 0.26535 7.9463  0.001 ***
  Residual            22   5.6967 0.73465                  
Total               23   7.7543 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1394 0.27871 8.5009  0.001 ***
  Residual            22   5.5367 0.72129                  
Total               23   7.6760 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_9 - Methanol control (7D)`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1   0.7241 0.1176 2.9321  0.008 **
  Residual            22   5.4333 0.8824                 
Total               23   6.1575 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`16 - Water control (14D)_vs_04 - Water acclimated`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6414 0.24986 6.9948  0.001 ***
  Residual            21   4.9278 0.75014                  
Total               22   6.5691 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`16 - Water control (14D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6250 0.25419 7.1574  0.001 ***
  Residual            21   4.7678 0.74581                  
Total               22   6.3928 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`16 - Water control (14D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.5634 0.25104 7.0389  0.001 ***
  Residual            21   4.6644 0.74896                  
Total               22   6.2279 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.6269 0.09485 2.3053  0.002 **
  Residual            22   5.9829 0.90515                 
Total               23   6.6099 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.3805 0.19015 5.1655  0.001 ***
  Residual            22   5.8796 0.80985                  
Total               23   7.2601 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_9 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.3970 0.19631 5.3736  0.001 ***
  Residual            22   5.7196 0.80369                  
Total               23   7.1167 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

attr(,"class")
[1] "pwadstrata" "list"  



##Bean interaction(permanova)


metadata<- read.csv("sample-metadata-bean-no-bean0.csv")

keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_bean <- prune_samples(keep.samples, phyloseq_merged) ##changed to _merged from _rarefy
phyloseq_bean



##bean without mesocosm 0

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 153 samples ]
sample_data() Sample Data:       [ 153 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]




set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_bean, method = "bray")
phyloseq_bray


# Adonis test
sample_df <- data.frame(sample_data(phyloseq_bean))
View(sample_df)

#set.seed(1)

#perm <- how(nperm = 999)
#setBlocks(perm) <- with(sample_df, Mesocosm)

set.seed(1)

adonis2(phyloseq_bray ~ Treatment*Timepoint,
        data = sample_df)

set.seed(1)
pairwise.adonis2(phyloseq_bray ~ Treatment_Timepoint,
         data = sample_df)





##Run adonis2 (no bean 0)
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ Treatment * Timepoint, data = sample_df)
Df SumOfSqs      R2       F Pr(>F)    
Treatment             4    8.284 0.18999 10.2754  0.001 ***
  Timepoint             2    2.671 0.06126  6.6265  0.001 ***
  Treatment:Timepoint   6    4.431 0.10163  3.6644  0.001 ***
  Residual            140   28.215 0.64713                   
Total               152   43.601 1.00000                   
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




###Run pairwiseadonis2 (no bean 0)
$parent_call
[1] "phyloseq_bray ~ Treatment_Timepoint , strata = Null , permutations 999"

$`04 - Water acclimated_vs_05 - Methanol control (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.6756 0.14756 3.8082  0.002 **
  Residual            22   3.9032 0.85244                 
Total               23   4.5788 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_06 - ABA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.6572 0.15228 3.9518  0.002 **
  Residual            22   3.6584 0.84772                 
Total               23   4.3155 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.7886 0.15268 3.9643  0.006 **
  Residual            22   4.3762 0.84732                 
Total               23   5.1648 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_09 - Methanol control (7D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.1158 0.2195 6.1871  0.001 ***
  Residual            22   3.9677 0.7805                  
Total               23   5.0835 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.9960 0.14891 3.8491  0.016 *
  Residual            22   5.6924 0.85109                
Total               23   6.6884 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.3056 0.20891 5.8097  0.001 ***
  Residual            22   4.9442 0.79109                  
Total               23   6.2498 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_12 - Water control (7D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.4004 0.27302 8.262  0.001 ***
  Residual            22   3.7289 0.72698                 
Total               23   5.1293 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.4865 0.23699 6.2119  0.001 ***
  Residual            20   4.7861 0.76301                  
Total               21   6.2726 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1694 0.20801 5.5155  0.001 ***
  Residual            21   4.4523 0.79199                  
Total               22   5.6217 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.5048 0.24421 7.1085  0.001 ***
  Residual            22   4.6572 0.75579                  
Total               23   6.1620 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7613 0.36497 12.644  0.001 ***
  Residual            22   3.0646 0.63503                  
Total               23   4.8258 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`04 - Water acclimated_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.6148 0.11164 2.7648  0.034 *
  Residual            22   4.8919 0.88836                
Total               23   5.5067 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_06 - ABA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.4698 0.12751 3.2153  0.002 **
  Residual            22   3.2147 0.87249                 
Total               23   3.6845 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.5442 0.12156 3.0444  0.006 **
  Residual            22   3.9325 0.87844                 
Total               23   4.4767 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_09 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.4045 0.10297 2.5254  0.003 **
  Residual            22   3.5240 0.89703                 
Total               23   3.9285 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.0996 0.17321 4.6091  0.003 **
  Residual            22   5.2488 0.82679                 
Total               23   6.3484 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.2285 0.21444 6.0054  0.001 ***
  Residual            22   4.5005 0.78556                  
Total               23   5.7290 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6104 0.32895 10.784  0.001 ***
  Residual            22   3.2853 0.67105                  
Total               23   4.8957 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1538 0.20993 5.3142  0.001 ***
  Residual            20   4.3424 0.79007                  
Total               21   5.4962 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1774 0.22704 6.1683  0.001 ***
  Residual            21   4.0086 0.77296                  
Total               22   5.1860 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1230 0.21043 5.8633  0.001 ***
  Residual            22   4.2135 0.78957                  
Total               23   5.3365 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1010 0.44495 17.636  0.001 ***
  Residual            22   2.6209 0.55505                  
Total               23   4.7219 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`05 - Methanol control (1D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   0.8078 0.15369 3.9952  0.001 ***
  Residual            22   4.4482 0.84631                  
Total               23   5.2560 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_07 - SA (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.5217 0.12393 3.1121   0.01 **
  Residual            22   3.6877 0.87607                 
Total               23   4.2094 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_09 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.0041 0.23442 6.7364  0.001 ***
  Residual            22   3.2792 0.76558                  
Total               23   4.2833 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_10 - ABA (7D)`
Df SumOfSqs      R2     F Pr(>F)   
Treatment_Timepoint  1   0.9187 0.15512 4.039  0.007 **
  Residual            22   5.0040 0.84488                
Total               23   5.9227 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_11 - SA (7D)`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1   0.9203 0.1778 4.7575  0.004 **
  Residual            22   4.2557 0.8222                 
Total               23   5.1760 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.6736 0.35502 12.11  0.001 ***
  Residual            22   3.0405 0.64498                 
Total               23   4.7141 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7564 0.30003 8.5728  0.001 ***
  Residual            20   4.0976 0.69997                  
Total               21   5.8540 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_14 - ABA (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.0272 0.2144 5.7312  0.001 ***
  Residual            21   3.7638 0.7856                  
Total               22   4.7910 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_15 - SA (14D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.6696 0.29611 9.255  0.001 ***
  Residual            22   3.9687 0.70389                 
Total               23   5.6383 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_16 - Water control (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   2.0479 0.4629 18.961  0.001 ***
  Residual            22   2.3761 0.5371                  
Total               23   4.4240 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`06 - ABA (1D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   0.8596 0.16978 4.4988  0.001 ***
  Residual            22   4.2034 0.83022                  
Total               23   5.0630 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_09 - Methanol control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   0.7874 0.16458 4.3341  0.001 ***
  Residual            22   3.9970 0.83542                  
Total               23   4.7845 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_10 - ABA (7D)`
Df SumOfSqs     R2      F Pr(>F)  
Treatment_Timepoint  1   0.8968 0.1355 3.4483  0.029 *
  Residual            22   5.7218 0.8645                
Total               23   6.6186 1.0000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.6083 0.10897 2.6906  0.045 *
  Residual            22   4.9736 0.89103                
Total               23   5.5818 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.8696 0.33221 10.944  0.001 ***
  Residual            22   3.7583 0.66779                  
Total               23   5.6279 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.2818 0.21023 5.3239  0.001 ***
  Residual            20   4.8154 0.78977                  
Total               21   6.0972 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2199 0.21396 5.7164  0.002 **
  Residual            21   4.4817 0.78604                 
Total               22   5.7016 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.1413 0.19584 5.3576  0.002 **
  Residual            22   4.6866 0.80416                 
Total               23   5.8279 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.4502 0.44195 17.423  0.001 ***
  Residual            22   3.0939 0.55805                  
Total               23   5.5442 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`07 - SA (1D)_vs_08 - Water control (1D)`
Df SumOfSqs     R2      F Pr(>F)  
Treatment_Timepoint  1   0.7668 0.1348 3.4277  0.031 *
  Residual            22   4.9212 0.8652                
Total               23   5.6880 1.0000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_10 - ABA (7D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2154 0.18616 5.0324  0.002 **
  Residual            22   5.3133 0.81384                 
Total               23   6.5286 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_11 - SA (7D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.5339 0.2515 7.3921  0.001 ***
  Residual            22   4.5650 0.7485                  
Total               23   6.0989 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6233 0.32642 10.661  0.001 ***
  Residual            22   3.3498 0.67358                  
Total               23   4.9731 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   0.7080 0.13842 3.2132  0.005 **
  Residual            20   4.4069 0.86158                 
Total               21   5.1149 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1460 0.21957 5.9083  0.001 ***
  Residual            21   4.0731 0.78043                  
Total               22   5.2191 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.1003 0.20459 5.6586  0.001 ***
  Residual            22   4.2780 0.79541                  
Total               23   5.3784 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.2718 0.45828 18.612  0.001 ***
  Residual            22   2.6854 0.54172                  
Total               23   4.9572 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`09 - Methanol control (7D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   0.9310 0.17103 4.5389  0.001 ***
  Residual            22   4.5127 0.82897                  
Total               23   5.4437 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_11 - SA (7D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.7562 0.10732 2.6449  0.031 *
  Residual            22   6.2898 0.89268                
Total               23   7.0460 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_12 - Water control (7D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.6692 0.24751 7.2364  0.001 ***
  Residual            22   5.0745 0.75249                  
Total               23   6.7437 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.0833 0.15014 3.5334  0.006 **
  Residual            20   6.1317 0.84986                 
Total               21   7.2149 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)
Treatment_Timepoint  1   0.4851 0.07721 1.7571  0.123
Residual            21   5.7979 0.92279              
Total               22   6.2830 1.00000              

$`10 - ABA (7D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.3101 0.17915 4.8016  0.005 **
  Residual            22   6.0028 0.82085                 
Total               23   7.3129 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.1787 0.33066 10.868  0.001 ***
  Residual            22   4.4102 0.66934                  
Total               23   6.5889 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`10 - ABA (7D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)  
Treatment_Timepoint  1   0.9995 0.13811 3.5254  0.012 *
  Residual            22   6.2375 0.86189                
Total               23   7.2370 1.00000                
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_12 - Water control (7D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   2.0228 0.3186 10.286  0.001 ***
  Residual            22   4.3263 0.6814                  
Total               23   6.3491 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.4611 0.21347 5.4281  0.001 ***
  Residual            20   5.3834 0.78653                  
Total               21   6.8445 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_14 - ABA (14D)`
Df SumOfSqs     R2      F Pr(>F)   
Treatment_Timepoint  1   1.0838 0.1767 4.5072  0.003 **
  Residual            21   5.0496 0.8233                 
Total               22   6.1334 1.0000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2559 0.19291 5.2584  0.002 **
  Residual            22   5.2545 0.80709                 
Total               23   6.5105 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.5340 0.40898 15.224  0.001 ***
  Residual            22   3.6619 0.59102                  
Total               23   6.1959 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`11 - SA (7D)_vs_08 - Water control (1D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.4707 0.2113 5.8942  0.001 ***
  Residual            22   5.4892 0.7887                  
Total               23   6.9599 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_13 - Methanol control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.4819 0.26227 7.1104  0.001 ***
  Residual            20   4.1681 0.73773                  
Total               21   5.6500 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.7109 0.30853 9.3699  0.001 ***
  Residual            21   3.8344 0.69147                  
Total               22   5.5452 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.9147 0.32158 10.428  0.001 ***
  Residual            22   4.0393 0.67842                  
Total               23   5.9540 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_16 - Water control (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1  0.43142 0.1499 3.8792  0.001 ***
  Residual            22  2.44666 0.8501                  
Total               23  2.87808 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`12 - Water control (7D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   0.9853 0.18735 5.0719  0.001 ***
  Residual            22   4.2740 0.81265                  
Total               23   5.2593 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_14 - ABA (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.2788 0.20725 4.9673  0.001 ***
  Residual            19   4.8915 0.79275                  
Total               20   6.1703 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_15 - SA (14D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.1511 0.18425 4.5174  0.004 **
  Residual            20   5.0964 0.81575                 
Total               21   6.2475 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_16 - Water control (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   1.9666 0.3595 11.226  0.001 ***
  Residual            20   3.5038 0.6405                  
Total               21   5.4704 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`13 - Methanol control (14D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)   
Treatment_Timepoint  1   1.2029 0.18409 4.5126  0.003 **
  Residual            20   5.3311 0.81591                 
Total               21   6.5339 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_15 - SA (14D)`
Df SumOfSqs      R2     F Pr(>F)    
Treatment_Timepoint  1   1.5773 0.24879 6.955  0.001 ***
  Residual            21   4.7626 0.75121                 
Total               22   6.3400 1.00000                 
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_16 - Water control (14D)`
Df SumOfSqs     R2      F Pr(>F)    
Treatment_Timepoint  1   2.1685 0.4062 14.366  0.001 ***
  Residual            21   3.1700 0.5938                  
Total               22   5.3385 1.0000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`14 - ABA (14D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.4104 0.22011 5.9269  0.001 ***
  Residual            21   4.9973 0.77989                  
Total               22   6.4077 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_16 - Water control (14D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   2.4703 0.42262 16.103  0.001 ***
  Residual            22   3.3749 0.57738                  
Total               23   5.8453 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`15 - SA (14D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.2942 0.19922 5.4731  0.001 ***
  Residual            22   5.2022 0.80078                  
Total               23   6.4964 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$`16 - Water control (14D)_vs_08 - Water control (1D)`
Df SumOfSqs      R2      F Pr(>F)    
Treatment_Timepoint  1   1.5598 0.30174 9.5068  0.001 ***
  Residual            22   3.6096 0.69826                  
Total               23   5.1694 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

attr(,"class")
[1] "pwadstrata" "list" 
     
 
##linked resource webpage: https://stat.ethz.ch/pipermail/r-sig-ecology/2018-November/005830.html

#################################
######## DESeq analysis #########
################################

#####Bean Drought DeSeq analysis 
otu = read.csv("activefinalnozero_method2_changing_to_dnaabun_final.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
tax = read.csv("taxonomy-rep-seqs-or-99-edited.csv", sep=",", row.names=1)
tax = as.matrix(tax)
#metadata = read.csv("sample-metadata-decontam-edited.csv", sep=",", row.names=1)

##changing to edaphic metadata for constrained ordinations
metadata = read.csv("sample-metadata-decontam-edited1.csv", sep=",", row.names=1)

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 324 samples ]
sample_data() Sample Data:       [ 324 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]


metadata<- read.csv("deseq_bean_sa.csv") 
View(metadata)
metadata <- metadata %>% filter(Treatment=="Water-control"|Treatment=="Water-acclimated")
View(metadata)
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) 
phyloseq_merged_final


##bean aba (aba/water-acclim)
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 24 samples ]
sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]


##sw aba
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 24 samples ]
sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]
> 

##bean SA
  phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 24 samples ]
sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]
  

##SW SA
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 24 samples ]
sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]


##methanol control and water acclim SW

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 23 samples ]
sample_data() Sample Data:       [ 23 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]


##water control and water acclimated SW ##says lenght of dimnames not equal to array extent
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 24 samples ]
sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]


##bean methanol control and water acclimated
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 24 samples ]
sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]


##bean water control and water acclimated
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 24 samples ]
sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]

  
head(sample_data(phyloseq_merged_final)$Timepoint, 24)

if (!requireNamespace("BiocManager", quietly = TRUE))
  
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
packageVersion("DESeq2")

##convert phyloseq format to DESEq dataset with dispersions estimated using the experimental design formula
diagdds = phyloseq_to_deseq2(phyloseq_merged_final, ~ Timepoint)
diagdds$Timepoint
diagdds

diagdds$Timepoint <-factor(diagdds$Timepoint, levels = c("Water-acclimated (PreHormone)", "1-Day"))
diagdds$Timepoint #make sure that Control is the first level in the treatment factor, so that the
#default log2 fold changes are calculated as treatment over control and not the other way around.

#diagdds$group <- factor(paste0(diagdds$Timepoint))

#design(diagdds) <- ~ group

#diagdds <- DESeq(diagdds)
#resultsNames(diagdds)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

##The default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function.

#The following results function call creates a table of the results of the tests. Very fast. 
#The hard work was already stored with the rest of the DESeq2-related data in our latest version 
#of the diagdds object (see above). I then order by the adjusted p-value, removing the entries 
#with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display.


res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_merged_final)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

#sigtab_positives<- sigtab %>% filter(!log2FoldChange<0)
#View(sigtab_positives) 


#Let's look at the OTUs that were positively enriched in the planted bean soil samples compared to the unplanted ones. The following makes a nice ggplot2 summary of the results.

library("ggplot2")
#install.packages("Polychrome")
library(Polychrome)
P52<- createPalette(52,  c("#010101", "#ff0000"))
names(P52) <- NULL
library(viridis)


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "P52", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

View(sigtab) 


write.csv(sigtab, "sigtab_bean_aba.csv") ##aba 1 day and water acclim
write.csv(sigtab, "sigtab_sw_aba.csv")##aba 1 day and water acclim
write.csv(sigtab, "sigtab_bean_sa.csv")##sa 1 day and water acclim
write.csv(sigtab, "sigtab_sw_sa.csv")##sa 1 day and water acclim

write.csv(sigtab, "sigtab_sw_methanol.csv")##methanol 1 day and water acclim
write.csv(sigtab, "sigtab_bean_methanol.csv")##methanol 1 day and water acclim
write.csv(sigtab, "sigtab_bean_watercontrol.csv")#water 1 day and water acclim

## Figure 5 E - Bean drought positively enriched OTUs in planted treatments compared to unplanted

plot_bean_aba<- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  #scale_x_reordered() +
  geom_point(size=3) + scale_color_manual(values=as.vector(P52))+ 
  ggtitle("Bean water control")+
  theme(plot.title = element_text(hjust = 0.5, size=18))+
  theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_text(size=18))+
  theme(axis.text.y=element_text(size=13), axis.title.y=element_text(size=18)) + 
  theme(legend.title = element_text(size=18))+ theme(legend.text = element_text(size=13)) +
  guides(col=guide_legend(ncol=1)) +theme(strip.text.x = element_text(size = 30)) +geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  ylab("Log2FoldChange")+coord_flip() 

plot_bean_aba

ggsave(filename = "Deseq2-bean-watercontrol.tiff", plot = plot_bean_aba, ##CHANGED NAMES AS NEEDED. 
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)





#########################################################
################Alpha Diversity #########################
#########################################################

#############################
####alpha diversity lm model 

##richness estimates, inverse simpson, observed OTUs
###############################


# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(phyloseq_merged)
trials = 100
richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(phyloseq_merged)
evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(phyloseq_merged)
# It is always important to set a seed when you subsample so your result is replicable, note of caution that this will yield different results of alpha diversity depending on the version of R being used. This is because set.seed function can vary across R versions. Here I am reporting results from R v.3.4.0 
set.seed(3)
for (i in 1:100) {
  # Subsample
  #r <- rarefy_even_depth(phyloseq_merged, sample.size = 1936, verbose = FALSE, replace = TRUE) #no rarefaction a second time.
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(phyloseq_merged, measures = "Observed"))) ##changed r to phyloseq_merged, no second rarefaction
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(phyloseq_merged, measures = "InvSimpson"))) ##changed r to phyloseq_merged, no second rarefaction
  evenness[ ,i] <- even
}
# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)
# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd) ##sd is always 0, why?
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)
alpha <- rbind(rich_stats, even_stats)
#s <- data.frame(sample_data(phyloseq_rarefy))
s <- read.csv("sample-metadata-decontam.csv")
alphadiv <- merge(alpha, s, by = "SampleID") 
write.csv(alphadiv,file = "alphadiv_withsingleton_nosecondrarefaction.csv")
#alphadiv <- order_dates(alphadiv)
alphadiv_sw<-alphadiv[alphadiv$Plant=="Switchgrass",]
dim(alphadiv_sw)
alphadiv_bean<-alphadiv[alphadiv$Plant=="Bean",]
dim(alphadiv_bean)
View(alphadiv_sw)
View(alphadiv_bean)
write.csv(alphadiv_bean, file="alphadiv_bean_withsingleton_nosecondrarefaction.csv")
write.csv(alphadiv_sw, file="alphadiv_switchgrass_withsingleton_nosecondrarefaction.csv")
#levels(alphadiv_WA$Year_Month)<-c("Spring2015","Fall2015",  "Spring2016", "Fall2016")
#levels(alphadiv_TN$Year_Month)<-c("Spring2015","Fall2015",  "Spring2016", "Fall2016")
group_sw<-alphadiv_sw%>%
  group_by(Treatment_Timepoint,measure) %>%
  summarise(mean_update = mean(mean))
write.csv(group_sw, file="alphadivmean_sw_withsingleton_nosecondrarefaction.csv")
group_bean<-alphadiv_bean%>%
  group_by(Treatment_Timepoint,measure) %>%
  summarise(mean_update = mean(mean))
write.csv(group_bean, file="alphadivmean_bean_withsingleton_nosecondrarefaction.csv")


library(ggplot2)
library(dplyr)
library(vegan)
library(reshape2)
library(scales)
library(grid)
library(readxl)


##combine bean and switchgrass richness and evenness together in panel


alphadiv_combined<- read.csv("alphadiv_edited_forplot.csv") ##rename harvest to days in csv file
alphadiv_richness<- alphadiv_combined %>% filter(measure=="Richness")
View(alphadiv_richness)
#alphadiv_richness<- alphadiv_combined %>% filter(measure=="Inverse Simpson")
#alphadiv_richness$group <- recode(alphadiv_richness$Timepoint, 
                          #    '1D' = '1-Day',
                            #  '7D' = '7-Day',
                           #   '14D' = '14-Day')
                             
alphadiv_richness$grid1 = factor(alphadiv_richness$grid, levels=c('Water control','Methanol control','Abscisic acid','Salicylic acid'))                             
#alphadiv_richness<- as.data.frame(alphadiv_richness)
#alphadiv_combined_plot<-ggplot(alphadiv_richness, aes(y=mean, x= Timepoint, fill=Plant,linetype=Plant, color=Plant)) +theme_classic()+
  alphadiv_combined_plot<-ggplot(alphadiv_richness, aes(y=mean, x= Timepoint, fill=Plant, color=Plant)) +theme_classic()+
    scale_fill_manual(values=c("#CC79A7", "#0072B2"))+
    scale_color_manual(values=c("#CC79A7", "#0072B2"))+
#scale_fill_manual(values=c("#39568CFF", "#55C667FF", "#440154FF", "#4DE725FF"))+
  #scale_color_manual(values=c("#39568CFF", "#55C667FF", "#440154FF", "#4DE725FF"))+ 
  #geom_point(aes(),size=2) + 
  #geom_smooth(method=lm, se=FALSE)+
  scale_x_discrete(limits= c("Field-Soil", "Predry", "Postdry", "Water-acclimated", "1-Day", "7-Day", "14-Day"))+
  facet_grid(~grid1) +  
  geom_boxplot()+
  #geom_smooth(method=lm, se=TRUE, aes(group=Plant), inherit.aes=TRUE) + 
  #xlab("")+
  labs(y = "Richness(Number of observed OTUs)", color='Plant')+ylab("Richness(Number of observed OTUs)")+
  #ggtitle("Alpha Diversity Estimates")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=12)) +
  theme(axis.title.x = element_blank())+
  #theme(axis.title.y = element_blank())+
  theme(axis.text=element_text(size=12))+
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=12)) +
  guides(linetype = guide_legend(override.aes= list(color = "black")))

alphadiv_combined_plot

ggsave(filename="alphadiv_combined_richness.TIFF", plot=alphadiv_combined_plot, width=10, height=8, units="in", dpi=300)
ggsave(filename="alphadiv_combined_richness_09.25.24.TIFF", plot=alphadiv_combined_plot, width=7, height=5, units="in", dpi=300)
ggsave(filename="alphadiv_combined_diversity_09.25.24.TIFF", plot=alphadiv_combined_plot, width=7, height=5, units="in", dpi=300)


####Alpha diversity statistics

library(ggplot2)
library(grid)
library(lattice)
library(multcompView)
library(tidyverse)
library(dplyr)
alphadiv <- read.csv("alphadiv_edited_forplot.csv") ## combined richness and evenness estimates of all samples
alphadiv<- alphadiv %>% filter(Treatment!="Field-Soil", Treatment!="Predry", Treatment!="Postdry")
alphadiv<- alphadiv %>% filter(Treatment!="Field-Soil", Treatment!="Predry", Treatment!="Postdry", Treatment!="Water-acclimated")
#"there are aliased coefficients in the model", possibly unbalanced data
pd<-position_dodge(0.7)

##with water acclimated samples
alphadiv_inv =  alphadiv %>% filter(measure=="Inverse Simpson")
alphadiv_rich =  alphadiv %>% filter(measure=="Richness")

##removing water acclimated samples
alphadiv_inv_treat =  alphadiv %>% filter(measure=="Inverse Simpson")
alphadiv_rich_treat =  alphadiv %>% filter(measure=="Richness")


##Conduct the following code for each of the datasets selected above
View(alphadiv)
options(contrasts = c("contr.sum", "contr.poly")) ### must set this before running the model for type III tests
model = lm(mean~ Plant+Treatment 
           + Plant:Treatment,
           data=alphadiv_inv_treat) ##change the inputs based on the datasets above

shapiro.test(resid(model)) #W statistic > 0.9 means normality okay #check normality to flag outliers 
ee= as.matrix(resid(model))
qqnorm (ee) #look for outliers (far away from curve) and test for normality 
ee

library(car)
#options(contrasts = c("contr.sum", "contr.poly"))
Anova(model, type=c("III"))



##removing field soil, predry and postdry inverse simpson
Shapiro-Wilk normality test

data:  resid(model)
W = 0.71256, p-value < 2.2e-16

Anova Table (Type III tests)

Response: mean
                 Sum Sq  Df  F value  Pr(>F)    
(Intercept)      633982   1 161.5143 < 2e-16 ***
Plant              8928   1   2.2745 0.13237    
Timepoint         35840   3   3.0436 0.02886 *  
Plant:Timepoint   10478   3   0.8898 0.44647    
Residuals       1456263 371                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##removing field soil, predry and postdry richness

Shapiro-Wilk normality test

data:  resid(model)
W = 0.9772, p-value = 1.092e-05

Response: mean
                 Sum Sq  Df   F value    Pr(>F)    
(Intercept)     6314730   1 4995.3393 < 2.2e-16 ***
Plant             46434   1   36.7322 3.333e-09 ***
Timepoint         31165   3    8.2179 2.620e-05 ***
Plant:Timepoint    2251   3    0.5936    0.6195    
Residuals        468990 371                        
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##removing field soil, predry, post dry and water acclim richness

Shapiro-Wilk normality test

data:  resid(model)
W = 0.98617, p-value = 0.008002

Anova Table (Type III tests)

Response: mean
                 Sum Sq  Df  F value    Pr(>F)    
(Intercept)     4419896   1 5249.814 < 2.2e-16 ***
Plant             29094   1   34.557 1.196e-08 ***
Treatment        146487   3   57.998 < 2.2e-16 ***
Plant:Treatment   32290   3   12.784 7.606e-08 ***
Residuals        231527 275                       
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##removing field soil, predry, post dry and water acclim inv simpson
Shapiro-Wilk normality test

data:  resid(model)
W = 0.69641, p-value < 2.2e-16

Anova Table (Type III tests)

Response: mean
                Sum Sq  Df  F value    Pr(>F)    
(Intercept)     366859   1 111.8846 < 2.2e-16 ***
Plant             1449   1   0.4418  0.506801    
Treatment        47335   3   4.8121  0.002778 ** 
Plant:Treatment  24613   3   2.5021  0.059673 .  
Residuals       901700 275                       
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1





###################################################################################################
##USING THE INDICATOR SPECIES HEATMAP CODE FROM NATCOMMS paper below 
###################################################################################################
Active_DNA<- read.csv("activefinalnozero_method2_changing_to_dnaabun_final_activity.csv", row.names=1)
View(Active_DNA)
#Active_new <- cbind(OTU = rownames(Active_DNA), Active_DNA)

#ind_sw_abund<- left_join(crop_ind_sw, Active_new, by="OTU")
#View(ind_sw_abund)
#rownames(ind_sw_abund) <- ind_sw_abund[,1]
#ind_sw_mat = ind_sw_abund[,16:ncol(ind_sw_abund)] ##chnage the column number as appropriate
#View(ind_sw_mat)
#ind_sw_50 <- ind_sw_mat[order(rowSums(ind_sw_mat), decreasing = TRUE),]
#View(ind_sw_50)

###moved this to later so that decreasing order is set after subsetting to switchgrass/bean samples
#ind_bean_50 <- ind_bean_mat[order(rowSums(ind_bean_mat), decreasing = TRUE),]
#ind_bean_50 <- ind_bean_50[c(1:50),] ##select top 50 most abundant Classes
#View(ind_bean_50)

#ind_bean_20 <- ind_bean_mat[order(rowSums(ind_bean_mat), decreasing = TRUE),]
#ind_bean_20 <- ind_bean_20[c(1:20),] ##select top 50 most abundant Classes
#View(ind_bean_20)

####start again from here
##subset to correct samples before computing max standardization
Active_DNA_long<- Active_DNA %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="Abundance", -ASV)#gather columns into key value pairs
View(Active_DNA_long)

metadata<-read.csv("sample-metadata-decontam-edited2.csv")

newdf = left_join(Active_DNA_long, metadata, by=c("SampleID"))
View(newdf)
#new_df_sw<- newdf %>% filter(crop=="bean", drought=="well-watered", planted=="planted") ##change to new_df_bean or keep using the same object name, but change the parameters as needed for the crop and within crop factors, should match the condition in line 3679 

##filter samples switchgrass
new_df_sw <- newdf %>% filter (Plant=="Switchgrass", Treatment != "Methanol-control", Treatment != "Water-control", Treatment != "Salicylic-acid") #ABA
new_df_sw<- newdf %>% filter (Plant=="Switchgrass", Treatment != "Methanol-control", Treatment != "Water-control", Treatment != "Abscisic-acid") #SA
#new_df_sw<- newdf %>% filter (Plant=="Switchgrass", Treatment != "Methanol-control", Treatment != "Salicylic-acid", Treatment != "Abscisic-acid") #WATER
new_df_sw<-  newdf %>% filter (Plant=="Switchgrass", Treatment != "Water-control", Treatment != "Salicylic-acid", Treatment != "Abscisic-acid") #METHANOL
View()

##filter samples bean
new_df_sw <- newdf %>% filter (Plant=="Bean", Treatment != "Methanol-control", Treatment != "Water-control", Treatment != "Salicylic-acid") #ABA
new_df_sw<- newdf %>% filter (Plant=="Bean", Treatment != "Methanol-control", Treatment != "Water-control", Treatment != "Abscisic-acid") #SA
#new_df_sw<- newdf %>% filter (Plant=="Bean", Treatment != "Methanol-control", Treatment != "Salicylic-acid", Treatment != "Abscisic-acid") #WATER
new_df_sw<-  newdf %>% filter (Plant=="Bean", Treatment != "Water-control", Treatment != "Salicylic-acid", Treatment != "Abscisic-acid") #METHANOL
View()


##ABA
new_df_sw$Treatment_Timepoint1 <- recode(new_df_sw$Treatment_Timepoint, 
                                                   '01 - Field Soil' = 'Field Soil',
                                                   '02 - Predry' = 'Predry',
                                                   '03 - Postdry' = 'Postdry',
                                                   '04 - Water acclimated' ='Water acclimated',
                                                   '06 - ABA (1D)' = 'ABA (1-Day)',
                                                   '10 - ABA (7D)' = 'ABA (7-Day)',
                                                   '14 - ABA (14D)' = 'ABA (14-Day)')

##SA
new_df_sw$Treatment_Timepoint1 <- recode(new_df_sw$Treatment_Timepoint, 
                                                   '01 - Field Soil' = 'Field Soil',
                                                   '02 - Predry' = 'Predry',
                                                   '03 - Postdry' = 'Postdry',
                                                   '04 - Water acclimated' ='Water acclimated',
                                                   '07 - SA (1D)' = 'SA (1-Day)',
                                                   '11 - SA (7D)' = 'SA (7-Day)',
                                                   '15 - SA (14D)' = 'SA (14-Day)')

#WATER CONTROL
new_df_sw$Treatment_Timepoint1 <- recode(new_df_sw$Treatment_Timepoint, 
                                                   '01 - Field Soil' = 'Field Soil',
                                                   '02 - Predry' = 'Predry',
                                                   '03 - Postdry' = 'Postdry',
                                                   '04 - Water acclimated' ='Water acclimated',
                                                   '08 - Water control (1D)' = 'Water control (1-Day)',
                                                   '12 - Water control (7D)' = 'Water control (7-Day)',
                                                   '16 - Water control (14D)' = 'Water control (14-Day)')

#METHANOL CONTROL
new_df_sw$Treatment_Timepoint1 <- recode(new_df_sw$Treatment_Timepoint, 
                                                   '01 - Field Soil' = 'Field Soil',
                                                   '02 - Predry' = 'Predry',
                                                   '03 - Postdry' = 'Postdry',
                                                   '04 - Water acclimated' ='Water acclimated',
                                                   '05 - Methanol control (1D)' = 'Methanol control (1-Day)',
                                                   '09 - Methanol control (7D)' = 'Methanol control (7-Day)',
                                                   '13 - Methanol control (14D)' = 'Methanol control (14-Day)')









#######skip this part for the new heatmaps
##pivot wider to plot averages, do the rowsums sorting after filtering to switchgrass? if i do it before then bean abundances may obscure the sw rowsum abudances (this is unlikely given that the OTUs picked are indicators for a specific crop and likely will have more abundances in that crop.)
new_df_sw_wide <- new_df_sw %>% select(ASV, Abundance, SampleID) %>%pivot_wider(names_from="SampleID", values_from="Abundance")  
View(new_df_sw_wide) 
rnames <-new_df_sw_wide$ASV

new_df_sw_wide<- new_df_sw_wide %>% select (!ASV)
rownames(new_df_sw_wide) <- rnames
View(new_df_sw_wide)

new_df_sw_wide<- as.matrix(new_df_sw_wide)
new_df_sw_wide <- new_df_sw_wide[order(rowSums(new_df_sw_wide), decreasing = TRUE),]
new_df_sw_wide <- new_df_sw_wide[c(1:50),] ##select top 50 most abundant Classes
View(new_df_sw_wide)
###########

##############
##to include the deseq2 otus in the heatmaps im copying the same code as the bubble plot for this section
###########

new_df_sw_wide<- new_df_sw %>% select(ASV, Abundance, SampleID) %>%pivot_wider(names_from="SampleID", values_from="Abundance")  
new_df_sw_wide_nodeseq <- new_df_sw_wide %>% filter(ASV!="c335e93a51cbc647ce7c8e354147e462" & ASV!="a78c7c528e0197b7e64605b4a076700a" & ASV!="2985946aa8517aed16c91ddc3afbbd59" & ASV!="76b7130f1162a090d1f340b4a29e632b")


rnames <-new_df_sw_wide_nodeseq$ASV
View(new_df_sw_wide_nodeseq) 
new_df_sw_wide_nodeseq<- new_df_sw_wide_nodeseq %>% select (!ASV)
rownames(new_df_sw_wide_nodeseq) <- rnames
View(new_df_sw_wide_nodeseq)

new_df_sw_wide_nodeseq<- as.matrix(new_df_sw_wide_nodeseq)
new_df_sw_wide_nodeseq <- new_df_sw_wide_nodeseq[order(rowSums(new_df_sw_wide_nodeseq), decreasing = TRUE),]
new_df_sw_wide_nodeseq <- new_df_sw_wide_nodeseq[c(1:50),] ##select top 50 most abundant Classes
View(new_df_sw_wide_nodeseq)

#add the deseq2 otus here 

#new_df_sw_deseq<- new_df_sw %>% filter(ASV==c("c335e93a51cbc647ce7c8e354147e462","a78c7c528e0197b7e64605b4a076700a","2985946aa8517aed16c91ddc3afbbd59","76b7130f1162a090d1f340b4a29e632b","a78c7c528e0197b7e64605b4a076700a"))
#new_df_sw_deseq_wide <- new_df_sw_deseq %>% select(ASV, Abundance, SampleID) %>%pivot_wider(names_from="SampleID", values_from="Abundance")  
 new_df_sw_wide <- new_df_sw %>% select(ASV, Abundance, SampleID) %>%pivot_wider(names_from="SampleID", values_from="Abundance")  
View(new_df_sw_wide) 

new_df_sw_wide_deseq<- new_df_sw_wide %>% filter(ASV=="c335e93a51cbc647ce7c8e354147e462" | ASV=="a78c7c528e0197b7e64605b4a076700a" | ASV=="2985946aa8517aed16c91ddc3afbbd59" | ASV=="76b7130f1162a090d1f340b4a29e632b")
rnames <-new_df_sw_wide_deseq$ASV
new_df_sw_wide_deseq<- new_df_sw_wide_deseq %>% select (!ASV)
rownames(new_df_sw_wide_deseq) <- rnames
View(new_df_sw_wide_deseq)

##rowbind the two above

new_df_sw_wide <- rbind(new_df_sw_wide_nodeseq, new_df_sw_wide_deseq)



##merging with dna cdna data to get NAs inserted in dataframe
##using rarefied DNA and cDNA copies with 15k reads.
dna<- read.csv("DNAcopy_final.csv", header=TRUE, row.names = 1)
cdna <- read.csv("RNAcopy_final.csv", header=TRUE, row.names = 1)

View(dna)
View(cdna)

dna_50taxa <- dna[rownames(dna)%in%rownames(new_df_sw_wide),]
View(dna_50taxa)

cdna_50taxa <- cdna[rownames(cdna)%in%rownames(new_df_sw_wide),]
View(cdna_50taxa)

dna1<- dna_50taxa %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="dna", -ASV)#gather columns into key value pairs
View(dna1)
cdna1<- cdna_50taxa %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="cdna", -ASV)
View(cdna1)
merge.df = full_join(dna1, cdna1, by=c("ASV", "SampleID"))
View(merge.df)

merge.df[is.na(merge.df)] <-0
View(merge.df)
metadata<-read.csv("sample-metadata-decontam-edited2.csv")

newdf = full_join(merge.df, metadata, by=c("SampleID"))
View(newdf)


finaldf = newdf %>%  mutate(ratioMethod1= ifelse(dna == 0 & cdna > 0, 100, cdna/dna)) %>%  
  mutate(dna2= ifelse(dna == 0, 1,dna))  %>%  mutate(ratioMethod2= cdna/dna2) %>% mutate(ratio_nophantoms = cdna/dna)
View(finaldf)
finaldf[is.na(finaldf)] <-0
finaldf$ratioMethod1[!is.finite(finaldf$ratioMethod1)] <- 0 ##using the code above for removing nas also works here as it changed the 0 denominator ratios from NA to 0, no infinite values only NaN
finaldf$ratioMethod2[!is.finite(finaldf$ratioMethod2)] <- 0 ##ratios computed by methods 1 and 2 will always be finite, because of the zero denominator being accounted for phantom taxa?
finaldf$ratio_nophantoms[!is.finite(finaldf$ratio_nophantoms)] <- 0
View(finaldf) 


##code taxa as inactive (0) or not detected in dna (NA)


finaldf_coded <- finaldf %>% #mutate(dna_code=ifelse(dna2==0, "not_detected", dna2)) 
  mutate(activity= ifelse(ratioMethod2 < 1, "inactive", "active"))
View(finaldf_coded) 

##Abundant_active dataset abundance combined with finaldf dataset 
Abundant_Active_new <- cbind(ASV = rownames(new_df_sw_wide), new_df_sw_wide)
View(Abundant_Active_new)
Abundant_Active_new<- as.data.frame(Abundant_Active_new)
Abundant_Active_long <-pivot_longer(Abundant_Active_new, !ASV, names_to="SampleID", values_to="AbundanceDNA_Active")
View(Abundant_Active_long)
merge_df_new<- merge(Abundant_Active_long, finaldf_coded, by=c("ASV", "SampleID"))
View(merge_df_new)
merge_df_new_NA<- merge_df_new %>% mutate(AbundanceDNA_Active_NA= ifelse(dna==0 & cdna==0, NA, AbundanceDNA_Active))
View(merge_df_new_NA)

##pivot wider
merge_df_new_NA<- merge_df_new_NA %>% mutate(Treat_Time_ID=paste(Treatment_Timepoint, SampleID))
merge_df_new_NA_wide <- merge_df_new_NA %>% select(ASV, AbundanceDNA_Active_NA, Treat_Time_ID) %>%pivot_wider(names_from="Treat_Time_ID", values_from="AbundanceDNA_Active_NA")  
View(merge_df_new_NA_wide) 
data<- as.data.frame(merge_df_new_NA_wide)

##create matrix


#option 2, both options work
rnames <-data[,1] 
rownames(data) <- rnames  # assign row names
View(data)
mat_data<- data %>% select(!ASV)
View(mat_data)


#REMOVE OTUS that cant be plotted in heatmap, note which ones are removed, if there are only numbers nothing is removed. if there are NAs some may be removed.

#mat_data<- data.matrix(mat_data)


sum(is.na(as.matrix(dist(mat_data))))

giveNAs = which(is.na(as.matrix(dist(mat_data))),arr.ind=TRUE)
head(giveNAs)
#mat_data[c(1,17),]

tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
  sum(is.na(as.matrix(dist(mat_data[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))] ##https://stackoverflow.com/questions/61469201/pheatmap-won-t-cluster-rows-na-nan-inf-in-foreign-function-call-arg-10 
rmv ## the intention with this code is not to remove rows but to see for which ASVs are we not able to calculate a distance for the heatmap

##use this line only if something needs to be removed
mat_data = mat_data[-as.numeric(rmv),] 
View(mat_data)


####merge taxonomy with OTU hash id

taxonomy<- read.csv("taxonomy-rep-seqs-or-99-edited.csv", row.names=1)
#View(taxonomy)
taxonomy_family<- taxonomy %>% select(c("Family", "Class"))
#View(taxonomy_family)

mat_data_merge <- merge(mat_data, taxonomy_family, by="row.names")
View(mat_data_merge)
mat_data_merge$OTUclass <- paste(mat_data_merge$Class, mat_data_merge$Row.names)
View(mat_data_merge)

rnames <-mat_data_merge$OTUclass

rownames(mat_data_merge) <- rnames  # assign row names
View(mat_data_merge)
mat_data_merge_otuclass<- mat_data_merge %>% select(!Row.names)%>% select (!Family) %>% select(!Class) %>% select(!OTUclass)
View(mat_data_merge_otuclass)

##renaming colnames to make it concise
# Number of characters to remove
num_characters_to_remove <- 5

# Rename columns using rename_with() and str_sub()
mat_data_merge_otuclass <- mat_data_merge_otuclass %>% 
  rename_with(~ str_sub(.x, num_characters_to_remove + 1))

# Strings to be replaced
old_string <- " X_"
new_string <- "_"

# Rename columns using rename_with() and str_replace()
mat_data_merge_otuclass<- mat_data_merge_otuclass %>% 
  rename_with(~ str_replace(.x, old_string, new_string))


##switchgrass
write.csv(mat_data_merge_otuclass, "mat_data_merge_otuclass_sw_aba_01.28.25.csv")
write.csv(mat_data_merge_otuclass, "mat_data_merge_otuclass_sw_sa_01.30.25.csv")
write.csv(mat_data_merge_otuclass, "mat_data_merge_otuclass_sw_met_con_01.30.25.csv")

#bean
write.csv(mat_data_merge_otuclass, "mat_data_merge_otuclass_bean_aba_01.30.25.csv")
write.csv(mat_data_merge_otuclass, "mat_data_merge_otuclass_bean_sa_01.30.25.csv")
write.csv(mat_data_merge_otuclass, "mat_data_merge_otuclass_bean_met_con_01.30.25.csv")


if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}
install_github("raivokolde/pheatmap")
library(pheatmap)

library(vegan)

mat_data_merge_otuclass<-read.csv("mat_data_merge_otuclass_sw_aba_01.28.25.csv", row.names=1, check.names=FALSE)

##CHANGING THE SAMPLE IDS in csv file before importing back in
#mat_data_merge_otuclass<- read.csv("mat_data_merge_otuclass_sw_aba.csv", row.names=1)

mat_data_merge_otuclass<- as.data.frame(mat_data_merge_otuclass)
#mat_data_merge_otuclass<- as_tibble(mat_data_merge_otuclass)

class(mat_data_merge_otuclass)
##replace NAs with zero
matrix<- as.matrix(mat_data_merge_otuclass)
matrix
View(matrix)
matrix[is.na(matrix)] <-0
View(matrix)

#matrix_new<-data.matrix(matrix)
#View(matrix_new)
##compute max standardization

matrix_max_standardize_bean <-decostand(matrix, method = "max", MARGIN = 1) ##from Jackson's phil trans paper
View(matrix_max_standardize_bean)

##put back NAs where they were before
matrix1<- as.matrix(mat_data_merge_otuclass)
matrix1
View(matrix1)
matrix1[is.na(matrix1)] <- -100
View(matrix1)
mask <- matrix1 == -100
mask

matrix_max_standardize_bean[mask] <- matrix1[mask]
View(matrix_max_standardize_bean)

matrix_max_standardize_bean[matrix_max_standardize_bean == -100]<- NA
View(matrix_max_standardize_bean) ##gives the standardized values with NA


###

library(viridis)

##pivot longer, merge with metadata
##convert to dataframe if doing this part

#sw_long<- matrix_max_standardize_sw %>%
#tibble::rownames_to_column(var="ASV") %>%
# tidyr::gather(key="SampleID", value="RA", -ASV)#gather columns into key value pairs
#View(sw_long)

#metadata<-read.csv("sample_metadata_dna_cdna_edited.csv")

#new_df_sw = left_join(sw_long, metadata, by=c("SampleID"))
#View(new_df_sw)

#compute averages by harvest and summarise new, might not work if there are NAs, might need the zeros . In that case dont put back the NAs after max standardization
#newdf_mean_sw<- new_df_sw %>% group_by(summarise, ASV) %>% summarise(mean=mean(RA))
#View(newdf_mean_sw)

##pivot wider to plot averages
#sw_ind_wide <- newdf_mean_sw %>% select(ASV, mean, summarise) %>%pivot_wider(names_from="summarise", values_from="mean")  
#View(sw_ind_wide) 

#class(sw_ind_wide)

#sw_ind_mat<- as.matrix(sw_ind_wide)
#class(sw_ind_mat)
#View(sw_ind_mat)

#rnames <-sw_ind_wide$ASV

#sw_ind_wide<- sw_ind_wide %>% select (!ASV)
#rownames(sw_ind_wide) <- rnames
#View(sw_ind_wide)


#plot

##bean planted drought, using default clustering "complete" linkage

#create the breaks
bk2 = unique(c(0.00000000, seq(0.00000001, 1.00000000, length=20)))

col1 = colorRampPalette(c("gray"))(1)

col2 =viridis(20)

colors2 <- c(col1, col2)

## SW ABA
# Define labels to color red
red_labels <- c("c__Myxococcia c335e93a51cbc647ce7c8e354147e462", "c__Actinobacteria a78c7c528e0197b7e64605b4a076700a", "c__Actinobacteria 2985946aa8517aed16c91ddc3afbbd59") 

# Create a vector of colors based on row names
colors <- ifelse(rownames(matrix_max_standardize_bean) %in% red_labels, "red", "black")


##SW_SA
red_labels<-c("c__Actinobacteria a78c7c528e0197b7e64605b4a076700a")
colors <- ifelse(rownames(matrix_max_standardize_bean) %in% red_labels, "red", "black")

##SW methanol
red_labels<-c("c__Myxococcia c335e93a51cbc647ce7c8e354147e462", "c__Actinobacteria a78c7c528e0197b7e64605b4a076700a", "c__Actinobacteria 2985946aa8517aed16c91ddc3afbbd59")
colors <- ifelse(rownames(matrix_max_standardize_bean) %in% red_labels, "red", "black")


##bean ABA
red_labels<-c()
colors <- ifelse(rownames(matrix_max_standardize_bean) %in% red_labels, "red", "black")


##bean SA  
red_labels<-c("c__Actinobacteria a78c7c528e0197b7e64605b4a076700a", "c__Actinobacteria 76b7130f1162a090d1f340b4a29e632b")
colors <- ifelse(rownames(matrix_max_standardize_bean) %in% red_labels, "red", "black")

##bean methanol control  
red_labels<-c("c__Actinobacteria a78c7c528e0197b7e64605b4a076700a", "c__Actinobacteria 76b7130f1162a090d1f340b4a29e632b")
colors <- ifelse(rownames(matrix_max_standardize_bean) %in% red_labels, "red", "black")

#matrix_max_standardize_bean<- as.data.frame(matrix_max_standardize_bean)
#matrix_max_standardize_bean_color<- matrix_max_standardize_bean %>% mutate(color=ifelse(row.names(matrix_max_standardize_bean)%in% deseq, "red", "black"))


#label_function <- function(x) {
 # paste0("<span style='color:", matrix_max_standardize_bean_color$color[match(x, row.names(matrix_max_standardize_bean_color))], ";'>", x, "</span>")
#}

library(pheatmap)


names(colors) <- rownames(matrix_max_standardize_bean) 

# Create annotation data frame

annotation_df <- data.frame(Color = colors)

rownames(annotation_df) <- rownames(matrix_max_standardize_bean)




plot_heatmap_sw <- pheatmap::pheatmap(matrix_max_standardize_bean, ##if plotting sample by sample and not the means then use matrix_max_standardize_sw
                                                        cluster_cols = F, border_color = "white", na_col = "white", 
                                                        main="Switchgrass ABA",
                                                         #labels_row=TRUE,
                                                        color = colors2,
                                                        breaks = bk2,
                                                        annotation_row = annotation_df, 
                                                        annotation_colors = list(Color = c(black = "black", red = "red")),clustering_distance_rows="euclidean",clustering_method = "complete"
                                      
                                                        
)



##new with deseq otus
ggsave(filename = "sw_aba_top50_heatmap_all_01.28.25.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)
ggsave(filename = "sw_sa_top50_heatmap_all_01.30.25.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)
ggsave(filename = "sw_met_con_top50_heatmap_all_01.30.25.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)

ggsave(filename = "bean_aba_top50_heatmap_all_01.30.25.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)
ggsave(filename = "bean_sa_top50_heatmap_all_01.30.25.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)
ggsave(filename = "bean_met_con_top50_heatmap_all_01.30.25.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)


ggsave(filename = "sw_aba_top50_heatmap_all_01.31.25.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)


##switchgrass

ggsave(filename = "sw_aba_top50_09.10.24.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)
ggsave(filename = "sw_sa_top50_09.20.24.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)
ggsave(filename = "sw_wat_con_top50_09.20.24.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)
ggsave(filename = "sw_met_con_top50_09.20.24.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)


##bean
ggsave(filename = "bean_aba_top50_09.24.24.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)

ggsave(filename = "bean_sa_top50_09.24.24.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)


ggsave(filename = "bean_wat_con_top50_09.20.24.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)


ggsave(filename = "bean_met_con_top50_09.24.24.tiff", plot = plot_heatmap_sw,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)


############
###bubble plots averaging abundance and including detection limit across replicates
#############

Active_DNA<- read.csv("activefinalnozero_method2_changing_to_dnaabun_final_activity.csv", row.names=1)
View(Active_DNA)


####start again from here
##subset to correct samples before computing max standardization
Active_DNA_long<- Active_DNA %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="Abundance", -ASV)#gather columns into key value pairs
View(Active_DNA_long)

metadata<-read.csv("sample-metadata-decontam-edited2.csv")

newdf = left_join(Active_DNA_long, metadata, by=c("SampleID"))
View(newdf)
#new_df_sw<- newdf %>% filter(crop=="bean", drought=="well-watered", planted=="planted") ##change to new_df_bean or keep using the same object name, but change the parameters as needed for the crop and within crop factors, should match the condition in line 3679 

##filter samples switchgrass
new_df_sw <- newdf %>% filter (Plant=="Switchgrass", Treatment != "Methanol-control", Treatment != "Water-control", Treatment != "Salicylic-acid") #ABA
new_df_sw<- newdf %>% filter (Plant=="Switchgrass", Treatment != "Methanol-control", Treatment != "Water-control", Treatment != "Abscisic-acid") #SA
#new_df_sw<- newdf %>% filter (Plant=="Switchgrass", Treatment != "Methanol-control", Treatment != "Salicylic-acid", Treatment != "Abscisic-acid") #WATER
new_df_sw<-  newdf %>% filter (Plant=="Switchgrass", Treatment != "Water-control", Treatment != "Salicylic-acid", Treatment != "Abscisic-acid") #METHANOL
View()

##filter samples bean
new_df_sw <- newdf %>% filter (Plant=="Bean", Treatment != "Methanol-control", Treatment != "Water-control", Treatment != "Salicylic-acid") #ABA
new_df_sw<- newdf %>% filter (Plant=="Bean", Treatment != "Methanol-control", Treatment != "Water-control", Treatment != "Abscisic-acid") #SA
#new_df_sw<- newdf %>% filter (Plant=="Bean", Treatment != "Methanol-control", Treatment != "Salicylic-acid", Treatment != "Abscisic-acid") #WATER
new_df_sw<-  newdf %>% filter (Plant=="Bean", Treatment != "Water-control", Treatment != "Salicylic-acid", Treatment != "Abscisic-acid") #METHANOL
View()


##ABA
new_df_sw$Treatment_Timepoint1 <- recode(new_df_sw$Treatment_Timepoint, 
                                                   '01 - Field Soil' = 'Field Soil',
                                                   '02 - Predry' = 'Predry',
                                                   '03 - Postdry' = 'Postdry',
                                                   '04 - Water acclimated' ='Water acclimated',
                                                   '06 - ABA (1D)' = 'ABA (1-Day)',
                                                   '10 - ABA (7D)' = 'ABA (7-Day)',
                                                   '14 - ABA (14D)' = 'ABA (14-Day)')

##SA
new_df_sw$Treatment_Timepoint1 <- recode(new_df_sw$Treatment_Timepoint, 
                                                   '01 - Field Soil' = 'Field Soil',
                                                   '02 - Predry' = 'Predry',
                                                   '03 - Postdry' = 'Postdry',
                                                   '04 - Water acclimated' ='Water acclimated',
                                                   '07 - SA (1D)' = 'SA (1-Day)',
                                                   '11 - SA (7D)' = 'SA (7-Day)',
                                                   '15 - SA (14D)' = 'SA (14-Day)')

#WATER CONTROL
new_df_sw$Treatment_Timepoint1 <- recode(new_df_sw$Treatment_Timepoint, 
                                                   '01 - Field Soil' = 'Field Soil',
                                                   '02 - Predry' = 'Predry',
                                                   '03 - Postdry' = 'Postdry',
                                                   '04 - Water acclimated' ='Water acclimated',
                                                   '08 - Water control (1D)' = 'Water control (1-Day)',
                                                   '12 - Water control (7D)' = 'Water control (7-Day)',
                                                   '16 - Water control (14D)' = 'Water control (14-Day)')

#METHANOL CONTROL
new_df_sw$Treatment_Timepoint1 <- recode(new_df_sw$Treatment_Timepoint, 
                                                   '01 - Field Soil' = 'Field Soil',
                                                   '02 - Predry' = 'Predry',
                                                   '03 - Postdry' = 'Postdry',
                                                   '04 - Water acclimated' ='Water acclimated',
                                                   '05 - Methanol control (1D)' = 'Methanol control (1-Day)',
                                                   '09 - Methanol control (7D)' = 'Methanol control (7-Day)',
                                                   '13 - Methanol control (14D)' = 'Methanol control (14-Day)')










##pivot wider to plot averages, do the rowsums sorting after filtering to switchgrass? if i do it before then bean abundances may obscure the sw rowsum abudances (this is unlikely given that the OTUs picked are indicators for a specific crop and likely will have more abundances in that crop.)

new_df_sw_wide<- new_df_sw %>% select(ASV, Abundance, SampleID) %>%pivot_wider(names_from="SampleID", values_from="Abundance")  
new_df_sw_wide_nodeseq <- new_df_sw_wide %>% filter(ASV!="c335e93a51cbc647ce7c8e354147e462" & ASV!="a78c7c528e0197b7e64605b4a076700a" & ASV!="2985946aa8517aed16c91ddc3afbbd59" & ASV!="76b7130f1162a090d1f340b4a29e632b")


rnames <-new_df_sw_wide_nodeseq$ASV
View(new_df_sw_wide_nodeseq) 
new_df_sw_wide_nodeseq<- new_df_sw_wide_nodeseq %>% select (!ASV)
rownames(new_df_sw_wide_nodeseq) <- rnames
View(new_df_sw_wide_nodeseq)

new_df_sw_wide_nodeseq<- as.matrix(new_df_sw_wide_nodeseq)
new_df_sw_wide_nodeseq <- new_df_sw_wide_nodeseq[order(rowSums(new_df_sw_wide_nodeseq), decreasing = TRUE),]
new_df_sw_wide_nodeseq <- new_df_sw_wide_nodeseq[c(1:50),] ##select top 50 most abundant Classes
View(new_df_sw_wide_nodeseq)

#add the deseq2 otus here 

#new_df_sw_deseq<- new_df_sw %>% filter(ASV==c("c335e93a51cbc647ce7c8e354147e462","a78c7c528e0197b7e64605b4a076700a","2985946aa8517aed16c91ddc3afbbd59","76b7130f1162a090d1f340b4a29e632b","a78c7c528e0197b7e64605b4a076700a"))
#new_df_sw_deseq_wide <- new_df_sw_deseq %>% select(ASV, Abundance, SampleID) %>%pivot_wider(names_from="SampleID", values_from="Abundance")  
 new_df_sw_wide <- new_df_sw %>% select(ASV, Abundance, SampleID) %>%pivot_wider(names_from="SampleID", values_from="Abundance")  
View(new_df_sw_wide) 

new_df_sw_wide_deseq<- new_df_sw_wide %>% filter(ASV=="c335e93a51cbc647ce7c8e354147e462" | ASV=="a78c7c528e0197b7e64605b4a076700a" | ASV=="2985946aa8517aed16c91ddc3afbbd59" | ASV=="76b7130f1162a090d1f340b4a29e632b")
rnames <-new_df_sw_wide_deseq$ASV
new_df_sw_wide_deseq<- new_df_sw_wide_deseq %>% select (!ASV)
rownames(new_df_sw_wide_deseq) <- rnames
View(new_df_sw_wide_deseq)

##rowbind the two above

new_df_sw_wide <- rbind(new_df_sw_wide_nodeseq, new_df_sw_wide_deseq)



##merging with dna cdna data to get NAs inserted in dataframe
##using rarefied DNA and cDNA copies with 15k reads.
dna<- read.csv("DNAcopy_final.csv", header=TRUE, row.names = 1)
cdna <- read.csv("RNAcopy_final.csv", header=TRUE, row.names = 1)

View(dna)
View(cdna)

dna_50taxa <- dna[rownames(dna)%in%rownames(new_df_sw_wide),]
View(dna_50taxa)

cdna_50taxa <- cdna[rownames(cdna)%in%rownames(new_df_sw_wide),]
View(cdna_50taxa)

dna1<- dna_50taxa %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="dna", -ASV)#gather columns into key value pairs
View(dna1)
cdna1<- cdna_50taxa %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="cdna", -ASV)
View(cdna1)
merge.df = full_join(dna1, cdna1, by=c("ASV", "SampleID"))
View(merge.df)

merge.df[is.na(merge.df)] <-0
View(merge.df)
metadata<-read.csv("sample-metadata-decontam-edited2.csv")

newdf = full_join(merge.df, metadata, by=c("SampleID"))
View(newdf)


finaldf = newdf %>%  mutate(ratioMethod1= ifelse(dna == 0 & cdna > 0, 100, cdna/dna)) %>%  
  mutate(dna2= ifelse(dna == 0, 1,dna))  %>%  mutate(ratioMethod2= cdna/dna2) %>% mutate(ratio_nophantoms = cdna/dna)
View(finaldf)
finaldf[is.na(finaldf)] <-0
finaldf$ratioMethod1[!is.finite(finaldf$ratioMethod1)] <- 0 ##using the code above for removing nas also works here as it changed the 0 denominator ratios from NA to 0, no infinite values only NaN
finaldf$ratioMethod2[!is.finite(finaldf$ratioMethod2)] <- 0 ##ratios computed by methods 1 and 2 will always be finite, because of the zero denominator being accounted for phantom taxa?
finaldf$ratio_nophantoms[!is.finite(finaldf$ratio_nophantoms)] <- 0
View(finaldf) 


##code taxa as inactive (0) or not detected in dna (NA)


finaldf_coded <- finaldf %>% #mutate(dna_code=ifelse(dna2==0, "not_detected", dna2)) 
  mutate(activity= ifelse(ratioMethod2 < 1, "inactive", "active"))
View(finaldf_coded) 

##Abundant_active dataset abundance combined with finaldf dataset 
Abundant_Active_new <- cbind(ASV = rownames(new_df_sw_wide), new_df_sw_wide)
View(Abundant_Active_new)
Abundant_Active_new<- as.data.frame(Abundant_Active_new)
Abundant_Active_long <-pivot_longer(Abundant_Active_new, !ASV, names_to="SampleID", values_to="AbundanceDNA_Active")
View(Abundant_Active_long)
merge_df_new<- left_join(Abundant_Active_long, finaldf_coded, by=c("ASV", "SampleID"), how='left')
View(merge_df_new)
merge_df_new_NA<- merge_df_new %>% mutate(AbundanceDNA_Active_NA= ifelse(dna==0 & cdna==0, NA, AbundanceDNA_Active))
View(merge_df_new_NA)

##pivot wider
merge_df_new_NA<- merge_df_new_NA %>% mutate(Treat_Time_ID=paste(Treatment_Timepoint, SampleID))
merge_df_new_NA_wide <- merge_df_new_NA %>% select(ASV, AbundanceDNA_Active_NA, Treat_Time_ID) %>%pivot_wider(names_from="Treat_Time_ID", values_from="AbundanceDNA_Active_NA")  
View(merge_df_new_NA_wide) 
data<- as.data.frame(merge_df_new_NA_wide)

##create matrix


#option 2, both options work
rnames <-data[,1] 
rownames(data) <- rnames  # assign row names
View(data)
mat_data<- data %>% select(!ASV)
View(mat_data)

####NOT DOING THIS FOR THE BUBBLE PLOT, ONLY FOR HEATMAP
#REMOVE OTUS that cant be plotted in heatmap, note which ones are removed, if there are only numbers nothing is removed. if there are NAs some may be removed.

#mat_data<- data.matrix(mat_data)


sum(is.na(as.matrix(dist(mat_data))))

giveNAs = which(is.na(as.matrix(dist(mat_data))),arr.ind=TRUE)
head(giveNAs)
#mat_data[c(1,17),]

tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
  sum(is.na(as.matrix(dist(mat_data[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))] ##https://stackoverflow.com/questions/61469201/pheatmap-won-t-cluster-rows-na-nan-inf-in-foreign-function-call-arg-10 
rmv ## the intention with this code is not to remove rows but to see for which ASVs are we not able to calculate a distance for the heatmap

##use this line only if something needs to be removed
mat_data = mat_data[-as.numeric(rmv),] 
View(mat_data)


##########################
###skip this step for bubble plot to preserve ASVs FOR MERGING LATER
####merge taxonomy with OTU hash id

taxonomy<- read.csv("taxonomy-rep-seqs-or-99-edited.csv", row.names=1)
View(taxonomy)
taxonomy_family<- taxonomy %>% select(c("Family", "Class"))
View(taxonomy_family)

mat_data_merge <- merge(mat_data, taxonomy_family, by="row.names")
View(mat_data_merge)
mat_data_merge$OTUclass <- paste(mat_data_merge$Class, mat_data_merge$Row.names)
View(mat_data_merge)

rnames <-mat_data_merge$OTUclass

rownames(mat_data_merge) <- rnames  # assign row names
View(mat_data_merge)
mat_data_merge_otuclass<- mat_data_merge %>% select(!Row.names)%>% select (!Family) %>% select(!Class) %>% select(!OTUclass)
View(mat_data_merge_otuclass)

####################
##write these files later for easy replotting preserving original sample names, remember to save the mat_data files here instead
##switchgrass
write.csv(mat_data, "mat_data_sw_aba_01.28.25.csv")
write.csv(mat_data, "mat_data_sw_sa_01.28.25.csv")
write.csv(mat_data, "mat_data_sw_met_con_01.28.25.csv")




#bean
write.csv(mat_data, "mat_data_bean_aba_01.28.25.csv")
write.csv(mat_data, "mat_data_bean_sa_01.28.25.csv")
write.csv(mat_data, "mat_data_bean_met_con_01.28.25.csv")
#########################

if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}
install_github("raivokolde/pheatmap")
library(pheatmap)

library(vegan)


##CHANGING THE SAMPLE IDS in csv file before importing back in

#######################
##importing here will not work for the bubble plot because the sample names need to be consistent with the mat_data object above and these were manually changed in the csv files
mat_data_merge_otuclass<- read.csv("mat_data_merge_otuclass_sw_aba.csv", row.names=1)
##########################



mat_data<- as.data.frame(mat_data)
#mat_data_merge_otuclass<- as_tibble(mat_data_merge_otuclass)

class(mat_data)
#mat_data_merge_otuclass<- as.numeric(mat_data_merge_otuclass)
##replace NAs with zero
matrix<- as.matrix(mat_data)
matrix
View(matrix)
matrix[is.na(matrix)] <-0
View(matrix)

mat_num<- matrix(as.numeric(matrix), ncol=ncol(matrix)) 
row.names(mat_num)<- row.names(mat_data)
colnames(mat_num)<- colnames(mat_data)


##compute max standardization

matrix_max_standardize_bean <-decostand(mat_num, method = "max", MARGIN = 1) ##from Jackson's phil trans paper
View(matrix_max_standardize_bean)

#####################################
##dont use code below for bubble plot
##put back NAs where they were before
matrix1<- as.matrix(mat_data_merge_otuclass)
matrix1
View(matrix1)
matrix1[is.na(matrix1)] <- -100
View(matrix1)
mask <- matrix1 == -100
mask

matrix_max_standardize_bean[mask] <- matrix1[mask]
View(matrix_max_standardize_bean)

matrix_max_standardize_bean[matrix_max_standardize_bean == -100]<- NA
View(matrix_max_standardize_bean) ##gives the standardized values with NA


#######################
####start from here for plotting
######################

###first convert this max standardize to long format and then merge back with merge_df_new_NA by sample_name and ASV
data<- as.data.frame(matrix_max_standardize_bean)
data$ASV<- row.names(data)
library(viridis)

data_long<- data %>% pivot_longer(!ASV, names_to="Treat_Time_ID", values_to="Abund_max_standard")

data_final<- left_join(merge_df_new_NA, data_long, by=c("ASV","Treat_Time_ID"))


##add the tazonomy here 
taxonomy<- read.csv("taxonomy-rep-seqs-or-99-edited.csv", row.names=1)
View(taxonomy)
taxonomy_family<- taxonomy %>% select(c("Family", "Class"))
View(taxonomy_family)
taxonomy_family$ASV<- row.names(taxonomy_family)


mat_data_merge <- left_join(data_final, taxonomy_family, by="ASV")
View(mat_data_merge)
mat_data_merge$OTUclass <- paste(mat_data_merge$Class, mat_data_merge$ASV)
View(mat_data_merge)

##SAVE THESE FILES TOO

write.csv(mat_data_merge, "mat_data_merge_sw_aba_01.28.25_withtax.csv")
write.csv(mat_data_merge, "mat_data_merge_sw_sa_01.28.25_withtax.csv")
write.csv(mat_data_merge, "mat_data_merge_sw_met_con_01.28.25_withtax.csv")

write.csv(mat_data_merge, "mat_data_merge_bean_aba_01.28.25_withtax.csv")
write.csv(mat_data_merge, "mat_data_merge_bean_sa_01.28.25_withtax.csv")
write.csv(mat_data_merge, "mat_data_merge_bean_met_con_01.28.25_withtax.csv")
########

#create the breaks
#bk2 = unique(c(0.00000000, seq(0.00000001, 1.00000000, length=20)))

col1 = colorRampPalette(c("gray"))(1)

col2 =viridis(20)

colors2 <- c(col1, col2)


#######not needed
plot_heatmap_sw <- pheatmap::pheatmap(matrix_max_standardize_bean, ##if plotting sample by sample and not the means then use matrix_max_standardize_sw
                                                        cluster_cols = F, border_color = "white", na_col = "white", clustering_distance_rows="euclidean",clustering_method = "complete",
                                                        #main="Bean indicator species (Unique to planted+drought condition)",
                                                        color = colors2,
                                                        breaks = bk2
)
##################
##replotting as per AS edits

sw_aba<-read.csv("mat_data_merge_sw_aba_01.28.25_withtax.csv")
sw_sa<- read.csv("mat_data_merge_sw_sa_01.28.25_withtax.csv")
bean_aba<- read.csv("mat_data_merge_bean_aba_01.28.25_withtax.csv")
bean_sa<- read.csv("mat_data_merge_bean_sa_01.28.25_withtax.csv")

#data_all<- rbind(sw_aba,sw_sa, bean_aba, bean_sa)

mat_data_merge_size_swaba<- sw_aba%>% group_by(ASV, Treatment_Timepoint) %>% mutate(count=sum(Abund_max_standard>0)/n()) %>% mutate(mean=mean(Abund_max_standard)) %>% mutate(inactive_NA=ifelse(all(dna==0 & cdna==0), NA, mean))
mat_data_merge_size_swsa<- sw_sa%>% group_by(ASV, Treatment_Timepoint) %>% mutate(count=sum(Abund_max_standard>0)/n()) %>% mutate(mean=mean(Abund_max_standard)) %>% mutate(inactive_NA=ifelse(all(dna==0 & cdna==0), NA, mean))
mat_data_merge_size_beanaba<- bean_aba%>% group_by(ASV, Treatment_Timepoint) %>% mutate(count=sum(Abund_max_standard>0)/n()) %>% mutate(mean=mean(Abund_max_standard)) %>% mutate(inactive_NA=ifelse(all(dna==0 & cdna==0), NA, mean))
mat_data_merge_size_beansa<- bean_sa%>% group_by(ASV, Treatment_Timepoint) %>% mutate(count=sum(Abund_max_standard>0)/n()) %>% mutate(mean=mean(Abund_max_standard)) %>% mutate(inactive_NA=ifelse(all(dna==0 & cdna==0), NA, mean))

mat_data_merge_size_swaba<- mat_data_merge_size_swaba %>% mutate(grid="Switchgrass ABA")
mat_data_merge_size_swsa<- mat_data_merge_size_swsa %>% mutate(grid="Switchgrass SA")
mat_data_merge_size_beanaba<- mat_data_merge_size_beanaba %>% mutate(grid="Bean ABA")
mat_data_merge_size_beansa<- mat_data_merge_size_beansa %>% mutate(grid="Bean SA")


bb<-c(0.820,0.600,0.400,0.200,0.00001) # define breaks.
lb<-c(1,0.8,0.6,0.4,0.2,0.08)
ll<-c("1","0.8","0.6","0.4","0.2","0.08") # labels.

#class(mat_data_merge_size$Treatment_Timepoint)

#SW ABA
mat_data_merge_size_swaba$Treatment_Timepoint_new<- dplyr::recode(mat_data_merge_size_swaba$Treatment_Timepoint, "01 - Field Soil" = "Field Soil",
                                                     "02 - Predry" ="Predry",
                                                     "03 - Postdry"="Postdry",
                                                     "04 - Water acclimated" = "Water acclimated",
                                                     "06 - ABA (1D)" = "ABA (1D)",
                                                     "10 - ABA (7D)" = "ABA (7D)",
                                                     "14 - ABA (14D)" = "ABA (14D)")
##SW SA
mat_data_merge_size_swsa$Treatment_Timepoint_new<- dplyr::recode(mat_data_merge_size_swsa$Treatment_Timepoint, "01 - Field Soil" = "Field Soil",
                                                     "02 - Predry" ="Predry",
                                                     "03 - Postdry"="Postdry",
                                                     "04 - Water acclimated" = "Water acclimated",
                                                     "07 - SA (1D)" = "SA (1D)",
                                                     "11 - SA (7D)" = "SA (7D)",
                                                     "15 - SA (14D)" = "SA (14D)")

#bean ABA
mat_data_merge_size_beanaba$Treatment_Timepoint_new<- dplyr::recode(mat_data_merge_size_beanaba$Treatment_Timepoint, "01 - Field Soil" = "Field Soil",
                                                     "02 - Predry" ="Predry",
                                                     "03 - Postdry"="Postdry",
                                                     "04 - Water acclimated" = "Water acclimated",
                                                     "06 - ABA (1D)" = "ABA (1D)",
                                                     "10 - ABA (7D)" = "ABA (7D)",
                                                     "14 - ABA (14D)" = "ABA (14D)")
##bean SA
mat_data_merge_size_beansa$Treatment_Timepoint_new<- dplyr::recode(mat_data_merge_size_beansa$Treatment_Timepoint, "01 - Field Soil" = "Field Soil",
                                                     "02 - Predry" ="Predry",
                                                     "03 - Postdry"="Postdry",
                                                     "04 - Water acclimated" = "Water acclimated",
                                                     "07 - SA (1D)" = "SA (1D)",
                                                     "11 - SA (7D)" = "SA (7D)",
                                                     "15 - SA (14D)" = "SA (14D)")



##SW met con (methanol specific code)
mat_data_merge_size$Treatment_Timepoint_new<- recode(mat_data_merge_size$Treatment_Timepoint, "01 - Field Soil" = "Field Soil",
                                                     "02 - Predry" ="Predry",
                                                     "03 - Postdry"="Postdry",
                                                     "04 - Water acclimated" = "Water acclimated",
                                                     "05 - Methanol control (1D)" = "Methanol control (1D)",
                                                     "09 - Methanol control (7D)" = "Methanol control (7D)",
                                                     "13 - Methanol control (14D)" = "Methanol control (14D)")

mat_data_merge_size$Treatment_Timepoint_new<- factor(mat_data_merge_size$Treatment_Timepoint_new, levels=c("Field Soil", "Predry", "Postdry", "Water acclimated", "Methanol control (1D)","Methanol control (7D)", "Methanol control (14D)"), ordered=TRUE)

#merging
data_all_size<- rbind(mat_data_merge_size_swaba, mat_data_merge_size_swsa, mat_data_merge_size_beanaba, mat_data_merge_size_beansa)

##ordering levels separately (did not do this)
mat_data_merge_size_swaba$Treatment_Timepoint_new<- factor(mat_data_merge_size_swaba$Treatment_Timepoint_new, levels=c("Field Soil", "Predry", "Postdry", "Water acclimated", "ABA (1D)","ABA (7D)", "ABA (14D)"), ordered=TRUE)
mat_data_merge_size_swsa$Treatment_Timepoint_new<- factor(mat_data_merge_size_swsa$Treatment_Timepoint_new, levels=c("Field Soil", "Predry", "Postdry", "Water acclimated", "SA (1D)","SA (7D)", "SA (14D)"), ordered=TRUE)

mat_data_merge_size_beanaba$Treatment_Timepoint_new<- factor(mat_data_merge_size_beanaba$Treatment_Timepoint_new, levels=c("Field Soil", "Predry", "Postdry", "Water acclimated", "ABA (1D)","ABA (7D)", "ABA (14D)"), ordered=TRUE)
mat_data_merge_size_beansa$Treatment_Timepoint_new<- factor(mat_data_merge_size_beansa$Treatment_Timepoint_new, levels=c("Field Soil", "Predry", "Postdry", "Water acclimated", "SA (1D)","SA (7D)", "SA (14D)"), ordered=TRUE)

##ordering levels together (this after rbind worked)
data_all_size$Treatment_Timepoint_new<- factor(data_all_size$Treatment_Timepoint_new, levels=c("Field Soil", "Predry", "Postdry", "Water acclimated", "ABA (1D)","ABA (7D)", "ABA (14D)","SA (1D)","SA (7D)", "SA (14D)"), ordered=TRUE)

#class(data_all_size$OTUclass)

#data_all_size$OTUclass<- factor(data_all_size$OTUclass, levels= unique(data_all_size$OTUclass))

##calling the deseq OTUS


## SW ABA
deseq<-c("c__Myxococcia c335e93a51cbc647ce7c8e354147e462", "c__Actinobacteria a78c7c528e0197b7e64605b4a076700a", "c__Actinobacteria 2985946aa8517aed16c91ddc3afbbd59")


##SW_SA
deseq<-c("c__Actinobacteria a78c7c528e0197b7e64605b4a076700a")

##SW methanol
deseq<-c("c__Myxococcia c335e93a51cbc647ce7c8e354147e462", "c__Actinobacteria a78c7c528e0197b7e64605b4a076700a", "c__Actinobacteria 2985946aa8517aed16c91ddc3afbbd59")


##bean SA  
deseq<-c("c__Actinobacteria a78c7c528e0197b7e64605b4a076700a", "c__Actinobacteria 76b7130f1162a090d1f340b4a29e632b")

##bean methanol control  
deseq<-c("c__Actinobacteria a78c7c528e0197b7e64605b4a076700a", "c__Actinobacteria 76b7130f1162a090d1f340b4a29e632b")


##calling deseq otus at once
deseq<-c("c__Myxococcia c335e93a51cbc647ce7c8e354147e462", "c__Actinobacteria a78c7c528e0197b7e64605b4a076700a", "c__Actinobacteria 2985946aa8517aed16c91ddc3afbbd59", "c__Actinobacteria 76b7130f1162a090d1f340b4a29e632b")


data_all_size_color<- data_all_size %>% mutate(color=ifelse(OTUclass %in% deseq, "red", "black"))
data_all_size_color$color<- factor(data_all_size_color$color, levels=c("red", "black"), ordered=TRUE)
data_all_size_color<- data_all_size_color %>% arrange(color)

# **CRITICAL:  Force the OTUclass factor order**
data_all_size_color$OTUclass <- factor(data_all_size_color$OTUclass, levels = unique(data_all_size_color$OTUclass), ordered = TRUE)


#data_ordered <- data_all_size_color[order(data_all_size_color$color),]

library(forcats)
#data_all_size_color$OTUclass<- fct_reorder(data_all_size_color$OTUclass, data_all_size_color$color)

#new label function
# Create the labels using a *lookup* – much more efficient and correct
label_function <- function(otu_class) {
  color_lookup <- setNames(data_all_size_color$color, data_all_size_color$OTUclass) # Named vector for lookup
  paste0("<span style='color:", color_lookup[otu_class], ";'>", otu_class, "</span>")
}



##old label function
label_function <- function(x) {
  paste0("<span style='color:", data_all_size_color$color[match(x, data_all_size_color$OTUclass)], ";'>", x, "</span>")
}

##bubble plot all

library(viridis)
library(ggtext)

plot<- ggplot(data_all_size_color,aes(x=Treatment_Timepoint_new,y=OTUclass,size=count,color=inactive_NA))+
  facet_grid(~grid, scales="free")+
  theme(strip.text.x = element_text(size=20, face="bold"))+
  geom_point()+
  geom_point(data=subset(data_all_size, is.na(inactive_NA)), aes(shape='NA'), size=3) +
  scale_shape_manual(values=c('NA'=3), name=NULL, labels="Not detected in DNA or cDNA") +
  #theme(strip.text.x=element_text(size=10))+
  #geom_text(aes(label = OTUclass, color = color), vjust = -0.5) +
  scale_colour_viridis(limits=c(0.00001,0.820),breaks=bb)+
  scale_size_continuous(limits=c(0.08,1),labels=ll,breaks=lb, name="Frequency of activity across \nexperimental replicates")+
  #ggtitle("Switchgrass ABA")+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  xlab("Treatment")+
  labs(color="Mean Relative\nAbundance")+
  scale_y_discrete(labels=label_function)+
   theme(axis.text.x = element_text(size=28, face="bold", angle = 90, hjust=1, vjust=0.5),
         axis.title.x = element_text(size=35, face="bold"))+
  theme(axis.text.y=element_markdown(size=15, face="bold"),
        axis.title.y=element_text(size=35, face="bold")) + 
  theme(legend.title = element_text(size=35))+ theme(legend.text = element_text(size=30))

plot
ggsave(filename = "bubble_plot_all_02.19.25.tiff", plot = plot,
       width = 70,
       height = 55, units = c("cm"),
       dpi = 300)

##bubble plot SW ABA
plot<- ggplot(mat_data_merge_size_color,aes(x=Treatment_Timepoint_new,y=OTUclass,size=count,color=inactive_NA))+
  geom_point()+
  geom_point(data=subset(mat_data_merge_size, is.na(inactive_NA)), aes(shape='NA'), size=3) +
  scale_shape_manual(values=c('NA'=3), name=NULL, labels="Not detected in DNA or cDNA") +
  #theme(strip.text.x=element_text(size=10))+
  #geom_text(aes(label = OTUclass, color = color), vjust = -0.5) +
  scale_colour_viridis(limits=c(0.00005,0.820),breaks=bb)+
  scale_size_continuous(limits=c(0.0001,1),labels=ll,breaks=lb, name="Frequency of activity across \nexperimental and technical replicates")+
  ggtitle("Switchgrass ABA")+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  xlab("Treatment")+
  labs(color="Mean Relative\nAbundance")+
  scale_y_discrete(labels=label_function)+
   theme(axis.text.x = element_text(size=12, face="bold", angle = 90, hjust=1, vjust=0.5),
         axis.title.x = element_text(size=13))+
  theme(axis.text.y=element_markdown(),
        axis.title.y=element_text(size=13)) + 
  theme(legend.title = element_text(size=14))+ theme(legend.text = element_text(size=13))
 # ylab("OTUClass") 


#SW SA


plot<- ggplot(mat_data_merge_size_color,aes(x=Treatment_Timepoint_new,y=OTUclass,size=count,color=inactive_NA))+
  geom_point()+
  geom_point(data=subset(mat_data_merge_size, is.na(inactive_NA)), aes(shape='NA'), size=3) +
  scale_shape_manual(values=c('NA'=3), name=NULL, labels="Not detected in DNA or cDNA") +
  #theme(strip.text.x=element_text(size=10))+
  #geom_text(aes(label = OTUclass, color = color), vjust = -0.5) +
  scale_colour_viridis(limits=c(0.00005,0.820),breaks=bb)+
  scale_size_continuous(limits=c(0.0001,1),labels=ll,breaks=lb, name="Frequency of activity across \nexperimental and technical replicates")+
  ggtitle("Switchgrass SA")+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  xlab("Treatment")+
  labs(color="Mean Relative Abundance")+
  scale_y_discrete(labels=label_function)+
   theme(axis.text.x = element_text(size=12, face="bold", angle = 90, hjust=1, vjust=0.5),
         axis.title.x = element_text(size=13))+
  theme(axis.text.y=element_markdown(),
        axis.title.y=element_text(size=13)) + 
  theme(legend.title = element_text(size=14))+ theme(legend.text = element_text(size=13))
 # ylab("OTUClass")

  
plot 

##SW methanol

plot<- ggplot(mat_data_merge_size_color,aes(x=Treatment_Timepoint_new,y=OTUclass,size=count,color=inactive_NA))+
  geom_point()+
  geom_point(data=subset(mat_data_merge_size, is.na(inactive_NA)), aes(shape='NA'), size=3) +
  scale_shape_manual(values=c('NA'=3), name=NULL, labels="Not detected in DNA or cDNA") +
  #theme(strip.text.x=element_text(size=10))+
  #geom_text(aes(label = OTUclass, color = color), vjust = -0.5) +
  scale_colour_viridis(limits=c(0.00005,0.820),breaks=bb)+
  scale_size_continuous(limits=c(0.0001,1),labels=ll,breaks=lb, name="Frequency of activity across \nexperimental and technical replicates")+
  ggtitle("Switchgrass Methanol Control")+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  xlab("Treatment")+
  labs(color="Mean Relative Abundance")+
  scale_y_discrete(labels=label_function)+
   theme(axis.text.x = element_text(size=12, face="bold", angle = 90, hjust=1, vjust=0.5),
         axis.title.x = element_text(size=13))+
  theme(axis.text.y=element_markdown(),
        axis.title.y=element_text(size=13)) + 
  theme(legend.title = element_text(size=14))+ theme(legend.text = element_text(size=13))
 # ylab("OTUClass")

  
plot 

##Bean ABA

plot<- ggplot(mat_data_merge_size,aes(x=Treatment_Timepoint_new,y=OTUclass,size=count,color=inactive_NA))+ ##NO DESEQ RESULTS FOR THIS SO USING THE mat_data_merge_size data frame
  geom_point()+
  geom_point(data=subset(mat_data_merge_size, is.na(inactive_NA)), aes(shape='NA'), size=3) +
  scale_shape_manual(values=c('NA'=3), name=NULL, labels="Not detected in DNA or cDNA") +
  #theme(strip.text.x=element_text(size=10))+
  #geom_text(aes(label = OTUclass, color = color), vjust = -0.5) +
  scale_colour_viridis(limits=c(0.00001,0.820),breaks=bb)+
  scale_size_continuous(limits=c(0.0001,1),labels=ll,breaks=lb, name="Frequency of activity across \nexperimental and technical replicates")+
  ggtitle("Bean ABA")+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  xlab("Treatment")+
  labs(color="Mean Relative Abundance")+
#scale_y_discrete(labels=label_function)+
   theme(axis.text.x = element_text(size=12, face="bold", angle = 90, hjust=1, vjust=0.5),
         axis.title.x = element_text(size=13))+
  theme(axis.text.y=element_markdown(),
        axis.title.y=element_text(size=13)) + 
  theme(legend.title = element_text(size=14))+ theme(legend.text = element_text(size=13))
 # ylab("OTUClass") 

plot


##BEAN SA
plot<- ggplot(mat_data_merge_size_color,aes(x=Treatment_Timepoint_new,y=OTUclass,size=count,color=inactive_NA))+
  geom_point()+
  geom_point(data=subset(mat_data_merge_size, is.na(inactive_NA)), aes(shape='NA'), size=3) +
  scale_shape_manual(values=c('NA'=3), name=NULL, labels="Not detected in DNA or cDNA") +
  #theme(strip.text.x=element_text(size=10))+
  #geom_text(aes(label = OTUclass, color = color), vjust = -0.5) +
  scale_colour_viridis(limits=c(0.00005,0.820),breaks=bb)+
  scale_size_continuous(limits=c(0.0001,1),labels=ll,breaks=lb, name="Frequency of activity across \nexperimental and technical replicates")+
  ggtitle("Bean SA")+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  xlab("Treatment")+
  labs(color="Mean Relative Abundance")+
  scale_y_discrete(labels=label_function)+
   theme(axis.text.x = element_text(size=12, face="bold", angle = 90, hjust=1, vjust=0.5),
         axis.title.x = element_text(size=13))+
  theme(axis.text.y=element_markdown(),
        axis.title.y=element_text(size=13)) + 
  theme(legend.title = element_text(size=14))+ theme(legend.text = element_text(size=13))
 # ylab("OTUClass")

  
plot 

##bean methanol
plot<- ggplot(mat_data_merge_size_color,aes(x=Treatment_Timepoint_new,y=OTUclass,size=count,color=inactive_NA))+
  geom_point()+
  geom_point(data=subset(mat_data_merge_size, is.na(inactive_NA)), aes(shape='NA'), size=3) +
  scale_shape_manual(values=c('NA'=3), name=NULL, labels="Not detected in DNA or cDNA") +
  #theme(strip.text.x=element_text(size=10))+
  #geom_text(aes(label = OTUclass, color = color), vjust = -0.5) +
  scale_colour_viridis(limits=c(0.00005,0.820),breaks=bb)+
  scale_size_continuous(limits=c(0.0001,1),labels=ll,breaks=lb, name="Frequency of activity across \nexperimental and technical replicates")+
  ggtitle("Bean Methanol Control")+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+
  xlab("Treatment")+
  labs(color="Mean Relative Abundance")+
  scale_y_discrete(labels=label_function)+
   theme(axis.text.x = element_text(size=12, face="bold", angle = 90, hjust=1, vjust=0.5),
         axis.title.x = element_text(size=13))+
  theme(axis.text.y=element_markdown(),
        axis.title.y=element_text(size=13)) + 
  theme(legend.title = element_text(size=14))+ theme(legend.text = element_text(size=13))
 # ylab("OTUClass")

  
plot 


##switchgrass

ggsave(filename = "sw_aba_top50_aggregated_01.28.25.tiff", plot = plot,
       width = 30,
       height = 27, units = c("cm"),
       dpi = 300)


ggsave(filename = "sw_sa_top50_aggregated_01.28.25.tiff", plot = plot,
       width = 30,
       height = 27, units = c("cm"),
       dpi = 300)

ggsave(filename = "sw_met_con_top50_aggregated_01.28.25.tiff", plot = plot,
       width = 30,
       height = 27, units = c("cm"),
       dpi = 300)

##bean 
ggsave(filename = "bean_aba_top50_aggregated_01.28.25.tiff", plot = plot,
       width = 30,
       height = 27, units = c("cm"),
       dpi = 300)

ggsave(filename = "bean_sa_top50_aggregated_01.28.25.tiff", plot = plot,
       width = 30,
       height = 27, units = c("cm"),
       dpi = 300)
ggsave(filename = "bean_met_con_top50_aggregated_01.28.25.tiff", plot = plot,
       width = 30,
       height = 27, units = c("cm"),
       dpi = 300)


#######################################
###CAP analysis with soil chemistry data
#######################################

##make a table for the soil characteristics post incubation and pre-incubation, check for acidity 

##making ph nitrate data join with metadata file 

data1 <- read.csv("sample-metadata-edaphic.csv")
data2<- read.csv("sample-metadata-decontam-edited.csv")
View(data1)
View(data2)
full_data<- full_join(data1, data2, by="SampleID")
View(full_data)
write.csv(full_data, "sample-metadata-edaphic-constrained.csv")

##sample names check with microbiome and soil chemistry data

data<- read.csv("sample-metadata-edaphic-constrained.csv", row.names=1)
cap_sw<- data %>% filter(Plant=="Switchgrass", Nucleic_Acid=="RNA", pH!="NA")
SampleID<- row.names(cap_sw)
cap_sw<- cap_sw %>% mutate(SampleID=SampleID)

#unique(cap_sw$SampleID) ## 51 SAMPLES

data<- read.csv("sample-metadata-edaphic-constrained.csv", row.names=1)
cap_bean<- data %>% filter(Plant=="Bean", Nucleic_Acid=="RNA", pH!="NA")
SampleID<- row.names(cap_bean)
cap_bean<- cap_bean %>% mutate(SampleID=SampleID)

#unique(cap_bean$SampleID) ## 51 SAMPLES


##reading in new phyloseq object with all samples in CB
otu_final<- read.csv("activefinalnozero_method2_changing_to_dnaabun_final.csv", row.names=1) ##change name to WLE or CB as intended
View(otu_final)
any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #20954 otus with no counts (wle), 15252 otus with no counts (cb)
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) #8887 otus total for bean and switchgrass...

#data_chem<- read.csv("sample-metadata-edaphic-constrained.csv")
#View(data_chem)
#metadata_cb_wle<- read.csv("sample-metadata-cb-wle-cap.csv")
#View(metadata_cb_wle)
#metadata_cap<- merge(data_chem, metadata_cb_wle, by="SampleID")

#View(metadata_cap)

#unique(metadata_cap$SampleID) #205 samples

#samples_lost<- anti_join(metadata_cb_wle, metadata_cap, by = "SampleID")
#View(samples_lost) ##negative controls lost only


##subset to CB samples

#cap_cb<- metadata_cap %>% filter(region.x=="Chesapeake Bay")

#View(cap_cb)
##subset to WLE samples
#cap_wle<- metadata_cap %>% filter(region.x=="Lake Erie ")


#View(cap_wle)


##melt to wide format 

#cap_sw_new<- cap_sw%>% dplyr::select(!analysis) %>% pivot_wider(names_from="name", values_from="value")  
#View(cap_sw_new) ##118 samples
##cap_bean_new<- cap_bean%>% dplyr::select(!analysis) %>% pivot_wider(names_from="name", values_from="value")  
#View(cap_bean_new) #87 samples

#write.csv(cap_cb_new, "cap_cb_new.csv") ##did not save this a second time

#write.csv(cap_wle_new, "cap_wle_new.csv") 

library(phyloseq)


tax<- read.csv("taxonomy-rep-seqs-or-99-edited.csv", sep=",", row.names=1)
#metadata <- read.csv("cap_cb_new.csv", sep=",", row.names=8) ##change as needed for cb and wle



##removing columns that are >50% NAs, adjust separately for CB and WLE
##SW
#metadata_new<- cap_sw %>% dplyr::select(!c(region.y, Sample.description, #Nitrate_meq100g, #Sulfate_meq100g, Phosphate_meq100g, #Bromide_meq100g, Ammonia_meq100g
                                         #    tree_number.y, transect.y, site.y, horizon.y, project, date, DNA.Well, DNA.Plate, PCR.Plate, PCR.Well, Index.ID, Barcode, Potassium_meq100g, dic_ugg, Sodium_meq100g, Fluoride_meq100g, Chloride_meq100g, Magnesium_meq100g, Calcium_meq100g,
                                            # percentOM, gwc_perc)) 
#Bean
#metadata_new<- cap_bean %>% dplyr::select(!c(region.y, Sample.description, #Nitrate_meq100g, #Sulfate_meq100g, #Bromide_meq100g, Ammonia_meq100g, Phosphate_meq100g,
                                         #    tree_number.y, transect.y, site.y, horizon.y, project, date, DNA.Well, DNA.Plate, PCR.Plate, PCR.Well, Index.ID, Barcode,Potassium_meq100g, dic_ugg, Sodium_meq100g, Fluoride_meq100g, Chloride_meq100g, Magnesium_meq100g, Calcium_meq100g,
                                         #   percentOM, gwc_perc)) 

metadata_new<- cap_sw
metadata_new<- cap_bean
View(metadata_new)

##removing rows that are NAs, need to do this because ordinate function does not work with NAs
#metadata_new=na.omit(metadata_new)

#metadata_new<- na.pass(metadata)
View(metadata_new)
tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata_new) ## 76 samples for wle because omit na removed some, 67 samples remain for CB

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged


  
  ##after NA are removed, for SW
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 50 samples ]
sample_data() Sample Data:       [ 50 samples by 14 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]

 
## after NA are removed for Bean
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8887 taxa and 47 samples ]
sample_data() Sample Data:       [ 47 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 8887 taxa by 7 taxonomic ranks ]

  
 
##weighted unifrac distance matrix needs to be computed for CAP plot to work below
set.seed(1)
phyloseq_bray <- phyloseq::distance(phyloseq_merged, method = "bray")
phyloseq_bray


##CAP plot

# First select the variables that significantly explain variation in community composition.
## To do this I'm using the ordistep function which steps through the variables removing ones that don't add to the significance. 
## There are better, more thurough explanations of that online.
### First set new names for all variables that fit nicer on the final figure.

##SW
var_goodnames = data.frame(labels = c("pH", "NO3.ppm", "NH4.ppm", "percent_OM"), ##chloride only present for may samples in CB
                           goodnames = c( "pH", "NO3.ppm", "NH4.ppm", "percent_OM"))

#Bean
var_goodnames = data.frame(labels = c("pH", "NO3.ppm", "NH4.ppm", "percent_OM"),
                           goodnames = c("pH", "NO3.ppm", "NH4.ppm", "percent_OM"))



### Full model with all variables
set.seed(4242)
cap_ord.full <- ordinate(physeq = phyloseq_merged, method = "CAP", distance=phyloseq_bray, #distance = phyloseq_wunifrac, 
                         formula = ~ pH + NO3.ppm + NH4.ppm + percent_OM)
### Null model with no variables
set.seed(4242)
cap_ord.null <- ordinate(physeq = phyloseq_merged, method = "CAP", distance=phyloseq_bray, #distance = phyloseq_wunifrac, 
                         formula = ~ 1)

class(cap_ord.null)

class(cap_ord.full)

### Model selection to get just significant variables
set.seed(4242)
library(vegan)
ordistep.res = ordistep(cap_ord.null, scope = formula(cap_ord.full), perm.max = 1000, trace=F)
goodform = ordistep.res$call$formula

# Get main CAP ordination
## This uses phyloseq package. The formula here is generated by that ordistep function.
set.seed(4242)
cap_ord <- ordinate(physeq = phyloseq_merged, method = "CAP", distance = phyloseq_bray, formula = goodform)

# CAP plot Switchgrass
## For this I use my own code rather than phyloseq plotting for more control. To do this I have to pull out the x (CAP1) and y (CAP2) axes from the ordination and add in my metadata
#key_match<- read.csv("key_match_cap.csv")
#View(key_match)
cap.ord.df1 <- data.frame(vegan::scores(cap_ord, display="sites")) %>%
  tibble::rownames_to_column(var="SampleID") %>%
  dplyr::select(SampleID,CAP1, CAP2) 
  #left_join(key_match, by="sample.id")
cap.ord.df = cap.ord.df1 %>%
  left_join(sample_data(phyloseq_merged), by = "SampleID")

# CAP plot Bean
## For this I use my own code rather than phyloseq plotting for more control. To do this I have to pull out the x (CAP1) and y (CAP2) axes from the ordination and add in my metadata

#key_match<- read.csv("key_match_cap.csv")
#View(key_match)
cap.ord.df1 <- data.frame(vegan::scores(cap_ord, display="sites")) %>%
tibble::rownames_to_column(var="SampleID") %>%
  dplyr::select(SampleID,CAP1, CAP2)
  #left_join(key_match, by="sample.id")
cap.ord.df = cap.ord.df1 %>%
  left_join(sample_data(phyloseq_merged), by = "SampleID")

View(cap.ord.df)
#cap.ord.df$transect.x <- gsub(fixed("wte"), "wetland", cap.ord.df$transect.x)
## Calculate eignevalues and fraction of variation explained by each CAP axis. I use this for the axis labels in the plot
eigvec = vegan::eigenvals(cap_ord)
fracvar = round(eigvec/sum(eigvec)*100, 2)

## Plot initial figure of points
library(ggplot2)

##14 day timepoint comparison
cap.ord.df.14<- cap.ord.df %>% filter(Timepoint=="14-Day")

#field vs lab comparison
cap.ord.df.all <- cap.ord.df


cap.ord.df$Treatment <- factor(cap.ord.df$Treatment,levels = c("Field-Soil (PreHormone)", "Salicylic-acid", "Abscisic-acid", "Methanol-control", "Water-control"), ordered=TRUE)



#cap.ord.df<- cap.ord.df %>% filter(!SampleID=="COMPASS.Dec2021.016")##filter out the outlier for wle COMPASS.Dec2021.016

cap_plot<- ggplot(data=cap.ord.df, aes(x=CAP1, y=CAP2)) +
  geom_point(aes(color=Treatment, shape=Timepoint), size=3)+
  scale_color_manual(values=c("black", "#E69F00",  "#0072B2", "#D55E00","#CC79A7"))+
  #scale_shape_manual(values=c("Abscisic-acid"=15, "Salicyclic-acid"=16, "Methanol-control"=17, "Water-control"=18))+ ##change site names for CB
  labs(x=paste("CAP1 (", fracvar[1], "%)", sep=""),y=paste("CAP2 (", fracvar[2], "%)", sep=""))+
  labs(shape="Timpeoint", fill="Treatment")+theme_classic()


#cap_plot$labels$site.x <- "site"
#cap_plot$labels$transect.x <- "transect"

cap_plot


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

## Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
 # mutate(labels = gsub("\\.", ":", labels))
colnames(arrowdf) = c("labels", "xend", "yend")
arrowdf = arrowdf %>%
  left_join(var_goodnames, by = "labels") %>%
  rename(old_labels = labels) %>%
  rename(labels = goodnames)

## Define the arrow aesthetic mapping
arrow_map <- aes(xend = xend, yend = yend, x = 0, y = 0, 
                 color = NULL)

label_map <- aes(x = xend + 0.02*xend/abs(xend), y = yend, 
                 color = NULL, label = labels)

arrowhead = arrow(length = unit(0.02, "npc"), type = "closed")

# Make a new graphic including labels and axes for variables.
cap.plot = cap_plot + 
  geom_segment(mapping = arrow_map, linewidth = 1.2, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_segment(mapping = arrow_map, linewidth = 0.5, data = arrowdf, color = "#44AA99", arrow = arrowhead) + 
  geom_label(mapping = label_map, data = filter(arrowdf, xend < 0), show.legend = FALSE, size=6*5/14, hjust=1, fill="#44AA99", color="black") +
  geom_label(mapping = label_map, data = filter(arrowdf, xend > 0), show.legend = FALSE, size=6*5/14, hjust=0, fill="#44AA99", color="black") +
  labs(title = paste("Variables explain ", round(100*RsquareAdj(cap_ord)$r.squared, 3), "% of the variation", sep="")) + xlim(-2, 2)+
  theme_classic() + labs(shape="Timepoint", color="Treatment")
 # theme(legend.position="bottom",
       # legend.direction = "vertical") +
  #guides(shape = guide_legend(order = 1, nrow=1),
      #   fill = guide_legend(order = 2, override.aes=list(shape=22), nrow=1))

cap.plot

ggsave(filename = "cap_sw_09.10.24_bray.TIFF", plot = cap.plot,
       width = 17,
       height = 13, units = c("cm"),
       dpi = 700)
ggsave(filename = "cap_bean_09.20.24_bray.TIFF", plot = cap.plot,
       width = 17,
       height = 13, units = c("cm"),
       dpi = 700)




