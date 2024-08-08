###-------------------- Bacterial community analyses of nifH amplicons in Racomitrium lanuginosum (Hedw.)Brid.--------------------###

# Code for moss nifH analyses. Dennis Escolástico-Ortiz. 2022

#devtools::install_github("benjjneb/dada2")

#----------------------------- 1. Checking Read quality -----------------------------------#

library("dada2")
packageVersion("dada2")

nif_sequences<-"~/moss_nifH_second_run"
list.files(nif_sequences)
nifFs <- sort(list.files(nif_sequences, pattern="_R1_001.fastq", full.names = TRUE))
nifRs <- sort(list.files(nif_sequences, pattern="_R2_001.fastq", full.names = TRUE))
nifH.names <- sapply(strsplit(basename(nifFs), "_"), `[`, 1)
plotQualityProfile(nifFs[73:81])
plotQualityProfile(nifRs[73:81])
# The expected length of the nifH amplicons using Ueda19F-R6 primers is 394 bp.

# Place filtered files in filtered/ subdirectory
filtnifFs <- file.path(nif_sequences, "filtered", paste0(nifH.names, "_filt_R1_001.fastq.gz"))
filtnifRs <- file.path(nif_sequences, "filtered", paste0(nifH.names, "_filt_R2_001.fastq.gz"))
names(filtnifFs) <- nifH.names 
names(filtnifRs) <- nifH.names 

out_nifH <- filterAndTrim(nifFs, filtnifFs, nifRs, filtnifRs, truncLen=c(290,250),
                     maxN=0, maxEE=c(2,2), truncQ=2,trimLeft=c(10,10), rm.phix=TRUE,
                     compress=TRUE, multithread=F) # On Windows set multithread=FALSE
head(out_nifH)

# Checking filter sequences
plotQualityProfile(filtnifFs[1:12])
plotQualityProfile(filtnifRs[1:12])


#----------------------------- 2. NifMAP - Analysis of nifH amplicons -----------------------------------#

# Continue to NifMAP pipeline 


#----------------------------- 3. OTU table and formatting to phyloseq -----------------------------------#

### Importing nifMAP results to phyloseq ###

# Installing Phyloseq
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.14")

#install.packages("devtools")

#BiocManager::install("phyloseq")
library("phyloseq"); packageVersion("phyloseq")

setwd("C:/Users/escol/OneDrive - Université Laval/6_Doctorat en Biologie/1_Thesis/Chapter3_Microbiome/Microbiome_bioinformatic_analyses")


OTU_nifH=read.table ("NifMAP_results/otu_table.txt",header=T, sep="")
OTU_nifH

# Extract row names
#otu_nifH[1]-> otu_names
#as.matrix(otu_names)->otu_names
#as.vector(otu_names,mode="any")->otu_names
#class(otu_names)

# Transform data.frame into matrix and add row names
data.matrix(OTU_nifH) -> OTU_nifH
class(OTU_nifH)
OTU_nifH

# Taxonomy table
taxmat=read.table("NifMAP_results/taxonomy_table.txt",header=T,sep="")
taxmat
as.matrix(taxmat)->taxmat
class(taxmat)

# Construct a phyloseq object

# OTU table and taxonomy information
OTU = otu_table(OTU_nifH,taxa_are_rows = TRUE)
TAX=tax_table(taxmat)

sampling_data=read.table("NifMAP_results/sample_data.txt",header=T,sep="")
VAR=sample_data(sampling_data)

OTU
TAX
phylo_nifH=phyloseq(OTU,TAX,VAR)
phylo_nifH

# Sequence data
dna_nifH <- readDNAStringSet("./NifMap_results/nifH_corr_nucl.fasta", format="fasta")
REF=refseq(dna_nifH)

# Sample data
phylo_nifH2=merge_phyloseq(phylo_nifH,VAR,REF)
phylo_nifH2

plot_bar(phylo_nifH2, fill = "Phylum")


#----------------------------- 4. Decontam - Remove contaminants -----------------------------------#

library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
theme_set(theme_bw())

# Check data
head(sample_data(phylo_nifH2))

# Inspect library sizes
df_pn <- as.data.frame(sample_data(phylo_nifH2)) 
df_pn$LibrarySize <- sample_sums(phylo_nifH2)
df_pn <- df_pn[order(df_pn$LibrarySize),]
df_pn
df_pn$Index <- seq(nrow(df_pn))
ggplot(data=df_pn, aes(x=Index, y=LibrarySize, color=Sample_status)) + geom_point() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Identify Contaminants - Frequency
contamdf.freq_nifH <- isContaminant(phylo_nifH2, method="frequency", conc="DNA_conc")
head(contamdf.freq_nifH)
table(contamdf.freq_nifH$contaminant)
head(which(contamdf.freq_nifH$contaminant))

# Identify Contaminants - Prevalence 
sample_data(phylo_nifH2)$is.neg <- sample_data(phylo_nifH2)$Sample_status == "Control Sample"
contamdf.prev_nifH <- isContaminant(phylo_nifH2, method="prevalence", neg="is.neg")
table(contamdf.prev_nifH$contaminant)
head(which(contamdf.prev_nifH$contaminant))

# Change the threshold for prevalence
contamdf.prev05_nifH <- isContaminant(phylo_nifH2, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05_nifH$contaminant)
head(which(contamdf.prev05_nifH$contaminant))
which(contamdf.prev05_nifH$contaminant)
contamdf.prev05_nifH[c(26 ,29,39,40 ,45 ,64 ,70 ,74,131,144),]


# Make phyloseq object of presence-absence in negative controls and true samples
mps_nifH.pa <- transform_sample_counts(phylo_nifH2, function(abund) 1*(abund>0))
mps_nifH.pa.neg <- prune_samples(sample_data(mps_nifH.pa)$Sample_status == "Control Sample", mps_nifH.pa)
mps_nifH.pa.pos <- prune_samples(sample_data(mps_nifH.pa)$Sample_status == "True Sample", mps_nifH.pa)

# Make data.frame of prevalence in positive and negative samples
df.mpa_nifH <- data.frame(pa.pos=taxa_sums(mps_nifH.pa.pos), pa.neg=taxa_sums(mps_nifH.pa.neg),
                     contaminant=contamdf.prev05_nifH$contaminant)
ggplot(data=df.mpa_nifH, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Remove ASVs according to prevalence

phylo_nifH2# 185 OTUs
p_nifH_noncontam <- prune_taxa(!contamdf.prev05_nifH$contaminant, phylo_nifH2)
p_nifH_noncontam # 175 OtUs


##----------------------- 5. Filter according to prevalence ----------------------- ##

# Before proceed we must filter out the negative controls
p_nifH_noncontam # 175 OtUs
p_nifH_noncontam@sam_data
negative_controls<-c("PCR46","PCR47","PCR48","PCR49","PCR58","PCR59")
p_nifH_filt <- prune_samples(!(sample_names(p_nifH_noncontam) %in% negative_controls), p_nifH_noncontam)
p_nifH_filt@sam_data
p_nifH_filt #175 OTUs and 75 samples

# Compute prevalence of each feature, store as data.frame
prevdf_nifH = apply(X = otu_table(p_nifH_filt),
               MARGIN = ifelse(taxa_are_rows(p_nifH_filt), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf_nifH  = data.frame(Prevalence = prevdf_nifH ,
                    TotalAbundance = taxa_sums(p_nifH_filt),
                    tax_table(p_nifH_filt))
plyr::ddply(prevdf_nifH, "Phylum", function(df1){cbind(mean(df1$Prevalence),
                                                  sum(df1$Prevalence))})

# Prevalence: number of samples in which a taxon appears at least once.
# Subset to the remaining phyla
prevdf1_nifH = subset(prevdf_nifH, Phylum %in% get_taxa_unique(p_nifH_filt, "Phylum"))
ggplot(prevdf1_nifH, aes(TotalAbundance, Prevalence / nsamples(p_nifH_filt),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Define prevalence threshold to keep ASV present in at least 5% of samples (Callahan et al. web):
prevalenceThreshold_nifH = 0.05 * nsamples(p_nifH_filt)
prevalenceThreshold_nifH #3.75

## Execute prevalence filter, using `prune_taxa()` function
keepTaxa_nifH = rownames(prevdf1_nifH)[(prevdf1_nifH$Prevalence >= prevalenceThreshold_nifH)]
p_nifH_2 = prune_taxa(keepTaxa_nifH,p_nifH_filt)
table(tax_table(p_nifH_2)[, "Phylum"], exclude = NULL)
p_nifH_2 # 105 OTUs 75 samples

# Check summary per habitat (number of samples, ASV, reads)
table(sample_data(p_nifH_2)$Habitat)
# Forest tundra= 42 samples; Shrub tundra = 33 samples
mean(sample_sums(p_nifH_2))
summary(sample_sums(p_nifH_2))
# Mean number of reads of 6460

#----------------------- 6. Relative abundance  -------------------------------------#

library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

write.table(tax_table(p_nifH_2),"tax_table_filtered.txt",sep = " ")

# Update taxonomic information only for filter OTUs based on TaxID
# Taxonomy table
taxtab=read.table("NifMAP_results/tax_table_update.txt",header=T,sep="")
taxtab
as.matrix(taxtab)->taxtab
class(taxtab)
TAXTABLE=tax_table(taxtab)
p_nifH_2=merge_phyloseq(p_nifH_2,TAXTABLE)


# Check total abundance per phylum
table(tax_table(p_nifH_2)[, "Phylum"], exclude = NULL)
moss_phyla_nifH <- tax_glom(p_nifH_2,taxrank = "Phylum")
taxa_sums(moss_phyla_nifH)
tax_table(moss_phyla_nifH)
write.table(taxa_sums(moss_phyla_nifH),"OTU_table_moss_phyla-nifH_abundance.txt", sep= "\t",  row.names = TRUE)

# Plot relative abundance (object my_phy2.ra) 
ps_nifH.ra <- transform_sample_counts(p_nifH_2, function(x){x / sum(x)})
ps_nifH.ra_table_phylum <- psmelt(ps_nifH.ra)

library("dplyr"); packageVersion("dplyr")

BarPlot_ps_nifH.ra_phylum <- ps_nifH.ra_table_phylum %>%
  ggplot(aes(x = Square, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position="stack") +
  labs(x = "Sample",
       y = "Relative Abundance",title="nifH gene") +
  facet_grid(~ Habitat, scales = "free", space='free') +
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
BarPlot_ps_nifH.ra_phylum

# Plot total abundance
# Abundances 
# Function to plot the abundance
# Do not modify this function, just run it and pass to the next argument. 
plot_abundance = function(physeq,title = "",
                          Facet = "Phylum", Color = "Phylum"){
  # Arbitrary subset, based on whatever you want, for plotting
  ps.f <- subset_taxa(physeq, Kingdom %in% c("Bacteria"))
  mphyseq <- psmelt(ps.f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Habitat",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")}


# Plot relative abundance (object my_phy2.ra) by taxonomy category (Phylum)
psOrd = subset_taxa(ps_nifH.ra, Kingdom %in% c("Bacteria"))
plot_abundance(psOrd, Facet = "Phylum", Color = 'Phylum')+
  labs(x = "Habitat",
       y = "Relative Abundance",
       title = "Bacteria phylum relative abundance for Racomitrim lanuginosum populations") + 
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Plot relative abundance by phylogenetic cluster
plot_bar(ps_nifH.ra,x="Square", fill="Phylo_clust_sub") + 
  geom_bar(aes(color = Phylo_clust_sub, fill = Phylo_clust_sub), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~Habitat, scales="free_x",nrow=1) +
  theme(panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

# Plot relative abundance by family
plot_bar(ps_nifH.ra,x="Square", fill="Family") + 
  geom_bar(aes(color = Family, fill = Family), stat="identity", position="stack") +
  labs(x = "Samples", y = "Relative Abundance\n") +
  facet_wrap(~Habitat, scales="free_x",nrow=1) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

table(tax_table(ps_nifH.ra)[, "Family"], exclude = NULL)


# Plot relative abundance by genus
plot_bar(ps_nifH.ra,x="Square", fill="Genus") + 
  geom_bar(aes(color = Genus, fill = Genus), stat="identity", position="stack") +
  labs(x = "Samples", y = "Relative Abundance\n") +
  facet_wrap(~Habitat, scales="free_x",nrow=1)

# Clustering Genera together to sum up the abundances. 
my_genus=tax_glom(ps_nifH.ra,taxrank = "Genus",NArm = FALSE)
write.table(my_genus@otu_table,"table_nifH_genus.txt", sep= "\t")
write.table(my_genus@tax_table,"tax_nifH_genus.txt", sep= "\t")
write.table(my_genus@sam_data,"data_nifH_genus.txt", sep= "\t")

plot_bar(my_genus, fill="Genus",) + 
  geom_bar(stat="identity", position="stack") +
  theme(panel.background = element_blank()) +
  facet_wrap(~Habitat, scales="free_x",nrow=1)

# Number of samples per habitat
nsamples(subset_samples(p_nifH_2, Habitat =="Forest tundra"))
nsamples(subset_samples(p_nifH_2, Habitat =="Shrub tundra"))

# Plot the 20 most abundant OTUs
nifH.top20 <- names(sort(taxa_sums(p_nifH_2), decreasing=TRUE))[1:20]
p_nifH_top <- prune_taxa(nifH.top20, p_nifH_2)
p_nifH_top <-transform_sample_counts(p_nifH_top, function(OTU) OTU/sum(OTU))

table(tax_table(p_nifH_top )[, "Family"], exclude = NULL)


plot_bar(p_nifH_top, x="Square", fill = "Family") + 
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") +
  labs(x = "Samples", y = "Relative Abundance\n") +
  facet_wrap(~Habitat, scales="free_x",nrow=1) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Clustering Genera together to sum up the abundances. 
my_genus2=tax_glom(p_nifH_top,taxrank = "Genus",NArm = FALSE)
table(tax_table(my_genus2)[, "Genus"], exclude = NULL)

write.table(my_genus2@otu_table,"table_nifH_genus2.txt", sep= "\t")
write.table(my_genus2@tax_table,"tax_nifH_genus2.txt", sep= "\t")
write.table(my_genus2@sam_data,"data_nifH_genus2.txt", sep= "\t")

plot_bar(my_genus2, fill="Genus",) + 
  geom_bar(stat="identity", position="stack") +
  theme(panel.background = element_blank()) +
  facet_wrap(~Habitat, scales="free_x",nrow=1)

library(knitr)

table(tax_table(p_nifH_top)[, "Family"], exclude = NULL)
table_nifH_top_20 <- psmelt(p_nifH_top)
write.table(table_nifH_top_20,"table_nifH_top_20.txt", sep= "\t")
kable(head(table_nifH_top_20))

# Heatmaps
plot_heatmap(ps_nifH.norm, taxa.label="Family",sample.label="Habitat",sample.order="Habitat",taxa.order="Family")
plot_heatmap(ps_nifH.norm, taxa.label="Family",sample.label="Habitat") # Ordination approach




#----------------------- 7. Non-normalized data - Alpha diversity   -------------------------------------#

theme_set(theme_bw())
# Alpha diversity

#plot_richness(my_phyloseq, x="Locality", measures=c("Shannon", "Simpson"), color="Habitat")
plot_richness(p_nifH_2, x="Habitat",color='Habitat', measures=c("Observed", "Shannon", "Simpson")) + geom_boxplot()


# Ordination analysis

# Transform data to proportions as appropriate for Bray-Curtis distances
ord.nmds.bray_nifH <- ordinate(ps_nifH.ra, method="NMDS", distance="bray")
plot_ordination(ps_nifH.ra, ord.nmds.bray_nifH, color="Habitat", title="Bray NMDS")
#                label="Collection")

#----------------------- 8. Normalize data using proportions  -------------------------------------#

# Normalization using proportions. The normalized data can be used for alpha diversity estimates.
library("vegan")
summary(sample_sums(p_nifH_2))
ps_nifH.ra
ps_nifH.norm <- transform_sample_counts(ps_nifH.ra, function(x) {x * 1000})
summary(sample_sums(ps_nifH.norm ))


#----------------------- 9. Alpha diversity comparison between habitats -------------------------------------#
library("vegan")
library(tidyr)
library(MASS)
library(dplyr)

# Diversity indexes
# !!! MARGIN is an option to transpose data. In this case our taxa are rows.
Shannon.nifH <- diversity(otu_table(ps_nifH.norm), MARGIN = 2, index = "shannon")
Simpson.nifH <- diversity(otu_table(ps_nifH.norm), MARGIN = 2, index = "simpson")
Pielou.nifH <- Shannon.nifH/log(specnumber(otu_table(ps_nifH.norm),MARGIN = 2))
variables_nifH<-sample_data(ps_nifH.norm)

# Merge results into a dataframe
Alpha_nifH <- data.frame(variables_nifH,Shannon.nifH,Simpson.nifH,Pielou.nifH)
names(Alpha_nifH)[9]<-"Shannon"
names(Alpha_nifH)[10]<-"Simpson"
names(Alpha_nifH)[11]<-"Pielou"

Alpha_nifH_means<- aggregate(Alpha.nifH, list(Alpha.nifH_F$Plot), FUN=mean)

# Comparison between habitats
Alpha.nifH_F<- Alpha_nifH [which(Alpha_nifH$Habitat=="Forest tundra"),]
Alpha.nifH_S<- Alpha_nifH [which(Alpha_nifH$Habitat=="Shrub tundra"),]

# Shannon means
Shannon_nifH_F <- aggregate(Alpha.nifH_F$Shannon, list(Alpha.nifH_F$Plot), FUN=mean)
Shannon_nifH_S <- aggregate(Alpha.nifH_S$Shannon, list(Alpha.nifH_S$Plot), FUN=mean)

truehist(Shannon_nifH_F[,2],prob=F)
truehist(Shannon_nifH_S[,2],prob=F)

# Normality tests
shapiro.test(Shannon_nifH_F[,2])
shapiro.test(Shannon_nifH_S[,2])

# Comparison using Mann-Whitney test
wilcox.test(Shannon_nifH_F[,2],Shannon_nifH_S[,2]) # Comparison using Mann-Whitney
# The Shannon index do not differ between habitats W = 59, p-value = 0.3445

# Simpson means
Simpson_nifH_F <- aggregate(Alpha.nifH_F$Simpson, list(Alpha.nifH_F$Plot), FUN=mean)
Simpson_nifH_S <- aggregate(Alpha.nifH_S$Simpson, list(Alpha.nifH_S$Plot), FUN=mean)

truehist(Simpson_nifH_F[,2],prob=F)
truehist(Simpson_nifH_S[,2],prob=F)

# Normality tests
shapiro.test(Simpson_nifH_F[,2])
shapiro.test(Simpson_nifH_S[,2]) # Not normally distributed

# Comparison using Mann-Whitney test
wilcox.test(Simpson_nifH_F[,2],Simpson_nifH_S[,2])# Comparison using Mann-Whitney
# The Simpson index do not differ between habitats W = 60, p-value = 0.373

# Pielou means
Pielou_nifH_F <- aggregate(Alpha.nifH_F$Pielou, list(Alpha.nifH_F$Plot), FUN=mean)
Pielou_nifH_S <- aggregate(Alpha.nifH_S$Pielou, list(Alpha.nifH_S$Plot), FUN=mean)

truehist(Pielou_nifH_F[,2],prob=F)
truehist(Pielou_nifH_S[,2],prob=F)

# Normality tests
shapiro.test(Pielou_nifH_F[,2]) # Normally distributed
shapiro.test(Pielou_nifH_S[,2]) # Not normally distributed

# Comparison between habitats using Mann-Whitney test
wilcox.test(Pielou_nifH_F[,2],Pielou_nifH_S[,2])# Comparison using Mann-Whitney
# The Pielou index do not differ between habitats W = 60, p-value = 0.373

#----------------------- 10. Phylogenetic tree   -------------------------------------#

# BiocManager::install("phangorn")
library(phangorn); packageVersion("phangorn")
library(Biostrings); packageVersion("Biostrings")
library("DECIPHER"); packageVersion("DECIPHER")

dna_nifH<- Biostrings::DNAStringSet(refseq(ps_nifH.norm))
names(dna_nifH) <- taxa_names(ps_nifH.norm)
aligned_nifH_seqs <- AlignSeqs(dna_nifH)
distances_nifH <- dist.ml(as.phyDat(as.DNAbin(aligned_nifH_seqs)))
nj_JC69_nifH <- NJ(distances_nifH)

# Tree inference
fit_nifH = pml(nj_JC69_nifH, data=as.phyDat(as.DNAbin(aligned_nifH_seqs)))
fitGTR_nifH <- optim.pml(fit_nifH, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
save(fitGTR_nifH, file="phylogenetic_tree_nifH.Rdata")
rooted_tree_nifH <- midpoint(fitGTR_nifH$tree)

#----------------------- 11. Faith's Phylogenetic diversity -------------------------------------#
#install.packages("picante")
library(picante); packageVersion("picante")
library(tidyr)
library(MASS)
library(ggplot2); packageVersion("ggplot2")

class(rooted_tree_nifH) # Check if the tree is rooted
getRoot(rooted_tree_nifH)
Faith1<- pd(t(otu_table(ps_nifH.norm)),rooted_tree_nifH,include.root=F)

Alpha_nifH <-data.frame(Alpha_nifH,Faith1)
write.table(Alpha_nifH,"Alpha_nifH.txt", sep= "\t",  row.names = TRUE)

# Comparison between habitats
Alpha.nifH_F<- Alpha_nifH [which(Alpha_nifH$Habitat=="Forest tundra"),]
Alpha.nifH_S<- Alpha_nifH [which(Alpha_nifH$Habitat=="Shrub tundra"),]

# Faith
Faith_nifH_F <- aggregate(Alpha.nifH_F$PD, list(Alpha.nifH_F$Plot), FUN=mean)
Faith_nifH_S <- aggregate(Alpha.nifH_S$PD, list(Alpha.nifH_S$Plot), FUN=mean)

truehist(Faith_nifH_F[,2],prob=F)
truehist(Faith_nifH_S[,2],prob=F)

# Normality tests
shapiro.test(Faith_nifH_F[,2])
shapiro.test(Faith_nifH_S[,2])

# Comparison using Mann-Whitney test
wilcox.test(Faith_nifH_F[,2],Faith_nifH_S[,2]) # Comparison using Mann-Whitney
# The Faith index do not differ between habitats W = 68, p-value = 0.6475

# Transform data for plotting
names(Alpha_nifH)[12]<-"Faith"
names(Alpha_nifH)[13]<-"Richness"

Alpha.nifH.tall<- Alpha_nifH %>% gather (key=Diversity, value=Value,Shannon:Faith)
Alpha.nifH.tall$Diversity<-as.factor(Alpha.nifH.tall$Diversity)
levels(Alpha.nifH.tall$Diversity)
Alpha.nifH.tall$Diversity<-factor(Alpha.nifH.tall$Diversity, levels=c("Shannon","Simpson","Pielou", "Faith"))

theme_set(theme_bw())
Alpha_nifH_plot<-ggplot(Alpha.nifH.tall,aes(x=Habitat,y=Value)) +
  geom_boxplot(aes(fill=Habitat)) +
  facet_wrap(~Diversity,scale="free")+
  labs(x = "nifH",
       y = "Alpha diversity measure")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

# To use vertical x-axis titles  guides(x = guide_axis(angle = 90))


# Plotting both 16S and nifH diversities estimates exporting a table
Alpha_div_mean_all<- read.delim("Alpha_div_mean_all.txt",sep ="")

Alpha.tall<- Alpha_div_mean_all %>% gather (key=Diversity, value=Value,Shannon:Faith)
Alpha.tall$Diversity<-as.factor(Alpha.tall$Diversity)
levels(Alpha.tall$Diversity)
Alpha.tall$Diversity<-factor(Alpha.tall$Diversity, levels=c("Shannon","Simpson","Pielou", "Faith"))

library("ggplot2")

theme_set(theme_bw())
Alpha_plot<-ggplot(Alpha.tall,aes(x=Gene,y=Value)) +
  geom_boxplot(aes(fill=Habitat)) +
  facet_wrap(~Diversity,scale="free")+
  labs(x = "",
       y = "Alpha diversity measure")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Alpha_plot
#----------------------- 11.2 Optional - Abundance in Cyanobacteria  -------------------------------------#
rooted_tree_nifH
ps_nifH.norm=merge_phyloseq(ps_nifH.norm, rooted_tree_nifH)
ps_nifH_Cyano <- subset_taxa(ps_nifH.norm, Phylum=="Cyanobacteria")
plot_tree(ps_nifH_Cyano, color="Habitat", shape="Family", label.tips="Genus", size="Abundance")


#----------------------- 12. GUnifrac distances  -------------------------------------#


#install.packages("GUniFrac")
library("GUniFrac")

#Create the GUnifrac distance matrix
unifracs_nifH<-GUniFrac(t(otu_table(ps_nifH.norm)),rooted_tree_nifH, alpha = 0.5)
Gunifrac_dist_nifH <- unifracs_nifH$unifracs[, , "d_0.5"]
class(Gunifrac_dist_nifH)
Gunifrac_dist_nifH<-as.dist(Gunifrac_dist_nifH)
class(Gunifrac_dist_nifH)

#----------------------- 13. Ordinations and PERMANOVA -------------------------------------#

library("vegan"); packageVersion("vegan")

theme_set(theme_bw(),theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#Example
# NMDS using Bray-Curtis
ord.nmds.nifH <- ordinate(ps_nifH.norm, method="NMDS", distance = "bray")
NMDS.nifH<- plot_ordination(ps_nifH.norm, ord.nmds.nifH, color="Habitat", title="moss nifH - Bray-Curtis NMDS")
NMDS.nifH + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw()

# NMDS using GUnifrac
ord.nmds.nifH <- ordinate(ps_nifH.norm, method="NMDS", distance = Gunifrac_dist_nifH)
NMDS.nifH<- plot_ordination(ps_nifH.norm, ord.nmds.nifH, color="Habitat", title="nifH - GUnifrac NMDS")
NMDS.nifH + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw()

# PCoA using GUnifrac
ord.pcoa.nifH <- ordinate(ps_nifH.norm, method="PCoA", distance = Gunifrac_dist_nifH, correction="lingoes")
# Correct the negative values !!!
plot_ordination(ps_nifH.norm, ord.pcoa.nifH, color="Habitat", title="nifH - GUnifrac PCoA")

ord.pcoa.nifH <- pcoa(Gunifrac_dist_nifH, correction="lingoes",row.names(Alpha_nifH))
PCoA_nifH <- plot_ordination(ps_nifH.norm,ord.pcoa.nifH  , color = "Habitat", title = "nifH - GUnifrac PCoA")
PCoA_nifH + 
  stat_ellipse(type = "norm", linetype = 2) +
  labs(x="Axis 1 (17.95%)",y="Axis 2 (11.01%)") +
  theme_bw()

# Plot both NMDS in the same plot

library("ggpubr")
ggarrange(NMDS.nifH + stat_ellipse(type = "norm", linetype = 2), PCoA_nifH + 
            stat_ellipse(type = "norm", linetype = 2) +
            labs(x="Axis 1 (17.95%)",y="Axis 2 (11.01%)"), 
          labels = c("A)", "B)"),
          ncol = 1, nrow = 2)

theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# PERMANOVA 

# Transform Plots to factors.
class(Alpha_nifH$Plot)
Alpha_nifH$Plot<-as.factor(Alpha_nifH$Plot)
levels(Alpha_nifH$Plot)

# Using adonis function
permanova_nifH0 <- adonis(Gunifrac_dist_nifH ~ Habitat, data = Alpha_nifH, permutations = 9999)
permanova_nifH0

permanova_nifH1 <- adonis(Gunifrac_dist_nifH ~ Habitat, strata = Alpha_nifH$Plot, data = Alpha_nifH, permutations = 9999)
permanova_nifH1

# Using adonis2 function with the perm and how command.
# Null hypothesis: Groups have the same centroid.
perm <- how(within = Within(type = "free"), plots = Plots(strata=Alpha_nifH$Plot, type = "free"),nperm = 9999)
#perm <- how(nperm = 10000)

permanova_nifH<- adonis2(Gunifrac_dist_nifH ~ Alpha_nifH$Habitat, permutations = perm)
permanova_nifH


# Check the dispersion in both habitats

beta_nifH <- betadisper(Gunifrac_dist_nifH,Alpha_nifH$Habitat)
permutest(beta_nifH)
# The test indicate that the there is not enough proofs 
# to confirm that the homogeneity of multivariate dispersion between habitats differ


#----------------------- 14. LEfSe -------------------------------------#

write.table(otu_table(ps_nifH.norm),"OTU_table_nifH_norm.txt", sep= "\t",  row.names = TRUE)
write.table(sample_data(ps_nifH.norm),file = "Data_nifH_norm.txt",sep= "\t",  row.names = TRUE)
write.table(otu_table(ps_nifH.ra),file = "OTU_table_nifH_ra.txt.txt",sep= "\t",  row.names = TRUE)
write.table(tax_table(ps_nifH.norm),file = "Tax_table_nifH_norm.txt",sep= "\t",  row.names = TRUE)

# UPDATE FORMATTING 08-06-2023
# Change according to nifH data (The script have the 16s objects)

library("data.table")
library("phyloseq")
library("ggplot2")
library("grid")
library("tidyverse")

# install.packages("vctrs")
# install.packages("cli")
# install.packages("tidyverse")

# You must have the script "lefse.R" to run the following analyses

source("lefse.R", local = TRUE)
#install.packages("magrittr") 
# tidyverse_update(recursive = FALSE, repos = getOption("repos"))
#install.packages("remotes")
# remotes::install_github("mikemc/speedyseq")
library(speedyseq)

taxa_are_rows(ps_nifH.norm)

tax2 = lefse_1name_obj(ps_nifH.norm, sample_data(ps_nifH.norm)$Habitat,subject = T)
lefse2 = lefse_obj(ps_nifH.norm)
lefse2 = rbind(tax2, lefse2)

# Replace unsupported chars with underscore
lefse2$name = gsub(" ","_",lefse2$name)
lefse2$name = gsub("-","_",lefse2$name)
lefse2$name = gsub("/","_",lefse2$name)

# Output the prepared LEfSe input to file
write.table(lefse2, file="lefse_table_nifH_UPDATED.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#----------------------- 15. Racomitrium lanuginosum core diazotrophic community -------------------------------------#

# Venn Diagram

#install.packages("eulerr") # If not installed
library(eulerr); packageVersion("eulerr") 
library(microbiome)
#devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities);packageVersion("microbiomeutilities")

# Transform to relative abundances using microbiome package
phylo_nifH.microbiome <- microbiome::transform(p_nifH_2, "compositional")

# Select the type of samples
Habitat_type <- unique(as.character(meta(phylo_nifH.microbiome)$Habitat))
list_nifH <- c() # an empty object to store information

# Write a for loop to go through each of the sample_type one by one and
# combine identified core taxa into a list.

for (n in Habitat_type){ # for each variable n in sample_type
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(phylo_nifH.microbiome, Habitat_type == n) # Choose sample from Sample_type by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from soil 
                         detection = 0.001, # 0.001 in at least 90% samples 
                         prevalence = 0.50)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each Sample_type.
  list_nifH[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}


mycols2 <- c(Forest_tundra="#ff9999", Shrub_tundra="#99ccff") 
plot(venn(list_nifH),
     fills = mycols2)
length(core_m) # Number of total core OTUs for Racomitrium lanuginosum

# Prevalence 0.75
# Forest tundra: "JC2961-nifH.687910".
# 

# Shrub tundra: "JC2803-nifH.183622" "JC2961-nifH.687910"
# 

# Prevalence 0.50
# Forest tundra: ""PCR49-A-nifH.921623"(Dolichospermum), "JC2944-nifH.475063(Trichormus)" JC2792-nifH.65992"(Mastigocladus),  "JC2788-nifH.4349 (Fischerella)" "JC3086-nifH.740049 (Nostoc)" 
# 

# Shrub tundra: ""JC2789-nifH.15630 (Anabaena)"  "JC3213-nifH.841527 (Nostoc)" "JC2958-nifH.647079 (Dolichospermum)"



#----------------------- 16. Racomitrium lanuginosum core diazotrophic community -------------------------------------#

# BiocManager::install("microbiome")
library(microbiome); packageVersion("microbiome") 
library(phyloseq)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(viridis)

phylo_nifH.microbiome <- prune_taxa(taxa_sums(phylo_nifH.microbiome) > 0, phylo_nifH.microbiome)
taxa_names(phylo_nifH.microbiome)[1:3]

# Analysis of core microbiome using a threshold of 0.75
core_75 <- core_members(phylo_nifH.microbiome, detection = 0.001, prevalence = 75/100, include.lowest = FALSE)
taxonomy <- as.data.frame(tax_table(phylo_nifH.microbiome))
# Subset the core microbiome
core_taxa_id75 <- subset(taxonomy, rownames(taxonomy) %in% core_75)
core.abundance.75 <- sample_sums(core(phylo_nifH.microbiome, detection = 0.001, prevalence = 75/100))
# core diazotrophic community at 0.75 = JC2961-nifH.687910 Nostoc

# Analysis of core microbiome using a threshold of 0.50
core_50 <- core_members(phylo_nifH.microbiome, detection = 0.001, prevalence = 50/100, include.lowest = FALSE)
#taxonomy <- as.data.frame(tax_table(phylo_nifH.microbiome))
# Subset the core microbiome
core_taxa_id50 <- subset(taxonomy, rownames(taxonomy) %in% core_50)
core.abundance.50 <- sample_sums(core(phylo_nifH.microbiome, detection = 0.001, prevalence = 50/100))

# Plotting results
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

plot_core(phylo_nifH.microbiome, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()

## Core microbiome with compositional (prevalence 75 -> min.prevalence = .75):
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)
detections <- round(detections,3)
## Also define gray color palette
gray <- gray(seq(0,1,length=5))
p1 <- plot_core(phylo_nifH.microbiome, 
                plot.type = "heatmap", 
                colours = gray,
                prevalences = prevalences, 
                detections = detections, min.prevalence = .50) +
  xlab("Detection Threshold (Relative Abundance (%))")

p1 <- p1 + theme_bw() + ylab("OTUs")
p1


print(p1 + scale_fill_viridis())


# Plot with taxonomic information
library(RColorBrewer)
library(knitr)

# Extract taxonomic information
df <- p1$data 
list <- df$Taxa 
tax <- as.data.frame(tax_table(phylo_nifH.microbiome))
tax$ASV <- rownames(tax)
tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 
tax.unit <- tidyr::unite(tax2, Taxa_level,c("Phylum", "Class", "Order", "Family", "Genus", "ASV"), sep = "_;", remove = TRUE)
tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)
df$Taxa <- tax.unit$Taxa_level
knitr::kable(head(df))
p1$data <- df

# Plot again with the updated taxonomic information
print(p1 + scale_fill_viridis())



