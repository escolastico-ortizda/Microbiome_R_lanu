###-------------------- Bacterial community analyses of Racomitrium lanuginosum 16S amplicons (moss analyses)--------------------###

# Code for moss 16S analyses. Dennis Escolástico-Ortiz. 2022

# install.packages("BiocManager")
# install.packages("devtools")
# BiocManager::install("phyloseq", force = TRUE )
# BiocManager::install("decontam")
# BiocManager::install("DECIPHER")
# BiocManager::install("dada2")

# Install DADA2

# Load packages devtools and DADA2
library("devtools")
library("dada2")
packageVersion("dada2")

##-------------------- 1. Analyses of 16S sequences from R. lanuginosum from the first sequencing run --------------------##

# Load the amplicons sequences
getwd()
Path16S_run1 <-"~/Microbiome_bioinformatic_analyses_working_directotry/16S_moss_data/16S_1run_moss_fq_files"
list.files(Path16S_run1)
##-------------------- 1.1 Filter and trim --------------------##

# Forward and reverse fastq filenames have the following format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs1 <- sort(list.files(Path16S_run1, pattern="_R1_001.fastq", full.names = TRUE))
fnRs1 <- sort(list.files(Path16S_run1, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs1), "_"), `[`, 1)

# Inspect read quality profiles
#plotQualityProfile(fnFs1[29:32]) #forward.Cut in 270,200,290,290,280,280,280,100,150,290,100,290,290,290,290,250,(JC2954)
#plotQualityProfile(fnRs1[29:32]) #reverse.Cut in 250,250,250,250,250,240,240,240(BadQ: JC2789,2931,2932,2934,2935,2954)

# Filter and trim. Place filtered files in filtered/ subdirectory.
filtFs1 <- file.path(Path16S_run1, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs1 <- file.path(Path16S_run1, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs1) <- sample.names
names(filtRs1) <- sample.names


# Trimming
out1 <- filterAndTrim(fnFs1, filtFs1, fnRs1, filtRs1, truncLen=c(280,230),
                     maxN=0, maxEE=c(2,2), truncQ=2, trimLeft=c(10,10), rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

head(out1)

# Re-check quality after trimming
#plotQualityProfile(filtFs1[1:12])
#plotQualityProfile(filtRs1[1:12])

##-------------------- 1.2 Infer amplicon sequence variants --------------------##

errF1 <- learnErrors(filtFs1, multithread=TRUE)
errR1 <- learnErrors(filtRs1, multithread=TRUE) # Two hours minimum.
#plotErrors(errF1, nominalQ=TRUE) # Report errors
#plotErrors(errR1, nominalQ=TRUE) # Report errors

# Sample inference
dadaFs1 <- dada(filtFs1, err=errF1, multithread=TRUE, pool="pseudo")
dadaRs1 <- dada(filtRs1, err=errR1, multithread=TRUE, pool="pseudo")
dadaFs1[[1]]
dadaRs1[[1]]

# Merge paired reads
mergers1 <- mergePairs(dadaFs1, filtFs1, dadaRs1, filtRs1, verbose=TRUE)
head(mergers1[[1]])

# Construct sequence table
seqtab1 <- makeSequenceTable(mergers1)
dim(seqtab1)

# Inspect distribution of sequence lengths 
table(nchar(getSequences(seqtab1))) 
# Check if length match the expected amplicon length.

# Removing of sequence with lower expected length (<400bp) or that did not merged well 
seqtab1b <- seqtab1[,nchar(colnames(seqtab1)) %in% seq(400,450)] #There is xx sequences with a length of xx bp
table(nchar(getSequences(seqtab1b)))


##-------------------- 2. Analyses of 16S sequences from R. lanuginosum from the second sequencing run --------------------##

# Load the amplicons sequences
getwd()
Path16S_run2 <-"~/Microbiome_bioinformatic_analyses_working_directotry/16S_moss_data/16S_2run_moss_fq_files"
list.files(Path16S_run2)

##-------------------- 2.1 Filter and trim --------------------##

# Forward and reverse fastq filenames have the following format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs2 <- sort(list.files(Path16S_run2, pattern="_R1_001.fastq", full.names = TRUE))
fnRs2 <- sort(list.files(Path16S_run2, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs2), "_"), `[`, 1)

# Inspect read quality profiles
#plotQualityProfile(fnFs2[50:56]) #forward.Cut in 290,280,290,280,280,280,290
#plotQualityProfile(fnRs2[41:50]) #reverse.Cut in 230,240,240,240,230

# Filter and trim. Place filtered files in filtered/ subdirectory.
filtFs2 <- file.path(Path16S_run2, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs2 <- file.path(Path16S_run2, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs2) <- sample.names
names(filtRs2) <- sample.names

# Trimming
out2 <- filterAndTrim(fnFs2, filtFs2, fnRs2, filtRs2, truncLen=c(280,230),
                     maxN=0, maxEE=c(2,2), truncQ=2, trimLeft=c(10,10), rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

head(out2)

# Re-check quality after trimming
#plotQualityProfile(filtFs2[1:4])
#plotQualityProfile(filtRs2[1:4])

##-------------------- 2.2 Infer amplicon sequence variants --------------------##

errF2 <- learnErrors(filtFs2, multithread=TRUE)
errR2 <- learnErrors(filtRs2, multithread=TRUE) # Two hours minimum.
plotErrors(errF2, nominalQ=TRUE) # Report errors
plotErrors(errR2, nominalQ=TRUE) # Report errors

# Sample inference
dadaFs2 <- dada(filtFs2, err=errF2, multithread=TRUE, pool="pseudo")
dadaRs2 <- dada(filtRs2, err=errR2, multithread=TRUE, pool="pseudo")
dadaFs2[[1]]
dadaRs2[[1]]

# Merge paired reads
mergers2 <- mergePairs(dadaFs2, filtFs2, dadaRs2, filtRs2, verbose=TRUE)
head(mergers2[[1]])

# Construct sequence table
seqtab2 <- makeSequenceTable(mergers2)
dim(seqtab2)

# Inspect distribution of sequence lengths 
table(nchar(getSequences(seqtab2))) 
# Check if length match the expected amplicon length.

# Removing of sequence with lower expected length (<400bp) or that did not merged well 
seqtab2b <- seqtab2[,nchar(colnames(seqtab2)) %in% seq(400,465)] #There is xx sequences with a length of xx bp
table(nchar(getSequences(seqtab2b)))



##-------------------- 3.1 Remove chimeras and merge sequences tables --------------------##

# Sequence table of moss 1st, 2nd run, and SOIL run.

# This data was saved from independent data processing.
load("./soil_16S_analyses.RData")
load("./16_data_analyses.RData")
dim(seqtab1b)
dim(seqtab2b)
dim(seqtab_soilb)

# Merge tabs from moss 16S first and second run
alltab0 <-mergeSequenceTables(seqtab1b,seqtab2b)
dim(alltab0)

# Merge tabs from moss 16S and soil
alltab <- mergeSequenceTables(alltab0,seqtab_soilb)
dim(alltab)

# Remove chimeras in all samples
seqtab.nochim.16all <- removeBimeraDenovo(alltab, method="consensus", multithread=TRUE, verbose=TRUE)
save(seqtab.nochim.16all,file="seqtab.nochim.16all.RData")
dim(seqtab.nochim.16all)
sum(seqtab.nochim.16all)/sum(alltab) # Check if this is the good object for comparison.

##-------------------- 3.2 Assign taxonomy for all samples --------------------##

set.seed(100) # Random number generator
taxa <- assignTaxonomy(seqtab.nochim.16all, "./silva_nr_v132_train_set.fa.gz", multithread=TRUE, minBoot = 80)
taxa_all<- taxa
dim(taxa_all) # 17402     6
unname(taxa_all)

#----------------------- 4. Handoff to phyloseq -------------------------------------#
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())
library("DECIPHER")
library("phangorn")

# Construct a phyloseq object with metadata
getwd()
#setwd("REPLACE_your_directory")
data_16S_all <- read.delim("./Microbiome_bioinformatic_analyses/Script_Bacterial_community/16S_all_data.txt", sep = "\t", header = TRUE, row.names = 1)
str(data_16S_all)
data_16S_all$Plot <-as.factor(data_16S_all$Plot)
str(data_16S_all)

phylo_16S_all <- phyloseq(otu_table(seqtab.nochim.16all, taxa_are_rows=FALSE),
                        sample_data(data_16S_all), 
                        tax_table(taxa_all))
dna_16S_all <- Biostrings::DNAStringSet(taxa_names(phylo_16S_all))
names(dna_16S_all) <- taxa_names(phylo_16S_all)
phylo_16S_all <- merge_phyloseq(phylo_16S_all, dna_16S_all)
taxa_names(phylo_16S_all) <- paste0("ASV", seq(ntaxa(phylo_16S_all)))

phylo_16S_all #17,402 ASVs
head(sample_data(phylo_16S_all), 20)
rank_names(phylo_16S_all)

# Exclude chloroplast and mitochondria
phylo_16S_all <- subset_taxa(phylo_16S_all, !is.na(Phylum) & !Order %in% c("Chloroplast") 
                          & !Family %in% c("Mitochondria"))
phylo_16S_all # 15,784 ASVs

table(tax_table(phylo_16S_all)[, "Phylum"], exclude = NULL)
phylo_16S_all <- subset_taxa(phylo_16S_all, !Phylum %in% c("Arthropoda","BRC1", "Ciliophora", "Deinococcus-Thermus","Entotheonellaeota","FCPU426","Latescibacteria", "Omnitrophicaeota","Phragmoplastophyta","Thaumarchaeota", "WS2")) #Less than 10
phylo_16S_all # 15,900 113 samples

#----------------------------- 5. Decontam - Remove contaminants -----------------------------------#

#BiocManager::install("decontam")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(dplyr)

# It is important to keep all the samples for these analyses (even low-read samples).
head(sample_data(phylo_16S_all))

# Inspect library sizes
df_16S_all <- as.data.frame(sample_data(phylo_16S_all)) 
df_16S_all$LibrarySize <- sample_sums(phylo_16S_all)
df_16S_all <- df_16S_all[order(df_16S_all$LibrarySize),]
df_16S_all # 113 samples

# Samples with no taxa and low read samples
# Use the df_mps object to select the true samples with low number of reads <50 reads.
Low_Lsize <- df_16S_all[which(df_16S_all$Sample_status =="TrueSample" & df_16S_all$LibrarySize < 50),]
Low_Lsize <- row.names.data.frame(Low_Lsize)
phylo_16S_all <- prune_samples(!(sample_names(phylo_16S_all) %in% Low_Lsize ), phylo_16S_all)
phylo_16S_all <- prune_samples(sample_sums(phylo_16S_all)>0,phylo_16S_all) # Remove Negative controls with 0 reads

phylo_16S_all # 15,750 ASVs & 97 samples

df_16S_all <- as.data.frame(sample_data(phylo_16S_all)) 
df_16S_all$LibrarySize <- sample_sums(phylo_16S_all)
df_16S_all <- df_16S_all[order(df_16S_all$LibrarySize),]
df_16S_all
df_16S_all$Index <- seq(nrow(df_16S_all))
df_16S_all

#Plotting library size per Sample status
ggplot(data=df_16S_all, aes(x=Index, y=LibrarySize, color=Sample_status)) + geom_point() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Identify Contaminants - Frequency
contamdf.freq <- isContaminant(phylo_16S_all, method="frequency", conc="DNA_conc")
# Samples with 0 counts will be removed for this analysis (7 samples in my case)
# As we did not specify the threshold, the default value of threshold = 0.1 was used, and $contaminant=TRUE if $p < 0.1.
head(contamdf.freq)

# This calculation has returned a data.frame with several columns, 
# the most important being $p which contains the probability that was used for 
# classifying contaminants, and $contaminant which contains TRUE/FALSE classification values
# with TRUE indicating that the statistical evidence that the associated sequence feature is
# a contaminant exceeds the user-settable threshold. 
table(contamdf.freq$contaminant)
which(contamdf.freq$contaminant)

plot_frequency(phylo_16S_all, taxa_names(phylo_16S_all)[c(1550,1599)], conc="DNA_conc") + 
  xlab("DNA Concentration (Nanodrop DNA concentrations)")

set.seed(100)
plot_frequency(phylo_16S_all, taxa_names(phylo_16S_all)[sample(which(contamdf.freq$contaminant),23)], conc="DNA_conc") + 
  xlab("DNA Concentration (Nanodrop DNA concentrations)")

# Remove contaminants
phylo_16S_all # 15,750
mps.noncontam <- prune_taxa(!contamdf.freq$contaminant, phylo_16S_all)
mps.noncontam # 15,713 ASVs 97 samples

# Identify Contaminants - Prevalence
sample_data(mps.noncontam)$is.neg <- sample_data(phylo_16S_all)$Sample_status == "ControlSample"
contamdf.prev <- isContaminant(mps.noncontam, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant)) # 3 contaminants

# Note that as before, the default threshold for a contaminant is that it reaches a 
# probability of 0.1 in the statistical test being performed. In the prevalence test 
# there is a special value worth knowing, threshold=0.5, that will identify as contaminants 
# all sequences that are more prevalent in negative controls than in positive samples.
# Let's try using this more aggressive classification threshold rather than the default.

contamdf.prev05 <- isContaminant(mps.noncontam, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant) # 37 contaminants
head(which(contamdf.prev05$contaminant))

mps.pa <- transform_sample_counts(mps.noncontam, function(abund) 1*(abund>0))
mps.pa.neg <- prune_samples(sample_data(mps.pa)$Sample_status == "ControlSample", mps.pa)
mps.pa.pos <- prune_samples(sample_data(mps.pa)$Sample_status == "TrueSample", mps.pa)

# Make data.frame of prevalence in positive and negative samples
df.mpa <- data.frame(pa.pos=taxa_sums(mps.pa.pos), pa.neg=taxa_sums(mps.pa.neg),
                     contaminant=contamdf.prev$contaminant)
ggplot(data=df.mpa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Remove ASVs according to prevalence
mps.noncontam # 15,713 ASVs
mps.noncontam2 <- prune_taxa(!contamdf.prev05$contaminant, mps.noncontam)
mps.noncontam2 # 15,676 ASVs

# Remove contaminants at the same time resulted from independent frequency and prevalence tests.
#contamdf.comb<- isContaminant(phylo_16S_all, conc="DNA_conc", method="combined", neg="is.neg", threshold=0.5)
#table(contamdf.comb$contaminant)
#head(which(contamdf.comb$contaminant))

##-----------------------  6. Filter according to prevalence  ----------------------- ##
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())


# Before proceed, filter out the negative controls
mps.noncontam2 # 15,676 ASVs 97 samples
mps.noncontam2@sam_data
mps.filtered <- prune_samples(sample_data(mps.noncontam2)$Sample_status == "TrueSample",mps.noncontam2)
mps.filtered@sam_data # 15,676 ASVs 89 samples


#Phyloseq object for moss samples
phylo_moss <-subset_samples(mps.filtered,Sample_type == "Moss")
phylo_moss 
table(tax_table(phylo_moss)[, "Phylum"], exclude = NULL)
phylo_moss <-prune_taxa(taxa_sums(phylo_moss) >=1,phylo_moss)
phylo_moss

# Compute prevalence of each feature, store as data.frame
prevdf2 = apply(X = otu_table(phylo_moss),
               MARGIN = ifelse(taxa_are_rows(phylo_moss), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf2 = data.frame(Prevalence = prevdf2,
                    TotalAbundance = taxa_sums(phylo_moss),
                    tax_table(phylo_moss))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),
                                                  sum(df1$Prevalence))})
# Prevalence: number of samples in which a taxon appears at least once.
# Subset to the remaining phyla
prevdf2.0 = subset(prevdf2, Phylum %in% get_taxa_unique(phylo_moss, "Phylum"))
ggplot(prevdf2.0, aes(TotalAbundance, Prevalence / nsamples(phylo_moss),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

  ## Filter out phylum with prevalence = 1.00000
## Define prevalence threshold to keep ASV present in at least 3% of samples (Callahan et al. web):
prevalenceThreshold2 = 0.05 * nsamples(phylo_moss)
prevalenceThreshold2 # 3.1 samples threshold for soil samples

## Execute prevalence filter, using `prune_taxa()` function
keepTaxa2 = rownames(prevdf2.0)[(prevdf2.0$Prevalence >= prevalenceThreshold2)]
phylo_moss2 = prune_taxa(keepTaxa2,phylo_moss)
table(tax_table(phylo_moss2)[, "Phylum"], exclude = NULL)
phylo_moss2 # 1,322 ASVs 62 samples

# Check summary per habitat (number of samples, ASV, reads)
table(sample_data(phylo_moss2)$Habitat)
# Forest tundra= 36 samples; Shrub tundra = 26 samples
mean(sample_sums(phylo_moss2))
summary(sample_sums(phylo_moss2))
# Mean number of reads of 7472

#----------------------- 7. Relative abundance  -------------------------------------#

library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

# Check total abundance per phylum (total relative abundance)
table(tax_table(phylo_moss2)[, "Phylum"], exclude = NULL)
moss_phyla <- tax_glom(phylo_moss2,taxrank = "Phylum")
taxa_sums(moss_phyla)
#write.table(taxa_sums(moss_phyla),"OTU_table_moss_phyla-abundance.txt", sep= "\t",  row.names = TRUE)
moss_phyla_ra <- transform_sample_counts(moss_phyla, function(x){x / sum(x)})
taxa_sums(moss_phyla_ra)
tax_table(moss_phyla_ra)


# Transform to relative abundance
p_moss2.ra <- transform_sample_counts(phylo_moss2, function(x){x / sum(x)})
taxa_sums(p_moss2.ra)

# To remove the undescore on Forest_tundra and Shrub_tundra
p_moss2.ra@sam_data
#write.table(sample_data(p_moss2.ra),"sam_data_p_moss2.ra.txt", sep= "\t",  row.names = TRUE)
#data_p_moss2_ra <- read.delim("./sam_data_p_moss2.ra.txt", sep = "\t", header = TRUE, row.names = 1)
#data_p_moss2_ra$Plot <-as.factor(data_p_moss2_ra$Plot)
#str(data_p_moss2_ra)
#sample_data(p_moss2.ra) <- data_p_moss2_ra

# Function to plot the abundance
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
psOrd = subset_taxa(p_moss2.ra , Kingdom %in% c("Bacteria"))
plot_abundance(psOrd, Facet = "Phylum", Color = 'Phylum')+
  labs(x = "Habitat",
       y = "Relative Abundance",
       title = "") + 
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# Plot relative abundance (object ps2.ra) using BarPlots and Boxplots
pmoss2.ra_table_phylum <- psmelt(p_moss2.ra)

# BarPlot 
# install.packages("dplyr")
library("dplyr"); packageVersion("dplyr")

# This the graphic for the article
BarPlot_pmosss2ra_phylum <- pmoss2.ra_table_phylum %>%
  ggplot(aes(x = Collection, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position="stack") +
  labs(x = "Sample",
       y = "Relative Abundance",title= "16S rRNA gene") +
  facet_grid(~ Habitat, scales = "free", space='free') +
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )
BarPlot_pmosss2ra_phylum

### If you want to plot both 16S and nifH abundance, load both datasets and run this command:
library("ggpubr")
ggarrange(BarPlot_pmosss2ra_phylum, BarPlot_ps_nifH.ra_phylum, 
        labels = c("a", "b"),
        ncol = 1, nrow = 2)
# Save as high resolution figures

### IMPORTANT!!! To save the barplots in tiff format use the "tiff" function. Otherwise, the plots will be save with weird lines.

tiff('Figure_2_Moss_ASV_OTU_Phylum_relative_abundance_174_210_300dpi.tiff', units="mm", width=174, height=210, res=300, compression = 'lzw')
# Run the command for the figure
dev.off()

#ggsave("Figure_2_Moss_ASV_OTU_Phylum_relative_abundance_174_210_300dpi.svg",units="mm", width=174, height=210, dpi=300)

# Barplot per phylum in two habitats
theme_set(theme_bw()) # Apply this before plotting, otherwise lines will not be shown.
plot_bar(p_moss2.ra, x="Collection", fill="Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(panel.background = element_blank()) +
  facet_wrap(~Habitat, scales="free_x",nrow=1) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


  # Barplot per order in two habitats
plot_bar(p_moss2.ra, x="Collection", fill="Order") + 
  geom_bar(aes(color = Order, fill = Order), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(panel.background = element_blank()) +
  facet_wrap(~Habitat, scales="free_x",nrow=1) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

  # Bar Plot - Top 50 ASVs
top_moss <- names(sort(taxa_sums(phylo_moss2), decreasing=TRUE))[1:50]

phylo_moss2.top2 <- prune_taxa(top_moss, phylo_moss2)
phylo_moss2.top2 <- transform_sample_counts(phylo_moss2.top2, function(OTU) OTU/sum(OTU))
table(tax_table(phylo_moss2.top2 )[, "Family"], exclude = NULL)

plot_bar(phylo_moss2.top2,x="Collection", fill="Family") + 
  geom_bar(aes(color = Family, fill = Family), stat="identity", position="stack") +
  labs(x = "Samples", y = "Relative Abundance\n") +
  facet_wrap(~Habitat, scales="free_x",nrow=1) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#----------------------- 8. Normalize data using proportions  -------------------------------------#

# Normalization using proportions. The normalized data can be used for alpha diversity estimates.
library("vegan")

summary(sample_sums(phylo_moss2))
summary(sample_sums(p_moss2.ra))
p_moss2.ra
ps_moss.norm <- transform_sample_counts(p_moss2.ra, function(x) {x * 1000})
summary(sample_sums(ps_moss.norm))

#----------------------- 9. Alpha diversity of moss between habitats -------------------------------------#

library("vegan")
library(tidyr)
library(MASS)
library(dplyr)

# Diversity indexes
Shannon.moss <- diversity(otu_table(ps_moss.norm), index = "shannon")
Simpson.moss <- diversity(otu_table(ps_moss.norm), index = "simpson")
Pielou.moss <- Shannon.moss/log(specnumber(otu_table(ps_moss.norm)))
variables_moss <- sample_data(ps_moss.norm)

# Merge into one dataframe
Alpha_moss <- data.frame(variables_moss,Shannon.moss,Simpson.moss,Pielou.moss)
names(Alpha_moss)[9]<-"Shannon"
names(Alpha_moss)[10]<-"Simpson"
names(Alpha_moss)[11]<-"Pielou"
Alpha_moss$Plot <- as.factor(Alpha_moss$Plot)

# Datsets per habitat
Alpha_moss_F<- Alpha_moss[which(Alpha_moss$Habitat=="Forest_tundra"),]
Alpha_moss_S<- Alpha_moss[which(Alpha_moss$Habitat=="Shrub_tundra"),]

## SHANNON ##
# Shannon means
Shannon_moss_F <- aggregate(Alpha_moss_F$Shannon, list(Alpha_moss_F$Plot), FUN=mean)
Shannon_moss_S <- aggregate(Alpha_moss_S$Shannon, list(Alpha_moss_S$Plot), FUN=mean)

truehist(Shannon_moss_F[,2],prob=F)
truehist(Shannon_moss_S[,2],prob=F)

# Normality tests
shapiro.test(Shannon_moss_F[,2]) # Not normal
shapiro.test(Shannon_moss_S[,2]) # Not normal

# Comparison using Mann-Whitney test
wilcox.test(Shannon_moss_F[,2],Shannon_moss_S[,2]) # Comparison using Mann-Whitney
# The Shannon index do not differ between habitats W = 98, p-value = 0.1083

## SIMPSON ##
# Simpson means
Simpson_moss_F <- aggregate(Alpha_moss_F$Simpson, list(Alpha_moss_F$Plot), FUN=mean)
Simpson_moss_S <- aggregate(Alpha_moss_S$Simpson, list(Alpha_moss_S$Plot), FUN=mean)

truehist(Simpson_moss_F[,2],prob=F)
truehist(Simpson_moss_S[,2],prob=F)

# Normality tests
shapiro.test(Simpson_moss_F[,2])
shapiro.test(Simpson_moss_S[,2]) # Not normally distributed

# Comparison using Mann-Whitney test
wilcox.test(Simpson_moss_F[,2],Simpson_moss_S[,2])# Comparison using Mann-Whitney
# The Simpson index do not differ between habitats W = 73, p-value = 0.8859

## PIELOU ##
# Pielou means
Pielou_moss_F <- aggregate(Alpha_moss_F$Pielou, list(Alpha_moss_F$Plot), FUN=mean)
Pielou_moss_S <- aggregate(Alpha_moss_S$Pielou, list(Alpha_moss_S$Plot), FUN=mean)

truehist(Pielou_moss_F[,2],prob=F)
truehist(Pielou_moss_S[,2],prob=F)

# Normality tests
shapiro.test(Pielou_moss_F[,2]) # Normally distributed
shapiro.test(Pielou_moss_S[,2]) # Normally distributed

# Comparison between habitats using Mann-Whitney test
wilcox.test(Pielou_moss_F[,2],Pielou_moss_S[,2])# Comparison using Mann-Whitney
# The Pielou index do not differ between habitats W = 59, p-value = 0.5458

#----------------------- 10. Phylogenetic tree   -------------------------------------#

# BiocManager::install("phangorn")
library(phangorn); packageVersion("phangorn")
library(Biostrings); packageVersion("Biostrings")
library("DECIPHER"); packageVersion("DECIPHER")

dna_moss<- Biostrings::DNAStringSet(refseq(ps_moss.norm))
names(dna_moss) <- taxa_names(ps_moss.norm)
aligned_moss_seqs <- AlignSeqs(dna_moss)
distances_moss <- dist.ml(as.phyDat(as.DNAbin(aligned_moss_seqs)))
nj_JC69_moss <- NJ(distances_moss)

# Tree inference
fit_moss = pml(nj_JC69_moss, data=as.phyDat(as.DNAbin(aligned_moss_seqs)))
fitGTR_moss <- optim.pml(fit_moss, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
save(fitGTR_moss, file="phylogenetic_tree_moss.Rdata")
rooted_tree_moss <- midpoint(fitGTR_moss$tree)

#----------------------- 11. Faith's Phylogenetic diversity -------------------------------------#

#install.packages("picante")
library(picante); packageVersion("picante")
library(tidyr)
library(MASS)
library(ggplot2); packageVersion("ggplot2")

class(rooted_tree_moss) # Check if the tree is rooted
getRoot(rooted_tree_moss)
Faith_moss<- pd(otu_table(ps_moss.norm),rooted_tree_moss,include.root=F)
Alpha_moss <-data.frame(Alpha_moss,Faith_moss)

# Comparison between habitats
Alpha_moss_F<- Alpha_moss [which(Alpha_moss$Habitat=="Forest_tundra"),]
Alpha_moss_S<- Alpha_moss [which(Alpha_moss$Habitat=="Shrub_tundra"),]

# Faith
Faith_moss_F <- aggregate(Alpha_moss_F$PD, list(Alpha_moss_F$Plot), FUN=mean)
Faith_moss_S <- aggregate(Alpha_moss_S$PD, list(Alpha_moss_S$Plot), FUN=mean)

truehist(Faith_moss_F[,2],prob=F)
truehist(Faith_moss_S[,2],prob=F)

# Normality tests
shapiro.test(Faith_moss_F[,2])
shapiro.test(Faith_moss_S[,2])

# Comparison using Mann-Whitney test
wilcox.test(Faith_moss_F[,2],Faith_moss_S[,2]) # Comparison using Mann-Whitney
# The Faith index differ between habitats W = 115, p-value = 0.007251

## Plotting ##
# Transform data for plotting
names(Alpha_moss)[9]<-"Shannon"
names(Alpha_moss)[10]<-"Simpson"
names(Alpha_moss)[12]<-"Faith"
names(Alpha_moss)[13]<-"Richness"
Alpha.moss.tall<- Alpha_moss %>% gather (key=Diversity, value=Value,Shannon:Faith)
Alpha.moss.tall$Diversity<-as.factor(Alpha.moss.tall$Diversity)
levels(Alpha.moss.tall$Diversity)
Alpha.moss.tall$Diversity<-factor(Alpha.moss.tall$Diversity, levels=c("Shannon","Simpson","Pielou","Faith"))

theme_set(theme_bw())
Alpha_moss_plot<-ggplot(Alpha.moss.tall,aes(x=Habitat,y=Value)) +
  geom_boxplot(aes(fill=Habitat)) +
  facet_wrap(~Diversity,scale="free")+
  labs(x = "moss",
       y = "Alpha diversity measure")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

Alpha_moss_plot
write.table(Alpha_moss,"Alpha_moss.txt",sep = "\t",  row.names = TRUE)

# Plotting both 16S and nifH diversities estimates exporting a table
# Select the working directory
# setwd("C:/Users/escol/OneDrive - Universit� Laval/6_Doctorat en Biologie/1_Thesis/Chapter3_Microbiome/Microbiome_bioinformatic_analyses/Script_Bacterial_community")
Alpha_div_mean_all<- read.delim("Alpha_div_mean_all.txt",sep ="")
Alpha.all.tall<- Alpha_div_mean_all %>% gather (key=Diversity, value=Value,Shannon:Faith)
Alpha.all.tall$Diversity<-as.factor(Alpha.all.tall$Diversity)
levels(Alpha.all.tall$Diversity)
Alpha.all.tall$Diversity<-factor(Alpha.all.tall$Diversity, levels=c("Shannon","Simpson","Pielou","Faith"))

theme_set(theme_bw())
Alpha_all_plot<-ggplot(Alpha.all.tall,aes(x=Gene,y=Value)) +
  geom_boxplot(aes(fill=Habitat)) +
  stat_summary(fun = mean, geom= "point", pch=21, col = "black", bg = "white", size = 3, aes(group=Habitat), position = position_dodge(0.6)) +
  facet_wrap(~Diversity,scale="free") +
  labs(x = "",
       y = "Alpha diversity measure") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) #+
#  guides(x = guide_axis(angle = 90))

Alpha_all_plot


#----------------------- 12. GUnifrac distances  -------------------------------------#

#install.packages("GUniFrac")
library("GUniFrac")

#Create the GUnifrac distance matrix
unifracs_moss <-GUniFrac(otu_table(ps_moss.norm),rooted_tree_moss, alpha =c(0.5,1))
Gunifrac_dist_moss <- unifracs_moss$unifracs[, , "d_0.5"]
class(Gunifrac_dist_moss)
Gunifrac_dist_moss<-as.dist(Gunifrac_dist_moss)
class(Gunifrac_dist_moss)

Wunifrac_dist_moss <- unifracs_moss$unifracs[, , "d_1"]
Wunifrac_dist_moss<-as.dist(Wunifrac_dist_moss)

#----------------------- 13. Ordinations and PERMANOVA -------------------------------------#

library("vegan"); packageVersion("vegan")

# Example
# NMDS using Bray-Curtis
ord.nmds.moss <- ordinate(ps_moss.norm, method="NMDS", distance = "bray")
NMDS.moss<- plot_ordination(ps_moss.norm, ord.nmds.moss, color="Habitat", title="moss 16S - Bray-Curtis NMDS")
NMDS.moss + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw()

# NMDS using GUnifrac
ord.nmds.moss <- ordinate(ps_moss.norm, method="NMDS", distance = Gunifrac_dist_moss)
NMDS.moss<- plot_ordination(ps_moss.norm, ord.nmds.moss, color="Habitat", title="moss 16S - GUnifrac NMDS")
NMDS.moss + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw()

library(phangorn); packageVersion("phangorn")
ord.pcoa.moss <- pcoa(Gunifrac_dist_moss, correction="lingoes", row.names(Alpha_moss))
PCoA_moss <- plot_ordination(ps_moss.norm,ord.pcoa.moss, color = "Habitat", title = "moss 16S - GUnifrac PCoA")
PCoA_moss  + 
  stat_ellipse(type = "norm", linetype = 2) +
  labs(x="Axis 1 (17.80%)",y="Axis 2 (14.80%)") +
  theme_bw()

ord.pcoa.moss_Wunifrac <- pcoa(Wunifrac_dist_moss, correction="lingoes", row.names(Alpha_moss))
PCoA_moss_Wunifrac <- plot_ordination(ps_moss.norm,ord.pcoa.moss_Wunifrac, color = "Habitat", title = "moss 16S - WUnifrac PCoA")
PCoA_moss_Wunifrac  + 
  stat_ellipse(type = "norm", linetype = 2) +
  labs(x="Axis 1 (17.80%)",y="Axis 2 (14.80%)") +
  theme_bw()

# Plot both NMDS in the same plot

library("ggpubr")
ggarrange(NMDS.moss + 
            stat_ellipse(type = "norm", linetype = 2)
          , PCoA_moss  + 
            stat_ellipse(type = "norm", linetype = 2) +
            labs(x="Axis 1 (17.80%)",y="Axis 2 (14.80%)"), 
          labels = c("A)", "B)"),
          ncol = 1, nrow = 2)

ggarrange(NMDS.moss + 
            stat_ellipse(type = "norm", linetype = 2)
          , PCoA_moss  + 
            stat_ellipse(type = "norm", linetype = 2) +
            labs(x="Axis 1 (17.80%)",y="Axis 2 (14.80%)"),
          NMDS.nifH + stat_ellipse(type = "norm", linetype = 2), PCoA_nifH + 
            stat_ellipse(type = "norm", linetype = 2) +
            labs(x="Axis 1 (17.95%)",y="Axis 2 (11.01%)"),
          NMDS.soil+ stat_ellipse(type = "norm", linetype = 2), PCoA_soil+ 
            stat_ellipse(type = "norm", linetype = 2) +
            labs(x="Axis 1 (23.23%)",y="Axis 2 (21.71%)"),
          labels = c("a", "b","c","d","e","f"),
          ncol = 2, nrow = 3)

#ggarrange(PCoA_moss, PCoA_nifH , 
#          labels = c("A", "B"),
#          ncol = 1, nrow = 2)

theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# PERMANOVA 

# Transform Plots to factors.
class(Alpha_moss$Plot)
Alpha_moss$Plot<-as.factor(Alpha_moss$Plot)
levels(Alpha_moss$Plot)

# Null hypothesis: Groups have the same centroid.

# Using adonis function
permanova_0 <- adonis(Gunifrac_dist_moss ~ Habitat, data = Alpha_moss, permutations = 9999)
permanova_0

permanova_1 <- adonis(Gunifrac_dist_moss ~ Habitat, strata = Alpha_moss$Plot, data = Alpha_moss, permutations = 9999)
permanova_1


# Using adonis2 function with the perm and how command.
#perm <- how(within = Within(type = "free"), plots = Plots(strata=Alpha_moss$Plot, type = "free"),nperm = 10000)
perm <- how(nperm = 9999)

permanova_moss<- adonis2(Gunifrac_dist_moss ~ Alpha_moss$Habitat, permutations = perm)
permanova_moss

# Check the dispersion in both habitats

beta_moss <- betadisper(Gunifrac_dist_moss,Alpha_moss$Habitat)
permutest(beta_moss)
# The test indicate that the there is not enough proofs 
# to confirm that the homogeneity of multivariate dispersion between habitats differ


#----------------------- 14. LEfSe -------------------------------------#

# Export data to formatting according to LEfSe analyses
#write.table(otu_table(ps_moss.norm),"OTU_table_moss_norm.txt", sep= "\t",  row.names = TRUE)
#write.table(sample_data(p_moss2.ra),file = "Data_moss_norm.txt",sep= "\t",  row.names = TRUE)
#write.table(t(otu_table(p_moss2.ra)),file = "ASV_table_moss_ra.txt.txt",sep= "\t",  row.names = TRUE)
#write.table(tax_table(p_moss2.ra),file = "Tax_table_moss_norm.txt",sep= "\t",  row.names = TRUE)

# UPDATE FORMATTING 08-06-2023

library("data.table")
library("phyloseq")
library("ggplot2")
library("grid")
library("tidyverse")

setwd("C:/Users/escol/OneDrive - Universit� Laval/6_Doctorat en Biologie/1_Thesis/Chapter3_Microbiome/Microbiome_bioinformatic_analyses/LEfSe")

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

# For general 16S
taxa_are_rows(ps_moss.norm)
y <- ps_moss.norm %>% orient_taxa(as = "rows")
taxa_are_rows(y)

# Script from : https://ycl6.github.io/16S-Demo/3_lefse_tutorial.html#Perform_LEfSe_analysis
# 16S rDNA amplicon sequencing analysis using R (Part 3: LEfSe & GraPhlAn)
# I-Hsuan Lin. University of Manchester. July 19, 2021

tax1 = lefse_1name_obj(y, sample_data(y)$Habitat,subject = T)
lefse1 = lefse_obj(y)
lefse1 = rbind(tax1, lefse1)

# Replace unsupported chars with underscore
lefse1$name = gsub(" ","_",lefse1$name)
lefse1$name = gsub("-","_",lefse1$name)
lefse1$name = gsub("/","_",lefse1$name)

# For **Silva v132** users, apply below to fix identifical phylum & class names (Actinobacteria and Deferribacteres)
lefse1 = fix_taxa_silva132(lefse1)

# Output the prepared LEfSe input to file
write.table(lefse1, file="lefse_table_16S_UPDATED.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#install.packages("reticulate")
library(reticulate)
conda_install ("r-reticulate", "lefse")

# Set the full-paths to LEfSe's python scripts if they are not in PATH
#Sourcing Python scripts - The source_python() function enables you to source a Python script the same way you would source() an R script (Python functions and objects defined within the script become directly available to the R session).
lefse_format_input = "/lefse.py"   # e.g. /home/ngs/lefse/lefse-format_input.py
run_lefse = "/lefse_run.py"         # e.g. /home/ngs/lefse/run_lefse.py

# Run Lefse
system2(lefse_format_input, 
        args = c("/lefse_table_16S_UPDATED.txt", "/lefse_table_16S_UPDATED.in", 
                 "-c 1", "-o 1000000"))

# set Kruskal-Wallis alpha (-a) to 1 to allow returning of all P-value to perform adjustment later
system2(run_lefse, 
        args = c("/lefse_table_16S_UPDATED.in", "/lefse_table_16S_UPDATED.res", 
                 "-b 100", "-a 1", "-l 1"), stdout = TRUE)





#----------------------- 15. Racomitrium lanuginosum core microbiome based on soil samples -------------------------------------#

# Venn Diagram

#install.packages("eulerr") # If not installed
library(eulerr); packageVersion("eulerr") 
library(microbiome)
#devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities);packageVersion("microbiomeutilities")

# Combine phyloseq objects - Moss and soil samples

phylo_16S_all_filt <- merge_phyloseq(phylo_moss2,phylo_soil2)
phylo_16S_all_filt@sam_data
summary(phylo_16S_all_filt@sam_data)
summary(sample_sums(phylo_16S_all_filt))


# Transform to relative abundances using microbiome package
phylo_16S_all.microbiome <- microbiome::transform(phylo_16S_all_filt, "compositional")

# Select the type of samples
Sample_types <- unique(as.character(meta(phylo_16S_all_filt)$Sample_type))
list_core <- c() # an empty object to store information

# Write a for loop to go through each of the sample_type one by one and
# combine identified core taxa into a list.

for (n in Sample_types){ # for each variable n in sample_type
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(phylo_16S_all.microbiome, Sample_type == n) # Choose sample from Sample_type by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from soil 
                         detection = 0.001, # 0.001 in at least 90% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each Sample_type.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

mycols <- c(Moss="#99ff99", SoilCRC="#ffcc99") 
plot(venn(list_core),
     fills = mycols)
length (core_m) # Number of total core ASVs for Racomitrium lanuginosum



### Venn diagram for moss core microbiome in two habitats (Using only R. lanuginosum samples) ####

# Select the type of samples
Habitat_types <- unique(as.character(meta(p_moss_micro)$Habitat))
list_core2 <- c() # an empty object to store information

# Write a for loop to go through each of the sample_type one by one and
# combine identified core taxa into a list.

for (n in Habitat_types){ # for each variable n in sample_type
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub2 <- subset_samples(p_moss_micro, Habitat_types == n) # Choose sample from DiseaseState by n
  
  core_m2 <- core_members(ps.sub2, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core2[[n]] <- core_m2 # add to a list core taxa for each group.
  #print(list_core)
}

mycols2 <- c(Forest_tundra="#ff9999", Shrub_tundra="#99ccff") 
plot(venn(list_core2),
     fills = mycols2)
length(core_m2) # Number of total core ASVs for Racomitrium lanuginosum

# Forest tundra: "ASV79", "ASV81", "ASV103", "ASV153", "ASV161".
# Acetobacteraceae, Acidiphilium, Conexibacter, Granulicella, Acetobacteraceae


# Shrub tundra: "ASV70", "ASV111", "ASV194" "ASV196"
# Conexibacter, WPS_2, Acidisphaera, Beijerinckiaceae-1174-901-12.

#----------------------- 16. Racomitrium lanuginosum core microbiome -------------------------------------#

# BiocManager::install("microbiome")
library(microbiome); packageVersion("microbiome") 
library(phyloseq)
library(RColorBrewer)
library(ggpubr)
library(dplyr)  

p_moss_micro <- microbiome::transform(phylo_moss2, "compositional") # relative abundance in microbiome
p_moss_micro<- prune_taxa(taxa_sums(p_moss_micro ) > 0, p_moss_micro )

taxa_names(p_moss_micro)[1:3]

# Analysis of core microbiome using a threshold of 0.75
core_75 <- core_members(p_moss_micro, detection = 0.001, prevalence = 75/100, include.lowest = FALSE)
taxonomy <- as.data.frame(tax_table(p_moss_micro))
# Subset the core microbiome
core_taxa_id75 <- subset(taxonomy, rownames(taxonomy) %in% core_75)
core.abundance.75 <- sample_sums(core(p_moss_micro, detection = 0.001, prevalence = 75/100))
View(core.abundance.75)


# Plotting results
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

plot_core(p_moss_micro, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()

## Core microbiome with compositional (prevalence 75 -> min.prevalence = .75):
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)
detections <- round(detections,3)
## Also define gray color palette
gray <- gray(seq(0,1,length=5))
p1 <- plot_core(p_moss_micro, 
                plot.type = "heatmap", 
                colours = gray,
                prevalences = prevalences, 
                detections = detections, min.prevalence = .75) +
  xlab("Detection Threshold (Relative Abundance (%))")

p1 <- p1 + theme_bw() + ylab("ASVs")
p1

library(viridis)
print(p1 + scale_fill_viridis())

# Plot with taxonomic information
library(RColorBrewer)
library(knitr)

# Extract taxonomic information
df <- p1$data 
list <- df$Taxa 
tax <- as.data.frame(tax_table(p_moss_micro))
tax$ASV <- rownames(tax)
tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 
tax.unit <- tidyr::unite(tax2, Taxa_level,c("Phylum", "Class", "Order", "Family", "Genus", "ASV"), sep = "_;", remove = TRUE)
tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)
df$Taxa <- tax.unit$Taxa_level
knitr::kable(head(df))
p1$data <- df

# Plot again with the updated taxonomic information
print(p1 + scale_fill_viridis())


# Analysis of core microbiome using a threshold of 0.90
core_90 <- core_members(p_moss_micro, detection = 0.0001, prevalence = 90/100, include.lowest = FALSE)
taxonomy <- as.data.frame(tax_table(p_moss_micro))
# Subset the core microbiome
core_taxa_id <- subset(taxonomy, rownames(taxonomy) %in% core_90)
core.abundance.90 <- sample_sums(core(p_moss_micro, detection = 0.0001, prevalence = 90/100))
View(core.abundance.90)

## Core microbiome with compositional (prevalence 90 -> min.prevalence = .90:
p2 <- plot_core(p_moss_micro, 
                plot.type = "heatmap", 
                colours = gray,
                prevalences = prevalences, 
                detections = detections, min.prevalence = .90) +
  xlab("Detection Threshold (Relative Abundance (%))")

p2 <- p2 + theme_bw() + ylab("ASVs")
p2
print(p2 + scale_fill_viridis())

# Extract taxonomic information
df2 <- p2$data 
list2 <- df2$Taxa 
tax <- as.data.frame(tax_table(p_moss_micro))
tax$ASV <- rownames(tax)
tax2 <- dplyr::filter(tax, rownames(tax) %in% list2) 
tax.unit <- tidyr::unite(tax2, Taxa_level,c("Phylum", "Class", "Order", "Family", "Genus", "ASV"), sep = "_;", remove = TRUE)
tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)
df2$Taxa <- tax.unit$Taxa_level
knitr::kable(head(df))
p2$data <- df2

# Plot again with the updated taxonomic information
print(p2 + scale_fill_viridis())

# Get the best classification of the ASVs directly. Change colors (optional plot)
p_moss_micro_format <- microbiomeutilities::format_to_besthit(p_moss_micro)

p.core <- plot_core(p_moss_micro_format , 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .75) + 
  xlab("Detection Threshold (Relative Abundance (%))")

p.core + theme(axis.text.y = element_text(face="italic"))

#----------------------- 17. Microbial interaction networks - CoNet -------------------------------------#

# Formatting files for CoNet analyses

# Subset data set according to habitat
# Try with row counts = phylo_moss2

# Check normalization
summary(sample_sums(phylo_moss2))

Forest_tundra_moss16S<- subset_samples(phylo_moss2, Habitat=="Forest_tundra")
table(tax_table(Forest_tundra_moss16S)[, "Phylum"], exclude = NULL)
Forest_tundra_moss16S <-prune_taxa(taxa_sums(Forest_tundra_moss16S) >=1,Forest_tundra_moss16S)
Forest_tundra_moss16S # 36 samples, 1303 ASVs


Shrub_tundra_moss16S<- subset_samples(phylo_moss2, Habitat=="Shrub_tundra")
table(tax_table(Shrub_tundra_moss16S)[, "Phylum"], exclude = NULL)
Shrub_tundra_moss16S <-prune_taxa(taxa_sums(Shrub_tundra_moss16S) >=1,Shrub_tundra_moss16S)
Shrub_tundra_moss16S # 26 samples, 1216 ASVs


# Change the order of taxa to taxa are rows. The phyloseq objects have "taxa are columns".
library(speedyseq)

taxa_are_rows(Forest_tundra_moss16S)
Forest_tundra_moss16S <- Forest_tundra_moss16S %>% orient_taxa(as = "rows")
taxa_are_rows(Forest_tundra_moss16S)

taxa_are_rows(Shrub_tundra_moss16S)
Shrub_tundra_moss16S <- Shrub_tundra_moss16S %>% orient_taxa(as = "rows")
taxa_are_rows(Shrub_tundra_moss16S)

# Export data into tables

#write.table(Forest_tundra_moss16S@otu_table, "ASV_table_Forest_tundra_moss16S.txt",row.names = T, sep = "\t")
#write.table(Shrub_tundra_moss16S@otu_table, "ASV_table_Shrub_tundra_moss16S.txt",row.names = T, sep = "\t")

#write.table(Forest_tundra_moss16S@tax_table, "Taxa_table_Forest_tundra_moss16S.txt",row.names = T, sep = "\t")
#write.table(Shrub_tundra_moss16S@tax_table, "Taxa_table_Shrub_tundra_moss16S.txt",row.names = T, sep = "\t")

# You have to edit those files for CoNet.
# Open the files and paste them into a new Excel file.
# For the ASV_table, add the "ASV_names" column and move the row with the sample name.
# For the Taxa_table, make sure that all the "" symbols are deleted when you paste them in the Excel sheet.
# Also, delete the first row with the classification names (Kingdom, Phylum, Class, etc.)


