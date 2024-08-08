###-------------------- Soil bacterial community analyses of 16S amplicons related to R. lanuginosum --------------------###

# Code for soil 16S analyses. Dennis Escolástico-Ortiz. 2022


# Load packages devtools and DADA2
library("devtools")
library("dada2")
packageVersion("dada2")

##-------------------- 1. Analyses of 16S sequences from soil --------------------##

# Load the amplicons sequences
#setwd("REPLACE_your_directory")
getwd()
Path16S_soil <-"./Amplicon_sequencing/Amplicon_Sequence_files/soil_16S"
list.files(Path16S_soil)

##-------------------- 1.1 Filter and trim --------------------##
# Forward and reverse fastq filenames have the following format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs_soil <- sort(list.files(Path16S_soil, pattern="_R1_001.fastq", full.names = TRUE))
fnRs_soil <- sort(list.files(Path16S_soil, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names_soil <- sapply(strsplit(basename(fnFs_soil), "_"), `[`, 1)

# Inspect read quality profiles

plotQualityProfile(fnFs_soil[28:36])
# Forward cut in : 270,270, 260

plotQualityProfile(fnRs_soil[28:30])
# Reverse cut in : 250,230,210, 220, 210, 200, 200

# Filter and trim. Place filtered files in filtered/ subdirectory.
filtfnFs_soil <- file.path(Path16S_soil, "filtered", paste0(sample.names_soil, "_F_filt.fastq.gz"))
filtfnRs_soil <- file.path(Path16S_soil, "filtered", paste0(sample.names_soil, "_R_filt.fastq.gz"))
names(filtfnFs_soil) <- sample.names_soil
names(filtfnRs_soil) <- sample.names_soil

# Trimming
out_soil <- filterAndTrim(fnFs_soil, filtfnFs_soil, fnRs_soil,filtfnRs_soil, truncLen=c(280,230),
             maxN=0, maxEE=c(2,2), truncQ=2, trimLeft=c(10,10), rm.phix=TRUE,
             compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

head(out_soil)

# Re-check quality after trimming
plotQualityProfile(filtfnFs_soil[28:30])
plotQualityProfile(filtfnRs_soil[28:30])

##-------------------- 1.2 Infer amplicon sequence variants from soil samples --------------------##

# Estimate errors. Two hours minimum.
errF_soil <- learnErrors(filtfnFs_soil, multithread = TRUE)
errR_soil <- learnErrors(filtfnRs_soil, multithread = TRUE)
plotErrors(errF_soil, nominalQ=TRUE) # Report errors
plotErrors(errR_soil, nominalQ=TRUE)

# Sample inference
dadaFs_soil <- dada(filtfnFs_soil, err=errF_soil, multithread = TRUE, pool = "pseudo")
dadaRs_soil <- dada(filtfnRs_soil, err=errF_soil, multithread = TRUE, pool = "pseudo")
dadaFs_soil[[1]]
dadaRs_soil[[1]]

# Merge paired reads
mergers_soil <- mergePairs(dadaFs_soil, filtfnFs_soil, dadaRs_soil, filtfnRs_soil, verbose=TRUE)
head(mergers_soil[[1]])

# Construct sequence table
seqtab_soil <- makeSequenceTable(mergers_soil)
dim(seqtab_soil)

# Inspect distribution of sequence lengths 
table(nchar(getSequences(seqtab_soil)))
# Check if length match the expected amplicon length.

# Removing of sequence with lower expected length (<400 bp) or that did not merged well
seqtab_soilb<- seqtab_soil[,nchar(colnames(seqtab_soil)) %in% seq(400,450)] 
table(nchar(getSequences(seqtab_soilb)))

# Remove chimeras
#seqtab.nochim_soil<- removeBimeraDenovo(seqtab_soilb, method = "consensus", multithread=TRUE, verbose = TRUE)
#dim(seqtab.nochim_soil)
#sum(seqtab.nochim_soil)/sum(seqtab_soilb)
# 63.40% of the merged reads were chimeras.

# Track reads through the pipeline
#getN <- function(x) sum(getUniques(x))
#track_soil <- cbind(out_soil, sapply(dadaFs_soil, getN), sapply(dadaRs_soil, getN), sapply(mergers_soil, getN), rowSums(seqtab.nochim_soil))
#colnames(track_soil) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
#rownames(track_soil) <- sample.names_soil
#head(track_soil)
#track_soil
#write.table(track_soil,"Processed_reads_soil_16S.txt", sep= "\t",
#            row.names = TRUE)

# Sequence table of moss 1st run and 2 run, and soil run.

load("./soil_16S_analyses.RData")
load("./16_data_analyses.RData")
dim(seqtab1b)
dim(seqtab2b)
dim(seqtab_soilb)

alltab0 <-mergeSequenceTables(seqtab1b,seqtab2b)
dim(alltab0)

alltab <- mergeSequenceTables(alltab0,seqtab_soilb)
dim(alltab)

# Remove chimeras in all samples
seqtab.nochim.16all <- removeBimeraDenovo(alltab, method="consensus", multithread=TRUE, verbose=TRUE)
save(seqtab.nochim.16all,file="seqtab.nochim.16all.RData")
dim(seqtab.nochim.16all)
sum(seqtab.nochim.16all)/sum(alltab) # Check if this is the good object for comparison.

##-------------------- 1.3 Assign taxonomy for all samples --------------------##

set.seed(100) # Random number generator
taxa <- assignTaxonomy(seqtab.nochim.16all, "./silva_nr_v132_train_set.fa.gz", multithread=TRUE, minBoot = 80)
taxa_all<- taxa
dim(taxa_all) # 17402     6
unname(taxa_all)

#----------------------- 2. Handoff to phyloseq -------------------------------------#
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())
library("DECIPHER")
library("phangorn")

# Construct a phyloseq object with metadata
getwd()
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
phylo_16S_all # 15,750 113 samples

#----------------------------- 3. Decontam - Remove contaminants -----------------------------------#

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

##-----------------------  4. Filter according to prevalence  ----------------------- ##
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())


# Before proceed, filter out the negative controls
mps.noncontam2 # 15,676 ASVs 97 samples
mps.noncontam2@sam_data
mps.filtered <- prune_samples(sample_data(mps.noncontam2)$Sample_status == "TrueSample",mps.noncontam2)
mps.filtered@sam_data # 15,676 ASVs 89 samples

# Phyloseq object only for soil samples
phylo_soil <-subset_samples(mps.filtered,Sample_type == "Soil")
# phylo_soil <- prune_samples(sample_data(mps.filtered)$Sample_type == "Soil", mps.filtered)
phylo_soil@sam_data
table(tax_table(phylo_soil)[, "Phylum"], exclude = NULL) # All the ASV are stored in this object even the
# moss ASVs but they have 0 counts. So, they could be remove with the prune function.
phylo_soil <- prune_taxa(taxa_sums(phylo_soil) >=1,phylo_soil)
phylo_soil 

#Phyloseq object for moss samples
phylo_moss <-subset_samples(mps.filtered,Sample_type == "Moss")
phylo_moss 
table(tax_table(phylo_moss)[, "Phylum"], exclude = NULL)
phylo_moss <-prune_taxa(taxa_sums(phylo_moss) >=1,phylo_moss)
phylo_moss

# Soil data
df_soil.filtered <- as.data.frame(sample_data(phylo_soil)) 
df_soil.filtered$LibrarySize <- sample_sums(phylo_soil)
df_soil.filtered <- df_soil.filtered[order(df_soil.filtered$LibrarySize),]
df_soil.filtered
summary(df_soil.filtered) # Library size (min: 12,486 max: 29,266)

# All 16S samples
df_mps.filtered <- as.data.frame(sample_data(mps.filtered)) 
df_mps.filtered$LibrarySize <- sample_sums(mps.filtered)
df_mps.filtered <- df_mps.filtered[order(df_mps.filtered$LibrarySize),]
df_mps.filtered
summary(df_mps.filtered) # LIbrary Size min: 1,236 max: 29,266

### Soil samples prevalence ###

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(phylo_soil),
               MARGIN = ifelse(taxa_are_rows(phylo_soil ), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phylo_soil),
                    tax_table(phylo_soil))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),
                                                  sum(df1$Prevalence))})

# Prevalence: number of samples in which a taxon appears at least once.
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(phylo_soil, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(phylo_soil),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## Filter out phylum with prevalence = 1.00000
## Define prevalence threshold to keep ASV present in at least 3% of samples (Callahan et al. web):
prevalenceThreshold = 0.05 * nsamples(phylo_soil)
prevalenceThreshold # 1.35 samples threshold for soil samples

## Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
phylo_soil2 = prune_taxa(keepTaxa,phylo_soil)
table(tax_table(phylo_soil2)[, "Phylum"], exclude = NULL)
phylo_soil2 # 3,777 ASVs 27 samples


### Moss samples prevalence ###

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

#----------------------- 5. Relative abundance  -------------------------------------#

library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

# Check total abundance per phylum
table(tax_table(phylo_soil2)[, "Phylum"], exclude = NULL)
soil_phyla <- tax_glom(phylo_soil2,taxrank = "Phylum")
taxa_sums(soil_phyla)
tax_table(soil_phyla)
write.table(taxa_sums(soil_phyla),"OTU_table_soil_phyla-abundance.txt", sep= "\t",  row.names = TRUE)

# Transform to relative abundance
p_soil2.ra <- transform_sample_counts(phylo_soil2, function(x){x / sum(x)})

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
psOrd = subset_taxa(p_soil2.ra , Kingdom %in% c("Bacteria"))
plot_abundance(psOrd, Facet = "Phylum", Color = 'Phylum')+
  labs(x = "Habitat",
       y = "Relative Abundance",
       title = "Bacteria phylum relative abundance for soil associated to Racomitrim lanuginosum populations") + 
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# Barplot per phylum in two habitats
plot_bar(p_soil2.ra, x="Collection", fill="Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(panel.background = element_blank()) +
  facet_wrap(~Habitat, scales="free_x",nrow=1)

# Barplot per order in two habitats
plot_bar(p_soil2.ra, x="Collection", fill="Order") + 
  geom_bar(aes(color = Order, fill = Order), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(panel.background = element_blank()) +
  facet_wrap(~Habitat, scales="free_x",nrow=1) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Bar Plot - Top 50 ASVs
top_soil <- names(sort(taxa_sums(phylo_soil2), decreasing=TRUE))[1:50]

phylo_soil2.top2 <- prune_taxa(top_soil, phylo_soil2)
phylo_soil2.top2 <- transform_sample_counts(phylo_soil2.top2, function(OTU) OTU/sum(OTU))
table(tax_table(phylo_soil2.top2)[, "Phylum"], exclude = NULL)
table(tax_table(phylo_soil2.top2 )[, "Order"], exclude = NULL)

plot_bar(phylo_soil2.top2,x="Collection", fill="Family") + 
  geom_bar(aes(color = Family, fill = Family), stat="identity", position="stack") +
  labs(x = "Samples", y = "Relative Abundance\n") +
  facet_wrap(~Habitat, scales="free_x",nrow=1) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#----------------------- 6. Normalize data using proportions  -------------------------------------#

# Normalization using proportions. The normalized data can be used for alpha diversity estimates.
library("vegan")

summary(sample_sums(phylo_soil2))
p_soil2.ra
ps_soil.norm <- transform_sample_counts(p_soil2.ra, function(x) {x * 1000})
summary(sample_sums(ps_soil.norm))


#----------------------- 7. Alpha diversity of soil between habitats -------------------------------------#

library("vegan")
library(tidyr)
library(MASS)
library(dplyr)

# Diversity indexes
Shannon.soil <- diversity(otu_table(ps_soil.norm), index = "shannon")
Simpson.soil <- diversity(otu_table(ps_soil.norm), index = "simpson")
Pielou.soil <- Shannon.soil/log(specnumber(otu_table(ps_soil.norm)))
variables_soil <- sample_data(ps_soil.norm)

# Merge results into a dataframe
Alpha_soil <- data.frame(variables_soil,Shannon.soil,Simpson.soil,Pielou.soil)

# Comparison between habitats
Alpha.soil_F <- Alpha_soil[which(Alpha_soil$Habitat=="Forest_tundra"),]
Alpha.soil_S <- Alpha_soil[which(Alpha_soil$Habitat=="Shrub_tundra"),]

# Shannon
truehist(Alpha.soil_F[,9],prob=F)
truehist(Alpha.soil_S[,9],prob=F)

# Normality tests
shapiro.test(Alpha.soil_F[,9]) # No evidence to say that it is not normally distributed
shapiro.test(Alpha.soil_S[,9]) # No evidence to say that it is not normally distributed

# Comparison using Mann-Whitney test
wilcox.test(Alpha.soil_F[,9],Alpha.soil_S[,9]) 
# The Shannon index for soil samples do not differ between habitats W =117, p-value = 0.1995

# Simpson
truehist(Alpha.soil_F[,10],prob=F)
truehist(Alpha.soil_S[,10],prob=F)

# Normality tests
shapiro.test(Alpha.soil_F[,10]) # It is not normally distributed
shapiro.test(Alpha.soil_S[,10]) # It is not normally distributed

# Comparison using Mann-Whitney test
wilcox.test(Alpha.soil_F[,10],Alpha.soil_S[,10]) 
# The Shannon index for soil samples do not differ between habitats W =126, p-value = 0.08322

# Pielou means
truehist(Alpha.soil_F[,11],prob=F)
truehist(Alpha.soil_S[,11],prob=F)

# Normality tests
shapiro.test(Alpha.soil_F[,11]) # Normally distributed
shapiro.test(Alpha.soil_S[,11]) # Normally distributed

# Comparison between habitats using Mann-Whitney test
wilcox.test(Alpha.soil_F[,11],Alpha.soil_S[,11])# Comparison using Mann-Whitney
# The Pielou index differ between habitats for soil bacterial community.


#----------------------- 8. Phylogenetic tree   -------------------------------------#

# BiocManager::install("phangorn")
library(phangorn); packageVersion("phangorn")
library(Biostrings); packageVersion("Biostrings")
library("DECIPHER"); packageVersion("DECIPHER")

dna_soil<- Biostrings::DNAStringSet(refseq(ps_soil.norm))
names(dna_soil) <- taxa_names(ps_soil.norm)
aligned_soil_seqs <- AlignSeqs(dna_soil)
distances_soil<- dist.ml(as.phyDat(as.DNAbin(aligned_soil_seqs)))
nj_JC69_soil <- NJ(distances_soil)

# Tree inference
fit_soil = pml(nj_JC69_soil, data=as.phyDat(as.DNAbin(aligned_soil_seqs)))
fitGTR_soil <- optim.pml(fit_soil, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
save(fitGTR_soil, file="phylogenetic_tree_soil.Rdata")
rooted_tree_soil<- midpoint(fitGTR_soil$tree)

#----------------------- 9. Faith's Phylogenetic diversity -------------------------------------#

#install.packages("picante")
library(picante); packageVersion("picante")
library(tidyr)
library(MASS)
library(ggplot2); packageVersion("ggplot2")

class(rooted_tree_soil) # Check if the tree is rooted
getRoot(rooted_tree_soil)
Faith_soil<- pd(otu_table(ps_soil.norm),rooted_tree_soil,include.root=F)
Alpha_soil <-data.frame(Alpha_soil,Faith_soil)

# Comparison between habitats
Alpha_soil_F <- Alpha_soil [which(Alpha_soil$Habitat=="Forest_tundra"),]
Alpha_soil_S <- Alpha_soil [which(Alpha_soil$Habitat=="Shrub_tundra"),]

# Faith
truehist(Alpha_soil_F[,11],prob=F)
truehist(Alpha_soil_S[,11],prob=F)

# Normality tests
shapiro.test(Alpha_soil_F[,11])
shapiro.test(Alpha_soil_S[,11])

# Comparison using Mann-Whitney test
wilcox.test(Alpha_soil_F[,11],Alpha_soil_S[,11]) # Comparison using Mann-Whitney
# The Faith index do not differ between habitats W = 106, p-value = 0.4559

# Transform data for plotting
names(Alpha_soil)[9]<-"Shannon"
names(Alpha_soil)[10]<-"Simpson"
names(Alpha_soil)[11]<-"Pielou"
names(Alpha_soil)[12]<-"Faith"
names(Alpha_soil)[13]<-"Richness"
Alpha.soil.tall<- Alpha_soil %>% gather (key=Diversity, value=Value,Shannon:Faith)
Alpha.soil.tall$Diversity<-as.factor(Alpha.soil.tall$Diversity)
levels(Alpha.soil.tall$Diversity)
Alpha.soil.tall$Diversity<-factor(Alpha.soil.tall$Diversity, levels=c("Shannon","Simpson","Pielou", "Faith"))

theme_set(theme_bw())
Alpha_soil_plot<-ggplot(Alpha.soil.tall,aes(x=Habitat,y=Value)) +
  geom_boxplot(aes(fill=Habitat)) +
  facet_wrap(~Diversity,scale="free")+
  labs(x = "soil",
       y = "Alpha diversity measure")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

Alpha_soil_plot # Plot diversity of soil samples
write.table(Alpha_soil,"Alpha_soil.txt", sep= "\t",  row.names = TRUE)

#----------------------- 10. GUnifrac distances  -------------------------------------#

#install.packages("GUniFrac")
library("GUniFrac")

#Create the GUnifrac distance matrix
unifracs_soil<-GUniFrac(otu_table(ps_soil.norm), rooted_tree_soil, alpha = 0.5)
Gunifrac_dist_soil <- unifracs_soil$unifracs[, , "d_0.5"]
class(Gunifrac_dist_soil)
Gunifrac_dist_soil<-as.dist(Gunifrac_dist_soil)
class(Gunifrac_dist_soil)

#----------------------- 11. PERMANOVA -------------------------------------#

library("vegan"); packageVersion("vegan")

# Null hypothesis: Groups have the same centroid.
perm <- how(nperm = 10000)
permanova_soil <- adonis2(Gunifrac_dist_soil ~ Alpha_soil$Habitat, permutations = perm)
permanova_soil

# NMDS using GUnifrac

ord.nmds.soil <- ordinate(ps_soil.norm, method="NMDS", distance = Gunifrac_dist_soil)
NMDS.soil<- plot_ordination(ps_soil.norm, ord.nmds.soil, color="Habitat", title="soil 16S - GUnifrac NMDS")
# shape="Sample_type" for separate moss and soil samples

theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
NMDS.soil

# PCoA using GUnifrac

ord.pcoa.soil <- pcoa(Gunifrac_dist_soil, correction="lingoes",row.names(Alpha_soil))
PCoA_soil <- plot_ordination(ps_soil.norm,ord.pcoa.soil  , color = "Habitat", title = "soil 16S - GUnifrac PCoA")
PCoA_soil+ 
  stat_ellipse(type = "norm", linetype = 2) +
  labs(x="Axis 1 (23.23%)",y="Axis 2 (21.71%)") +
  theme_bw()

# Plot both NMDS and PCoA for soil in the same plot
library("ggpubr")
ggarrange(NMDS.soil+ stat_ellipse(type = "norm", linetype = 2), PCoA_soil+ 
            stat_ellipse(type = "norm", linetype = 2) +
            labs(x="Axis 1 (23.23%)",y="Axis 2 (21.71%)"), 
          labels = c("A)", "B)"),
          ncol = 1, nrow = 2)

theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#----------------------- 12. LEfSe (Optional) -------------------------------------#

# write.table(otu_table(ps_soil.norm),"OTU_table_soil_norm.txt", sep= "\t",  row.names = TRUE)
# write.table(sample_data(ps_soil.norm),file = "Data_soil_norm.txt",sep= "\t",  row.names = TRUE)
# write.table(otu_table(ps_soil.ra),file = "OTU_table_soil_ra.txt.txt",sep= "\t",  row.names = TRUE)
# write.table(tax_table(ps_soil.norm),file = "Tax_table_soil_norm.txt",sep= "\t",  row.names = TRUE)
