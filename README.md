# Microbial analyses and associated N<sup>2</sup>-fixation rates of Racomitrium lanuginosum 

## 1. Overview
This repository contains the scripts and dataset related to the project “Differentially abundant bacteria drive the N<sup>2</sup>-fixation of a widespread moss in the forest-tundra transition zone ”. The main goal of the code is to provide detailed information on the project/article analyses (reproducibility).

> **Abstract**: Bryophytes maintain symbiosis with bacteria influencing the local nutrient budget. Moss bacterial communities are composed of a core microbiome and bacteria recruited from environmental sources. Notably, symbiotic N<sup>2</sup>-fixing bacteria contribute to the N budget in northern ecosystems through biological nitrogen fixation. This process may be affected by the abundance of diazotrophs and moss nutrient content. We used the abundant moss *Racomitrium lanuginosum* in a forest tundra and shrub tundra in Northern Quebec, Canada, to investigate the bacterial and diazotrophic communities associated with habitat type using amplicon sequencing of the bacterial 16S rRNA and nifH genes and test whether the moss core microbiome has recruitment from the soil bacteria community. The nifH amplicons and element analysis were used to test the effect of diazotrophic abundance and moss nutrient content on N<sup>2</sup>-fixation activity estimated by acetylene reduction assays. Moss microbial communities between tundra types hosted similar bacterial diversity but differentially abundant groups and characteristic microbial interaction patterns. The core microbiome of *R. lanuginosum* is composed of bacteria strongly associated with northern mosses with no significant recruitment from the soil. The relative abundances of dominant diazotrophs are significantly correlated with acetylene reduction rates. In contrast, the moss nutrient content did not significantly drive N<sup>2</sup>-fixation. The proteobacterial genera *Azorhizobium* and *Rhodomicrobium* represent newly reported bacteria associated with N<sup>2</sup>-fixation rates in the tundra. We identified critical bacterial groups related to moss-bacterial symbiosis and N<sup>2</sup>-fixation in the forest-tundra transition zone, a changing environment susceptible to climate warming.
> **Keywords**: *Azorhizobium* · Biological nitrogen fixation · Core microbiome · Moss symbiosis · *Racomitrium lanuginosum* · *Rhodomicrobium*.

## 2. Software information
The code was run on the bioinformatic platform of the [Institut de Biologie Intégrative et des Systèmes (IBIS)](https://www.ibis.ulaval.ca/en/services-2/bioinformatics/documentation-servers/) at Laval University and on the statistic environment R. Check each script for details.

The following software and R packages were employed for the analyses:

**Software**
- LefSE
- CoNet
- NifMAP
- usearch/9.2.64_i86linux32
- hmmer/3.3
- FrameBot
- ncbiblast/2.11.0
- mafft/7.453
- RAxML/8.2.9
- seqmagick/0.8.4

**R/v.4.1.3 packages**
- dada2
- phyloseq
- Biostrings
- ggplot2
- phangorn
- DECIPHER
- decontam
- dplyr
- vegan
- tidyr
- MASS
- picante
- GUniFrac
- data.table
- grid
- vctrs
- cli
- tidyverse
- magrittr
- speedyseq
- eulerr
- microbiome
- microbiomeutilities
- RColorBrewer
- ggpubr
- lme4
- lmerTest
- car
- glmmTMB

## 3. Data availability
Files used for analyses are provided in the [Datasets](/Datasets) folder.

The moss vouchers are deposited at the QFAherbarium corresponding to catalogue numbers QFA-637608 to QFA-637686. The raw sequences of this study are deposited in the NCBI Sequence Read Archive (SRA) database under the [BioProject PRJNA893897](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA893897/). The moss and soil samples are associated with Bio-Samples SAMN31436785–SAMN31436859 and SAMN31439125–SAMN31439151, respectively. The moss 16S sequences correspond to SRA accessions SRR22028313–SRR22028356, moss nifHsequences to SRR22032306–SRR22032359, and soil 16S sequences to SRR22031903– SRR22031923. For detailed information, see Supplementary Information Table S1.

Differentially abundant bacteria drive the N<sup>2</sup>-fixation of a widespread moss in the forest-tundra transition zone. Available from: [ResearchGate - Dennis Escolástico.](https://www.researchgate.net/publication/372956582_Differentially_abundant_bacteria_drive_the_N2-fixation_of_a_widespread_moss_in_the_forest-tundra_transition_zone)

## 4. Citation
If you use part or the entire code in your work, please cite it using the following reference:

Escolástico-Ortiz, D.A., Blasi, C., Bellenger, JP., Derome, N. & Villarreal-A, J.C. 2023. Differentially abundant bacteria drive the N<sup>2</sup>-fixation of a widespread moss in the forest-tundra transition zone. Symbiosis 90, 193–211. https://doi.org/10.1007/s13199-023-00930-y

> [!NOTE]
Be aware that the code may contain issues related to updated software or packages and IS NOT intended to serve as an optimized pipeline but as supplementary information for the research.
