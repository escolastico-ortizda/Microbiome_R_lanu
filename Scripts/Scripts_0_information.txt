Bioinformatic scripts of N2-fixation associated to the northern distributed moss Racomitrium lanuginosum.

1. Overview
Bryophytes maintain symbiosis with bacteria influencing the local nutrient budget. Moss bacterial communities are com-posed of a core microbiome and bacteria recruited from environmental sources. Notably, symbiotic N2-fixing bacteriacontribute to the N budget in northern ecosystems through biological nitrogen fixation. This process may be affected bythe abundance of diazotrophs and moss nutrient content. We used the abundant moss Racomitrium lanuginosum in a foresttundra and shrub tundra in Northern Quebec, Canada, to investigate the bacterial and diazotrophic communities associatedwith habitat type using amplicon sequencing of the bacterial 16S rRNA and nifH genes and test whether the moss coremicrobiome has recruitment from the soil bacteria community. The nifH amplicons and element analysis were used to testthe effect of diazotrophic abundance and moss nutrient content on N2-fixation activity estimated by acetylene reductionassays. Moss microbial communities between tundra types hosted similar bacterial diversity but differentially abundantgroups and characteristic microbial interaction patterns. The core microbiome of R. lanuginosum is composed of bacteriastrongly associated with northern mosses with no significant recruitment from the soil. The relative abundances of domi-nant diazotrophs are significantly correlated with acetylene reduction rates. In contrast, the moss nutrient content did notsignificantly drive N2-fixation. The proteobacterial genera Azorhizobium and Rhodomicrobium represent newly reportedbacteria associated with N2-fixation rates in the tundra. We identified critical bacterial groups related to moss-bacterialsymbiosis and N 2-fixation in the forest-tundra transition zone, a changing environment susceptible to climate warming.

Keywords Azorhizobium · Biological nitrogen fixation · Core microbiome · Moss symbiosis · Racomitrium lanuginosum · Rhodomicrobium.

2. Software information

Software
LefSE
CoNet
NifMAP
usearch/9.2.64_i86linux32
hmmer/3.3
FrameBot
ncbiblast/2.11.0
mafft/7.453
RAxML/8.2.9
seqmagick/0.8.4

R/v.4.1.3 packages
dada2
phyloseq
Biostrings
ggplot2
phangorn
DECIPHER
decontam
dplyr
vegan
tidyr
MASS
picante
GUniFrac
data.table
grid
vctrs
cli
tidyverse
magrittr
speedyseq
eulerr
microbiome
microbiomeutilities
RColorBrewer
ggpubr
lme4
lmerTest
car
glmmTMB

3. Data availability
 
The raw sequences of this study are deposited in the NCBI Sequence Read Archive (SRA) database under the BioProject PRJNA893897. The moss and soil samples are associated with BioSamples SAMN31436785–SAMN31436859 and SAMN31439125– SAMN31439151, respectively. The moss 16S sequences correspond to SRA accessions SRR22028313–SRR22028356, moss nifH sequences to 136 SRR22032306–SRR22032359, and soil 16S sequences to SRR22031903– SRR22031923. For detailed information, see Supplementary Information Table S1-1 in the published article.

Files used for analyses are provided in the dataset folder.

4. Citation
If you use part or the entire code in your work, please cite it using the following referece:

Escolástico-Ortiz D.A., Blasi C., Bellenger J.P., Derome N., Villarreal-A. J.C. 2023. Differentially abundant bacteria drive the N2-fixation of a widespread moss in the forest-tundra transition zone. Symbiosis, 90(2),193–211. https://doi.org/10.1007/s13199-023-00930-y  


1. Moss 16S data processing (DADA2, etc.). R environment
2. Moss nifH data processing (DADA2, NifMap results, etc.). R environment
3.1 - 3.4 NifMAP pipeline. 
4. Soil 16S data processing. R environment
5. ARA & nutrient content
6. Acetylene reduction regression analyses. R environment
