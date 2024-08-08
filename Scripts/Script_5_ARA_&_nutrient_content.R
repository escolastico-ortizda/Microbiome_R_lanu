## ----------- Acetylene reduction rates of Racomitrium lanuginosum microbiome in two type of tundras ----------------- ##

# Code for statistical analyses of acetylene reduction rates. Dennis Escolástico-Ortiz. 2022

library(ggplot2); packageVersion("ggplot2")
library(tidyr)

# Acetylene reduction assay results

ARA<- read.delim("C:/Users/escol/OneDrive - Université Laval/6_Doctorat en Biologie/1_Thesis/Chapter3_Microbiome/ARA/Mean_ARA_results.txt",header=T)

theme_set(theme_bw())

ggplot(ARA, aes(x=Habitat,y=Mean)) + 
  geom_boxplot(aes(fill=Habitat),show.legend = FALSE) +
  stat_summary(fun = mean, geom= "point", pch=21, col = "black", bg = "white", size = 3, aes(group=Habitat), position = position_dodge(0.6)) +
  ylab(expression(paste("Acetylene reduction " ~ "(nmol. ", h^{-1},".moss ", g^{-1}, ")")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_discrete(name="")

tiff('Figure_6_Acetylene_reduction_rates_R_lanuginosum_two_tundras.tiff', units="mm", width=120, height=100, res=300, compression = 'lzw')
# Run the command for the figure
dev.off()
ggsave("Figure_6_Acetylene_reduction_rates_R_lanuginosum_two_tundras.svg",units="mm", width=120, height=100, dpi=300)



#labs(x = "",
#     y = "Acetylene reduction (nmol.h-1.moss g -1)",
#     title = "") +

# Nitrogen fixation comparison between two types of tundra

library(MASS)

ARA_F<- ARA[which(ARA$Habitat=="Forest tundra"),]
ARA_S<- ARA[which(ARA$Habitat=="Shrub tundra"),]

truehist(ARA_F[,4],prob=F)
truehist(ARA_S[,4],prob=F)

# Shapiro-Wilk test of normality for one variable
# Null hypothesis: The population is normally distributed
# From the output, the p-value > 0.05 implying that the distribution of the data are significantly different from normal distribution.

shapiro.test(ARA_F[,4])
mean (ARA_F[,4])
shapiro.test(ARA_S[,4])
mean (ARA_S[,4], na.rm=TRUE)

# Mann-Whitney Wilcoxon test for independent populations without normal distributions.
# If both x and y are given and paired is FALSE, a Wilcoxon rank sum test (equivalent to the Mann-Whitney test) is carried out.
# In this case, the null hypothesis is that the distributions of x and y differ by a location shift of mu
# and the alternative is that they differ by some other location shift (and the one-sided alternative "greater" is that x is shifted to the right of y).

wilcox.test(ARA_F[,4],ARA_S[,4])

# There is not enough evidence to say that the means are not similar
# W = 57, p-value = 0.169
# Using the whole dataset
# W = 467, p-value = 0.01198

## ------------------------- Metal content of Racomitrium lanuginosum in two types of tundra  ------------------------- ##

# Metal content results

Metals<- read.delim("C:/Users/escol/OneDrive - Université Laval/6_Doctorat en Biologie/1_Thesis/Chapter3_Microbiome/Element_concentration_moss/Mean_element_concentrations.txt",header=T)
Metals<-Metals[,-3]
# Metal concentration comparison between two types of tundra

library(MASS)

Metals_F<- Metals[which(Metals$Habitat=="Forest tundra"),]
Metals_S<- Metals[which(Metals$Habitat=="Shrub tundra"),]

# P concentrations
truehist(Metals_F[,2],prob=F) # P Forest tundra
truehist(Metals_S[,2],prob=F) # P Shrub tundra

shapiro.test(Metals_F[,2])
shapiro.test(Metals_S[,2])

wilcox.test(Metals_F[,2],Metals_S[,2]) # Comparison using Mann-Whitney
# The mean P concentrations do not differ between tundras.W = 99, p-value = 0.1174

# Fe concentrations
truehist(Metals_F[,5],prob=F) # Fe Forest tundra
truehist(Metals_S[,5],prob=F) # Fe Shrub tundra

shapiro.test(Metals_F[,5])
shapiro.test(Metals_S[,5])

wilcox.test(Metals_F[,5],Metals_S[,5]) # Comparison using Mann-Whitney
# The mean Fe concentrations do not differ between tundras.W = 96, p-value = 0.1639 

# Mo concentrations
truehist(Metals_F[,6],prob=F) # Mo Forest tundra
truehist(Metals_S[,6],prob=F) # Mo Shrub tundra

shapiro.test(Metals_F[,6])
shapiro.test(Metals_S[,6])

wilcox.test(Metals_F[,6],Metals_S[,6]) # Comparison using Mann-Whitney
# The mean Mo concentrations do not differ between tundras.W = 70, p-value = 0.9538

# V concentrations
truehist(Metals_F[,4],prob=F) # V Forest tundra
truehist(Metals_S[,4],prob=F) # V Shrub tundra

shapiro.test(Metals_F[,4])
shapiro.test(Metals_S[,4])

wilcox.test(Metals_F[,4],Metals_S[,4]) # Comparison using Mann-Whitney
# The mean Mo concentrations do not differ between tundras.W = 103, p-value = 0.07224

# Transform data to plot in a desired order.
Metals.tall<- Metals %>% gather (key=Element, value=Value,P:Mo)
Metals.tall

# Log transform
Metals.tall$logValue=log(Metals.tall$Value+1)
Metals.tall

# Change order of factor levels
Metals.tall$Element=as.factor(Metals.tall$Element)
Metals.tall$Element <- ordered(Metals.tall$Element,levels=c("P","Fe","V","Mo"))
str(Metals.tall)

# Plot
ggplot(Metals.tall,aes(x=Element,y=logValue)) + 
  geom_boxplot(aes(fill=Habitat)) +
  stat_summary(fun = mean, geom= "point", pch=21, col = "black", bg = "white", size = 3, aes(group=Habitat), position = position_dodge(0.6)) +
  labs(x = "Nutrient",
       y = "Concentration in ppm log(Y+1)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# Save as high resolution figures
tiff("Figure_7_Moss_nutrient_content_700_500.tiff", units="mm", width=174, height=120, res=300, compression = 'lzw')
# Run command for figure
dev.off()
ggsave("Figure_7_Moss_nutrient_content_700_500.svg",units="mm", width=174, height=120, dpi=300)
