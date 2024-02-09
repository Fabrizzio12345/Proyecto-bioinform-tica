rm(list=ls())

# Verificar si hay dispositivos gráficos abiertos antes de intentar cerrarlos
if (length(dev.list()) > 0) {
  dev.off()
}


BiocManager::install("DESeq2")
BiocManager::install("vegan")
BiocManager::install("forcats")
BiocManager::install("PMCMR")
BiocManager::install("ggplot2")

library("phyloseq")
library("DESeq2")
library("ggplot2")
library("vegan")
library ("forcats")
library ("PMCMR")

chosen_directory <- choose.dir()

# Verifica si la selección fue exitosa
if (!identical(chosen_directory, "")) {
  setwd(chosen_directory)
}
getwd()

#Pre-processing
######################

#Import the Phyloseq Object
JH02_JH16_phyloseq <- readRDS("JH02_NT_JH16_rare_25K_phyloseq_ASVs.rds")
JH02_JH16_phyloseq 

#subset for JH16 samples 
JH16_sample <- subset_samples(JH02_JH16_phyloseq, LibraryID == "JH16")
JH16_sample
sample_sums(JH16_sample)



#subset for rhizosphere samples
JH16_rhizo <- subset_samples(JH16_sample, Microhabitat == "Rhizosphere")
#extract dryweight information
design_rhizo <- as.data.frame(as.matrix(sample_data(JH16_rhizo)))
#data distribution
hist(as.numeric(design_rhizo$Dryweight))
shapiro.test(as.numeric(design_rhizo$Dryweight))
#stat
stat <-aov(as.numeric(Dryweight) ~ Description, data = design_rhizo)
summary(stat)
TukeyHSD(stat)
#observed
design_rhizo$Description <- ordered(design_rhizo$Description, levels=c("Modern.N","Modern.A", "B1K.N","B1K.A"))
p <- ggplot(design_rhizo, aes(x=Description, y=as.numeric(Dryweight), fill=Description)) + geom_boxplot() 
p + geom_jitter( size=4,shape=22, position=position_jitter(0.2))+ scale_fill_manual(values = c("#CC79A7","#0072B2","#E69F00", "#56B4E9"))


########################
#Figure 7C
########################

#CAP analysis for plotting
JH16_CAP <- ordinate(JH16_sample, "CAP", "bray", ~ Soil * Microhabitat)
plot_ordination(JH16_sample, JH16_CAP, color = "Description", shape ="Treatment")

#assign shapes to Soil and color to Ecotype
p=plot_ordination(JH16_sample, JH16_CAP, color = "Description", shape = "Description")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_shape_manual(values = c(15, 15, 15, 15, 15, 15))
p = p + scale_colour_manual(values = c("#56B4E9","black","#E69F00","#0072B2", "white","#CC79A7"))
p + ggtitle("CAP 16S data, Bray distance")

#Permanova calculation rhizosphere samples
BC <- phyloseq::distance(JH16_rhizo, "bray")
adonis(BC ~ Soil * Treatment, data= design_rhizo, permutations = 5000)

#Permutation: free
#Number of permutations: 5000

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Soil            1    0.1845  0.1845  1.0751 0.01446 0.3223    
#Treatment       1    3.6569  3.6569 21.3142 0.28664 0.0002 ***
#  Soil:Treatment  1    0.1664  0.1664  0.9701 0.01305 0.4033    
#Residuals      51    8.7501  0.1716         0.68586           
#Total          54   12.7579                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Figure 7D
#######################

#create a DESeq object
#extract count data 
JH16_counts_integer <- otu_table(JH16_sample)
countData = as.data.frame(JH16_counts_integer)
colnames(JH16_counts_integer)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(JH16_sample)[colnames(JH16_counts_integer), ]))
rownames(colData)
class(colData)

#construct a DESeq dataset combining count data and sample information t
JH16_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Description)

#execute the differential count analysis with the function DESeq 
JH16_cds_test <- DESeq(JH16_cds, fitType="local", betaPrior=FALSE) 

#Morex conditioned calculation
##################




