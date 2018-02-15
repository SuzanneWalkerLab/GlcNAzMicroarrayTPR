##########################
#Analyses of single array replicates for comparison of detection methods 

#######################
###needed libraries and functions
###if any are missing, install them
#######################

library(limma)
library(MASS)
library(ggplot2)
library(reshape2)
source('Functions.R')

#######################
###Read in arrays and obtain mean values of test protein spots on a per-Accession basis
#######################
targets <- readTargets('TargetsMethodComparison.txt', row.names = 'ArrayID')
Data <- ReadInArrays(targets)
#correct for HA20301 having some NA descriptions, to allow full analysis of data...
Data$HA20301$genes <- readGAL(galfile='HA20301.gal', path='./lotspecific/HA20301/' )
#separate out control spots
Data$HA20358$genes$Description <- modifyDescription(Data$HA20358$genes)
Data$HA20301$genes$Description <- modifyDescription(Data$HA20301$genes)
#find only proteins that overlap between batches based upon gene list in .gpr file
Data.Overlap <- ProteinOverlap(Data)
#obtain the mean intensity value of each protein spot
Data.Means <- Collapse2Means(Data.Overlap$HA20358)
Data.Means2 <- Collapse2Means(Data.Overlap$HA20301)
Data.Means <- cbind(Data.Means, Data.Means2[,3:dim(Data.Means2)[2]])
#add HGNC and whether these are known O-GlcNAc proteins
ArrayContent <- read.csv('ArrayContent.csv', stringsAsFactors = FALSE) 
#above based upon data mining, provided with code
Data.Means <- addHGNC(Data.Means, ArrayContent, indivDB = TRUE)
#remove noncoding RNAs (long and antisense) as well as those lacking an HGNC or controls
Data.Means.Prot <- Data.Means[!grepl('LINC',Data.Means$HGNC) & !grepl('-AS',Data.Means$HGNC) &
                                !(Data.Means$HGNC=='') & !(is.na(Data.Means$HGNC)) &
                                !(grepl('Empty',Data.Means$Description)),]
#######################
###Data Analysis
#######################

#Convert to intensity ratios
Intensity.Ratio <- data.frame('CTD' = Data.Means.Prot$CTDwt1/
                                Data.Means.Prot$CTDcont1,
                              'Biotin' = Data.Means.Prot$Biotin_wt1/
                                Data.Means.Prot$Biotin_cont1,
                              'Description' = Data.Means.Prot$Description,
                              'HGNC' = Data.Means.Prot$HGNC,
                              'Known.GlcNAc' = Data.Means.Prot$Known.GlcNAc,
                              'Accession'=Data.Means.Prot$Accession)
#standardize the data using mean and standard deviation, so same cutoff can be applied
Intensity.Ratio$CTD.std <- 2^((log2(Intensity.Ratio$CTD) - mean(log2(Intensity.Ratio$CTD)))/sd(log2(Intensity.Ratio$CTD)))
Intensity.Ratio$Biotin.std <- 2^((log2(Intensity.Ratio$Biotin) - mean(log2(Intensity.Ratio$Biotin)))/
                                   sd(log2(Intensity.Ratio$Biotin)))
#hit calling
Intensity.Ratio$Hit.Biotin <- Intensity.Ratio$Biotin.std > 4
Intensity.Ratio$Hit.CTD <- Intensity.Ratio$CTD.std > 4

#######################
###Figures S2 and 1C,D 
#######################

#S2A: waterfall plots showing each staining method and relative signal intensity on each analyzed protein 
CTD.ordered <- Data.Means.Prot[]
CTD.ordered$Accession <- factor(Data.Means.Prot$Accession, levels = (
  Data.Means.Prot$Accession[order(
    Intensity.Ratio$CTD
  )]))
CTD.ordered <- melt(CTD.ordered[,c(1,5,6)])
ggplot(CTD.ordered, aes(x = Accession, y = value, color = variable )) + 
  geom_point(alpha = 0.2, show.legend = FALSE) +
  scale_y_log10() + theme_bw() + ylab('') + xlab('') +
  scale_x_discrete(breaks = c('NM_016530.1',#note- these values were chosen based upon looking at#sorted data and manually annotating
                              'BC028424.1',
                              'BC031997.1',
                              'BC008438.1',
                              'BC057811.1',
                              'NM_173822.1',
                              'NM_014667.1'),
                   labels = c('0.033',#0.033
                              '0.75',#1.3
                              '1.0',
                              '1.25',
                              '1.5',
                              '4.0',
                              '997')
  )

Biotin.ordered <- Data.Means.Prot[]
Biotin.ordered$Accession <- factor(Data.Means.Prot$Accession, levels = (
  Data.Means.Prot$Accession[order(
    Intensity.Ratio$Biotin
  )]))
Biotin.ordered <- melt(Biotin.ordered[,c(1,4,3)])
ggplot(Biotin.ordered, aes(x = Accession, y = value, color = variable )) + 
  geom_point(alpha = 0.2, show.legend = FALSE) +
  theme_bw()  + scale_y_log10()  + ylab('') + xlab('') +
  scale_x_discrete(breaks = c('NM_002710.1',#note- these values were chosen based upon looking at
                              'NM_182970.2',#sorted data and manually annotating
                              'BC059393.1',
                              'NM_012121.4',
                              'BC012815.1',
                              'BC024002.2',
                              'BC013171.1',
                              'NM_001280.1',
                              'NM_032786.1'),
                   labels = c('0.04',
                              '0.5',
                              '0.75',
                              '1.0',
                              '1.25',
                              '1.5',
                              '2.0',
                              '4.0',
                              '188')
  )
#S2B: density plots of un-standardized and standardized data
IRorder.un <- melt(Intensity.Ratio[,c('CTD','Biotin')])
ggplot(IRorder.un, aes(x=value, color = variable, fill=variable)) + 
  geom_histogram(alpha = 0.3, position = 'identity', show.legend = FALSE, bins = 100) + 
  scale_x_continuous(trans ='log2')+xlab('')+ylab('')+theme_bw() + 
  scale_x_log10()+xlab('')+ylab('')+theme_bw()

IRorder.std <- melt(Intensity.Ratio[,c('CTD.std','Biotin.std')])
ggplot(IRorder.std, aes(x=value, color = variable, fill=variable)) + 
  geom_histogram(alpha = 0.3, position = 'identity', show.legend = FALSE, bins = 100) + 
  scale_x_continuous(trans ='log2')+xlab('')+ylab('')+theme_bw() + geom_vline(xintercept = 4, linetype = 'dashed')

#S2C: qq plots of intensity ratio in log2 support a cutoff of 2.
#CTD
qqnorm(log2(Intensity.Ratio$CTD))
qqline(log2(Intensity.Ratio$CTD))
abline(v= 2, lty = 2, col = 'red')


#Biotin
qqnorm(log2(Intensity.Ratio$Biotin))
qqline(log2(Intensity.Ratio$Biotin))
abline(v= 2, lty = 2, col = 'red')

#Figure 1C: scatterplot of both detection methods
Intensity.Ratio.ordered.b <- Intensity.Ratio
Intensity.Ratio.ordered.b$Accession <- factor(Intensity.Ratio$Accession, levels =
                                                Intensity.Ratio$Accession[order(
                                                  Intensity.Ratio$Biotin.std)])
Intensity.Ratio.ordered.b <- melt(Intensity.Ratio.ordered.b[,c(6,7,8)])
ggplot(Intensity.Ratio, aes(x = Biotin.std, y = CTD.std)) + 
  geom_point(alpha = 0.5, show.legend = FALSE, size = 4) + theme_bw() +
  scale_y_log10() + scale_x_log10() + 
  ylab('CTD110.6 Signal Ratio') + xlab('Biotin Signal Ratio') +
  geom_hline(yintercept = 4, linetype = 'dashed',color = 'red', size = 2) + 
  geom_vline(xintercept = 4, linetype = 'dashed', color = 'red', size = 2)

#Data for Figure 1D
byHGNC.compare <-as.data.frame(group_by(Intensity.Ratio, HGNC) %>% summarize(Hit.CTD = max(Hit.CTD),
                                                                             Hit.Biotin = max(Hit.Biotin),
                                                                             Known.GlcNAc = max(Known.GlcNAc)))
vennDiagram(byHGNC.compare[,c('Known.GlcNAc','Hit.CTD','Hit.Biotin')])
###hit numbers and data used to make circle diameters for Venn diagram
#number of hits CTD
sum(byHGNC$Hit.CTD)#204, sqrt(204/pi)*2= diameter circle = 16.11
#number of hits Biotin
sum(byHGNC$Hit.Biotin)#251, sqrt(251/pi)*2= diameter circle = 17.87
#number of hits GlcNAc
sum(byHGNC$Known.GlcNAc)#275, sqrt(275/pi)*2 =  18.712
(22+29+9)/275 #60 proteins, 21.8%

#######################
###Tables S1 - S3 
#######################
TableS1 <- data.frame('cDNA Description' = Intensity.Ratio$Description,
                      'cDNA NCBI Accession' = Intensity.Ratio$Accession,
                      'HGNC Symbol' = Intensity.Ratio$HGNC,
                      'In PhosphoSitePlus?' = Data.Means.Prot$PhosphoSitePlus,
                      'In Wang et al.?' = Data.Means.Prot$Wang,
                      'CTD Signal Intensity Control' = Data.Means.Prot$X52641.CTD110.control.May.31.2014,
                      'CTD Signal Intensity OGT Treated' = Data.Means.Prot$X52630.CTD110.Exp.May.31.2014,
                      'CTD Signal Intensity Ratio (unstandardized)' = Intensity.Ratio$CTD,
                      'CTD Signal Intensity Ratio (standardized)' = Intensity.Ratio$CTD.std,
                      'CTD Hit?' = Intensity.Ratio$Hit.CTD,
                      'Biotin Signal Intensity Control' = Data.Means.Prot$MO125_121914_control,
                      'Biotin Signal Intensity OGT Treated' = Data.Means.Prot$MO125_121914_experiment,
                      'Biotin Signal Intensity Ratio (unstandardized)' = Intensity.Ratio$Biotin,
                      'Biotin Signal Intensity Ratio (standardized)' = Intensity.Ratio$Biotin.std,
                      'Biotin Hit?' = Intensity.Ratio$Hit.Biotin)
write.csv(TableS1, 'TableS1.csv') #Names modified to as written above and sorted to form final .xlsx

#Table S2
IVTOrtizMeoz <- c('OTX2','HGS','NR3C1','E2F8','C20orf112','BAIAP2','EPB49','SSBP3',
                  'MAPK1IP1L','MEF2D','SSBP2','MEF2A')
IVTOrtizMeoz.HGNC <- alias2Symbol(IVTOrtizMeoz)
#how many of these are in the analyzed data?
sum(IVTOrtizMeoz.HGNC %in% Intensity.Ratio$HGNC)
IVTOrtizMeoz.HGNC <- IVTOrtizMeoz.HGNC[IVTOrtizMeoz.HGNC %in% Intensity.Ratio$HGNC]
BiotinOverlap <- IVTOrtizMeoz.HGNC[IVTOrtizMeoz.HGNC %in% Intensity.Ratio$HGNC[Intensity.Ratio$Hit.Biotin]]
CTDOverlap <- IVTOrtizMeoz.HGNC[IVTOrtizMeoz.HGNC %in% Intensity.Ratio$HGNC[Intensity.Ratio$Hit.CTD]]
fullDataIVTMeoz <- Intensity.Ratio[Intensity.Ratio$HGNC %in% IVTOrtizMeoz.HGNC,]
#Table S3
CFProteins <- c('ARNT2','SKA3', 'HDAC4','HDAC7','MBNL3', 'RUNX1T1', 
                'SCEL', 'TLE3') #all are already HGNC
CF.IR <- Intensity.Ratio[Intensity.Ratio$HGNC %in% CFProteins,]