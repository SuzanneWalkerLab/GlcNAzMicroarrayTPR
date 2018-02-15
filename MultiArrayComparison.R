##########################
#Analyses of multiple array replicates for comparison of wild-type and 5N5A mutant 

#######################
###needed libraries and functions
###if any are missing, install them
#######################
library(limma)
library(qvalue)
library(ggplot2)
library(MASS)
library(reshape2)
source('Functions.R')
targets <- readTargets('Targets.txt', row.names = 'ArrayID')

#######################
###Read in arrays 
#######################

Raw.Data <- list()
spottypes <- readSpotTypes()
for(i in unique(targets$Batch)){
  Raw.Data[[i]] <- read.maimages(targets$FileName[targets$Batch==i],
                                 source='genepix.median', 
                                 annotation = c('Block', 'Row', 'Column', 'ID', 'Name',
                                                'Description'),
                                 wt.fun = wtflags(weight=0,cutoff=-50))
  Raw.Data[[i]]$genes$Status <- controlStatus(spottypes, Raw.Data[[i]])
}

#correct background on all arrays using normal-exponential model
Data <- list()
for(i in unique(targets$Batch)){
  Data[[i]] <- backgroundCorrect(Raw.Data[[i]], method = 'normexp')
}

#correct for HA20301 having some missing descriptions
Data$HA20301$genes <- readGAL(galfile='HA20301.gal', path='./lotspecific/HA20301/' )

#find control files and use them to add concentrations for V5-tagged biotinylated controls
ControlFiles <- list()
for(i in names(Data)){
  ControlFiles[[i]]<-paste('./lotspecific/',as.character(i),'/',as.character(i),'_Control.txt',sep='')
  ControlFiles[[i]]<-read.table(ControlFiles[[i]], header = T, sep='\t', row.names = NULL)
  colnames(ControlFiles[[i]]) <-
    c(colnames(ControlFiles[[i]])[2:length(colnames(ControlFiles[[i]]))], 'X')
  ControlFiles[[i]] <- ControlFiles[[i]][,colnames(ControlFiles[[i]])!='X']#manipulations needed because there is a trailing \t in the data files...
}

for(i in names(Data)){
  Data[[i]]$genes <- addV5Conc(Data[[i]]$genes,ControlFiles[[i]])
}
#pull together GlcNAz data into a single list rather than separating by batches
Data.GlcNAz <- Data$HA20358
Data.GlcNAz$R <- cbind(Data.GlcNAz$R,Data$HA20446$R)
Data.GlcNAz$G <- cbind(Data.GlcNAz$G,Data$HA20446$G)
Data.GlcNAz$weights <- cbind(Data.GlcNAz$weights,Data$HA20446$weights)
Data.GlcNAz$targets <- c(Data.GlcNAz$targets,Data$HA20446$targets)
Data.GlcNAz$genes$Description <- modifyDescription(Data.GlcNAz$genes)
Data.GlcNAz$genes$Name <- modifyAccession(Data.GlcNAz$genes)

#######################
###Normalize arrays and get averages
#######################
#normalize for staining
Data.GlcNAz.norm <- normStain(Data.GlcNAz)

#normalize for ERalpha activity

colFactors <- c(rep(1,times=8),rep(0,6), 2, 0, rep(2,2))

Data.GlcNAz.ERa.Norm <- NormERalpha(Data.GlcNAz.norm,colFactors,relativeActivity = c(1,0.434))

MeansArray <- Collapse2Means(Data.GlcNAz.ERa.Norm)
ArrayContent <- read.csv('ArrayContent.csv',
                         stringsAsFactors = FALSE)
MeansArray <- addHGNC(MeansArray, ArrayContent, indivDB = TRUE)
MeansArray.Prot <- MeansArray[!grepl('LINC',MeansArray$HGNC) & !grepl('-AS',MeansArray$HGNC) &
                                !(MeansArray$HGNC=='') & !(is.na(MeansArray$HGNC))&
                                !(grepl('Empty',MeansArray$Description)),]

#######################
###Permutation test
#######################

eBayesResult <- eBayesArrays(
  MeansArray.Prot, 
  c(rep('wt',times=8),rep('control',6), 'nc5N5A', 'control', rep('nc5N5A',2))
)
eBayesResultCorr <- eBayesArrays(MeansArray.Prot, 
                                 c(rep('wt',times=8),rep('control',6), 'nc5N5A', 'control',
                                   rep('nc5N5A',2)), TRUE)
wtcontPermuted <- permuteBackgroundwtcont(
  eBayesResultCorr$p.value[,1],
  eBayesResultCorr$t[,1],
  MeansArray.Prot[,!(colnames(MeansArray.Prot) %in%
                       c('Accession','Description','HGNC','Known.GlcNAc',
                         'PhosphoSitePlus','Wang'))],
  c(rep('wt',times=8),rep('control',6), 'nc5N5A',
    'control', rep('nc5N5A',2)),
  c('wt','control'), alpha = 0.1,oneway = TRUE, N = 100
)
nc5N5AcontPermuted <- permuteBackground5N5Acont(
  eBayesResultCorr$p.value[,2],
  eBayesResultCorr$t[,2],
  MeansArray.Prot[,!(colnames(MeansArray.Prot) %in%
                       c('Accession','Description','HGNC','Known.GlcNAc',
                         'PhosphoSitePlus','Wang'))],
  c(rep('wt',times=8),rep('control',6), 'nc5N5A',
    'control', rep('nc5N5A',2)),
  c('nc5N5A','control'), alpha = 0.1,oneway = TRUE, N = 100
)
wtv5N5AcontPermuted <- permuteBackgroundwtv5N5A(
  eBayesResult$p.value[,3],
  eBayesResult$t[,3],
  MeansArray.Prot[,!(colnames(MeansArray.Prot) %in%
                       c('Accession','Description','HGNC','Known.GlcNAc',
                         'PhosphoSitePlus','Wang'))],
  c(rep('wt',times=8),rep('control',6), 'nc5N5A',
    'control', rep('nc5N5A',2)),
  c('wt', 'nc5N5A'), alpha = 0.1,oneway = FALSE, N = 100
)
wtcont.qval <- qvalue(wtcontPermuted$pvals[,1],pi0.method = 'bootstrap')
nc5N5A.qval <- qvalue(nc5N5AcontPermuted$pvals[,1],pi0.method = 'bootstrap')
wtv5N5A.qval <- qvalue(wtv5N5AcontPermuted$pvals[,1],pi0.method = 'bootstrap')

#make table of final results
FinalResults <- MeansArray.Prot[,c(1,2,21,22,23,24,11:16,18,3:10,17,19,20)]
#add Median-subtracted data
FinalResults.medCorrect <- medianCorrect(FinalResults, 7:13, list('wt' = 14:21, '5N5A' = 22:24))
colnames(FinalResults.medCorrect) <- paste(colnames(FinalResults.medCorrect), 'MedianCorrected')
FinalResults <- cbind(FinalResults, FinalResults.medCorrect[,7:24])
FinalResultsMeds <- medianSignal(FinalResults, list('control' = 7:13, 
                                                    'wt'= 14:21,
                                                    '5N5A' = 22:24,
                                                    'controlMedCor' = 25:31,
                                                    'wtMedCor' = 32:39,
                                                    '5N5AMedCor' = 40:42))
FinalResults.sum <- cbind(FinalResults[,1:6],FinalResultsMeds[,c(1:3,5,6)]) #note- this excludes individual arrays. Can be obtained from NCBI GEO and
#use of analysis scripts.
FinalResults.sum[,c('wt glycosylation qvalue','5N5A glycosylation qvalue','wt vs. 5N5A qvalue')] <- cbind(wtcont.qval$qvalues,
                                                                                                          nc5N5A.qval$qvalues,
                                                                                                          wtv5N5A.qval$qvalues)
FinalResults.sum[,c('wt fold signal','5N5A fold signal','wt fold signal (Median Corrected)','5N5A fold signal (Median Corrected)',
                    'wt vs. 5N5A fold signal')] <- cbind(
                      FinalResults.sum[,'wt Median']/FinalResults.sum[,'control Median'],
                      FinalResults.sum[,'5N5A Median']/FinalResults.sum[,'control Median'],
                      FinalResults.sum[,'wtMedCor Median']/FinalResults.sum[,'control Median'],
                      FinalResults.sum[,'5N5AMedCor Median']/FinalResults.sum[,'control Median'],
                      FinalResults.sum[,'wt Median']/FinalResults.sum[,'5N5A Median']
                    )
FinalResults.sum[,c('wt Hit','5N5A Hit')] <- cbind(#865 wt hits, 35 5N5A hits
  (FinalResults.sum$`wt glycosylation qvalue` < 0.05) & (FinalResults.sum$`wt fold signal (Median Corrected)` > 2),
  (FinalResults.sum$`5N5A glycosylation qvalue` < 0.05) & (FinalResults.sum$`5N5A fold signal (Median Corrected)`> 2)
)
FinalResults.sum[,'wt v. 5N5A Hit'] <- (FinalResults.sum$`wt vs. 5N5A qvalue` < 0.05) & (FinalResults.sum$`wt vs. 5N5A fold signal` > 2) &
  (FinalResults.sum$`wt Hit` | FinalResults.sum$`5N5A Hit`)#519 5N5A hits. Note I looked into less than 0.5 signal- none show up.

byHGNC.all <-as.data.frame(group_by(FinalResults.sum, HGNC) %>% summarize(`wt v. 5N5A Hit` = max(`wt v. 5N5A Hit`),
                                                                          `wt Hit` = max(`wt Hit`),
                                                                          `5N5A Hit` = max(`5N5A Hit`),
                                                                          Known.GlcNAc = max(Known.GlcNAc)))

#######################
###Figures S6-S10, 3A
#######################

#S6A: boxplots of individual arrays
arrays <- colnames(Data.GlcNAz$R)
Unstained <- melt(Data.GlcNAz$R[!grepl('Control', Data.GlcNAz$genes$Description),])#don't plot controls
Unstained$OGT <- ''
Unstained$OGT[Unstained$Var2 %in% arrays[c(1:8)]] <- 'wt'
Unstained$OGT[Unstained$Var2 %in% arrays[c(9:14,16)]] <- 'control'
Unstained$OGT[Unstained$Var2 %in% arrays[c(15,17,18)]] <- '5N5A'
Unstained$Var2 <- factor(Unstained$Var2, levels = arrays[c(9:14,16,1:8,15,17,18)])
ggplot(Unstained, aes(fill = OGT, x = Var2, y = value)) + 
  geom_boxplot(outlier.color = alpha('black', 0.1), show.legend = FALSE) + 
  scale_y_log10() + theme_bw() + theme(axis.text.x = element_blank(),
                                       panel.grid.major.x = element_blank(),
                                       panel.grid.minor.x = element_blank(),
                                       axis.ticks.x = element_blank()
  ) + xlab('') + ylab('')
#S6B linear fits of select array blocks (one each type, show 3 blocks 1 array per treatment),
#then aftercorrection
BlockUnstain <- cbind(Data.GlcNAz$genes[,c('Block','Description','V5conc')],Data.GlcNAz$R[,c(1,10,17)])
#note- arrays for this comparison arbitrarily chosen
BlockUnstain$Block <- as.character(BlockUnstain$Block)
colnames(BlockUnstain) <- c(colnames(BlockUnstain)[1:3],'wt','control','nc5N5A')
BlockUnstain <- BlockUnstain[BlockUnstain$Block %in% c('1','20'),]
V5controls.un <- BlockUnstain[grepl('Control: V5control', BlockUnstain$Description),]
V5controls.un <- melt(V5controls.un, id.vars = c('Block','Description','V5conc'),
                      measure.vars = c('wt','control','nc5N5A'))
ggplot(V5controls.un, aes(x = V5conc, y = value, color = interaction(variable, Block))) + 
  geom_point(alpha = 0.5, show.legend = FALSE) + geom_smooth(method = 'glm', alpha = 0.5,show.legend = FALSE) + 
  scale_color_manual(values = c(rep(hue_pal()(3)[c(3,2,1)],2))) + theme_bw()
#corrected
BlockStain <- cbind(Data.GlcNAz.norm$genes[,c('Block','Description','V5conc')],Data.GlcNAz.norm$R[,c(1,10,17)])
BlockStain$Block <- as.character(BlockStain$Block)
colnames(BlockStain) <- c(colnames(BlockStain)[1:3],'wt','control','nc5N5A')
BlockStain <- BlockStain[BlockStain$Block %in% c('1','20'),]
V5controls <- BlockStain[grepl('Control: V5control', BlockStain$Description),]
V5controls <- melt(V5controls, id.vars = c('Block','Description','V5conc'),
                   measure.vars = c('wt','control','nc5N5A'))
ggplot(V5controls, aes(x = V5conc, y = value, color = interaction(variable, Block))) + 
  geom_point(alpha = 0.5,show.legend = FALSE) + geom_smooth(method = 'glm', show.legend = FALSE) + 
  scale_color_manual(values = c(rep(hue_pal()(3)[c(3,2,1)],2))) + theme_bw()

#S6C show control boxplots before and after correction
Stained <- melt(Data.GlcNAz.norm$R[!grepl('Control', Data.GlcNAz.norm$genes$Description),])#don't plot controls
Stained$OGT <- ''
Stained$OGT[Stained$Var2 %in% arrays[c(1:8)]] <- 'wt'
Stained$OGT[Stained$Var2 %in% arrays[c(9:14,16)]] <- 'control'
Stained$OGT[Stained$Var2 %in% arrays[c(15,17,18)]] <- '5N5A'
Stained$Var2 <- factor(Stained$Var2, levels = arrays[c(9:14,16,1:8,15,17,18)])

#after correction
ggplot(Stained[Stained$OGT=='control',], aes(x = Var2, y = value)) + 
  geom_boxplot(outlier.color = alpha('black', 0.1), show.legend = FALSE, fill = '#00BA38') + 
  scale_y_log10() + theme_bw() + theme(axis.text.x = element_blank(),
                                       panel.grid.major.x = element_blank(),
                                       panel.grid.minor.x = element_blank(),
                                       axis.ticks.x = element_blank()
                                       ) + xlab('') + ylab('')
#before correction
ggplot(Unstained[Unstained$OGT=='control',], aes(x = Var2, y = value)) + 
  geom_boxplot(outlier.color = alpha('black', 0.1), show.legend = FALSE, fill = '#00BA38') + 
  scale_y_log10() + theme_bw() + theme(axis.text.x = element_blank(),
                                       panel.grid.major.x = element_blank(),
                                       panel.grid.minor.x = element_blank(),
                                       axis.ticks.x = element_blank()
  ) + xlab('') + ylab('')
#S6D
#after correction
ggplot(Stained[Stained$OGT!='control',], aes(fill=OGT, x = Var2, y = value)) + 
  geom_boxplot(outlier.color = alpha('black', 0.1), show.legend = FALSE) + 
  scale_y_log10() + theme_bw() + theme(axis.text.x = element_blank(),
                                       panel.grid.major.x = element_blank(),
                                       panel.grid.minor.x = element_blank(),
                                       axis.ticks.x = element_blank()
  ) + xlab('') + ylab('') + scale_fill_manual(values = c('#F8766D','#619CFF'))
#before correction
ggplot(Unstained[Unstained$OGT!='control',], aes(fill=OGT, x = Var2, y = value)) + 
  geom_boxplot(outlier.color = alpha('black', 0.1), show.legend = FALSE) + 
  scale_y_log10() + theme_bw() + theme(axis.text.x = element_blank(),
                                       panel.grid.major.x = element_blank(),
                                       panel.grid.minor.x = element_blank(),
                                       axis.ticks.x = element_blank()
  ) + xlab('') + ylab('')+ scale_fill_manual(values = c('#F8766D','#619CFF'))

#Figure S7A

ERkinetics <- read.csv(file = 'ERa_timecourseRLU.csv',
                       stringsAsFactors = FALSE)
ERk <- melt(ERkinetics, id.vars = 'time')
ERk$nc5N5A <- 0
ERk$nc5N5A[ERk$variable=='X5N5A']<- ERk$time[ERk$variable=='X5N5A']
ERk$wt <- 0
ERk$wt[ERk$variable=='wt']<- ERk$time[ERk$variable=='wt']
x <- lm(value ~ nc5N5A + wt, data = ERk[ERk$time<11,])
summary(x) #y = 14263±785 + time(wt) * 5337.4±128.6 + time(5N5A) * 2299.2±128.6
ERk$variable <- factor(ERk$variable, levels = c('X5N5A','wt'))
predictTime <- seq(2,10,0.01)
predictFake <- data.frame('wt' = predictTime)
predictFake$nc5N5A <- 0
predictFake$variable <- 'wt'
predictFake2 <- data.frame('nc5N5A' = predictTime)
predictFake2$wt <- 0
predictFake2$variable <- 'X5N5A'
predictFake <- rbind(predictFake,predictFake2)
y <- predict(x, predictFake, interval = 'predict')
y <- as.data.frame(y)
y$wt <- predictFake$wt
y$X5N5A <- predictFake$nc5N5A
ggplot(ERk[ERk$time<11,], aes(x = time, y = value, color = variable)) + geom_point(show.legend = FALSE) + 
  scale_color_manual(values =c('#F8766D','#619CFF')) +
  geom_abline(intercept = 14263, slope = 5337.4, color = '#619CFF') + #3color blued
  geom_abline(intercept = 14263, slope = 2299.2, color = '#F8766D') + # 3 color red
  geom_line(data = y[predictFake$variable == 'wt',], aes(x = wt, y = lwr), color = '#619CFF',
            linetype = 'dashed')+
  geom_line(data = y[predictFake$variable == 'wt',], aes(x = wt, y = upr), color = '#619CFF',
            linetype = 'dashed')+
  geom_line(data = y[predictFake$variable == 'X5N5A',], aes(x = X5N5A, y = lwr), color = '#F8766D',
            linetype = 'dashed')+
  geom_line(data = y[predictFake$variable == 'X5N5A',], aes(x = X5N5A, y = upr), color = '#F8766D',
            linetype = 'dashed')+theme_bw()

#Figure S7B
aa <- melt(Data.GlcNAz.ERa.Norm$R)
aa$OGT <- ''
aa$OGT[aa$variable %in% colnames(Data.GlcNAz.ERa.Norm$R)[c(1:8)]] <- 'wt'
aa$OGT[aa$variable %in% colnames(Data.GlcNAz.ERa.Norm$R)[c(9:14,16)]] <- 'control'
aa$OGT[aa$variable %in% colnames(Data.GlcNAz.ERa.Norm$R)[c(15,17,18)]] <- '5N5A'
aa$variable <- factor(aa$variable, levels = colnames(Data.GlcNAz.ERa.Norm$R)[c(9:14,16,1:8,15,17,18)])
ggplot(aa, aes(x = variable, fill = OGT, y = value)) + geom_boxplot(outlier.color = alpha('black', 0.1),
                                                                    show.legend = FALSE) +
  scale_y_log10()+theme_bw() +
  theme(axis.text.x = element_blank(),panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank()
  )+ xlab('') + ylab('')

#S8: median normalization- use density plots

aa <- melt(MeansArray.Prot, id.vars = c("Accession", "Description", "HGNC", "Known.GlcNAc", 
                                        "PhosphoSitePlus", "Wang"))
aa$OGT <- ''
aa$OGT[aa$variable %in% colnames(MeansArray.Prot)[c(3:10)]] <- 'wt'
aa$OGT[aa$variable %in% colnames(MeansArray.Prot)[c(11:16,18)]] <- 'control'
aa$OGT[aa$variable %in% colnames(MeansArray.Prot)[c(17,19,20)]] <- '5N5A'
aa$OGT <- factor(aa$OGT, levels = c('5N5A','control','wt'))
aa$variable <- factor(aa$variable, levels = colnames(MeansArray.Prot)[c(11:16,18,3:10,17,19,20)])


#S8A: violin plot
ggplot(aa, aes(y=value, x = variable, fill = OGT))+geom_violin()+scale_y_log10() + 
  geom_boxplot(fill = 'white', alpha = 0.8, width = 0.3, outlier.alpha = 0)+
  xlab('')+ylab('')+theme_bw() + geom_hline(yintercept = median(aa$value[aa$OGT=='control']),
                                            color = 'darkgreen', linetype = 'dashed')+
  geom_hline(yintercept = median(aa$value[aa$OGT=='wt']), color = 'blue', linetype = 'dashed') +
  geom_hline(yintercept = median(aa$value[aa$OGT=='5N5A']), color = 'red', linetype = 'dashed') +
  theme(axis.text.x = element_blank(),panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank())

#S8B: make medians histogram pre-median normalization

bb <- melt(FinalResults.sum[,c('wt Median','5N5A Median','control Median','Accession')],
           id.vars = 'Accession')
bb$variable <- factor(bb$variable, levels = c('5N5A Median', 'control Median','wt Median'))
ggplot(bb, aes(x = value, fill = variable)) + geom_histogram(alpha = 0.5, position = 'identity', bins = 100) + 
  scale_x_log10()+ theme_bw()

#S8C
qqnorm(log10(FinalResults.sum$`wt Median`/FinalResults.sum$`control Median`), axes = FALSE,frame.plot = TRUE)
qqline(log10(FinalResults.sum$`wt Median`/FinalResults.sum$`control Median`))
axis(2, at = c(-1,0,1,2,log10(30),log10(3),log10(0.3)),labels = c(0.1,1,10,100,30,3,0.3))
sdwt <- sd(log10(FinalResults.sum$`wt Median`/FinalResults.sum$`control Median`)) #0.2836255
meanwt <- mean(log10(FinalResults.sum$`wt Median`/FinalResults.sum$`control Median`)) #0.499797
breakpoints <- c((1-meanwt)/sdwt,(2-meanwt)/sdwt,(0-meanwt)/sdwt,
                 (-1-meanwt)/sdwt,(log10(30)-meanwt)/sdwt,(log10(3)-meanwt)/sdwt,
                 (log10(0.3)-meanwt)/sdwt)
axis(1, at = breakpoints, labels = c(10,100,1,0.1,30,3,0.3))
qqnorm(scale(log10(FinalResults.sum$`wt Median`)))
qqline(scale(log10(FinalResults.sum$`wt Median`)))
abline(v = scale(log10(FinalResults.sum$`wt Median`))[FinalResults.sum$`wt Hit`],lty = 2, col = 'red')

#S8D
bb <- melt(FinalResults.medCorrect, id.vars = c("Accession MedianCorrected",
                                                "Description MedianCorrected",
                                                "HGNC MedianCorrected", "Known.GlcNAc MedianCorrected", 
                                                "PhosphoSitePlus MedianCorrected", "Wang MedianCorrected"))
bb$OGT <- ''
bb$OGT[bb$variable %in% colnames(FinalResults.medCorrect)[c(14:21)]] <- 'wt'
bb$OGT[bb$variable %in% colnames(FinalResults.medCorrect)[c(7:13)]] <- 'control'
bb$OGT[bb$variable %in% colnames(FinalResults.medCorrect)[c(22:24)]] <- '5N5A'
bb$OGT <- factor(bb$OGT, levels = c('5N5A','wt','control'))
ggplot(bb, aes(y=value, x = variable, fill = OGT))+geom_violin()+scale_y_log10() + 
  geom_boxplot(fill = 'white', alpha = 0.8, width = 0.3, outlier.alpha = 0)+
  xlab('')+ylab('')+theme_bw() + geom_hline(yintercept = median(bb$value[aa$OGT=='control']),
                                            color = 'darkgreen', linetype = 'dashed')+
  geom_hline(yintercept = median(bb$value[aa$OGT=='wt']), color = 'blue', linetype = 'dashed') +
  geom_hline(yintercept = median(bb$value[aa$OGT=='5N5A']), color = 'red', linetype = 'dashed') +
  theme(axis.text.x = element_blank(),panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank())

#S9


#S9A
vennDiagram(byHGNC.all[,c('wt Hit', 'wt v. 5N5A Hit', 'Known.GlcNAc')])
sqrt(157+16+66+500)*2/pi#17.30622, diameter wt Hit
sqrt(222+16+66)*2/pi #11.09985 diameter known GlcNAc
sqrt(66+500)*2/pi #15.14566

#S9B
aa <- melt(FinalResults.sum[,c('Accession','5N5A fold signal','wt fold signal')], id.vars = 'Accession')
aa$Accession <- factor(aa$Accession, levels = FinalResults.sum$Accession[order(FinalResults.sum[,'wt fold signal'])])
aa$variable <- factor(aa$variable, levels = c('5N5A fold signal', 'wt fold signal'))
ggplot(aa, aes(x = Accession, y = value, color = variable)) + 
  geom_point(alpha = 0.2, show.legend = FALSE) +
  theme_bw()  + scale_y_log10()  + ylab('') + xlab('') + 
  theme(axis.text.x = element_blank(),panel.grid.major.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank()) # waterfall plot
ggplot(aa, aes(x = value, fill = variable, color = variable)) + 
  geom_density(alpha = 0.5, show.legend = FALSE) +
  theme_bw()  + scale_x_log10()  + ylab('') + xlab('') 

#S9C

median_correctedVolcano <- FinalResults.sum[,c('Accession','Known.GlcNAc','HGNC',
                                               '5N5A fold signal (Median Corrected)',
                                               '5N5A glycosylation qvalue')]
colnames(median_correctedVolcano)[4:5] <- c('FoldSignal','qvalue')
median_correctedVolcano$OGT <- '5N5A'
median_correctedVolcano2 <- FinalResults.sum[,c('Accession','Known.GlcNAc','HGNC',
                                                'wt fold signal (Median Corrected)',
                                                'wt glycosylation qvalue')]
colnames(median_correctedVolcano2)[4:5] <- c('FoldSignal','qvalue')
median_correctedVolcano2$OGT <- 'wt'
median_correctedVolcano <- rbind(median_correctedVolcano,median_correctedVolcano2)
median_correctedVolcano$pQval <- -log10(median_correctedVolcano$qvalue)
ggplot(median_correctedVolcano, aes(x = FoldSignal, y = pQval, color = OGT)) + 
  geom_point(alpha = 0.5)+theme_bw()+scale_x_log10()+
  geom_hline(yintercept = -log10(0.05), col = 'black', linetype= 'dashed')+
  geom_vline(xintercept = 2, col = 'black', linetype= 'dashed')

#S9D
fortHist <- data.frame('wtSignal' = eBayesResultCorr$t[,1],
                       '5N5ASignal' = eBayesResultCorr$t[,2],
                       'wtv5N5A' = eBayesResult$t[,3])
cc <- melt(fortHist)
ggplot(cc, aes(fill = variable, color = variable, x = value)) + 
  geom_histogram(alpha = 0.3, bins = 50, position = 'identity')+
  theme_bw()

#S9E
Volcano_wtv5N5A <- FinalResults.sum[,c('Accession','Known.GlcNAc','HGNC',
                                       'wt vs. 5N5A fold signal',
                                       'wt vs. 5N5A qvalue')]
colnames(Volcano_wtv5N5A)[4:5] <- c('FoldSignal','qvalue')
Volcano_wtv5N5A$pQval <- -log10(Volcano_wtv5N5A$qvalue)
ggplot(Volcano_wtv5N5A, aes(x = FoldSignal, y = pQval)) + 
  geom_point(alpha = 0.3)+theme_bw()+scale_x_log10()+
  geom_hline(yintercept = -log10(0.05), col = 'red', linetype= 'dashed')+
  geom_vline(xintercept = c(0.5,2), col = 'red', linetype= 'dashed') 

#########
#S10

ggplot(FinalResults.sum[FinalResults.sum$`wt Hit`,], aes(x = `wt vs. 5N5A fold signal`)) +
  geom_histogram(bins = 50)+scale_x_log10()+
  geom_vline(xintercept = FinalResults.sum[FinalResults.sum$HGNC=='TAB1','wt vs. 5N5A fold signal'][1],
             color='red') + theme_bw()
geom_vline(xintercept = median(FinalResults.sum[FinalResults.sum$`wt Hit`,'wt vs. 5N5A fold signal']), 
           color='green', linetype = 'dashed') +theme_bw()


####################
#Table S4
####################
write.csv(FinalResults.sum, row.names = FALSE, file = 'TableS4.csv')