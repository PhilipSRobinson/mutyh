# Robinson et al (2021) "Inherited MUTYH mutations cause elevated somatic mutation rates and distinctive mutational signatures in normal human cells"

#############################
# CODE TO RECREATE FIGURE 1 #
#############################

require(RColorBrewer)
require(dplyr)
require(readxl)
options(stringsAsFactors = F)

dft <- read_xlsx("Extended_Data_Table1.xlsx",sheet=1)
## Select microdissected samples
df <- dft[dft$sequencing=="low-input/LCM",]
## Read in mutation burdens from normal health controls Lee-Six et al 2019
norm1 <- read.delim("nomal_control_sbsrates.txt",header=T)
dfncomb <- rbind(df[df$type=="normal",c("sample", "patient", "sbstotal_corr","germline_mutation", "age")],norm1[,c("sample", "patient", "sbstotal_corr","germline_mutation", "age")])
## Arrange by germline mutation and order by age
dfncomb<- dfncomb[order(dfncomb$germline_mutation, rank(-dfncomb$age), decreasing=T),]
#Convert patient to factor to retain order of the boxplots
dfncomb$patient <- factor(dfncomb$patient , levels=unique(dfncomb$patient))
## Arrange by germline mutation and order by age
df<- df[order(df$germline_mutation, rank(-df$age), decreasing=T),]
#Convert patient to factor to retain order of the boxplots
df$patient <- factor(df$patient , levels=unique(df$patient))
#Set up colour palette
coul <- colorRampPalette(brewer.pal(8, "Spectral"))(13)
coul <- coul[c(3,6,10,13)]
df$colour[df$germline_mutation == "Y179C +/- G396D +/-"] = coul[1]
df$colour[df$germline_mutation == "Y179C -/-"] = coul[2]
df$colour[df$germline_mutation == "Y104* -/-"] = coul[3]
df$colour[df$germline_mutation == "G286E -/-"] = coul[4]

## Normal ##
dfn <- df[df$type == "normal",]
## Polyp ##
dfp <- df[df$type == "adenoma",]

meta <- unique(df[,c("patient","gender","age", "colour")])
meta$gender = ifelse(meta$gender == "Female", "F", "M")
labels = paste0(meta$patient,"\n", meta$age, meta$gender)
boxcol = meta$colour
dotcol = c(rep(coul[1],4),rep(coul[2],3),rep(coul[3],2),rep(coul[4],1)) 
atvec = c(1.5, 2.5, 3.5, 4.5,6, 7,8,9.5,10.5,12,13.5)
length(atvec)=length(boxcol)=length(dotcol)
boxwidth=rep(0.7,5)

dfn$sbsrate = dfn$sbstotal_corr / dfn$age

######## FIGURE 1A ########

pdf("plots/sbstotal_boxplot.pdf",width=10,height=7)
## COMPOUND HETEROZYGOTE
boxplot(sbstotal_corr ~ patient,data=dfncomb[dfncomb$patient %in% c("PD44890","PD44887","PD50744","PD44891"),], col=coul[1],names=labels,las=1, outline=F,ylim = c(-2500, 35000),xlab="",ylab="", frame.plot=F, xaxt='n', at=atvec, boxwex=boxwidth)
stripchart(sbstotal_corr ~ patient, vertical = TRUE, data = dfncomb[dfncomb$patient %in% c("PD44887","PD44890","PD50744","PD44891"),], 
           method = "jitter", add = TRUE, pch = 21, bg = coul[1], col="black", at=atvec)
## HOMOZYGOUS Y179C
par(new=TRUE) # Ensures the plots are overlaid - circumvents issues with colour of dot in stripcharts
boxplot(sbstotal_corr ~ patient,data=dfncomb[dfncomb$patient %in% c("PD50746","PD50747","PD50745"),], col=coul[2],las=1, outline=F, names=labels[5:7],ylim = c(-2500, 35000),xlab="",ylab="",frame.plot=F,xaxt='n',yaxt='n', at=atvec, boxwex=boxwidth)
stripchart(sbstotal_corr ~ patient, vertical = TRUE, data = dfncomb[dfncomb$patient %in% c("PD50746","PD50747","PD50745"),], 
           method = "jitter", add = TRUE, pch = 21, bg = coul[2], col="black", at=atvec)
abline(h=0, lty=1, col="black", lwd=1)
## HOMOZYGOUS TRUNCATING 
par(new=TRUE) # Ensures the plots are overlaid - circumvents issues with colour of dot in stripcharts
boxplot(sbstotal_corr ~ patient,data=dfncomb[dfncomb$patient %in% c("PD44888","PD44889"),], col=coul[3],las=1, outline=F, names=labels[8:9],ylim = c(-2500, 35000),xlab="",ylab="",frame.plot=F,xaxt='n',yaxt='n', at=atvec, boxwex=boxwidth)
stripchart(sbstotal_corr ~ patient, vertical = TRUE, data = dfncomb[dfncomb$patient %in% c("PD44888","PD44889"),], 
           method = "jitter", add = TRUE, pch = 21, bg = coul[3], col="black", at=atvec)
abline(h=0, lty=1, col="black", lwd=1)
## HOMOZYGOUS MISSENSE
par(new=TRUE) # Ensures the plots are overlaid - circumvents issues with colour of dot in stripcharts
boxplot(sbstotal_corr ~ patient,data=dfncomb[dfncomb$patient %in% c("PD50743"),], col=coul[4],las=1, outline=F, names=labels[10],ylim = c(-2500, 35000),xlab="",ylab="",frame.plot=F,xaxt='n',yaxt='n', at=atvec, boxwex=boxwidth)
stripchart(sbstotal_corr ~ patient, vertical = TRUE, data = dfncomb[dfncomb$patient %in% c("PD50743"),], 
           method = "jitter", add = TRUE, pch = 21, bg = coul[4], col="black", at=atvec)
abline(h=0, lty=1, col="black", lwd=1)
text(x = atvec,y = -2500,labels = labels,cex = 0.8)
## WILD TYPE CONTROLS > 50 YRS AGE
par(new=TRUE) # Ensures the plots are overlaid - circumvents issues with colour of dot in stripcharts
boxplot(sbstotal_corr ~ patient,data=dfncomb[dfncomb$patient %in% c("control"),], col="lightgrey",las=1, outline=F,ylim = c(-2500, 35000),xlab="",ylab="",names=labels,frame.plot=F,xaxt='n',yaxt='n', at=atvec, boxwex=boxwidth)
stripchart(sbstotal_corr ~ patient,data=dfncomb[dfncomb$patient %in% c("control"),], vertical = TRUE,
           method = "jitter", add = TRUE, pch = 21, bg = "lightgrey", col="black", at=atvec)
abline(h=0, lty=1, col="black", lwd=1)
text(x = atvec,y = -2500,labels = labels,cex = 0.8)
dev.off()

lwd=2

######## FIGURE 1B ######## 
foldint <- read.delim("fold_changes_model_20211122.txt",header=T)
par(oma=c(1,1,1,1))
foldvalues=c(1,2,3,4,5,30,45)
scale=log(foldvalues,100)
lines=log(c(2,3,4,5,30,45),100)
normal_pos=c(0.5,1.5,2.5,3.5,4.5)
lwd=1
pdf("plots/SBS_rate_fold_dotwhisker_genotype.pdf", h=5.5,w=4)
plot(normal_pos,log(foldint[3:7,"fold.change"],100), xlab="", ylab="",
     col=boxcol,lwd=lwd, frame.plot = F, cex=1.5, pch=21,bg=boxcol[c(1,2,5,8,10)], ylim=c(0,1),xlim=c(0,5),yaxt='n', xaxt='n')
axis(2,at=scale,las=2,lwd=lwd, labels=foldvalues,cex=0.8)
segments(normal_pos,log(foldint[3:7,"fold.lower"],100),normal_pos,log(foldint[3:7,"fold.upper"],100),col = boxcol[c(1,2,5,8,10)], lwd=3,cex=5)
abline(h=0,lwd=lwd)
abline(h=lines,lwd=lwd,lty=5,col="lightgrey")
par(new=TRUE)
plot(normal_pos,log(foldint[3:7,"fold.change"],100), xlab="", ylab="",
     col=boxcol[c(1,2,5,8,10)],lwd=1, frame.plot = F, cex=1.5 , pch=21,bg=boxcol[c(1,2,5,8,10)], ylim=c(0,1),xlim=c(0,5),yaxt='n', xaxt='n')
axis(2,at=scale,las=2,lwd=lwd, labels=foldvalues,cex=0.8)
segments(normal_pos,log(foldint[3:7,"fold.lower"],100),normal_pos,log(foldint[3:7,"fold.upper"],100),col = boxcol[c(1,2,5,8,10)], lwd=3,cex=5)
dev.off()

## Adenomas - absolute burden and fold change
dfp$sbstotal_corr_rate=dfp$sbstotal_corr/dfp$age
dfp$patient= as.character(dfp$patient)
psbs <- data.frame(
  patient = aggregate(dfp$sbstotal_corr, by = list(dfp$patient), max)[,1],
  max = aggregate(dfp$sbstotal_corr, by = list(dfp$patient), max)[,2],
  mean = aggregate(dfp$sbstotal_corr, by = list(dfp$patient), mean)[,2],
  min = aggregate(dfp$sbstotal_corr, by = list(dfp$patient), min)[,2])
rownames(psbs)=psbs$patient
psbs <- psbs[c("PD44890","PD44887", "PD50747", "PD44889", "PD44888"),]
lwd=2

dfn$patient= as.character(dfn$patient)
sbsburden <- data.frame(
  patient = aggregate(dfn$sbstotal_corr, by = list(dfn$patient), max)[,1],
  max = aggregate(dfn$sbstotal_corr, by = list(dfn$patient), max)[,2],
  mean = aggregate(dfn$sbstotal_corr, by = list(dfn$patient), mean)[,2],
  min = aggregate(dfn$sbstotal_corr, by = list(dfn$patient), min)[,2])
rownames(sbsburden)=sbsburden$patient
sbsburden <- sbsburden[c("PD44890","PD44887","PD50747","PD44889", "PD44888"),]
lwd=1

######## FIGURE 1C #########

pdf("Adenoma_SBS_abs_dotwhisker.pdf")
normal_cols = "darkgrey"
normal_pos = c(1, 3,5.5, 8, 10)
adenoma_pos = c(1.5, 3.5, 6,8.5 ,10.5)
#boxcolad= boxcol[c(1,2,4:5)]
boxcolad=c("#F57748","#F57748","#FDDB87","#99D6A4","#99D6A4")
plot(adenoma_pos,psbs[,"mean"],xlab="", ylab="", 
     col=boxcolad,lwd=2, frame.plot = F, cex=2, pch=21,bg=boxcolad,axes=F, ylim=c(0,36000), xlim=c(0,11))
axis(2,at=c("0", "5000", "15000", "25000", "35000"),las=2, lwd=lwd)
abline(h=c("5000", "15000", "25000", "35000"),lwd=lwd,lty=5,col="lightgrey")
par(new=TRUE)
plot(adenoma_pos,psbs[,"mean"],xlab="", ylab="", 
     col=boxcolad,lwd=2, frame.plot = F, cex=2, pch=21,bg=boxcolad,axes=F, ylim=c(0,36000), xlim=c(0,11))
axis(2,at=c("0", "5000", "15000", "25000", "35000"),las=2, lwd=lwd)
segments(adenoma_pos,psbs[,"min"],adenoma_pos,psbs[,"max"],col = boxcolad, lwd=5,cex=5)
par(new=TRUE) # Ensures the plots are overlaid - circumvents issues with colour of dot in stripcharts
plot(normal_pos,sbsburden[c("PD44890","PD44887","PD50747","PD44889", "PD44888"),"mean"],xlab="", ylab="", 
     col=normal_cols,lwd=2, frame.plot = F, cex=2, pch=21,bg=normal_cols,axes=F, ylim=c(0,36000),xlim=c(0,11))
segments(normal_pos,sbsburden[c("PD44890","PD44887","PD50747","PD44889", "PD44888"),"min"],normal_pos,sbsburden[c("PD44890","PD44887","PD50747","PD44889", "PD44888"),"max"],col = normal_cols, lwd=5,cex=5)
abline(h=0, lwd=lwd)
dev.off()

