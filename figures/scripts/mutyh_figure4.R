# Robinson et al (2021) "Inherited MUTYH mutations cause elevated somatic mutation rates and distinctive mutational signatures in normal human cells"

#############################
# CODE TO RECREATE FIGURE 4 #
#############################

options(stringsAsFactors = F)
suppressWarnings(library(dplyr))
suppressWarnings(library(RColorBrewer))
suppressWarnings(library(data.table))
suppressWarnings(library(readxl))
suppressWarnings(library(nlme))
suppressWarnings(library(MuMIn))
options(stringsAsFactors = F)

#Read in intestinal crypt mutation burden
a <- read_xlsx("Extended_Data_Table2.xlsx",sheet=5)
b <- a[a$type=="normal",] #Select normal crypts
c <- b %>% dplyr::group_by(patient, germline_mutation, age) %>% dplyr::summarise(mean = mean(sbstotal_corr),
                                                                                 upper = max(sbstotal_corr),
                                                                                 lower = min(sbstotal_corr))
#Read in nanoseq data
df <- read.delim("nanoseq_burdens_all.txt",header=T)
d <- read.delim("nanoseqblood_intestine_burdens.txt",header=T)
df$cohort=ifelse(df$germline_mutation=="wt","wt","mutyh")


####### FIGURE 4A #########
# Load Abascal et al 2021 - normal granulocyte NanoSeq mutation burdens 
genomesize = 5722652910 # callable genomesize
nsd <- as.data.frame(read_xlsx("41586_2021_3477_MOESM3_ESM.xlsx",sheet=4))
nsd1 <- nsd[nsd$cell_type=="grans" &grepl("grans",nsd$sampleID),c("sampleID","cell_type","burden", "age")]
nsd1$sbs_percell = nsd1$burden*genomesize
nsd1$burden=NULL
colnames(nsd1)=c("sample","type","age","sbs_percell")
dfcomb <- rbind(df[df$type=="blood",c("sample","type","sbs_percell","age")],nsd1[,c("sample","type","sbs_percell","age")])
dfcomb$patient=substr(df$sample,1,7)
dfcomb$cohort=ifelse(grepl("ds",dfcomb$sample),"mutyh", "wt")

# Mixed effects model to assess rate of SBS mutations in white blood cells from peripheral bulk blood
model.age.blood <- lme(fixed = sbs_percell ~ age:cohort, #Interaction term to fit separate slopes for 1.wt and 2.mutyh
                       random = list(patient = pdSymm(form = cohort ~ 1)), 
                       weights = varIdent(form= ~ 1 | cohort),
                       data = dfcomb[!dfcomb$patient=="PD44890",], method="ML") #PD44890 excluded as extreme outlier
r.squaredGLMM(model.age.blood)
intervals(model.age.blood, which="fixed")
#Precise p-value estimates
s$tTable[,"p-value"]

lwd=2
pdf("SBS_bloodonly_dotwhisker.pdf", height=4,width=4)
plot(df$age[df$type=="blood"],df$sbs_percell[df$type=="blood"],xlab="", ylab="", 
     col=df$colour[df$type=="blood"],lwd=2, frame.plot = F, cex=1.5, pch=21,bg=df$colour[df$type=="blood"],axes=F, ylim=c(0,5200), xlim=c(0,80))
segments(df$age[df$type=="blood"],df$sbs_percell_lci[df$type=="blood"],df$age[df$type=="blood"],df$sbs_percell_uci[df$type=="blood"],col = df$colour[df$type=="blood"], lwd=5,cex=5)
par(new=T)
axis(2,at=c("0","1000", "2000", "3000", "4000", "5000"),las=2, lwd=lwd)
axis(1,at=c("0","20", "40", "60", "80"),las=1, lwd=lwd)
abline(a=model.age.lymph$coefficients$fixed[1], b=model.age.lymph$coefficients$fixed[3],lwd=lwd,lty=5)
abline(a=model.age.lymph$coefficients$fixed[1], b=model.age.lymph$coefficients$fixed[2],lwd=lwd,lty=3)
dev.off()

####### FIGURE 4B #########
df2 <- read.delim("NanoSeq_blood_colon_mutyhsigsburden.txt",header=T)

# Linear model to assess relationship between blood and colon MUTYH-related mutation burdens 
fit <- lm(colonsigs~bloodsigs-1,data = df2[!df2$patient=="PD44890",])
coint <- confint(fit)

# Plot of MUTYH mutation burdens (SBS18 and SBS36) in blood vs. colon
pdf("nanoseq_blood_colon_mutyhsigs_all_fitmiunusD44890.pdf",height=4,width=4)
plot(df2$bloodsigs,df2$colonsigs, pch=19,col=df2$colour,ylim=c(0,25000), xlim=c(0,5000), axes=F, frame.plot = F, xlab = "",ylab="",cex=1.5)
axis(1,at=c("0","1000", "2000","3000","4000","5000"), lwd=lwd)
axis(2,at=c("0","5000", "10000", "15000", "20000", "25000"),las=2, lwd=lwd)
abline(fit,lty=1)
abline(a=1,b=coint[,1],lty=3)    
abline(a=1,b=coint[,2],lty=3)
dev.off()


####### FIGURE 4C #########

#Read in NanoSeq trinucleotide mutation counts
counts <- read.delim("NanoSeq_trinuc_counts.txt",header=T)
#Read in nanoseq normal control granulocytes from Abascal et al 2021
ns  <- as.data.frame(read_xlsx("41586_2021_3477_MOESM3_ESM.xlsx",sheet=7))
ns1 <- ns[grepl("_grans_nanoseqv2",ns$SampleID) | ns$SampleID == "PD43976_grans_nanoseq",]
rownames(ns1)=ns1$SampleID
ns1$SampleID=NULL
colnames(ns1)=colnames(counts)
counts <- rbind(counts,ns1)

#Setup Colours
sig_order=1:5
names(sig_order)=c("SBS1","SBS5","SBS18","SBS36","SBS9")
all_cols=c("#9E0142","#E95D47","#54AEAC","#BEE5A0", "lightblue")
names(all_cols)=names(sig_order)

#Load signatures - reference from cosmic but SBS36 replaced with extracted version due to reduced contamination.
final_sigs_reference_new <- read.table("sigs_used_treefit.txt")

blood=c("PD43979b_grans_nanoseqv2",
        "PD43980b_grans_nanoseqv2", 
        "PD43987b_grans_nanoseqv2",
        "PD43981bR2_grans_nanoseqv2", 
        "PD43976_grans_nanoseq",
        "PD43982bR2_grans_nanoseqv2", 
        "PD43988b_grans_nanoseqv2", 
        "PD43983b_grans_nanoseqv2",
        "PD43984b_grans_nanoseqv2",
        "PD44890i_ds0001",
        "PD44887i_ds0001",
        "PD50744c_ds0001",
        "PD44891c_ds0001",
        "PD50746c_ds0001",
        "PD50745d_ds0001",
        "PD50747b_ds0001", 
        "PD44889e_ds0001",
        "PD44888e_ds0001",
        "PD50743c_ds0001") 
sigs_select=c("SBS1", "SBS5", "SBS18", "SBS36")
fit=fit_signatures(counts=counts[blood,], 
                   signatures = final_sigs_reference_new[sigs_select,],
                   iter = 20000,
                   warmup = 10000,
                   model="poisson",
                   chains = 2)

pars=retrieve_pars(fit, 
                   par = "exposures", 
                   hpd_prob = 0.95)

data_percentage <- apply(t(pars$mean), 2, function(x){x*100/sum(x,na.rm=T)})
data_percentage <- data_percentage/100
colnames(data_percentage) = blood
dp_reorder <- data_percentage[,blood]

pdf("nanoseq_batch_signature_contribution.pdf",height=3.5,width=8.5)
barplot(dp_reorder, col=all_cols[c("SBS1", "SBS5", "SBS18", "SBS36")] ,xaxt='n', yaxt='n', space = c(rep(0.1,9),0.5,rep(0.1,9)))
axis(2,at=c(0,0.2,0.4,0.6,0.8, 1),las=1,lwd=2)
dev.off()

# Test the significance of difference in proportion of MUTYH signatures proportion in blood from wild-type vs MAP
wilcox.test(colSums(dp_reorder[c("SBS18","SBS36"),1:9]),colSums(dp_reorder[c("SBS18","SBS36"),10:19]))



####### FIGURE 4D #########

# Load Abascal et al data
nsd <- as.data.frame(read_xlsx("41586_2021_3477_MOESM3_ESM.xlsx",sheet=4))
nsd1 <- nsd[nsd$cell_type=="grans" &grepl("grans",nsd$sampleID),c("sampleID","cell_type","burden", "age")]
nsd1$sbs_percell = nsd1$burden*genomesize
nsd1$burden=NULL
colnames(nsd1)=c("sample","type","age","sbs_percell")
dfcomb <- rbind(df[df$type=="blood",c("sample","type","sbs_percell","age")],nsd1[,c("sample","type","sbs_percell","age")])
dfcomb$patient=substr(df$sample,1,7)
dfcomb$cohort=ifelse(grepl("ds",dfcomb$sample),"mutyh", "wt")


# BLOOD LINEAR MIXED EFFECT MODEL
model.age.lymph <- lme(fixed = sbs_percell ~ age:cohort, # Once again interaction term used to reflect the continuous accumualtion of mutations
                       random = list(patient = pdSymm(form = ~ 1)),
                       weights = varIdent(form= ~ 1 | cohort),
                       data = df[df$type=="lymphocytes"& !df$patient=="PD44890" ,], method="ML")
r.squaredGLMM(model.age.lymph)
intervals(model.age.lymph, which="fixed")

#Precise p-value estimates
s$tTable[,"p-value"]

pdf("SBS_lymphocytes_minusPD44890_dotwhisker.pdf", height=4.5,width=4)
plot(df$age[df$type=="lymphocytes"],df$sbs_percell[df$type=="lymphocytes"],xlab="", ylab="", 
     col=df$colour[df$type=="lymphocytes"],lwd=2, frame.plot = F, cex=1.5, pch=21,bg=df$colour[df$type=="lymphocytes"],axes=F, ylim=c(0,6000), xlim=c(0,80))
segments(df$age[df$type=="lymphocytes"],df$sbs_percell_lci[df$type=="lymphocytes"],df$age[df$type=="lymphocytes"],df$sbs_percell_uci[df$type=="lymphocytes"],col = df$colour[df$type=="lymphocytes"], lwd=5,cex=5)
axis(2,at=c("0","1000", "2000", "3000", "4000", "5000", "6000"),las=2, lwd=lwd)
axis(1,at=c("0","20", "40", "60", "80"),las=1, lwd=lwd)
abline(a=model.age.lymph$coefficients$fixed[1], b=model.age.lymph$coefficients$fixed[2], lwd=lwd, lty=3)
abline(a=model.age.lymph$coefficients$fixed[1], b=model.age.lymph$coefficients$fixed[3], lwd=lwd, lty=5)
dev.off()

####### FIGURE 4E ########

lymphocytes=c("PD45782c_ds0001","PD42834c_ds0001","PD44886b_ds0001","PD42835c_ds0001","PD44881b_ds0001","PD40840c_ds0001", "PD41851b_ds0001", "PD44890e_ds0001", "PD44891b_ds0001", "PD44891b_ds0002","PD44889b_ds0001","PD44888d_ds0001","PD44888d_ds0002") 
sigs_select=c("SBS1", "SBS5", "SBS18", "SBS36", "SBS9")
fit=fit_signatures(counts=counts[lymphocytes,], 
                   signatures = final_sigs_reference_new[sigs_select,],
                   iter = 20000, 
                   warmup = 10000,
                   model="poisson",
                   chains = 2)

pars=retrieve_pars(fit, par = "exposures", hpd_prob = 0.95)

data_percentage <- apply(t(pars$mean), 2, function(x){x*100/sum(x,na.rm=T)})
data_percentage <- data_percentage/100
colnames(data_percentage) = lymphocytes
dp_reorder_l <- data_percentage[,lymphocytes]

pdf("nanoseq_batch_signature_contribution_lymphocytes_wPD44890.pdf",height=3.5,width=6)
barplot(dp_reorder_l, col=all_cols ,xaxt='n', yaxt='n', space = c(rep(0.1,7),0.5,rep(0.1,5)))
axis(2,at=c(0,0.2,0.4,0.6,0.8, 1),las=1,lwd=2)
dev.off()

# Test the significance of difference in proportion of MUTYH signatures in blood from wild-type vs MAP
wilcox.test(colSums(dp_reorder_l[c("SBS18","SBS36"),1:6]),colSums(dp_reorder_l[c("SBS18","SBS36"),7:12]))
