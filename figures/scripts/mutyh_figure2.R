# Robinson et al (2021) "Inherited MUTYH mutations cause elevated somatic mutation rates and distinctive mutational signatures in normal human cells"

#############################
# CODE TO RECREATE FIGURE 2 #
#############################

options(stringsAsFactors = FALSE)
library(data.table)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(deconstructSigs)

###### FIGURE 2a #######
## Plot COSMIC reference spectra
# COSMIC referencd signatures (accessed october 2020)
sigs <- read.delim("COSMIC_referencesigs_v3.1.txt",header=T)

## Plot SBSOGG1 reference spectrum
#Supplementary File from Zou et al 2021 Nature Cancer
zou <- as.data.frame(readxl::read_xlsx("43018_2021_200_MOESM5_ESM.xlsx",sheet=1))
#Normalise mutatio counts to give a proportion
zou_norm <- zou$OGG1/sum(zou$OGG1)
sum(zou_norm)==1
plot_standard_SNV_spectrums("ZOU_OGG1", as.data.frame(round(zou_norm*10000)))

#Initialise plotting function
plot_mini_spectrum = function(vector,plotname,height,width){
  col_vec_num <- rep(16,6)
  sig_cat = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sig_cat,each=16),rep(ctx_vec,times=6),sep=",")
  snv_context = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  col_vec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  
  maxy=max(vector)
  
  pdf(paste0(plotname, "_sbs_mini_profile.pdf"),height=height,width=width)
  b=barplot(vector,col=col_vec,ylim=c(0,maxy*1.7),yaxt='n')
  
  for (j in 1:length(sig_cat)) {
    xpos = b[c(sum(col_vec_num[1:j])-col_vec_num[j]+1,sum(col_vec_num[1:j]))]
    rect(xpos[1]-0.5, maxy*1.4, xpos[2]+0.5, maxy*1.5, border=NA, col=unique(col_vec)[j])
    #text(x=mean(xpos), pos=3, y=maxy*1.55, label=sig_cat[j], cex = 1)
  } 
  dev.off()
}

plot_mini_spectrum(zou_norm,"SBS_ZOU",height=3,width=6)
plot_mini_spectrum(sigs[trinuc,"SBS18"],"COSMIC_SBS18",height=3,width=6)
plot_mini_spectrum(sigs[trinuc,"SBS36"],"COSMIC_SBS36",height=3,width=6)

#HDP OGG1/MUTYH
sigcomps <- read.delim("HDP_signatures_2021-05-11.txt", header=T)

# One plot per HDP signature component N1:N3
plot_mini_spectrum(sigcomps[,"N1"],"HDP_N1",height=3,width=6)
plot_mini_spectrum(sigcomps[,"N2"],"HDP_N2",height=3,width=6)
plot_mini_spectrum(sigcomps[,"N3"],"HDP_N3",height=3,width=6)


###### FIGURE 2b ####### 
#Load mutation matrix - 96 trinucleotide contexts - each row corresponds to one sample
mutations_file = "/lustre/scratch116/casm/cgp/users/pr10/sigprofiler/vcfs/map/matrix/20210810_persample_normals/output/SBS/pass_vcf.SBS96.all"
mut_count.raw <- read.delim(mutations_file,header = T,sep = "\t")

#Reformat the order of the sigprofiler generated matrix
trinuc_pos_three <- rep(c("A","C","G","T"), times = 24)
trinuc_pos_two.2 <- rep(c("C","T"), each = 48)
trinuc_pos_two.1 <- rep(c("A","G", "T","A", "C", "G"), each = 16)
trinuc_pos_one <- rep(rep(c("A","C","G","T"), each = 4), times = 6)
neworder <-  as.factor(paste0(trinuc_pos_one,"[",trinuc_pos_two.2,">",trinuc_pos_two.1,"]",trinuc_pos_three))
rownames(mut_count.raw) <- mut_count.raw[,1]
mutations <- t(mut_count.raw[neworder,][,-1])
rownames(mutations)=gsub("_normal_sample_snp","",rownames(mutations))

#Aggregate by germline mutation
hom104=colSums(mutations[c("PD44888","PD44889"),])
hom286=mutations["PD50743",]
hom179=colSums(mutations[c("PD50745","PD50746","PD50747"),])
chet=colSums(mutations[c("PD44887","PD50744","PD44891"),])
PD44890=mutations["PD44890",]

mutagg <- rbind(hom104,hom286,hom179,chet,PD44890)

# Plotting function - per-patient aggregated spectra
plot_mini_spectrum = function(vector,plotname,height,width){
  col_vec_num <- rep(16,6)
  sig_cat = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sig_cat,each=16),rep(ctx_vec,times=6),sep=",")
  snv_context = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  col_vec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  
  maxy=max(vector)
  
  pdf(paste0(plotname, "_sbs_mini_profile.pdf"),height=height,width=width)
  b=barplot(vector,col=col_vec,ylim=c(0,maxy*1.7),yaxt='n',names.arg="")
  
  for (j in 1:length(sig_cat)) {
    xpos = b[c(sum(col_vec_num[1:j])-col_vec_num[j]+1,sum(col_vec_num[1:j]))]
    rect(xpos[1]-0.5, maxy*1.4, xpos[2]+0.5, maxy*1.5, border=NA, col=unique(col_vec)[j])
    #text(x=mean(xpos), pos=3, y=maxy*1.55, label=sig_cat[j], cex = 1)
  } 
  dev.off()
}
for (type in rownames(mutagg)){
  print(type)
  plot_mini_spectrum(mutagg[type,]/sum(mutagg[type,]),paste0(type,"_NORM_AGGSPECT_REL"),height=3,width=6)
}

# Plotting function - per-patient aggregated spectra
plot_mini_spectrum_scale = function(vector,plotname,height,width){
  col_vec_num <- rep(16,6)
  sig_cat = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sig_cat,each=16),rep(ctx_vec,times=6),sep=",")
  snv_context = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  col_vec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  
  maxy=max(vector)
  
  pdf(paste0(plotname, "_sbs_mini_profile.pdf"),height=height,width=width)
  b=barplot(vector,col=col_vec,ylim=c(0,maxy*1.7),yaxt='n',names.arg="")
  
  for (j in 1:length(sig_cat)) {
    xpos = b[c(sum(col_vec_num[1:j])-col_vec_num[j]+1,sum(col_vec_num[1:j]))]
    rect(xpos[1]-0.5, maxy*1.4, xpos[2]+0.5, maxy*1.5, border=NA, col=unique(col_vec)[j])
    #text(x=mean(xpos), pos=3, y=maxy*1.55, label=sig_cat[j], cex = 1)
  }   
  top <- max(round(mutagg[type,],-3))
  axis(2,at = seq(0,top,top/2),cex.axis=0.8)
  dev.off()
}
for (type in rownames(mutagg)){
  print(type)
  plot_mini_spectrum_scale(mutagg[type,],paste0(type,"_NORM_AGGSPECT_ABS"),height=3,width=6)
}
