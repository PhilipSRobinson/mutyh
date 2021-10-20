# mutyh

Repository to accompany manuscript Robinson et al. 2021 "Inherited MUTYH mutations cause elevated somatic mutation rates and distinctive mutational signatures in normal human cells".

This repository contains the source data and code to recreate figures in the manuscript and underage statistical modelling. 

DNA sequencing data are deposited in the European Genome-Phenome Archive (EGA) with accession code: EGAD00001007958 and EGAD00001007997.

### Mutation Calling ###
Mutation calling was undertaken using CaVEMan, Pindel, ASCAT, GRIDSS and TraFiC-mem. Further detail is provided in the manuscript. They were run as part of the Sanger pipeline and are all openly available here: https://github.com/cancerit/

### Filtering ###
Mutation filtering for single-base substitution (SBS) and insertion and deletion (ID) mutations intestinal crypts was performed by intitally running the mutation calling unmatched against a synthetic normal BAM file. An allele counter, CgpVAF, was run to create a pile-up of read counts for wach mutation in each fsample per-individual. The first stage of filtering involed running filters written by Mathijs Sanders designed to remove artefacts that are specific to low-input DNA sequencing available here:  https://github.com/MathijsSanders/SangerLCMFiltering 

Thereafter filtering of mutations was performed using binomial and beta-binomial filters to remove germline and artefactual mutations, code available here: https://github.com/TimCoorens/Unmatched_NormSeq

### Tree building ###
Phylogenetic trees were created using mpboot http://www.iqtree.org/mpboot/

Mutations were assigned to branches of the trees using treemut https://github.com/NickWilliamsSanger/treemut

### Signature analysis ###
Software for mutational signature analysis is available from the following urls:

HDP https://github.com/nicolaroberts/hdp 
SigFit https://github.com/kgori/sigfit 
SigProfiler https://github.com/AlexandrovLab. 

### Analysis of nanoseq data ###
Software for analysis of raw duplex / NanoSeq data is provided at https://github.com/cancerit/NanoSeq. 

