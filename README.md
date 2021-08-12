# mutyh

Repository to accompany manuscript Robinson et al. 2021 Inherited MUTYH mutations cause elevated somatic mutation rates and distinctive mutational signatures in normal human cells.

This repository contains the source data and code to recreate figures in the manuscript and underage statistical modelling. 

Raw sequencing data can be accessed from: 

### Mutation Calling ###
Mutation calling was undertaken using CaVEMan, Pindel, ASCAT, GRIDSS and TraFiC-mem. Further detail is provided int he manuscript. They were run as part of the Sanger pipeline and are all openly available here: https://github.com/cancerit/

### Filtering ###
Mutation filtering for single-base substitution (SBS) and insertion and deletion (ID) mutations intestinal crypts was performed by intitally running the mutation calling unmatched against a synthetic normal BAM file. The first stage of filtering involed running filters written by Mathijs Sanders designed to remove artefacts that are specific to low-input DNA sequencing available here:  https://github.com/MathijsSanders/SangerLCMFiltering 

Thereafter filtering of mutations was performed using binomial and beta-binomial filters to remove germline and artefactual mutations, code available here: https://github.com/TimCoorens/Unmatched_NormSeq

### Tree building ###
Phylogenetic trees were created and mutations were assigned using scripts in this repository.

### Signature analysis ###


