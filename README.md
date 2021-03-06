# Genomic Prediction of Tocochromanols in Exotic Maize

This repository contains code and data used in the manuscript "Genomic prediction of tocochromanols in exotic maize" (manuscript in preparation). For more details, please contact ltibbs@iastate.edu.

## Optimal Training Population Design 

The script `PAM_FURS_MaxCD.R` contains code used to apply the optimal training population design methods PAM, FURS, and MaxCD to inbreds (originally applied to hybrids in https://doi.org/10.1016/j.molp.2018.12.022).

## Consensus Sequence

The script `R_consensus.R` is based on the publicly-available Python script https://github.com/mdzievit/Genomic_Prediction/tree/master/ThePlantGenome_Version/Consensus_Script and is used to create a consensus sequence for accessions genotyped multiple times (vcf format, original data from `ZeaGBSv27_publicSamples_imputedV5_AGPv4-181023.vcf` available from Panzea).

