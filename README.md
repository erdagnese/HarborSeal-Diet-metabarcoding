# HarborSeal-Diet-metabarcoding


This is a repository for the scripts and example inputs for running the harbor seal diet metabarcoding analysis post sequencing

It includes the scripts for:

1. Demultiplexing using cutadapt to separate the MiFish and 18S reads for the same samples, while also removing the primers.
2. Then guidelines for setting up running the tourmaline pipeline for each marker's sequences.
3. The post trourmaline decontamination of reads
4. The determination of final taxonomic assignment of decontaminated reads
5. The scripts needed to run the model to correct amplification bias
6. The scripts needed to correct the tissue DNA density bias
7. The scripts to correct the digestion bias using controlled feeding study scats
8. Finally a comparison of the proportion of diet changes throughout each step of bias correction
9. Scripts for visualizing data
