# This is the qiime2 script for trimming primers using Cutapadt
cd /data 

conda activate qiime2-2023.5 

# this imports the Mifish data, make sure to be in the proper directory for each marker
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path  ComBay_run3_18S_fastq \
  --input-format  CasavaOneEightLanelessPerSampleDirFmt \
  --output-path ComBay-run3-18S_demux-paired-end.qza

qiime cutadapt trim-paired \
--p-cores 8 \
--i-demultiplexed-sequences ComBay_run2_demux-paired-end.qza \
--p-adapter-f GTCGGTAAAACTCGTGCCAGC...CAAACTGGGATTAGATACCCCACTATG \
--p-adapter-r CATAGTGGGGTATCTAATCCCAGTTTG...GCTGGCACGAGTTTTACCGAC \
--p-discard-untrimmed \
--o-trimmed-sequences MiFish_ComBay2_trimmed.qza \
--p-error-rate 0.2 \
--verbose &> MiFishprimerCB2_trimed.log 

qiime demux summarize \
  --i-data MiFish_ComBay2_trimmed.qza  \
  --o-visualization MiFish_ComBay2_trimmed.qzv

# forward trim at 180, reverse trim at 180


# go look at BP quality before you go modify the configure.yaml for tourmaline

# 18S now

qiime cutadapt trim-paired \
--p-cores 8 \
--i-demultiplexed-sequences ComBay-run3-18S_demux-paired-end.qza \
--p-adapter-f GGTCTGTGATGCCCTTAGATG...CCCTGCCCTTTGTACACACC \
--p-adapter-r GGTGTGTACAAAGGGCAGGG...CATCTAAGGGCATCACAGACC \
--p-discard-untrimmed \
--o-trimmed-sequences 18s_ComBay3Trimmed.qza \
--p-error-rate 0.2 \
--verbose &> 18Sprimer_CB3_trimed.log 

qiime demux summarize \
  --i-data 18s_ComBay3Trimmed.qza \
  --o-visualization 18s_ComBay3Trimmed.qzv

