# Tourmaline pipeline mofified from tutorial here: https://github.com/aomlomics/tourmaline

# run the docker image
docker pull aomlomics/tourmaline
docker run -v $HOME:/mnt/d/Projects/2022_wdfw_SPSdiet -it aomlomics/tourmaline


# EVERY NEW DATASET (OR NEW PARAMS ON A DATASET) CLONE THE TOURMALINE REPOSITORY AND NAME IT FOR YOUR PROJECT
cd /data # your repository name
git clone https://github.com/aomlomics/tourmaline.git

# rename it for your project (tourmaline-myproject-YYYYMMDD)
mv tourmaline tourmaline-18S_CB3_17122024 # this is an example

cd tourmaline-18S_CB3_17122024

# Store the pre-imported reference FASTA and taxonomy .qza files as 01-imported/refseqs.qza and 01-imported/reftax.qza.
conda activate snakemake

# note: files into tourmaline must be demultiplexed first and also primers trimmed with cutadapt in qiime2 
# Store your pre-imported and trimmed FASTQ .qza files as 01-imported/fastq_pe.qza (paired-end)

# Edit the configuration file config.yaml to set DADA2 parameters for the run

# now run snakemake
snakemake --use-conda dada2_pe_denoise --cores 24

# Pausing after the denoise step allows you to make changes before proceeding:

# Check the table summaries and representative sequence lengths to determine if DADA2 parameters need to be modified. If so, you can rename or delete the output directories and then rerun the denoise rule.
# View the table visualization to decide an appropriate subsampling (rarefaction) depth. Then modify the parameters "alpha_max_depth" and "core_sampling_depth" in config.yaml.
# Decide whether to filter your feature table and representative sequences by taxonomy or feature ID. After the taxonomy step, you can examine the taxonomy summary and bar plot to aid your decision. 
# If you do filter your data, all output from that point on will go in a separate folder so you can compare output with and without filtering

snakemake --use-conda dada2_pe_taxonomy_unfiltered --cores 24

#if there are too many singletons in some datasets it may cause the taxonomy to fail so let's try the filtering options
snakemake --use-conda dada2_pe_taxonomy_filtered --cores 30
