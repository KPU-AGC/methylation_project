# Methylation Project Data Analysis Protocol
## Background
The NGS pipeline for processing targeted amplicon bisulfite sequencing data is based on `Bismark`, a commonly used methylation mapping and calling package written in Perl. This protocol assumes that you are familiar with `Python`, `Conda`, and Unix command-line (`Bash`). Ensure that you have set up the `conda` environment necessary for running the pipeline. The scripts and documentation for the pipeline can be found in the lab [Github](https://github.com/KPU-AGC/methylation_project). 

## Outline of Protocol
This outline is laid out in roughly the same order that the each procedure should be performed in. 

### Setting up environment and retrieving necessary data files
	
1.  Set up `conda` environment
2.  Downloading and generating the reference genome files for Bismark mapper
    * Downloading reference genome from NCBI
    * Generating variant call file for genome masking
    * Generating a masked genome using `bedtools maskfasta`
    * Generating a bisulfite-converted reference genome from standard reference genome or masked genome
3. Generating a primer.bed file for analysis scripts
4. Retrieving run files from the Illumina MiSeq

### Processing FASTQ files using the NGS pipeline
1. Setting up the directory tree for analysis
2. Running `ngs_process.sh`
    * Running individual scripts instead of the entire pipeline
    * `SNPSplit` for separating reads by genotyping data

### QC and processing Bismark outputs
1. Generating QC results for the entire run
    * `process-bismark-reports.py`
    * `process-bismark-incomplete-reports.py`
2. `process-bismark-coverage.py` - generating summary files and human readable .csv files
3. `cgvisualizer.py` - plot methylation calls for each target site for each sample
4. `assess-by-primer.py` - assess coverage for target sites across all samples

## Setting up environment and retrieving necessary data files
### Set up conda environment
For more information, refer to the [quick conda guide](link). A `conda` environment file can be found in the repository and used to generate the environment for all data analysis related to this project. Run the following command: 

```
conda env create --file <methylation_project.yml>
```

To activate the the environment, use `conda activate methylation_project`. Whenever a shell script or Python script is run, this environment must be active.  

### Downloading and generating the reference genome files for Bismark mapper
#### Downloading a reference genome from NCBI
The reference genomes for your organism can be downloaded from NCBI. There are several methods of downloading the genomes, but the most important thing is that the genomes are in the `multi-FASTA` file format. Methods available are: 

1. Using [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/)
2. Accessing the NCBI [FTP](https://ftp.ncbi.nlm.nih.gov/) site

Detailed instructions on using either of those resources can be found [here](link).

#### Generating variant call file for genome masking
TBD - unnecessary atm
#### Generating a masked genome using bedtools maskfasta
TBD - unnecesssary atm
#### Generating a bisulfite-converted reference genome
The method for generating a bisulfite-converted reference genome is the same for the masked and umasked genomes. Ensure that you're using the appropriate genome, and that the compressed (`gzip`) genome FASTA file is in a correctly labeled folder (e.g. arsucd1.2 for the ARS-UCD1.2 Bovine genome). `Bismark` has a tool called `bismark_genome_preparation` that will do an in-silico bisulfite conversion of the supplied reference genome, allowing for mapping to all four of the resulting converted sequences. The instructions for preparing the genome are adapted from [here](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html#i-running-bismark-genome-preparation).

To prepare the genome, run the following command: 

```
bismark_genome_preparation --verbose --bowtie2 </path/to/genome/folder>
```

A folder in the same directory as the genome FASTA file named `Bisulfite_Genome` will be generated containing the converted sequences. 

```
genome_name/
├─ Bisulfite_Genome/
├─ genome.fna.gz
```
When using Bismark, pass `genome_name/` as the argument for the genome path. 

### Generating a primer.bed file for analysis scripts
This file is generated using [`bedtools getfasta`](blah). Refer to the [UCSC page](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) for information about the BED file format. 

The BSPCR-short.bed file used for input into the assess-by-primer.py is slightly modified: 

| chrom | start coordinate | end coordinate | primer ID | sequence |
| ----- | ---------------- | -------------- | --------- | -------- |
| 3	| 101832510	| 101832748	| DMAP-BP1X	| TTCTTCA....TTCACCTA |

Normally, the fifth field of the BED file is a track score. However, bedtools getfasta can output a sequence to the fifth field, which we need for later analysis scripts.

#### Protocol

1. Generate the source tab file for bedtools, which contains the first four columns of the result BED file for each target site. 

2. Run getfasta with the source tab file as input to obtain the sequences for each primer site. 

    ```
    bedtools getfasta -fi </path/to/ref/genome> -bed <primers.tab> -bedout > <output.bed>
    ```

The resulting file can now be used for the analysis scripts. 

Note the following: 
* Bisearch coordinates are 1-based, whereas BED files used a half-open 0-based counting system. Therefore, ensure that the coordinates you are using are properly converted.
* bedtools is used to extract the fasta files. Once again, ensure that the coordinates are properly converted into half-open 0-based counting or else the first base of the target site will be missing. 
* There is a helper script called `generate-primer-bed.sh` that does this process automatically (assuming you have specified the correct genome and input files). To use this script, create the standard directory tree containing the shell scripts and Python helper script (`split-fasta.py`), then run the following command: 
    ```
    ./generate-primer-bed.sh primers.tab /path/to/genome
    ```
    This will generate your BED file, as well as FASTA files corresponding to each target site containined within primer.tab, allowing you to double check that the correct sequences have been obtained. 
    
    Directory structure: 
    ```
    2023-03-08_Analysis/
    ├─ bin/
    │  ├─ generate-primer-bed.sh
    │  ├─ split-fasta.py
    ├─ data/
    │  ├─ primers.tab
    ├─ results/
    ```

### Retrieving run files from the Illumina MiSeq
The Illumina MiSeq is directly connected to the Z drive, the lab NAS. After a run, the entire run folder should be transferred to the Z drive into the NGS folder (`Raw-data/NGS`). The run folder can be found under `TO-BE-FILLED`. 

For analysis, only the FASTQ files in `run_folder/Alignment_1/date_time/Fastq/` are needed.

## Processing FASTQ files using NGS pipeline
The entire pipeline can be run using one script. Alternatively, each step of the pipeline can be run individually. This is not recommended as it is significantly more complex, but can aid with troubleshooting if the script is crashing at a specific step.
### Setting up the correct directory tree
The directory tree for the analysis should be organized into the following structure:
```
2023-03-08_analysis/
├─ data/
│  ├─ .fastq
├─ bin/
│  ├─ ngs_process.sh
├─ results/
```
This directory structure is standard for analyses across other projects in the AGC.The script is written to detect the data and results directories in this format. Note that the `FASTQ` files are uncompressed. To uncompress the `FASTQ` files, run the following command in the directory containing the `FASTQ` files: 

```
pigz -d *.fastq.gz
```
`pigz` is an implementation of `gzip` that takes advantage of multiple cores. 
### Running ngs_process.sh and outputs
Ensure that the conda environment has been activated, then run the following script: 
```
./ngs_process.sh </path/to/ref/genome>
```
All steps of the analysis pipeline should occur automatically and output to the correct subdirectories in `results`. 

```
results/
├─ fastqc/
|   ├─ pre-trim/
|   ├─ post-trim/
├─ bismark/
|   ├─ indexed_bam/
|   |  ├─ all_sequences/
|   |  ├─ nonCG_filtered/
|   ├─ mapped_bam/
|   ├─ incomplete_conv_bam/
|   ├─ methylation_extraction/
|   │  ├─ all_sequences/
|   │  ├─ nonCG_filtered/
├─ trimmed_fastq/
```
| folder | description |
| ------ | ----------- |
| `fastqc/pre-trim/` | Contains `fastqc` and `multiqc` results for the `FASTQ` files directly from the MiSeq.
| `fastqc/post-trim/` | Contains `fastqc` and `multiqc` results for the `FASTQ` files after trimming with Trimmomatic.
| `trimmed_fastq/` | Contains `FASTQ` files after trimming with Trimmomatic. 
| `bismark/mapped_bam/` | Contains `BAM` files after mapping with `bowtie2` against the bisulfite-converted reference genome. Also contains `report.txt` files that describe the mapping results for each sample.|
| `bismark/incomplete_conv_bam/` | Contains `BAM` files with reads containing >3 unconverted non-CpG cytosines removed. Also contains `non-conversion_filtering.txt` files that describe the filtering results for each sample. Removed reads are also stored in a separate `BAM` file named `nonCG_removed_seqs.bam`. |
| `bismark/indexed_bam/` | Contains `BAM` files that have been sorted and indexed using `samtools`, allowing for viewing in the `Integrative Genome Viewer`. Subdirectories `all_sequences` and `nonCG_filtered` refer to unfiltered and filtered datasets. |
| `bismark/methylation_extraction/` | Contains compressed versions of the coverage and `bedGraph` outputs from the methylation caller. Additionally, the m-bias and splitting reports can be found in here, but these are only relevant for directional whole-genome bisulfite sequencing. 

## QC and processing Bismark outputs
### Generating QC results for the entire run
#### Review the `fastqc` and `multiqc` results
This protocol will not go over how to interpret `fastqc` results. A basic tutorial can be found [here](https://rtsf.natsci.msu.edu/genomics/tech-notes/fastqc-tutorial-and-faq/). `multiqc` is a useful tool for summarizing all of the fastqc results from a single run. The results (`.html` files) for both tools can be viewed in a web browser. 
#### process-bismark-reports.py and process-bismark-incomplete-reports.py
`process-bismark-reports.py` generates a .csv file containing summary statistics of the read mapping for each sample. It crawls through the `report.txt` files and extracts the relevant information. To use this script, run the following command: 
```
python ./process-bismark-reports.py <bismark/mapped_bam> -o </path/to/results/directory>
```
`process-bismark-incomplete-reports.py` generates a .csv file containing the results from filtering incompletely-converted reads. It crawls through the `non-conversion_filtering.txt` files and extracts the relevant information. To use this script, run the following command: 
```
python ./process-bismark-incomplete-reports.py <bismark/incomplete_conv_bam> -o </path/to/results/directory>
```
In a future revision of this pipeline, both scripts will be combined. For now, the data from both .csv files has to be manually combined to yield the final run results .csv file. The percentage of usable reads can now be determined. 

### process-bismark-coverage.py
The purpose of the script is to generate sample summary files and human-readable .csv files of the methylation call data. The outputs of this script are used in downstream analysis scripts. To run the script, run the following command: 

```
python ./process-bismark-coverage.py </path/to/coverage-files> </path/to/primer.bed> -o <results/>
```

The output files are organized in this following directory structure: 

```
results/
├─ individual_primers/
│  ├─ sample-name_primer_processed.csv
├─ summaries/
│  ├─ sample-name_summary.csv
```
#### `sample-name_primer_processed.csv`
A `sample-name_primer_processed.csv` file is generated for each target site that has mapped reads to it. 

| chromosome | start | end | methylation_percentage | count_methylated | count_unmethylated | seq_pos | coverage | basecall |
| ---------- | ----- | --- | -----------------------| ---------------- | ------------------ | ------- | -------- | -------- |
| 3 | 101832538 | 101832538 | 64.4225701 | 37894 | 20927 | 27 | 58821 | c |

Each row corresponds to a CG position in the target site. the `count_methylated` and `count_unmethylated` refer to the number of reads where at that `CG` position, a `C` or `T` was observed, respectively. `coverage` is the sum of the two counts. `seq_pos` is the relative position of the `CG` in the primer site sequence. `basecall` is the basecall at that position if it was methylated. For the positive strand, it would be `C`, whereas the negative stand would be `G`.

#### `sample-name_summary.csv`
A `sample-name_summary.csv` file is generated for each sample, and summarizes the methylation state of each target site that data was available for. 

| primer | #CG_total | #CG_covered | mean_methylation | std_methylation | mean_coverage | std_coverage | #CG<mean |
| ------ | --------- | ----------- | ---------------- | --------------- | ------------- | ------------ | -------- |
| DMAP1-BP1X | 15 | 15 | 62.13961 | 19.5795 | 536.8667 | 10.45307 | 6 |

Each row corresponds to a target site that there is data for. `#CG_total` refers to the total number of `CG` sites in the target region. `#CG_covered` refers to the number of `CG` sites in the region that have data for. `mean_methylation` and `std_methylation` are the summary statistics for the methylation status of the entire target site. `mean_coverage` and `std_coverage` are summary statistics for the number of reads that contribute to the calls of the region. 

### cgvisualizer.py
This script plots the methylation calls for each target site for each sample. The input files for this script are the coverage files and the primer.bed file. To run the script, run the following command: 
```
python ./cgvisualizer.py </path/to/coverage-files> </path/to/primer.bed> -o <results/>
```
Each sample gets its own `.png` file containing the methylation plots for all 21 target sites. Target sites that do not have data are labeled as 'no data'. There is also a line plot that corresponds to the coverage at each `CG` position overlayed on top of the methylation plot. 

### assess-by-primer.py
This script generates box plots of the site coverage for each primer set. The input files are the `sample-name_summary.csv` files from `process-bismark-coverage.py` for the group of samples that you want to analyze together, as well as a `sample-sequences.csv` file containing two columns: the sample names, and the final number of usable reads (refer to QC steps to get this value). To run the script, run the following command: 
```
python ./assess-by-primer.py </path/to/sample-name_summary.csv> -n </path/to/sample-sequences.csv> -t <tag> -o <results/>
```
Two outputs are generated: `tag_primer_summary.csv` files and a `tag_primer_boxplot_summary.png` file that contains the boxplot figure. `tag_primer_summary.csv` contains the data that was used to generate the boxplot for each target site. 