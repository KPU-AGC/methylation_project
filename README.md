# Methylation Project
This repository contains the scripts and pipelines for analyzing the targeted bisulfite amplicon sequencing data produced during the methylation project. 

## Project Outline
### Identifying target sites
21 genes of interest were identified by our partners. Some of these genes had CpG islands that were previously worked on by our partners and others, and so these CpG islands were designed as our targets. For other genes, we identified putative CpG islands near the site of interest. Detailed information on these sites can be found in the Teams OneDrive. 

### Primer design for target sites
After choosing each target site, two sets of primers were designed: 

1. Genomic sequencing primers : this set of primers targeted the unconverted genomic sequence of the target site. The purpose of this set was to facilitate genotyping of the target site. 
2. Bisulfite PCR primers : this set of primers targeted the bisulfite-converted sequence of the target site. The amplicon produced by this set would be used in the NGS run to get the methylation calls in the site. 

The genomic sequencing primers were designed using [PrimerBLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/), whereas the bisulfite PCR primers were designed using [Bisearch](http://bisearch.enzim.hu/). Due to the non-complementary nature of bisulfite-converted ssDNA, a primer set targeting the positive and negative converted strands was created where possible. 

The genomic sequencing primers were optimized with DreamTaq polymerase. The bisulfite PCR primers were optimized with Phusion U polymerase. The final optimized conditions can be found in the Teams OneDrive summary documents. 

### Sequencing and genotyping each target site for putative semen donors
The semen donors were genotyped at each target site. The genotyping calls allow us to determine the parental origin of our methylation calls, which is sigificant for some of our imprinted genes. 

The sequencing instrument was the SeqStudio, using BigDye Terminator v3.1 chemistry. ab1 files were quality controlled and filtered using [sanger_qc.py](https://github.com/KPU-AGC/general-resources/blob/main/sanger-processing/sanger_qc.py) and [sanger-sequence-trim.py](https://github.com/KPU-AGC/general-resources/blob/main/sanger-processing/sanger_sequence_trim.py). 

[tracy](https://www.gear-genomics.com/docs/tracy/) was the command-line program that was used to generate genotyping calls for each ab1 file. There is a [software-assisted genotyping program](https://github.com/KPU-AGC/batch-genotyping) that we created to QC tracy's genotype calls. 

### Bisulfite PCR amplification of samples 
Bisulfite PCR amplification is done for each target site for each sample. Both a singleplex and multiplex set up is currently being evaluated. Success of the PCRs was evaluated using agarose gel electrophoresis. Each agarose gel image of all of the PCRs is processed using the [GelAnalyzer](http://www.gelanalyzer.com/) software. The relative quantity of each target amplicon was determined by assessing its band's thickness and intensity in the gel image. The target amplicons were normalized and pooled together based on their relative quantity. A gel excision cleanup was used to separate the target size range from the incorrectly sized non-specific product. 

### NGS library preparation and MiSeq sequencing
NGS library preparation was done using the [NEB Ultra II DNA library prep kit](https://www.neb.ca/e7645) with the [NEBNext Multiplex Oligo for Illumina](https://www.neb.ca/E7600) indices. The [MiSeq Reagent Kit v3 (600-cycle)](https://www.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/miseq-reagent-kit-v3.html#:~:text=The%20MiSeq%20v3%20kits%20can,format%20that%20enables%20counting%20applications.) was used. A detailed protocol for how the library prep is organized and conducted can be found in the Teams OneDrive. 

### Processing pipeline summary
The FASTQ files generated by sequencer are first processed with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [multiqc](https://multiqc.info/) to assess the read quality and parameters of the run. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) v0.4.1 is used to trim the adaptor and index sequences. fastqc and multiqc is run again to get the post-trim QC data. 

The trimmed FASTQs are then processed by the [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/#:~:text=Bismark%20is%20a%20program%20to,of%20their%20samples%20straight%20away.) mapper, using the [bowtie2](https://www.google.com/search?q=bowtie2+aligner&rlz=1C1VDKB_enCA1033CA1033&oq=bowtie2+aligner&aqs=chrome..69i57.1429j0j7&sourceid=chrome&ie=UTF-8) aligner to generate BAM files. Sorting and indexing of the BAM files using samtools is done to allow for viewing of the reads in [IGV](https://software.broadinstitute.org/software/igv/home).

The Bismark methylation extractor is used to generate a coverage file and bedgraph file for the mapped reads, allowing us to asses the methylation calls for each CG of our target sites. Additional downstream scripts generate human-readable figures. 

## Scripts
### Bisulfite PCR primer design
[bisearch-primer-design.py](https://github.com/KPU-AGC/methylation_project/blob/main/misc/bisearch_primer_design.py) is available to batch design bisulfite PCR primers using [Bisearch](http://bisearch.enzim.hu/).

### Processing .fastq data with Bismark
[ngs-process.sh](https://github.com/KPU-AGC/methylation_project/blob/main/ngs_shell_scripts/ngs-process.sh) is the main script that runs the entire pipeline for processing .fastq data. 

The [individual_processes](https://github.com/KPU-AGC/methylation_project/tree/main/ngs_shell_scripts/individual_processes) directory contains scripts for each step of the process that can be executed independently.

### Processing Bismark output
There are two scripts that are important for assessing the run quality: 
* [process-bismark-reports.py]() : extracts the number of mapped reads for all of the samples and outputs it into a .csv. 
* [process-bismark-incomplete-reports.py]() : extracts the number of filtered reads for all samples and outputs it into a .csv. 

[process-bismark-coverage.py](https://github.com/KPU-AGC/methylation_project/blob/main/ngs_process_bismark/process-bismark-coverage.py) converts the Bismark coverage files into several .csv files containing methylation calls for each CG and site for further processing, as well as gathering some basic statistical data.

The output for [process-bismark-coverage.py](https://github.com/KPU-AGC/methylation_project/blob/main/ngs_process_bismark/process-bismark-coverage.py) can then be fed to two additional scripts: 
* [cgvisualizer.py](https://github.com/KPU-AGC/methylation_project/blob/main/ngs_process_bismark/cgvisualizer.py) : produces a figure for each sample containing the methylation and coverage plot for each target site. 
* [assess-by-primer.py](https://github.com/KPU-AGC/methylation_project/blob/main/ngs_process_bismark/assess-by-primer.py) : produces a figure containing a box plot for each target site that details its normalized coverage distribution. 
