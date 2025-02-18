## methyl-sim
A repository for simulating methylated reads with some experimental weights.


## Background
DNA methylation is a crucial epigenetic modification that plays a significant role in gene regulation, development, and cellular differentiation. The patterns of DNA methylation can vary widely between different cell types, developmental stages, and even within a single cell cycle.

CpG methylation patterns exhibit distinct tissue-specific profiles. For instance, non-CpG methylation is predominantly present in pluripotent cell types and decreases upon differentiation, showing near-complete absence in various somatic cell types. Similarly, specific profiles of CpG island hypermethylation exist for different tumor types, allowing classification within hierarchical clusters according to the originating tissue.


## Objective
To simulate methylated reads for stress-testing analysis pipelines in targeted-amplicon bisulfite-sequencing analysis, especially to apply experimental weights and methylation patterns.


## Requirements
The recommended way to install dependencies is via `conda` using the `environment.yml` provided in the root directory of this repository.

The following programs and packages should be installed.
- [`art`](https://manpages.debian.org/testing/art-nextgen-simulation-tools/art_illumina.1.en.html)
- [`pandas`](https://pandas.pydata.org/)

## Usage
### 1A. Generating simulated reads of a given region with all-or-nothing methylation.
This is currently surprisingly involved, so I'll probably adjust this usage. But the purpose here is to take a template, identify CpG positions, and then render them as completely un-/methylated.

The following will generate a paired-end set of completely methylated reads for "DNMT3A-BP1X".
```shell
python src/methyl-sim.py \
  --weights-path pkl/new_weights.pkl \
  --region "DNMT3A-BP1X" \
  --manual-weights 100 \
  --out-prefix "SIMULATED_SAMPLE"
```

The following will generate a paired-end set of half completely-methylated and half completely unmethylated reads for "PLAGL1-BN1X". Note: the proportions have to add up to 100.
```shell
python src/methyl-sim.py \
  --weights-path pkl/new_weights.pkl \
  --region "PLAGL1-BN1X" \
  --manual-weights 100:0.5,0:0.5 \
  --out-prefix "SIMULATED_SAMPLE"
```

I suppose, if you really wanted to, you could also generate reads that are variably methylated in different proportion. The program will try its best to hit the methylation value--this works better if the target region has lots of CpG sites.
```shell
python src/methyl-sim.py \
  --weights-path pkl/new_weights.pkl \
  --region "PLAGL1-BN1X" \
  --manual-weights 80:0.25,75:0.25,60:25 \
  --out-prefix "SIMULATED_SAMPLE"
```


### 1B. Generating simulated reads of a given region with experimental weights.
I think this implementation is much more elegant and true to biology--that methylation exists different states across a population of cells. In any case, this uses experimental weights to generate methylation patterns typical of a given region.

The following will generate a set of reads for "DNMT3A-BP1X" that display the characteristic W-pattern.
```shell
python src/methyl-sim.py \
  --weights-path pkl/new_weights.pkl \
  --region "DNMT3A-BP1X" \
  --out-prefix "SIMULATED_SAMPLE"
```

The following will generate a set of reads for "KCNQ1-BP1X" that is typical of samples with large offspring syndrome (LOS).
```shell
python src/methyl-sim.py \
  --weights-path pkl/new_weights.pkl \
  --region "KCNQ1-BP1X" \
  --metadata pkl/metadata.pkl \
  --apply-filter "large_calf=TRUE" \
  --out-prefix "SIMULATED_SAMPLE"
```

### 2. Debugging and investigating output.
I've included a script here for debugging which produces [`bismark`](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)-like `.cov` output. The intent here is that we should be able to easily integrate these simulated files into our analysis. This provides, at least, a cursory glance without dealing with bisulfite-converted alignment.

This relies on alignment and will require some quick pre-processing of files. The following should be installed:
- [`biopython`](https://biopython.org/)
- [`flash`](https://ccb.jhu.edu/software/FLASH/)

The following will stitch the forward and reverse reads to make alignment easier and produce a `.cov` file-like output for further parsing.
```shell
flash {R1} {R2} -m 250 -o {OUTPUT}
python src/debug.py {OUTPUT}.extendedFrags.fastq
```


### Note
This is where the R&D budget goes.

### TODO:
- [ ] include parameters for introducing parental SNPs.

## References