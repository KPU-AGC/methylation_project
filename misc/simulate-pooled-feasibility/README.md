## Overview
This repository contains a Python script to simulate bootstrap-based hierarchical read depth allocation for single-cell sequencing data. The script estimates read distributions across cells and regions, computes detection probabilities, and generates detailed statistics for various cell pooling configurations.

## Features
- Simulates read depth allocation across multiple pools of cells.
- Supports Gamma distribution modeling of read variability.    
- Implements Negative Binomial and Poisson models for read count distributions.
- Computes overall and per-region detection probabilities.
- Supports bootstrap resampling for statistical accuracy.
- Outputs detailed statistics to a CSV file.
- Configurable via command-line arguments or a JSON configuration file.

## Installation
This script requires Python 3 and the following dependencies:

```
pip install numpy pandas scipy argparse
```

## Usage

### Running the Simulation
Run the script with default parameters:

```
python3 script.py
```

Or specify custom parameters:
```
python3 script.py --new_total_reads 50000000 --cells_per_pool "80,80,80,80,80" \
                  --min_read_depth 15 --n_bootstrap 10000 --min_read_depth_frac 1.0 \
                  --sample_stats_file "pkl/reads-per-sample.pkl" \
                  --region_props_file "pkl/reads-per-region.pkl" \
                  --trim_frac 0.05 --gamma_shape 2.0 --gamma_scale 1.0 \
                  --random_seed 42 --debug --output_csv "simulation_detailed_stats.csv"
```

### Using a JSON Configuration File
You can specify parameters in a JSON config file:
```
{
    "new_total_reads": 50000000,
    "cells_per_pool": [80, 80, 80, 80, 80],
    "min_read_depth": 15,
    "n_bootstrap": 10000,
    "min_read_depth_frac": 1.0,
    "sample_stats_file": "pkl/reads-per-sample.pkl",
    "region_props_file": "pkl/reads-per-region.pkl",
    "trim_frac": 0.05,
    "gamma_shape": 2.0,
    "gamma_scale": 1.0,
    "random_seed": 42,
    "debug": true,
    "output_csv": "simulation_detailed_stats.csv"
}
```

Run the script with:

```
python3 script.py --config config.json
```

## Output
The script generates:

1. **Log messages** showing simulation progress and estimated probabilities/    
2. **A CSV file** with detailed statistics for each pool size and region.

## Example Log Output
```
INFO - Starting bootstrap simulation for read depth allocation...
INFO - Final estimated overall probability: 0.9234
INFO - Final estimated probability for Region 1: 0.8789
INFO - Final estimated probability for Region 2: 0.9012
INFO - Detailed read depth statistics saved to simulation_detailed_stats.csv
```

## Contributing
Feel free to fork and contribute! Open issues or submit pull requests for improvements.

## License
MIT License