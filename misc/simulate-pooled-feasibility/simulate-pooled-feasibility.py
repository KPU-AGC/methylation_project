#!/usr/bin/env python3
import numpy as np
import pandas as pd
from scipy.stats import nbinom, poisson, gamma
import argparse
import logging
import json
import sys
from typing import Tuple, Optional, List, Dict, Any

def setup_logging(log_level: int = logging.INFO, log_file: Optional[str] = None) -> logging.Logger:
    """Setup logging to file and console with specified log level."""
    logger = logging.getLogger()
    logger.setLevel(log_level)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    # Clear any existing handlers.
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    # File handler (if log_file provided)
    if log_file:
        try:
            fh = logging.FileHandler(log_file)
            fh.setFormatter(formatter)
            logger.addHandler(fh)
        except Exception as e:
            logger.error(f"Failed to set up file logging: {e}")
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger

def fit_gamma_params(cell_depth_data: np.ndarray) -> Tuple[float, float]:
    """
    Estimate Gamma distribution parameters (shape, scale) from per-cell read depth data.
    
    Parameters:
        cell_depth_data (np.ndarray): Array of read depth values per cell.
    
    Returns:
        Tuple[float, float]: Estimated shape and scale parameters.
    """
    shape, loc, scale = gamma.fit(cell_depth_data, floc=0)
    return shape, scale

def simulate_sample_read_counts(sample_stats_file: str, trim_frac: float, debug: bool = False) -> Tuple[float, float, bool, Optional[float], Optional[float]]:
    """
    Simulate sample-level read counts based on observed statistics.
    
    Parameters:
        sample_stats_file (str): Path to the sample statistics file.
        trim_frac (float): Fraction of outliers to trim from each end.
        debug (bool): Flag to enable debug logging.
    
    Returns:
        Tuple containing:
            - m (float): Mean of trimmed counts.
            - v (float): Variance of trimmed counts.
            - use_nb (bool): Whether to use Negative Binomial distribution.
            - r (Optional[float]): Negative Binomial parameter r (if applicable).
            - p_nb (Optional[float]): Negative Binomial parameter p (if applicable).
    """
    try:
        df = pd.read_csv(sample_stats_file, sep="\t", header=None, names=["file", "count"])
    except Exception as e:
        logging.error(f"Error reading sample stats file: {e}")
        sys.exit(1)
    df = df[~df["file"].str.contains("Undetermined")]
    counts = df["count"].values
    lower = np.percentile(counts, trim_frac * 100)
    upper = np.percentile(counts, 100 - trim_frac * 100)
    trimmed = counts[(counts >= lower) & (counts <= upper)]
    m = np.mean(trimmed)
    v = np.var(trimmed, ddof=1)
    use_nb = v > m
    if use_nb:
        r = m**2 / (v - m) if v > m else 1e6  # fallback to near-Poisson if needed
        p_nb = r / (r + m)
    else:
        r, p_nb = None, None
    if debug:
        logging.debug(f"Sample read counts: mean={m:.2f}, variance={v:.2f}, using {'Negative Binomial' if use_nb else 'Poisson'}")
    return m, v, use_nb, r, p_nb

def simulate_cell_allocation(T_i: int, n_cells: int, gamma_params: Tuple[float, float], debug: bool = False, pool_index: int = 0) -> np.ndarray:
    """
    Allocate total reads T_i among n_cells using Gamma variability.
    
    Parameters:
        T_i (int): Total reads allocated to the pool.
        n_cells (int): Number of cells in the pool.
        gamma_params (Tuple[float, float]): Parameters (shape, scale) for the Gamma distribution.
        debug (bool): Flag to enable debug logging.
        pool_index (int): Index of the pool (for logging purposes).
    
    Returns:
        np.ndarray: Array of per-cell read totals.
    """
    cell_factors = gamma.rvs(a=gamma_params[0], scale=gamma_params[1], size=n_cells)
    cell_probs = cell_factors / cell_factors.sum()
    cell_totals = np.random.multinomial(T_i, cell_probs)
    if debug:
        logging.debug(f"Pool {pool_index}: T_i={T_i}, n_cells={n_cells}")
        logging.debug(f"Pool {pool_index}: Gamma factors: {cell_factors}")
        logging.debug(f"Pool {pool_index}: Cell totals: {cell_totals}")
    return cell_totals

def simulate_region_allocation(cell_totals: np.ndarray, region_props: Optional[pd.DataFrame], pool_index: int, debug: bool = False) -> List[np.ndarray]:
    """
    For a given pool, allocate each cell's reads across regions using Dirichlet probabilities.
    
    Parameters:
        cell_totals (np.ndarray): Array of per-cell read totals.
        region_props (Optional[pd.DataFrame]): DataFrame containing region properties.
        pool_index (int): Index of the pool (for logging purposes).
        debug (bool): Flag to enable debug logging.
    
    Returns:
        List[np.ndarray]: List of arrays representing region counts for each cell.
    """
    n_cells = len(cell_totals)
    # Determine base region weights
    if region_props is not None:
        if region_props.shape[0] > pool_index:
            base_region = region_props.iloc[pool_index].values.astype(float)
        else:
            base_region = region_props.sample(n=1).values.flatten().astype(float)
    else:
        base_region = np.array([1.0])
    effective_alpha = np.maximum(base_region, 1e-6)
    if debug:
        logging.debug(f"Pool {pool_index}: Base region weights: {base_region}")
        logging.debug(f"Pool {pool_index}: Effective alpha: {effective_alpha}")
    # Generate Dirichlet probabilities for each cell
    p_regions = np.random.dirichlet(effective_alpha, size=n_cells)
    if debug:
        logging.debug(f"Pool {pool_index}: Dirichlet probabilities (first 3 cells): {p_regions[:3]}")
    # Allocate reads across regions for each cell using list comprehension.
    cell_region_counts = [np.random.multinomial(ct, p) for ct, p in zip(cell_totals, p_regions)]
    return cell_region_counts

def _simulate_iteration(sample_reads: np.ndarray, cells_per_pool: List[int], gamma_params: Tuple[float, float], region_props: Optional[pd.DataFrame], debug: bool = False) -> List[np.ndarray]:
    """
    Helper function to simulate a single bootstrap iteration for all pools.
    
    Parameters:
        sample_reads (np.ndarray): Array of allocated read counts for each pool.
        cells_per_pool (List[int]): List containing the number of cells per pool.
        gamma_params (Tuple[float, float]): Gamma distribution parameters.
        region_props (Optional[pd.DataFrame]): Region properties DataFrame.
        debug (bool): Flag to enable debug logging.
    
    Returns:
        List[np.ndarray]: List of arrays for each cell's region read counts across all pools.
    """
    all_cell_region_counts = []
    for i, T_i in enumerate(sample_reads):
        n_cells = cells_per_pool[i]
        cell_totals = simulate_cell_allocation(T_i, n_cells, gamma_params, debug, pool_index=i)
        cell_region_counts = simulate_region_allocation(cell_totals, region_props, pool_index=i, debug=debug)
        all_cell_region_counts.extend(cell_region_counts)
    return all_cell_region_counts

def simulate_bootstrap_min_depth_hierarchical(new_total_reads: int, cells_per_pool: List[int], min_read_depth: int,
                                              n_bootstrap: int, min_read_depth_frac: float, sample_stats_file: str,
                                              region_props_file: Optional[str], trim_frac: float, gamma_params: Tuple[float, float],
                                              random_seed: Optional[int] = None, debug: bool = False) -> Tuple[float, np.ndarray]:
    """
    Simulation that computes overall and per-region detection probabilities.
    
    Parameters:
        new_total_reads (int): Total number of reads to simulate.
        cells_per_pool (List[int]): List containing the number of cells per pool.
        min_read_depth (int): Minimum read depth threshold for detection.
        n_bootstrap (int): Number of bootstrap iterations.
        min_read_depth_frac (float): Fraction of cells required to meet min_read_depth.
        sample_stats_file (str): Path to sample statistics file.
        region_props_file (Optional[str]): Path to region properties file.
        trim_frac (float): Fraction to trim outliers in sample statistics.
        gamma_params (Tuple[float, float]): Gamma distribution parameters (shape, scale).
        random_seed (Optional[int]): Random seed for reproducibility.
        debug (bool): Flag to enable debug logging.
    
    Returns:
        Tuple containing overall detection probability and per-region detection probabilities.
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    n_samples = len(cells_per_pool)
    
    m, v, use_nb, r, p_nb = simulate_sample_read_counts(sample_stats_file, trim_frac, debug)
    
    # Load region properties if provided.
    region_props = None
    if region_props_file:
        try:
            region_props = pd.read_csv(region_props_file, sep="\t", index_col=0)
        except Exception as e:
            logging.error(f"Error reading region properties file: {e}")
    
    overall_success_iterations = 0
    region_success_counts = None
    cells_array = np.array(cells_per_pool)
    
    for b in range(n_bootstrap):
        # Simulate per-pool raw read counts.
        if use_nb:
            raw = nbinom.rvs(r, p_nb, size=n_samples)
        else:
            raw = poisson.rvs(m, size=n_samples)
        weighted = raw * cells_array
        total_weight = weighted.sum()
        if total_weight == 0:
            weighted = cells_array
            total_weight = weighted.sum()
        sample_reads = np.round(weighted / total_weight * new_total_reads).astype(int)
        sample_reads = sample_reads // 2 * 2  # ensure even number
        
        if debug and (b < 5 or b % 1000 == 0):
            logging.debug(f"Iteration {b}: Raw counts: {raw}, Weighted: {weighted}, Total weight: {total_weight}, Sample reads: {sample_reads}")
        
        all_cell_region_counts = _simulate_iteration(sample_reads, cells_per_pool, gamma_params, region_props, debug)
        all_cell_region_counts = np.array(all_cell_region_counts)
        region_detection = np.mean(all_cell_region_counts >= min_read_depth, axis=0)
        
        if region_success_counts is None:
            region_success_counts = np.zeros_like(region_detection)
        
        if debug and (b < 5 or b % 1000 == 0):
            logging.debug(f"Iteration {b}: Region detection fractions: {region_detection}")
        
        if np.any(region_detection >= min_read_depth_frac):
            overall_success_iterations += 1
        
        region_success_counts += (region_detection >= min_read_depth_frac)
    
    overall_probability = overall_success_iterations / n_bootstrap
    per_region_probabilities = region_success_counts / n_bootstrap
    
    if debug:
        logging.info(f"Total successful iterations (any region): {overall_success_iterations} out of {n_bootstrap}")
    return overall_probability, per_region_probabilities

def simulate_bootstrap_detailed_stats(new_total_reads: int, cells_per_pool: List[int], min_read_depth: int,
                                      n_bootstrap: int, min_read_depth_frac: float, sample_stats_file: str,
                                      region_props_file: Optional[str], trim_frac: float, gamma_params: Tuple[float, float],
                                      random_seed: Optional[int] = None, debug: bool = False) -> pd.DataFrame:
    """
    Extended simulation that accumulates detailed read depth statistics per pooled cell group.
    
    For each unique pool size, across all bootstrap iterations, it computes:
      - The average read depth per cell per region,
      - The median, minimum, and maximum read depth observed per cell per region,
      - The probability that the pooled group meets the specified read depth threshold,
      - The minimum, median, mean, and maximum proportion of cells that meet the read count threshold.
    
    Parameters:
        new_total_reads (int): Total number of reads to simulate.
        cells_per_pool (List[int]): List of cell counts per pool.
        min_read_depth (int): Minimum read depth threshold for detection.
        n_bootstrap (int): Number of bootstrap iterations.
        min_read_depth_frac (float): Fraction of cells required to meet min_read_depth.
        sample_stats_file (str): Path to sample statistics file.
        region_props_file (Optional[str]): Path to region properties file.
        trim_frac (float): Fraction to trim outliers in sample statistics.
        gamma_params (Tuple[float, float]): Gamma distribution parameters.
        random_seed (Optional[int]): Random seed for reproducibility.
        debug (bool): Flag to enable debug logging.
    
    Returns:
        pd.DataFrame: DataFrame containing detailed statistics per pooled cell group and region.
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    n_samples = len(cells_per_pool)
    
    m, v, use_nb, r, p_nb = simulate_sample_read_counts(sample_stats_file, trim_frac, debug)
    
    # Load region properties if provided.
    region_props = None
    if region_props_file:
        try:
            region_props = pd.read_csv(region_props_file, sep="\t", index_col=0)
        except Exception as e:
            logging.error(f"Error reading region properties file: {e}")
    
    # Count unique pool sizes.
    unique_pool_sizes: Dict[int, int] = {}
    for pool_size in cells_per_pool:
        unique_pool_sizes[pool_size] = unique_pool_sizes.get(pool_size, 0) + 1

    # Get number of regions using a dummy simulation.
    T_i_dummy = 1000
    dummy_cell_totals = simulate_cell_allocation(T_i_dummy, cells_per_pool[0], gamma_params, debug)
    dummy_region_counts = simulate_region_allocation(dummy_cell_totals, region_props, pool_index=0, debug=debug)
    n_regions = len(dummy_region_counts[0])
    
    # Determine region names.
    if region_props is not None:
        region_names = list(region_props.columns)
        if len(region_names) != n_regions:
            region_names = [f"Region {i}" for i in range(n_regions)]
    else:
        region_names = [f"Region {i}" for i in range(n_regions)]
    
    # Initialize aggregates for each pool group.
    group_aggregates: Dict[int, Dict[str, Any]] = {}
    for pool_size in unique_pool_sizes.keys():
        group_aggregates[pool_size] = {
            'total_sum': np.zeros(n_regions),
            'total_cells': 0,
            'min': np.full(n_regions, np.inf),
            'max': np.full(n_regions, -np.inf),
            'success_count': np.zeros(n_regions),
            'all_values': [ [] for _ in range(n_regions) ],
            'prop_list': [ [] for _ in range(n_regions) ]
        }
    
    for b in range(n_bootstrap):
        if use_nb:
            raw = nbinom.rvs(r, p_nb, size=n_samples)
        else:
            raw = poisson.rvs(m, size=n_samples)
        weighted = raw * np.array(cells_per_pool)
        total_weight = weighted.sum()
        if total_weight == 0:
            weighted = np.array(cells_per_pool)
            total_weight = weighted.sum()
        sample_reads = np.round(weighted / total_weight * new_total_reads).astype(int)
        sample_reads = sample_reads // 2 * 2
        
        if debug and (b < 5 or b % 1000 == 0):
            logging.debug(f"Iteration {b}: Raw counts: {raw}, Weighted: {weighted}, Total weight: {total_weight}, Sample reads: {sample_reads}")
        
        pool_group_data: Dict[int, List[np.ndarray]] = {}
        for i in range(n_samples):
            T_i = sample_reads[i]
            n_cells = cells_per_pool[i]
            cell_totals = simulate_cell_allocation(T_i, n_cells, gamma_params, debug, pool_index=i)
            cell_region_counts = simulate_region_allocation(cell_totals, region_props, pool_index=i, debug=debug)
            cell_region_counts = np.array(cell_region_counts)
            group_key = n_cells
            if group_key not in pool_group_data:
                pool_group_data[group_key] = []
            pool_group_data[group_key].append(cell_region_counts)
        
        for pool_size, list_of_arrays in pool_group_data.items():
            group_data = np.vstack(list_of_arrays)
            group_aggregates[pool_size]['total_sum'] += np.sum(group_data, axis=0)
            group_aggregates[pool_size]['total_cells'] += group_data.shape[0]
            group_aggregates[pool_size]['min'] = np.minimum(group_aggregates[pool_size]['min'], np.min(group_data, axis=0))
            group_aggregates[pool_size]['max'] = np.maximum(group_aggregates[pool_size]['max'], np.max(group_data, axis=0))
            # Store all cell read counts per region for median calculation
            for region_idx in range(n_regions):
                group_aggregates[pool_size]['all_values'][region_idx].append(group_data[:, region_idx])
            # Compute and store the proportion of cells meeting the threshold in this iteration
            group_prop = np.mean(group_data >= min_read_depth, axis=0)
            for region_idx in range(n_regions):
                group_aggregates[pool_size]['prop_list'][region_idx].append(group_prop[region_idx])
            success_indicator = (group_prop >= min_read_depth_frac)
            group_aggregates[pool_size]['success_count'] += success_indicator.astype(int)
    
    results = []
    for pool_size, agg in group_aggregates.items():
        avg = agg['total_sum'] / agg['total_cells']
        success_prob = agg['success_count'] / n_bootstrap
        # Compute median per region for read depths
        medians = []
        for region_idx in range(n_regions):
            region_values = np.concatenate(agg['all_values'][region_idx])
            median_val = np.median(region_values)
            medians.append(median_val)
        # Compute summary statistics for the proportion of cells meeting the threshold
        for region_idx in range(n_regions):
            prop_arr = np.array(agg['prop_list'][region_idx])
            min_prop = np.min(prop_arr)
            median_prop = np.median(prop_arr)
            mean_prop = np.mean(prop_arr)
            max_prop = np.max(prop_arr)
            num_pools = unique_pool_sizes[pool_size]
            pool_label = f"{num_pools}x{pool_size} cells"
            results.append({
                "Pool_Group": pool_label,
                "Region": region_names[region_idx],
                "Avg_Read_Depth": avg[region_idx],
                "Median_Read_Depth": medians[region_idx],
                "Min_Read_Depth": agg['min'][region_idx],
                "Max_Read_Depth": agg['max'][region_idx],
                "Success_Probability": success_prob[region_idx],
                "Min_Prop": min_prop,
                "Median_Prop": median_prop,
                "Mean_Prop": mean_prop,
                "Max_Prop": max_prop
            })
    results_df = pd.DataFrame(results)
    return results_df

def load_config(config_file: str) -> dict:
    """
    Load simulation parameters from a JSON configuration file.
    
    Parameters:
        config_file (str): Path to the configuration file.
    
    Returns:
        dict: Configuration parameters.
    """
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
    except Exception as e:
        logging.error(f"Error reading config file: {e}")
        sys.exit(1)
    return config

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Bootstrap-based hierarchical read depth allocation simulation for single-cell sequencing data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--new_total_reads", type=int, default=50000000, help="Total number of reads to simulate.")
    parser.add_argument("--cells_per_pool", type=str, default="80,80,80,80,80", help="Comma-separated list of cell counts per pool.")
    parser.add_argument("--min_read_depth", type=int, default=15, help="Minimum read depth per cell for detection.")
    parser.add_argument("--n_bootstrap", type=int, default=10000, help="Number of bootstrap iterations.")
    parser.add_argument("--min_read_depth_frac", type=float, default=1.0, help="Fraction of cells that must meet min_read_depth.")
    parser.add_argument("--sample_stats_file", type=str, default="pkl/reads-per-sample.pkl", help="Path to sample stats file.")
    parser.add_argument("--region_props_file", type=str, default="pklreads-per-region.pkl", help="Path to region properties file.")
    parser.add_argument("--trim_frac", type=float, default=0.05, help="Fraction to trim outliers in sample stats.")
    parser.add_argument("--gamma_shape", type=float, default=2.0, help="Shape parameter for the Gamma distribution.")
    parser.add_argument("--gamma_scale", type=float, default=1.0, help="Scale parameter for the Gamma distribution.")
    parser.add_argument("--random_seed", type=int, default=42, help="Random seed for reproducibility.")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging.")
    parser.add_argument("--log_file", type=str, default=None, help="File to write logs to.")
    parser.add_argument("--config", type=str, default=None, help="Path to a JSON configuration file.")
    parser.add_argument("--output_csv", type=str, default="simulation_detailed_stats.csv", help="Output CSV file for detailed stats.")
    args = parser.parse_args()

    if args.config:
        config = load_config(args.config)
        for key, value in config.items():
            setattr(args, key, value)
    args.cells_per_pool = [int(x.strip()) for x in args.cells_per_pool.split(",")]
    return args

def main() -> None:
    """Main function to run the simulation."""
    args = parse_args()
    
    log_level = logging.DEBUG if args.debug else logging.INFO
    setup_logging(log_level=log_level, log_file=args.log_file)
    
    logging.info("Starting bootstrap simulation for read depth allocation...")
    gamma_params = (args.gamma_shape, args.gamma_scale)
    
    overall_probability, per_region_probabilities = simulate_bootstrap_min_depth_hierarchical(
        new_total_reads=args.new_total_reads,
        cells_per_pool=args.cells_per_pool,
        min_read_depth=args.min_read_depth,
        n_bootstrap=args.n_bootstrap,
        min_read_depth_frac=args.min_read_depth_frac,
        sample_stats_file=args.sample_stats_file,
        region_props_file=args.region_props_file,
        trim_frac=args.trim_frac,
        gamma_params=gamma_params,
        random_seed=args.random_seed,
        debug=args.debug
    )
    
    logging.info(f"Final estimated overall probability: {overall_probability:.4f}")
    for idx, prob in enumerate(per_region_probabilities):
        region_names = None
        if args.region_props_file:
            try:
                region_props = pd.read_csv(args.region_props_file, sep="\t", index_col=0)
                region_names = list(region_props.columns)
            except Exception as e:
                logging.error(f"Error reading region properties file: {e}")

        region_name = region_names[idx] if region_names and idx < len(region_names) else f"Region {idx}"
        logging.info(f"Final estimated probability for {region_name}: {prob:.4f}")

    detailed_df = simulate_bootstrap_detailed_stats(
        new_total_reads=args.new_total_reads,
        cells_per_pool=args.cells_per_pool,
        min_read_depth=args.min_read_depth,
        n_bootstrap=args.n_bootstrap,
        min_read_depth_frac=args.min_read_depth_frac,
        sample_stats_file=args.sample_stats_file,
        region_props_file=args.region_props_file,
        trim_frac=args.trim_frac,
        gamma_params=gamma_params,
        random_seed=args.random_seed,
        debug=args.debug
    )
    
    try:
        detailed_df.to_csv(args.output_csv, index=False)
        logging.info(f"Detailed read depth statistics saved to {args.output_csv}")
    except Exception as e:
        logging.error(f"Failed to save detailed statistics to CSV: {e}")
    
    # Additionally, print out the proportion summary for each pool group and region.
    for index, row in detailed_df.iterrows():
        logging.info(f"Pool_Group: {row['Pool_Group']}, Region: {row['Region']}, "
                     f"Min_Prop: {row['Min_Prop']:.4f}, Median_Prop: {row['Median_Prop']:.4f}, "
                     f"Mean_Prop: {row['Mean_Prop']:.4f}, Max_Prop: {row['Max_Prop']:.4f}")

if __name__ == "__main__":
    main()
