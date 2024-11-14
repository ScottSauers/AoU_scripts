import pandas as pd
import os
import re
from tqdm import tqdm

def read_bim_file(bim_file, target_chrom='22'):
    """
    Reads the BIM file and returns a set of (chrom, pos) tuples and a DataFrame.
    Only considers the target chromosome.
    """
    print("Reading BIM file...")
    try:
        # BIM file columns: chrom, variant, ignore1, ignore2, allele1, allele2
        bim_df = pd.read_csv(
            bim_file,
            delim_whitespace=True,
            header=None,
            usecols=[0, 1, 4, 5],
            names=['chrom', 'variant', 'allele1', 'allele2'],
            dtype={'chrom': str, 'variant': str, 'allele1': str, 'allele2': str}
        )
        # Remove 'chr' prefix from chromosome
        bim_df['chrom'] = bim_df['chrom'].str.replace('chr', '', case=False).str.strip()
        
        # Extract position from 'variant' (assuming format 'chrom:pos:allele1:allele2')
        bim_split = bim_df['variant'].str.split(':', expand=True)
        bim_df['chrom_variant'] = bim_split[0].str.replace('chr', '', case=False).str.strip()
        bim_df['pos_variant'] = bim_split[1].str.strip()
        
        # Filter to target chromosome
        bim_df = bim_df[bim_df['chrom_variant'] == target_chrom]
        
        # Create set of (chrom, pos)
        bim_positions = set(zip(bim_df['chrom_variant'], bim_df['pos_variant']))
        print(f"Total positions in BIM file for chr{target_chrom}: {len(bim_positions)}\n")
        return bim_positions, bim_df
    except Exception as e:
        print(f"Error reading BIM file: {e}")
        return set(), pd.DataFrame()

def read_weights_file(weights_file, target_chrom='22'):
    """
    Reads the weights CSV, filters to target_chrom, and returns a DataFrame and set of positions.
    """
    print("Reading weights file...")
    try:
        # Read weights CSV with only necessary columns
        weights_df = pd.read_csv(
            weights_file,
            dtype={'chr': str, 'pos': str, 'effect_allele': str, 'weight': float, 'id': str},
            usecols=['chr', 'pos', 'effect_allele', 'weight', 'id']
        )
        # Remove 'chr' prefix and filter to target_chrom
        weights_df['chr'] = weights_df['chr'].str.replace('chr', '', case=False).str.strip()
        weights_df['pos'] = weights_df['pos'].str.strip()
        weights_df = weights_df[weights_df['chr'] == target_chrom]
        
        # Get unique positions
        weights_positions = set(zip(weights_df['chr'], weights_df['pos']))
        total_weights_positions = len(weights_positions)
        print(f"Total positions in weights file (chr{target_chrom}): {total_weights_positions}\n")
        
        # Display first 5
        print("First 5 parsed variants from weights file:")
        for idx, row in weights_df.head(5).iterrows():
            print(f"  Variant {idx+1}: Chromosome: {row['chr']}, Position: {row['pos']}, "
                  f"Effect Allele: {row['effect_allele']}, Weight: {row['weight']}, ID: {row['id']}")
        print("\n")
        return weights_df, weights_positions
    except Exception as e:
        print(f"Error reading weights file: {e}")
        return pd.DataFrame(), set()

def find_matches(weights_positions, bim_positions):
    """
    Finds matched positions and calculates percentage.
    Returns matched positions set.
    """
    print("Finding matched positions...")
    matched_positions = weights_positions.intersection(bim_positions)
    matched_count = len(matched_positions)
    total = len(weights_positions)
    percentage = (matched_count / total) * 100 if total > 0 else 0
    print(f"Total matched positions: {matched_count}")
    print(f"Percentage of matched positions: {percentage:.2f}% out of {total} total positions in weights file\n")
    return matched_positions, matched_count, percentage

def collect_examples(weights_df, matched_positions, bim_df, weights_positions, example_limit=5):
    """
    Collects example matches and non-matches.
    Returns lists of example matches and non-matches.
    """
    print("Collecting example matches and non-matches...")
    example_matches = []
    example_non_matches = []
    
    # For matches
    matched_df = weights_df[weights_df.set_index(['chr', 'pos']).index.isin(matched_positions)].copy()
    
    # Merge with BIM DataFrame
    if not bim_df.empty and not matched_df.empty:
        # Set index for faster lookup
        bim_df_indexed = bim_df.set_index(['chrom_variant', 'pos_variant'])
        matched_df_indexed = matched_df.set_index(['chr', 'pos'])
        
        # Rename bim_df_indexed index to match matched_df_indexed
        bim_df_indexed.index.names = ['chr', 'pos']
        
        # Perform join
        merged_df = matched_df_indexed.join(bim_df_indexed, how='left', lsuffix='_weights', rsuffix='_bim')
        
        # Reset index to access 'chr' and 'pos' as columns
        merged_df = merged_df.reset_index()
        
        # Replace alleles with 'X' using regex
        merged_df['allele1_repl'] = merged_df['allele1'].str.replace('[ACTG]', 'X', regex=True)
        merged_df['allele2_repl'] = merged_df['allele2'].str.replace('[ACTG]', 'X', regex=True)
        
        # Create variant string with X
        merged_df['variant_repl'] = merged_df['chr'] + ':' + merged_df['pos'] + ':' + merged_df['allele1_repl'] + ':' + merged_df['allele2_repl']
        
        # Take first example_limit rows
        example_matches = merged_df.head(example_limit).to_dict('records')
    
    # For non-matches
    non_matched_positions = weights_positions - matched_positions
    non_matched_df = weights_df[weights_df.set_index(['chr', 'pos']).index.isin(non_matched_positions)]
    example_non_matches = non_matched_df.head(example_limit).to_dict('records')
    
    return example_matches, example_non_matches

def main():
    weights_file = "pgs003725_processed_weights.csv"
    bim_file = "acaf_threshold.chr22.bim"
    target_chrom = '22'
    
    # Check if files exist
    if not os.path.isfile(weights_file):
        print(f"Error: Weights file '{weights_file}' does not exist.")
        return
    if not os.path.isfile(bim_file):
        print(f"Error: BIM file '{bim_file}' does not exist.")
        return
    
    # Read BIM file
    bim_positions, bim_df = read_bim_file(bim_file, target_chrom=target_chrom)
    
    # Read weights file
    weights_df, weights_positions = read_weights_file(weights_file, target_chrom=target_chrom)
    
    if not weights_positions:
        print(f"No positions found in weights file for chromosome {target_chrom}. Exiting.")
        return
    
    # Find matches
    matched_positions, matched_count, percentage = find_matches(weights_positions, bim_positions)
    
    # Collect examples
    example_matches, example_non_matches = collect_examples(weights_df, matched_positions, bim_df, weights_positions, example_limit=5)
    
    # Print results
    print("Processing complete.")
    print(f"Total matched positions: {matched_count}")
    print(f"Percentage of matched positions: {percentage:.2f}% out of {len(weights_positions)} total positions in weights file\n")
    
    # Example matches
    if example_matches:
        print("Example Matches:")
        for match in example_matches:
            print(f"  Weights File - Chromosome: {match['chr']}, Position: {match['pos']}, "
                  f"Effect Allele: {match['effect_allele']}, Weight: {match['weight']}, ID: {match['id']}")
            print(f"  BIM File - Chromosome: {match['chr']}, Position: {match['pos']}, Variant: {match['variant_repl']}\n")
    else:
        print("No example matches found.\n")
    
    # Example non-matches
    if example_non_matches:
        print("Example Non-Matches (from Weights File):")
        for non_match in example_non_matches:
            print(f"  Chromosome: {non_match['chr']}, Position: {non_match['pos']}, "
                  f"Effect Allele: {non_match['effect_allele']}, Weight: {non_match['weight']}, "
                  f"ID: {non_match['id']}, Match Status: No Match\n")
    else:
        print("No example non-matches found.\n")

if __name__ == "__main__":
    main()
