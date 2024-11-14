import pandas as pd
import os
from tqdm import tqdm

def read_bim_file(bim_file):
    """
    Reads the BIM file and returns a set of (chrom, pos) tuples and a DataFrame.
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
        # Extract position from 'variant' (assuming format 'chr:pos:allele1:allele2')
        split_variants = bim_df['variant'].str.split(':', expand=True)
        bim_df['chrom_variant'] = split_variants[0].str.replace('chr', '', case=False).str.strip()
        bim_df['pos_variant'] = split_variants[1].str.strip()
        # Create set of (chrom, pos)
        bim_positions = set(zip(bim_df['chrom_variant'], bim_df['pos_variant']))
        print(f"Total positions in BIM file: {len(bim_positions)}\n")
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

def collect_examples(weights_df, matched_positions, bim_df, example_limit=5):
    """
    Collects example matches and non-matches.
    Returns lists of example matches and non-matches.
    """
    print("Collecting example matches and non-matches...")
    example_matches = []
    example_non_matches = []
    
    # For matches
    if matched_positions:
        matched_list = list(matched_positions)[:example_limit]
        for variant in matched_list:
            chr_, pos = variant
            # Fetch corresponding row from weights_df
            weight_row = weights_df[(weights_df['chr'] == chr_) & (weights_df['pos'] == pos)].iloc[0]
            # Fetch corresponding row from bim_df
            bim_row = bim_df[(bim_df['chrom_variant'] == chr_) & (bim_df['pos_variant'] == pos)].iloc[0]
            # Replace alleles with 'X' for privacy
            allele1_repl = bim_row['allele1'].replace('[ACTG]', 'X', regex=True)
            allele2_repl = bim_row['allele2'].replace('[ACTG]', 'X', regex=True)
            variant_repl = f"{bim_row['chrom_variant']}:{bim_row['pos_variant']}:{allele1_repl}:{allele2_repl}"
            example_matches.append({
                'weights_chr': weight_row['chr'],
                'weights_pos': weight_row['pos'],
                'weights_effect_allele': weight_row['effect_allele'],
                'weights_weight': weight_row['weight'],
                'weights_id': weight_row['id'],
                'bim_chrom': bim_row['chrom'],
                'bim_pos': bim_row['pos_variant'],
                'bim_variant': variant_repl
            })
    
    # For non-matches
    non_matched_positions = weights_positions - matched_positions
    non_matched_list = list(non_matched_positions)[:example_limit]
    for variant in non_matched_list:
        chr_, pos = variant
        # Fetch corresponding row from weights_df
        weight_row = weights_df[(weights_df['chr'] == chr_) & (weights_df['pos'] == pos)].iloc[0]
        example_non_matches.append({
            'weights_chr': weight_row['chr'],
            'weights_pos': weight_row['pos'],
            'weights_effect_allele': weight_row['effect_allele'],
            'weights_weight': weight_row['weight'],
            'weights_id': weight_row['id'],
            'match_status': 'No Match'
        })
    
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
    bim_positions, bim_df = read_bim_file(bim_file)
    
    # Read weights file
    weights_df, weights_positions = read_weights_file(weights_file, target_chrom=target_chrom)
    
    if not weights_positions:
        print("No positions found in weights file for chromosome 22. Exiting.")
        return
    
    # Find matches
    matched_positions, matched_count, percentage = find_matches(weights_positions, bim_positions)
    
    # Collect examples
    example_matches, example_non_matches = collect_examples(weights_df, matched_positions, bim_df, example_limit=5)
    
    # Print results
    print("Processing complete.")
    print(f"Total matched positions: {matched_count}")
    print(f"Percentage of matched positions: {percentage:.2f}% out of {len(weights_positions)} total positions in weights file\n")
    
    # Example matches
    if example_matches:
        print("Example Matches:")
        for match in example_matches:
            print(f"  Weights File - Chromosome: {match['weights_chr']}, Position: {match['weights_pos']}, "
                  f"Effect Allele: {match['weights_effect_allele']}, Weight: {match['weights_weight']}, ID: {match['weights_id']}")
            print(f"  BIM File - Chromosome: {match['bim_chrom']}, Position: {match['bim_pos']}, Variant: {match['bim_variant']}\n")
    else:
        print("No example matches found.\n")
    
    # Example non-matches
    if example_non_matches:
        print("Example Non-Matches (from Weights File):")
        for non_match in example_non_matches:
            print(f"  Chromosome: {non_match['weights_chr']}, Position: {non_match['weights_pos']}, "
                  f"Effect Allele: {non_match['weights_effect_allele']}, Weight: {non_match['weights_weight']}, "
                  f"ID: {non_match['weights_id']}, Match Status: {non_match['match_status']}\n")
    else:
        print("No example non-matches found.\n")

if __name__ == "__main__":
    main()
