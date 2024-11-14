import pandas as pd
import psutil
import os
import csv
from tqdm import tqdm  # For advanced progress bars

def get_available_memory_fraction(memory_fraction=0.1):
    """Calculate available memory fraction for processing."""
    available_memory = psutil.virtual_memory().available
    return int(available_memory * memory_fraction)

def parse_weights_file(weights_file):
    """
    Parses the weights CSV file and returns a set of (chr, pos) tuples.
    Assumes the file has a header and one variant per line.
    """
    print("Parsing weights file...")
    weights_set = set()
    total_variants = 0
    example_variants = []

    with open(weights_file, 'r') as f:
        reader = csv.DictReader(f)
        for line_number, row in enumerate(reader, start=2):  # Start at 2 to account for header
            try:
                chr_ = row['chr'].strip()
                pos = row['pos'].strip()
                weights_set.add((chr_, pos))
                total_variants += 1

                # Store first 5 variants for example
                if len(example_variants) < 5:
                    example_variants.append((chr_, pos))
            except KeyError as e:
                print(f"Error: Missing expected column {e} in weights file at line {line_number}.")
                continue
            except Exception as e:
                print(f"Error parsing line {line_number}: {e}")
                continue

    print(f"Total variants parsed from weights file: {total_variants}\n")
    print("First 5 parsed variants from weights file:")
    for variant in example_variants:
        print(f"Chromosome: {variant[0]}, Position: {variant[1]}")
    
    return weights_set, total_variants

def count_common_variants(weights_set, bim_file, total_weights_variants):
    """
    Counts the number of common variants between weights_set and bim_file.
    Processes the BIM file in chunks for efficiency.
    """
    print("Starting to process BIM file in chunks...")
    common_count = 0
    example_bim_variants = []

    # Determine total number of lines in BIM file for progress tracking
    print("Calculating total number of lines in BIM file...")
    with open(bim_file, 'r') as f:
        total_lines = sum(1 for _ in f)
    print(f"Total lines in BIM file: {total_lines}\n")

    # Define chunk size based on available memory
    chunk_size = 100000  # Adjust based on system capabilities
    print(f"Processing BIM file with chunk size: {chunk_size}\n")

    # Initialize progress bar
    pbar = tqdm(total=total_lines, desc="Processing BIM File", unit="lines")

    # Read BIM file in chunks
    bim_columns = ["chrom", "variant", "ignore1", "ignore2", "allele1", "allele2"]
    try:
        bim_iter = pd.read_csv(
            bim_file,
            delim_whitespace=True,
            header=None,
            names=bim_columns,
            usecols=["chrom", "variant"],
            chunksize=chunk_size
        )
    except Exception as e:
        print(f"Error reading BIM file: {e}")
        return

    for i, chunk in enumerate(bim_iter, start=1):
        # Extract chromosome and position from the 'variant' column
        # Assuming 'variant' is in the format 'chr:pos:allele1:allele2'
        split_variants = chunk['variant'].str.split(':', expand=True)
        if split_variants.shape[1] < 2:
            print(f"Warning: Unexpected variant format in chunk {i}. Skipping chunk.")
            pbar.update(len(chunk))
            continue

        chunk['chrom'] = chunk['chrom'].str.replace('chr', '').str.strip()
        chunk['pos'] = split_variants[1].str.strip()

        # Store first 5 parsed BIM variants for example
        if i == 1:
            for idx, row in chunk.head(5).iterrows():
                allele_replaced = row['variant'].replace('A', 'X').replace('C', 'X').replace('T', 'X').replace('G', 'X')
                print(f"Example BIM Variant {idx+1}: Chromosome: {row['chrom']}, Position: {row['pos']}, Variant: {allele_replaced}")

        # Create a set of (chrom, pos) tuples for the current chunk
        bim_variants_set = set(zip(chunk['chrom'], chunk['pos']))

        # Calculate the intersection with weights_set
        common_in_chunk = bim_variants_set.intersection(weights_set)
        common_count += len(common_in_chunk)

        # Update progress bar
        pbar.update(len(chunk))

        # Detailed progress report every 10 chunks or at the end
        if i % 10 == 0 or i == (total_lines // chunk_size) + 1:
            progress = (pbar.n / total_lines) * 100
            print(f"Chunk {i}: {progress:.2f}% complete, found {len(common_in_chunk)} common variants in this chunk, cumulative common variants: {common_count}")

    pbar.close()
    print("\nProcessing complete.")
    print(f"Total common variants: {common_count}")
    common_percentage = (common_count / total_weights_variants) * 100
    print(f"Percentage of common variants: {common_percentage:.2f}% out of {total_weights_variants} total variants in weights file")

def main():
    weights_file = "pgs003725_processed_weights.csv"
    bim_file = "acaf_threshold.chr22.bim"

    # Check if files exist
    if not os.path.isfile(weights_file):
        print(f"Error: Weights file '{weights_file}' does not exist.")
        return
    if not os.path.isfile(bim_file):
        print(f"Error: BIM file '{bim_file}' does not exist.")
        return

    # Parse weights file
    weights_set, total_weights_variants = parse_weights_file(weights_file)
    if total_weights_variants == 0:
        print("No variants found in weights file. Exiting.")
        return

    # Count common variants
    count_common_variants(weights_set, bim_file, total_weights_variants)

if __name__ == "__main__":
    main()
