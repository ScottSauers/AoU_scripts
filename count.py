import pandas as pd
import psutil
import os

def get_optimal_chunksize(file_path, memory_fraction=0.1):
    available_memory = psutil.virtual_memory().available
    file_size = os.path.getsize(file_path)
    estimated_row_size = file_size / sum(1 for _ in open(file_path))
    return int((available_memory * memory_fraction) / estimated_row_size)

def count_common_variants(weights_file, bim_file):
    weights_data = pd.read_csv(weights_file, usecols=[0, 1], names=["chr", "pos"], skiprows=1)
    weights_data["chr"] = weights_data["chr"].astype(str).str.replace("chr", "")
    weights_data["pos"] = weights_data["pos"].astype(str).str.strip()
    total_weights_variants = len(weights_data)

    total_lines = sum(1 for _ in open(bim_file))
    chunksize = get_optimal_chunksize(bim_file)
    common_count = 0
    processed_lines = 0

    print("Starting to process BIM file in chunks...")

    for i, chunk in enumerate(pd.read_csv(bim_file, delim_whitespace=True, header=None, usecols=[0, 3], names=["chrom", "pos"], chunksize=chunksize)):
        chunk["chrom"] = chunk["chrom"].astype(str).str.replace("chr", "").str.strip()
        chunk["pos"] = chunk["pos"].astype(str).str.strip()
        merged = pd.merge(weights_data, chunk, left_on=["chr", "pos"], right_on=["chrom", "pos"])
        chunk_common_count = len(merged)
        common_count += chunk_common_count
        processed_lines += len(chunk)
        progress = (processed_lines / total_lines) * 100
        print(f"Chunk {i+1}: {progress:.2f}% complete, found {chunk_common_count} common variants in this chunk, cumulative common variants: {common_count}")

    common_percentage = (common_count / total_weights_variants) * 100
    print("\nProcessing complete.")
    print("Total common variants:", common_count)
    print(f"Percentage of common variants: {common_percentage:.2f}% out of {total_weights_variants} total variants in score file")

count_common_variants("pgs003725_processed_weights.csv", "acaf_threshold.chr22.bim")
