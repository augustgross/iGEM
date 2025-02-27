import pandas as pd
import os
import glob

# Load the main reuse table and excluded samples
reuse_table_path = 'CRyPTIC_reuse_table_20231208.csv'
excluded_samples_path = 'CRyPTIC_excluded_samples_20220607.tsv'
df = pd.read_csv(reuse_table_path)
excluded_samples = pd.read_csv(excluded_samples_path, sep='\t')

# Filter the DataFrame
df = df[~df['UNIQUEID'].isin(excluded_samples['UNIQUEID'])]
df = df.dropna(subset=['INH_BINARY_PHENOTYPE'])

# Function to read k-mer files for a given sample ID and merge them into a single DataFrame
def read_and_merge_kmer_files(sample_id, kmer_suffix, kmer_base_path='outputs'):
    all_kmers = pd.DataFrame()
    file_pattern = os.path.join(kmer_base_path, f"{sample_id}*", f"{sample_id}.kmer*{kmer_suffix}.csv")
    file_paths = glob.glob(file_pattern)
    
    for file_path in file_paths:
        if os.path.isfile(file_path):
            kmer_data = pd.read_csv(file_path, header=None)
            kmer_data.columns = ['kmer', 'count']
            kmer_data['kmer_type'] = os.path.basename(file_path).split('.')[1]  # Extract kmer type from file name
            kmer_data = kmer_data.pivot_table(index='kmer_type', columns='kmer', values='count', aggfunc='sum').fillna(0)
            if all_kmers.empty:
                all_kmers = kmer_data
            else:
                all_kmers = pd.concat([all_kmers, kmer_data])

    all_kmers = all_kmers.groupby(all_kmers.index).sum()
    return all_kmers

def create_feature_matrix(sample_ids, kmer_base_path='outputs'):
    feature_matrix = pd.DataFrame()

    for sample_id in sample_ids:
        # Process k-mer files with suffix 'A'
        sample_kmers_A = read_and_merge_kmer_files(sample_id, 'A', kmer_base_path)
        sample_kmers_A['sample_id'] = f"{sample_id}A"
        
        # Process k-mer files with suffix 'B'
        sample_kmers_B = read_and_merge_kmer_files(sample_id, 'B', kmer_base_path)
        sample_kmers_B['sample_id'] = f"{sample_id}B"

        # Append both A and B to the feature matrix
        feature_matrix = pd.concat([feature_matrix, sample_kmers_A, sample_kmers_B], ignore_index=True)

    return feature_matrix

# List of unique sample IDs and the path to k-mer count files
sample_ids = df['UNIQUEID'].tolist()
kmer_base_path = "outputs"

# Generate the combined feature matrix
feature_matrix = create_feature_matrix(sample_ids, kmer_base_path)

# Output the feature matrix to check
print(feature_matrix.head())

# Save data to csv
feature_matrix.to_csv('features.csv', index=False)
