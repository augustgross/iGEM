import pandas as pd
import os
import glob

# Assuming reuse_table_path and excluded_samples_path are already set and loaded as df
reuse_table_path = 'CRyPTIC_reuse_table_20231208.csv'
excluded_samples_path = 'CRyPTIC_excluded_samples_20220607.tsv'

# Load the main reuse table and excluded samples
df = pd.read_csv(reuse_table_path)
excluded_samples = pd.read_csv(excluded_samples_path, sep='\t')

# Filter the DataFrame
df = df[~df['UNIQUEID'].isin(excluded_samples['UNIQUEID'])]
df = df.dropna(subset=['INH_BINARY_PHENOTYPE'])

# Function to read k-mer files for a given sample ID and merge them into a single DataFrame
def read_and_merge_kmer_files(sample_id, kmer_base_path='outputs'):
    all_kmers = pd.DataFrame()
    # List all k-mer files for the sample
    file_pattern = os.path.join(kmer_base_path, f"{sample_id}*", f"{sample_id}.kmer*A.csv")
    file_paths = glob.glob(file_pattern)
    
    for file_path in file_paths:
        if os.path.isfile(file_path):
            kmer_data = pd.read_csv(file_path, header=None)
            kmer_data.columns = ['kmer', 'count']
            kmer_data['kmer_type'] = os.path.basename(file_path).split('.')[1]  # Extract kmer type from file name
            # Pivot k-mer counts into columns
            kmer_data = kmer_data.pivot_table(index='kmer_type', columns='kmer', values='count', aggfunc='sum').fillna(0)
            if all_kmers.empty:
                all_kmers = kmer_data
            else:
                all_kmers = pd.concat([all_kmers, kmer_data])

    # Sum counts across k-mer types for each sample
    all_kmers = all_kmers.groupby(all_kmers.index).sum()
    return all_kmers

def create_feature_matrix(sample_ids, kmer_base_path='outputs'):
    feature_matrix = pd.DataFrame()

    for sample_id in sample_ids:
        sample_kmers = read_and_merge_kmer_files(sample_id, kmer_base_path)
        # Flatten DataFrame so k-mers are columns, append sample_id as a column
        sample_features = sample_kmers.reset_index(drop=True)
        sample_features['sample_id'] = sample_id  # add sample_id for reference
        feature_matrix = pd.concat([feature_matrix, sample_features], ignore_index=True)

    return feature_matrix

# List of unique sample IDs and the path to k-mer count files
sample_ids = df['UNIQUEID'].tolist()
kmer_base_path = "outputs"

# Generate the feature matrix
feature_matrix = create_feature_matrix(sample_ids, kmer_base_path)

# Output the feature matrix to check
print(feature_matrix.head())

#print data to csv
feature_matrix.to_csv('feature_matrix.csv', index=False)
