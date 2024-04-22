import pandas as pd
import wget
import gzip
import shutil
import os

# Path to the dataset files
reuse_table_path = 'CRyPTIC_reuse_table_20231208.csv'
excluded_samples_path = 'CRyPTIC_excluded_samples_20220607.tsv'
# label_columns = [
#     'AMI_BINARY_PHENOTYPE', 'BDQ_BINARY_PHENOTYPE', 'CFZ_BINARY_PHENOTYPE', 
#     'DLM_BINARY_PHENOTYPE', 'EMB_BINARY_PHENOTYPE', 'ETH_BINARY_PHENOTYPE', 
#     'INH_BINARY_PHENOTYPE', 'KAN_BINARY_PHENOTYPE', 'LEV_BINARY_PHENOTYPE', 
#     'LZD_BINARY_PHENOTYPE', 'MXF_BINARY_PHENOTYPE', 'RIF_BINARY_PHENOTYPE', 
#     'RFB_BINARY_PHENOTYPE'
# ]
label_column = 'AMI_BINARY_PHENOTYPE'
# Load the main reuse table
df = pd.read_csv(reuse_table_path)

# Load the excluded samples list
excluded_samples = pd.read_csv(excluded_samples_path, sep='\t')
df = df[~df['UNIQUEID'].isin(excluded_samples['UNIQUEID'])]
df = df.dropna(subset=[label_column])

# Reset the DataFrame index to ensure it's sequential after filtering
df.reset_index(drop=True, inplace=True)

base = "https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/reproducibility/"

# Ensure the directory exists before attempting to write files into it
os.makedirs('VCF', exist_ok=True)

# Download and extract the VCF files
for i in range(min(10, len(df))):
    sample = df['VCF'][i]
    file_path = df['UNIQUEID'][i] + '.vcf.gz'
    wget.download(base + sample, file_path)
    
    # Unzip the file directly into the intended directory
    with gzip.open(file_path, 'rb') as f_in:
        with open(os.path.join('VCF', df['UNIQUEID'][i] + '.vcf'), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file_path)
