import pandas as pd
import wget
import gzip
import shutil
import os

# Path to the dataset files
reuse_table_path = 'CRyPTIC_reuse_table_20231208.csv'
excluded_samples_path = 'CRyPTIC_excluded_samples_20220607.tsv'

# Load the main reuse table
df = pd.read_csv(reuse_table_path)

# Load the excluded samples list
excluded_samples = pd.read_csv(excluded_samples_path, sep='\t')
df = df[~df['UNIQUEID'].isin(excluded_samples['UNIQUEID'])]

base = "https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/reproducibility/"

for i in range(10):
    sample = df['VCF'][i]
    wget.download(base + sample, df['UNIQUEID'][i] + '.vcf.gz')
    
#unzip
#move the files to a new directory
os.mkdir('VCF')

for i in range(10):
    with gzip.open(df['UNIQUEID'][i] + '.vcf.gz', 'rb') as f_in:
        with open('VCF/' + df['UNIQUEID'][i] + '.vcf', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(df['UNIQUEID'][i] + '.vcf.gz')