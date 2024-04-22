import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
from Bio import SeqIO
import numpy as np
import pandas as pd

# Define the CNN model
class DNAConvNet(nn.Module):
    def __init__(self, sequence_length, num_classes):
        super(DNAConvNet, self).__init__()
        self.conv1 = nn.Conv1d(1, 16, kernel_size=5, padding=2)
        self.pool = nn.MaxPool1d(2)
        self.conv2 = nn.Conv1d(16, 32, kernel_size=5, padding=2)
        self.fc1 = nn.Linear(32 * (sequence_length // 2 // 2), 128)
        self.fc2 = nn.Linear(128, num_classes)

    def forward(self, x):
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        x = x.view(-1, 32 * (sequence_length // 2 // 2))
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return x

# Function to generate normalized k-mer frequencies
def generate_normalized_kmer_frequencies(sequence, k_range):
    kmer_counts = {}
    total_kmers = 0
    for k in k_range:
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
            total_kmers += 1
    for kmer in kmer_counts:
        kmer_counts[kmer] /= total_kmers
    return kmer_counts

# Function to build a global index of k-mers
def build_global_kmer_index(sequences, k_range):
    global_kmer_index = {}
    for sequence in sequences:
        for k in k_range:
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                if kmer not in global_kmer_index:
                    global_kmer_index[kmer] = len(global_kmer_index)
    return global_kmer_index

# Function to vectorize sequences based on k-mer frequencies
def vectorize_sequence(sequence, global_kmer_index, k_range):
    vector = np.zeros(len(global_kmer_index))
    kmer_frequencies = generate_normalized_kmer_frequencies(sequence, k_range)
    for kmer, frequency in kmer_frequencies.items():
        if kmer in global_kmer_index:
            vector[global_kmer_index[kmer]] = frequency
    return vector

# Function to process FASTA files and vectorize
def process_fasta_files(fasta_files):
    sequences = []
    for fasta_file in fasta_files:
        sequences.extend([str(record.seq).upper() for record in SeqIO.parse(fasta_file, "fasta")])
    return sequences

# Main function integrating data loading, model training, and evaluation
def main(fasta_files, reuse_table_path, excluded_samples_path, label_columns):
    # Load labels and sequence data
    df = pd.read_csv(reuse_table_path)
    excluded_samples = pd.read_csv(excluded_samples_path, sep='\t')
    df = df[~df['UNIQUEID'].isin(excluded_samples['UNIQUEID'])]
    df = df.dropna(subset=label_columns)
    # Map 'S' to 0, 'R' to 1
    label_mapping = {'S': 0, 'R': 1}
    for col in label_columns:
        df[col] = df[col].map(label_mapping)
    

    labels = df[label_columns].values.astype(float)
    
    sequences = process_fasta_files(fasta_files)
    if len(sequences) != len(labels):
        raise ValueError("The number of sequences and labels must be the same.")

    k_range = range(1, 7)
    global_kmer_index = build_global_kmer_index(sequences, k_range)
    vectors = np.array([vectorize_sequence(seq, global_kmer_index, k_range) for seq in sequences])
    vectors = torch.tensor(vectors, dtype=torch.float32).unsqueeze(1)

    if vectors.shape[0] != labels.shape[0]:
        raise ValueError(f"Mismatch in dataset sizes, vectors: {vectors.shape[0]}, labels: {labels.shape[0]}")

    labels = torch.tensor(labels, dtype=torch.float32)
    dataset = TensorDataset(vectors, labels)
    dataloader = DataLoader(dataset, batch_size=32, shuffle=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    sequence_length = vectors.shape[2]  # This assumes all sequences are the same length
    num_classes = labels.shape[1]
    model = DNAConvNet(sequence_length, num_classes).to(device)
    
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    # Training loop
    num_epochs = 10
    for epoch in range(num_epochs):
        model.train()
        for inputs, targets in dataloader:
            inputs, targets = inputs.to(device), targets.to(device)
            outputs = model(inputs)
            loss = criterion(outputs, targets)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')

    # Evaluation (example setup)
    model.eval()
    correct = 0
    total = 0
    with torch.no_grad():
        for inputs, labels in dataloader:
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)
            predicted = torch.sigmoid(outputs) > 0.5
            correct += (predicted == labels.byte()).all(1).sum().item()
            total += labels.size(0)
    print(f'Accuracy of the model on test data: {100 * correct / total}%')

# Define paths and label columns
reuse_table_path = 'CRyPTIC_reuse_table_20231208.csv'
excluded_samples_path = 'CRyPTIC_excluded_samples_20220607.tsv'
label_columns = ['AMI_BINARY_PHENOTYPE', 'BDQ_BINARY_PHENOTYPE', 'CFZ_BINARY_PHENOTYPE', 
                 'DLM_BINARY_PHENOTYPE', 'EMB_BINARY_PHENOTYPE', 'ETH_BINARY_PHENOTYPE', 
                 'INH_BINARY_PHENOTYPE', 'KAN_BINARY_PHENOTYPE', 'LEV_BINARY_PHENOTYPE', 
                 'LZD_BINARY_PHENOTYPE', 'MXF_BINARY_PHENOTYPE', 'RIF_BINARY_PHENOTYPE', 
                 'RFB_BINARY_PHENOTYPE']
fasta_files = [
    "fasta/site.02.subj.0002.lab.2014222005.iso.1.fasta",
    "fasta/site.02.subj.0005.lab.2014222011.iso.1.fasta",
    "fasta/site.02.subj.0006.lab.2014222013.iso.1.fasta",
    "fasta/site.02.subj.0007.lab.2014222016.iso.1.fasta",
    "fasta/site.02.subj.0008.lab.2014222017.iso.1.fasta",
    "fasta/site.02.subj.0009.lab.2014222037.iso.1.fasta",
    "fasta/site.02.subj.0010.lab.2014222040.iso.1.fasta",
    "fasta/site.02.subj.0011.lab.2014222046.iso.1.fasta",
    "fasta/site.02.subj.0012.lab.2014222053.iso.1.fasta",
    "fasta/site.02.subj.0013.lab.2014222055.iso.1.fasta"
]

main(fasta_files, reuse_table_path, excluded_samples_path, label_columns)