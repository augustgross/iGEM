from Bio import SeqIO
import numpy as np
import csv

def generate_normalized_kmer_frequencies(sequence, k_range):
    kmer_counts = {}
    total_kmers = 0
    
    # Generate and count k-mers
    for k in k_range:
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
            total_kmers += 1
    
    # Normalize k-mer frequencies
    for kmer in kmer_counts:
        kmer_counts[kmer] /= total_kmers
    
    return kmer_counts

def build_global_kmer_index(sequences, k_range):
    global_kmer_index = {}
    for sequence in sequences:
        for k in k_range:
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                if kmer not in global_kmer_index:
                    global_kmer_index[kmer] = len(global_kmer_index)
    return global_kmer_index

def vectorize_sequence(sequence, global_kmer_index, k_range):
    vector = np.zeros(len(global_kmer_index))
    kmer_frequencies = generate_normalized_kmer_frequencies(sequence, k_range)
    
    for kmer, frequency in kmer_frequencies.items():
        if kmer in global_kmer_index:
            vector[global_kmer_index[kmer]] = frequency
            
    return vector

def process_fasta_files(fasta_files):
    sequences = []
    for fasta_file in fasta_files:
        sequences.extend([str(record.seq).upper() for record in SeqIO.parse(fasta_file, "fasta")])
    return sequences

def main(fasta_files):
    k_range = range(1, 7)  # K-mers from 1 to 6
    sequences = process_fasta_files(fasta_files)
    global_kmer_index = build_global_kmer_index(sequences, k_range)
    
    # Vectorize sequences
    vectors = np.array([vectorize_sequence(seq, global_kmer_index, k_range) for seq in sequences])
    
    # Reshape for CNN input
    vectors = vectors.reshape(vectors.shape[0], vectors.shape[1], 1)
    
    return vectors

fasta_files = [
    "fasta/site.02.subj.0001.lab.2014222001.iso.1.fasta",
    "fasta/site.02.subj.0002.lab.2014222005.iso.1.fasta",
    "fasta/site.02.subj.0004.lab.2014222010.iso.1.fasta",
    "fasta/site.02.subj.0005.lab.2014222011.iso.1.fasta",
    "fasta/site.02.subj.0006.lab.2014222013.iso.1.fasta",
    "fasta/site.02.subj.0007.lab.2014222016.iso.1.fasta",
    "fasta/site.02.subj.0008.lab.2014222017.iso.1.fasta",
    "fasta/site.02.subj.0009.lab.2014222037.iso.1.fasta",
    "fasta/site.02.subj.0010.lab.2014222040.iso.1.fasta",
    "fasta/site.02.subj.0011.lab.2014222046.iso.1.fasta"
]


# Use the main function with your FASTA file
sequence_vectors = main(fasta_files)

print(sequence_vectors.shape)

# Save to CSV
with open('output.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(sequence_vectors.reshape(sequence_vectors.shape[0], -1))  # Reshape back for CSV writing