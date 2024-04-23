from flask import Flask, request, jsonify
import joblib
from Bio import SeqIO
from io import StringIO
from collections import Counter
import pandas as pd
import numpy as np

app = Flask(__name__)

# Load pre-trained model and scaler
model = joblib.load('model.pkl')
scaler = joblib.load('scaler.pkl')

# Load feature order from feature_matrix.csv, excluding the last column if it's a non-feature
df_features = pd.read_csv('feature_matrix.csv')
sorted_kmers = df_features.columns.tolist()[:-1]  # Excludes the last column, assumed non-feature

def compute_kmers(sequence, k_size_range=(1, 7)):
    kmers = []
    for k_size in range(*k_size_range):
        kmers.extend([str(sequence[i:i+k_size]) for i in range(len(sequence) - k_size + 1)])
    return kmers

@app.route('/upload', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'GET':
        return '''
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <title>Upload FASTA File</title>
        </head>
        <body>
            <h1>Upload FASTA File for Antibiotic Resistance Prediction</h1>
            <form action="/upload" method="post" enctype="multipart/form-data">
                <input type="file" name="file" required>
                <button type="submit">Upload and Predict</button>
            </form>
        </body>
        </html>
        '''
    elif request.method == 'POST':
        file = request.files['file']
        file_stream = StringIO(file.read().decode('utf-8'))
        sequence = next(SeqIO.parse(file_stream, "fasta")).seq
        kmers = compute_kmers(sequence)
        kmer_counts = Counter(kmers)

        # Convert kmer_counts to a DataFrame with feature names for scaling
        features = np.array([kmer_counts.get(kmer, 0) for kmer in sorted_kmers])
        features_df = pd.DataFrame([features], columns=sorted_kmers)
        features_scaled = scaler.transform(features_df)

        # Define label columns for output
        label_columns = [
            'Amikacin', 'Bedaquiline', 'Clofazimine', 
            'Delamanid', 'Ethambutol', 'Ethionamide', 
            'Isoniazid', 'Kanamycin', 'Levofloxacin', 
            'Linezolid', 'Moxifloxacin', 'Rifampicin', 
            'Rifabutin'
        ]
        
        # Predict probabilities for each antibiotic resistance
        predictions = model.predict_proba(features_scaled)
        results = {label: float(pred[0, 1]) for label, pred in zip(label_columns, predictions)}

        # Return predictions as JSON
        return jsonify(results)

if __name__ == '__main__':
    app.run(debug=True)




