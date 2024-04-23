import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix, roc_curve, roc_auc_score
from sklearn.preprocessing import StandardScaler, LabelBinarizer
from sklearn.multioutput import MultiOutputClassifier
import joblib

# Paths to data files
reuse_table_path = 'CRyPTIC_reuse_table_20231208.csv'
excluded_samples_path = 'CRyPTIC_excluded_samples_20220607.tsv'
feature_matrix_path = 'features.csv'

# Load the data
df = pd.read_csv(reuse_table_path)
excluded_samples = pd.read_csv(excluded_samples_path, sep='\t')
df_features = pd.read_csv(feature_matrix_path)

# Filter out excluded samples
df = df[~df['UNIQUEID'].isin(excluded_samples['UNIQUEID'])]

# Specify phenotype columns
label_columns = [
    'AMI_BINARY_PHENOTYPE', 'BDQ_BINARY_PHENOTYPE', 'CFZ_BINARY_PHENOTYPE', 
    'DLM_BINARY_PHENOTYPE', 'EMB_BINARY_PHENOTYPE', 'ETH_BINARY_PHENOTYPE', 
    'INH_BINARY_PHENOTYPE', 'KAN_BINARY_PHENOTYPE', 'LEV_BINARY_PHENOTYPE', 
    'LZD_BINARY_PHENOTYPE', 'MXF_BINARY_PHENOTYPE', 'RIF_BINARY_PHENOTYPE', 
    'RFB_BINARY_PHENOTYPE'
]

# Clean and binarize phenotypic data
df = df.dropna(subset=label_columns)
for column in label_columns:
    df = df[df[column].isin(['S', 'R'])]
    label_binarizer = LabelBinarizer()
    df[column] = label_binarizer.fit_transform(df[column].values)

# Prepare labels DataFrame
df_labels = df[['UNIQUEID'] + label_columns]

# Adjust feature_matrix sample IDs to match the reuse_table
df_features['sample_id_base'] = df_features['sample_id'].str[:-1]  # Remove the last character (A or B)
df_features = df_features[df_features['sample_id_base'].isin(df['UNIQUEID'])]
# Merge features and labels on modified sample_id_base to UNIQUEID
df_merged = pd.merge(df_features, df_labels, left_on='sample_id_base', right_on='UNIQUEID', how='left')

# Prepare feature matrix and labels
X = df_merged.drop(['sample_id', 'sample_id_base', 'UNIQUEID'] + label_columns, axis=1)
y = df_merged[label_columns]

# Split and scale data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Multi-label Logistic Regression
model = MultiOutputClassifier(LogisticRegression(max_iter=2000))
model.fit(X_train_scaled, y_train)

# Predictions and evaluation
predictions = model.predict(X_test_scaled)
print(classification_report(y_test, predictions, target_names=label_columns))

probs = [estimator.predict_proba(X_test_scaled)[:, 1] for estimator in model.estimators_]

# Confusion matrices and ROC curves
for idx, label in enumerate(label_columns):
    print(f"--- {label} ---")
    print(confusion_matrix(y_test.iloc[:, idx], predictions[:, idx]))
    print(classification_report(y_test.iloc[:, idx], predictions[:, idx]))

    # ROC Curve
    fpr, tpr, thresholds = roc_curve(y_test.iloc[:, idx], probs[idx])
    auc = roc_auc_score(y_test.iloc[:, idx], probs[idx])
    plt.figure()
    plt.plot(fpr, tpr, label=f'ROC curve (area = {auc:.2f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'Receiver operating characteristic: {label}')
    plt.legend(loc="lower right")
    plt.show()

# Save model and scaler
joblib.dump(model, 'model.pkl')
joblib.dump(scaler, 'scaler.pkl')