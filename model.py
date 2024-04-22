import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFECV
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

# Assuming 'UNIQUEID' is your index or unique identifier common to both datasets.
df_features = pd.read_csv('feature_matrix.csv')  # Features including 'sample_id'
df_labels = df[['UNIQUEID', 'INH_BINARY_PHENOTYPE']]  # Labels must have the same UNIQUEID

# Merge features and labels on 'sample_id' or 'UNIQUEID'
df_merged = pd.merge(df_features, df_labels, left_on='sample_id', right_on='UNIQUEID')

# Check if any rows have null values which may indicate mismatched data
print(df_merged.isnull().sum())

# Prepare your feature matrix and labels
X = df_merged.drop(['sample_id', 'UNIQUEID', 'INH_BINARY_PHENOTYPE'], axis=1)
y = df_merged['INH_BINARY_PHENOTYPE']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

model = LogisticRegression(max_iter=1000)
model.fit(X_train_scaled, y_train)

predictions = model.predict(X_test_scaled)
print(confusion_matrix(y_test, predictions))
print(classification_report(y_test, predictions))

# Recursive Feature Elimination with Cross-Validation
selector = RFECV(model, step=1, cv=5)
selector.fit(X_train_scaled, y_train)

# Find features considered important by the selector
important_features = X.columns[selector.support_]
print("Important features:", important_features)

# Redefine X to only use important features
X_train_important = selector.transform(X_train_scaled)
X_test_important = selector.transform(X_test_scaled)

# Re-train model
model.fit(X_train_important, y_train)
new_predictions = model.predict(X_test_important)
print(classification_report(y_test, new_predictions))

import joblib
joblib.dump(model, 'logistic_regression_model.pkl')
joblib.dump(scaler, 'scaler.pkl')
