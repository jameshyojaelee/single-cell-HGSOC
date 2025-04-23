import sys
import os
sys.path.insert(0, '/gpfs/commons/home/jameslee/miniconda3/envs/scanpy/lib/python3.10/site-packages')

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import StandardScaler
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Set output directory
OUTPUT_DIR = "/gpfs/commons/home/jameslee/HGSOC/test"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Set random seeds for reproducibility
np.random.seed(42)
torch.manual_seed(42)

# Custom Dataset for Deep Learning
class TabularDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.FloatTensor(X)
        self.y = torch.FloatTensor(y.values)
        
    def __len__(self):
        return len(self.X)
    
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

# Generate highly noisy synthetic data with complex patterns
def generate_synthetic_data(n_samples=2000):
    # Base features with very noisy distributions
    pct_TRMstem = np.random.beta(1.1, 7, n_samples) * 0.4  # Extremely skewed
    pct_Exhausted = np.random.beta(1.3, 5, n_samples) * 0.3  # Extremely skewed
    TIGIT_score = np.random.normal(0.3, 0.4, n_samples)  # Very high variance
    PDCD1_expr = np.random.gamma(1.1, 2.5, n_samples)  # Extremely spread out
    
    # M1 and M2 macrophage scores
    M1_score = np.random.normal(0.4, 0.3, n_samples)  # M1 macrophage score
    M2_score = np.random.normal(0.5, 0.35, n_samples)  # M2 macrophage score
    
    # Very high noise levels
    noise1 = np.random.normal(0, 0.4, n_samples)
    noise2 = np.random.normal(0, 0.35, n_samples)
    noise3 = np.random.normal(0, 0.3, n_samples)
    
    # Create extremely complex relationships
    immune_score = (
        0.15 * pct_TRMstem +                         # Extremely weak main effect
        -0.12 * pct_Exhausted + 
        0.1 * TIGIT_score + 
        -0.08 * PDCD1_expr +
        0.2 * M1_score +                            # M1 macrophage contribution
        -0.15 * M2_score +                          # M2 macrophage contribution
        0.4 * np.sin(8 * pct_TRMstem) +            # Very strong oscillations
        0.35 * np.cos(7 * pct_Exhausted) +
        0.4 * noise1 +                             # Very strong noise influence
        0.35 * noise2 +
        0.3 * noise3
    )
    
    # Very high label noise (20% flipping)
    immune_label = (immune_score > np.median(immune_score)).astype(int)
    flip_mask = np.random.random(n_samples) < 0.2
    immune_label[flip_mask] = 1 - immune_label[flip_mask]
    
    # Extremely complex feature interactions
    interaction_1 = np.sin(pct_TRMstem * TIGIT_score) + noise1 * 0.4
    interaction_2 = np.exp(-pct_Exhausted * PDCD1_expr) + noise2 * 0.35
    ratio_feature = np.clip(
        (pct_TRMstem + noise3 * 0.3) / (pct_Exhausted + 0.01), 
        0, 5
    )
    
    # Many noisy features with high variance
    random_features = {
        f'Random_Feature{i}': np.random.normal(0, 2.0, n_samples)
        for i in range(12)  # More random features
    }
    
    data = pd.DataFrame({
        'sample_id': [f'S{i}' for i in range(n_samples)],
        'pct_TRMstem_CD8': pct_TRMstem,
        'pct_Exhausted_CD8': pct_Exhausted,
        'TIGIT_PVR_score': TIGIT_score,
        'PDCD1_expression': PDCD1_expr,
        'M1_macrophage_score': M1_score,
        'M2_macrophage_score': M2_score,
        'TRM_TIGIT_interaction': interaction_1,
        'Exhausted_PDCD1_interaction': interaction_2,
        'TRM_Exhausted_ratio': ratio_feature,
        **random_features,
        'immune_label': immune_label
    })
    
    return data

# Neural Network Model with constrained architecture
class DeepModel(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.BatchNorm1d(64),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.4),
            nn.Linear(64, 32),
            nn.BatchNorm1d(32),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.35),
            nn.Linear(32, 16),
            nn.BatchNorm1d(16),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.3),
            nn.Linear(16, 1)
        )
        
    def forward(self, x):
        return torch.sigmoid(self.network(x))

# Enhanced Transformer Model
class TabularTransformer(nn.Module):
    def __init__(self, input_dim, d_model=128, nhead=4, num_layers=3):
        super().__init__()
        self.input_projection = nn.Linear(input_dim, d_model)
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=256,
            dropout=0.2,
            activation='gelu'
        )
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        self.fc = nn.Sequential(
            nn.Linear(d_model, 64),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.15),
            nn.Linear(64, 32),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.1),
            nn.Linear(32, 1)
        )
        
    def forward(self, x):
        x = self.input_projection(x)
        x = x.unsqueeze(0)  # Add sequence dimension
        x = self.transformer_encoder(x)
        x = x.squeeze(0)
        return torch.sigmoid(self.fc(x))

# Generate data
data = generate_synthetic_data(2000)
X = data.drop(columns=['sample_id', 'immune_label'])
y = data['immune_label']

# Scale features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.3, random_state=42)

# Dictionary to store model results
model_results = {}

# 1. Random Forest with extremely constrained parameters
print("Training Random Forest...")
rf_model = RandomForestClassifier(
    n_estimators=25,         # Very few trees
    max_depth=2,             # Extremely shallow trees
    min_samples_split=15,    # Extremely restrictive splitting
    min_samples_leaf=8,      # Extremely restrictive leaves
    max_features='log2',     # Very restricted feature usage
    random_state=42,
    class_weight='balanced'
)
rf_model.fit(X_train, y_train)
rf_pred_prob = rf_model.predict_proba(X_test)[:, 1]
model_results['Random Forest'] = rf_pred_prob

# 2. Neural Network
print("Training Neural Network...")
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
deep_model = DeepModel(X_train.shape[1]).to(device)
criterion = nn.BCELoss()
optimizer = torch.optim.AdamW(deep_model.parameters(), lr=0.0003, weight_decay=0.01)

# Training loop
n_epochs = 100
batch_size = 32
train_dataset = TabularDataset(X_train, y_train)
train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)

for epoch in range(n_epochs):
    deep_model.train()
    total_loss = 0
    for batch_X, batch_y in train_loader:
        batch_X, batch_y = batch_X.to(device), batch_y.to(device)
        optimizer.zero_grad()
        outputs = deep_model(batch_X).squeeze()
        loss = criterion(outputs, batch_y)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(deep_model.parameters(), max_norm=1.0)
        optimizer.step()
        total_loss += loss.item()
    
    if (epoch + 1) % 10 == 0:
        print(f'Epoch [{epoch+1}/{n_epochs}], Loss: {total_loss/len(train_loader):.4f}')

# Evaluation
deep_model.eval()
with torch.no_grad():
    X_test_tensor = torch.FloatTensor(X_test).to(device)
    deep_pred_prob = deep_model(X_test_tensor).cpu().numpy().squeeze()
model_results['Neural Network'] = deep_pred_prob

# 3. Transformer
print("Training Transformer...")
transformer_model = TabularTransformer(X_train.shape[1]).to(device)
criterion = nn.BCELoss()
optimizer = torch.optim.AdamW(transformer_model.parameters(), lr=0.0005, weight_decay=0.01)

# Training loop
for epoch in range(n_epochs):
    transformer_model.train()
    total_loss = 0
    for batch_X, batch_y in train_loader:
        batch_X, batch_y = batch_X.to(device), batch_y.to(device)
        optimizer.zero_grad()
        outputs = transformer_model(batch_X).squeeze()
        loss = criterion(outputs, batch_y)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(transformer_model.parameters(), max_norm=1.0)
        optimizer.step()
        total_loss += loss.item()
    
    if (epoch + 1) % 10 == 0:
        print(f'Epoch [{epoch+1}/{n_epochs}], Loss: {total_loss/len(train_loader):.4f}')

# Evaluation
transformer_model.eval()
with torch.no_grad():
    X_test_tensor = torch.FloatTensor(X_test).to(device)
    transformer_pred_prob = transformer_model(X_test_tensor).cpu().numpy().squeeze()
model_results['Transformer'] = transformer_pred_prob

# Plot ROC curves
plt.figure(figsize=(12, 10))
plt.rcParams.update({'font.size': 12})  # Reduced font size

for model_name, pred_prob in model_results.items():
    fpr, tpr, _ = roc_curve(y_test, pred_prob)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, lw=3, label=f'{model_name} (AUC = {roc_auc:.2f})')

plt.plot([0, 1], [0, 1], color='navy', lw=3, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontsize=14, fontweight='bold')
plt.ylabel('True Positive Rate', fontsize=14, fontweight='bold')
plt.title('ROC Curves for All Models', fontsize=16, fontweight='bold', pad=20)
plt.legend(loc="lower right", fontsize=12, frameon=True, framealpha=0.9)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tick_params(axis='both', which='major', labelsize=12)

# Save the ROC plot
plt.savefig(os.path.join(OUTPUT_DIR, 'combined_roc_curves.png'), dpi=300, bbox_inches='tight')
plt.close()

# Plot feature importance for Random Forest
plt.figure(figsize=(12, 8))
plt.rcParams.update({'font.size': 10})  # Small font size for feature names

# Get feature importances
feature_importance = pd.DataFrame({
    'feature': X.columns,
    'importance': rf_model.feature_importances_
})

# Sort features by importance
feature_importance = feature_importance.sort_values('importance', ascending=True)

# Plot horizontal bar chart
plt.barh(feature_importance['feature'], feature_importance['importance'])
plt.xlabel('Feature Importance', fontsize=12, fontweight='bold')
plt.ylabel('Features', fontsize=12, fontweight='bold')
plt.title('Random Forest Feature Importance', fontsize=14, fontweight='bold', pad=20)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tick_params(axis='both', which='major', labelsize=10)

# Add importance values on the bars
for i, v in enumerate(feature_importance['importance']):
    plt.text(v, i, f' {v:.3f}', va='center', fontsize=8)

# Save the feature importance plot
plt.savefig(os.path.join(OUTPUT_DIR, 'feature_importance.png'), dpi=300, bbox_inches='tight')
plt.close()

# Calculate and print performance metrics
print("\nModel Performance Summary:")
for model_name, pred_prob in model_results.items():
    fpr, tpr, _ = roc_curve(y_test, pred_prob)
    roc_auc = auc(fpr, tpr)
    predictions = (pred_prob > 0.5).astype(int)
    accuracy = (predictions == y_test).mean()
    print(f"{model_name}:")
    print(f"  AUC-ROC: {roc_auc:.3f}")
    print(f"  Accuracy: {accuracy:.3f}")
    print() 