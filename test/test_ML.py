import sys
import os
sys.path.insert(0, '/gpfs/commons/home/jameslee/miniconda3/envs/scanpy/lib/python3.10/site-packages')

# Import numpy first
import numpy as np

# Import other packages
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import roc_curve, auc, confusion_matrix, classification_report
from sklearn.preprocessing import StandardScaler
# import xgboost as xgb  # Commented out
# import lightgbm as lgb  # Commented out
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

# Import matplotlib with a workaround
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns

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
        self.y = torch.FloatTensor(y.values)  # Convert pandas Series to numpy array first
        
    def __len__(self):
        return len(self.X)
    
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

# Improved Neural Network Model
class ImprovedNN(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.batch_norm1 = nn.BatchNorm1d(input_dim)
        self.dense1 = nn.Linear(input_dim, 64)
        self.batch_norm2 = nn.BatchNorm1d(64)
        self.dense2 = nn.Linear(64, 32)
        self.batch_norm3 = nn.BatchNorm1d(32)
        self.dense3 = nn.Linear(32, 16)
        self.dense4 = nn.Linear(16, 1)
        self.dropout = nn.Dropout(0.3)
        self.activation = nn.LeakyReLU(0.1)
        
    def forward(self, x):
        x = self.batch_norm1(x)
        x = self.dense1(x)
        x = self.activation(x)
        x = self.batch_norm2(x)
        x = self.dropout(x)
        
        x = self.dense2(x)
        x = self.activation(x)
        x = self.batch_norm3(x)
        x = self.dropout(x)
        
        x = self.dense3(x)
        x = self.activation(x)
        x = self.dropout(x)
        
        x = self.dense4(x)
        return torch.sigmoid(x)

# Generate more meaningful synthetic data with clear patterns
def generate_synthetic_data(n_samples=1000):
    # Base features with more complex distributions
    pct_TRMstem = np.random.beta(1.5, 5, n_samples) * 0.6  # More skewed
    pct_Exhausted = np.random.beta(2, 3, n_samples) * 0.5  # More skewed
    TIGIT_score = np.random.normal(0.5, 0.25, n_samples)  # More variance
    PDCD1_expr = np.random.gamma(1.5, 1.5, n_samples)  # More spread out
    
    # Increased noise
    noise1 = np.random.normal(0, 0.2, n_samples)  # More noise
    noise2 = np.random.normal(0, 0.18, n_samples)  # More noise
    noise3 = np.random.normal(0, 0.15, n_samples)  # More noise
    
    # Create more complex relationships
    immune_score = (
        0.35 * pct_TRMstem +                         # Reduced main effect
        -0.3 * pct_Exhausted + 
        0.2 * TIGIT_score + 
        -0.15 * PDCD1_expr +
        0.2 * np.sin(5 * pct_TRMstem) +            # More oscillations
        0.15 * np.cos(4 * pct_Exhausted) +
        0.2 * noise1 +                             # More noise influence
        0.15 * noise2 +
        0.1 * noise3
    )
    
    # Add more random flipping of labels (8% of labels)
    immune_label = (immune_score > np.median(immune_score)).astype(int)
    flip_mask = np.random.random(n_samples) < 0.08  # More label noise
    immune_label[flip_mask] = 1 - immune_label[flip_mask]
    
    # Create more complex feature interactions
    interaction_1 = np.sin(pct_TRMstem * TIGIT_score) + noise1 * 0.2
    interaction_2 = np.exp(-pct_Exhausted * PDCD1_expr) + noise2 * 0.15
    ratio_feature = np.clip(
        (pct_TRMstem + noise3 * 0.15) / (pct_Exhausted + 0.01), 
        0, 5
    )
    
    # Add more noisy features
    random_feature1 = np.random.normal(0, 1.0, n_samples)
    random_feature2 = np.random.uniform(-2, 2, n_samples)
    random_feature3 = np.random.exponential(1.5, n_samples)
    random_feature4 = np.random.gamma(1.5, 2, n_samples)
    
    data = pd.DataFrame({
        'sample_id': [f'S{i}' for i in range(n_samples)],
        'pct_TRMstem_CD8': pct_TRMstem,
        'pct_Exhausted_CD8': pct_Exhausted,
        'TIGIT_PVR_score': TIGIT_score,
        'PDCD1_expression': PDCD1_expr,
        'TRM_TIGIT_interaction': interaction_1,
        'Exhausted_PDCD1_interaction': interaction_2,
        'TRM_Exhausted_ratio': ratio_feature,
        'Random_Feature1': random_feature1,
        'Random_Feature2': random_feature2,
        'Random_Feature3': random_feature3,
        'Random_Feature4': random_feature4,
        'immune_label': immune_label
    })
    
    return data

# Generate enhanced synthetic data with more samples
data = generate_synthetic_data(3000)  # Increased sample size for more complexity

# Features and labels
X = data.drop(columns=['sample_id', 'immune_label'])
y = data['immune_label']

# Scale the features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Train-test split with more test data
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.3, random_state=42)

# Dictionary to store model results
model_results = {}

# 1. Random Forest with constrained parameters
print("Training Random Forest...")
rf_model = RandomForestClassifier(
    n_estimators=100,        # Fewer trees
    max_depth=5,             # Shallower trees
    min_samples_split=8,     # More samples required to split
    min_samples_leaf=4,      # More samples required in leaves
    max_features='log2',     # More restricted feature usage
    random_state=42,
    class_weight='balanced'  # Handle class imbalance
)
rf_model.fit(X_train, y_train)
rf_pred_prob = rf_model.predict_proba(X_test)[:, 1]
model_results['Random Forest'] = rf_pred_prob

# # 2. XGBoost (Commented out)
# xgb_model = xgb.XGBClassifier(random_state=42)
# xgb_model.fit(X_train, y_train)
# xgb_pred_prob = xgb_model.predict_proba(X_test)[:, 1]
# model_results['XGBoost'] = xgb_pred_prob

# # 3. LightGBM (Commented out)
# lgb_model = lgb.LGBMClassifier(random_state=42)
# lgb_model.fit(X_train, y_train)
# lgb_pred_prob = lgb_model.predict_proba(X_test)[:, 1]
# model_results['LightGBM'] = lgb_pred_prob

# 2. Deep Learning Model with enhanced architecture
print("Training Deep Learning Model...")
class DeepModel(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.BatchNorm1d(256),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.2),
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.2),
            nn.Linear(128, 64),
            nn.BatchNorm1d(64),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.15),
            nn.Linear(64, 32),
            nn.BatchNorm1d(32),
            nn.LeakyReLU(0.1),
            nn.Dropout(0.1),
            nn.Linear(32, 1)
        )
        
    def forward(self, x):
        return torch.sigmoid(self.network(x))

# Training setup
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
deep_model = DeepModel(X_train.shape[1]).to(device)
criterion = nn.BCELoss()
optimizer = torch.optim.AdamW(deep_model.parameters(), lr=0.001, weight_decay=0.01)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=5)

# Convert data to PyTorch tensors
X_train_tensor = torch.FloatTensor(X_train).to(device)
y_train_tensor = torch.FloatTensor(y_train.values).to(device)
X_test_tensor = torch.FloatTensor(X_test).to(device)

# Training loop
n_epochs = 200
batch_size = 64
n_batches = len(X_train) // batch_size

for epoch in range(n_epochs):
    deep_model.train()
    total_loss = 0
    
    # Shuffle data
    indices = torch.randperm(len(X_train))
    
    for i in range(n_batches):
        start_idx = i * batch_size
        end_idx = start_idx + batch_size
        batch_X = X_train_tensor[indices[start_idx:end_idx]]
        batch_y = y_train_tensor[indices[start_idx:end_idx]]
        
        optimizer.zero_grad()
        outputs = deep_model(batch_X).squeeze()
        loss = criterion(outputs, batch_y)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(deep_model.parameters(), max_norm=1.0)
        optimizer.step()
        
        total_loss += loss.item()
    
    avg_loss = total_loss / n_batches
    scheduler.step(avg_loss)
    
    if (epoch + 1) % 10 == 0:
        print(f'Epoch [{epoch+1}/{n_epochs}], Loss: {avg_loss:.4f}')

# Evaluation
deep_model.eval()
with torch.no_grad():
    deep_pred_prob = deep_model(X_test_tensor).cpu().numpy().squeeze()
model_results['Neural Network'] = deep_pred_prob

# Plot ROC curves
plt.figure(figsize=(12, 10))  # Increased figure size
plt.rcParams.update({'font.size': 16})  # Set base font size

for model_name, pred_prob in model_results.items():
    fpr, tpr, _ = roc_curve(y_test, pred_prob)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, lw=3, label=f'{model_name} (AUC = {roc_auc:.2f})')

plt.plot([0, 1], [0, 1], color='navy', lw=3, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontsize=18, fontweight='bold')
plt.ylabel('True Positive Rate', fontsize=18, fontweight='bold')
plt.title('ROC Curves for Different Models', fontsize=20, fontweight='bold', pad=20)
plt.legend(loc="lower right", fontsize=16, frameon=True, framealpha=0.9)

# Add grid for better readability
plt.grid(True, linestyle='--', alpha=0.7)

# Adjust tick label sizes
plt.tick_params(axis='both', which='major', labelsize=16)

# Save the plot with high DPI for better quality
plt.savefig(os.path.join(OUTPUT_DIR, 'roc_curves.png'), dpi=300, bbox_inches='tight')
plt.close()

# Feature importance for Random Forest only
print("Plotting feature importance...")
feature_importance_df = pd.DataFrame({
    'Feature': X.columns,
    'Random Forest': rf_model.feature_importances_
})

# Plot feature importance
plt.figure(figsize=(12, 6))
feature_importance_df.set_index('Feature').plot(kind='bar')
plt.title('Feature Importance (Random Forest)')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'feature_importance.png'))
plt.close()

# Model Performance Summary
print("Calculating performance metrics...")
performance_summary = pd.DataFrame(columns=['Model', 'AUC-ROC', 'Accuracy'])
for model_name, pred_prob in model_results.items():
    fpr, tpr, _ = roc_curve(y_test, pred_prob)
    roc_auc = auc(fpr, tpr)
    predictions = (pred_prob > 0.5).astype(int)
    accuracy = (predictions == y_test).mean()
    performance_summary = pd.concat([performance_summary, pd.DataFrame({
        'Model': [model_name],
        'AUC-ROC': [roc_auc],
        'Accuracy': [accuracy]
    })], ignore_index=True)

# Save performance summary
performance_summary.to_csv(os.path.join(OUTPUT_DIR, 'model_performance_summary.csv'), index=False)

print(f"Analysis complete! Check the output directory: {OUTPUT_DIR}")
print("\nModel Performance Summary:")
print(performance_summary) 