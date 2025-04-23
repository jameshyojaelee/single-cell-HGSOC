import sys
import os
sys.path.insert(0, '/gpfs/commons/home/jameslee/miniconda3/envs/scanpy/lib/python3.10/site-packages')

# Import numpy first
import numpy as np

# Import other packages
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, confusion_matrix, classification_report
from sklearn.preprocessing import StandardScaler
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

# Positional Encoding for Transformer
class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=5000):
        super().__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-np.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        self.register_buffer('pe', pe)

    def forward(self, x):
        return x + self.pe[:, :x.size(1)]

# Enhanced Transformer Model
class TabularTransformer(nn.Module):
    def __init__(self, input_dim, d_model=256, nhead=16, num_layers=4, dim_feedforward=512, dropout=0.1):
        super().__init__()
        
        # Input projection with residual connection
        self.input_proj = nn.Linear(input_dim, d_model)
        self.batch_norm1 = nn.BatchNorm1d(d_model)
        self.dropout1 = nn.Dropout(dropout)
        
        # Positional encoding
        self.pos_encoder = PositionalEncoding(d_model)
        
        # Transformer encoder with more layers and heads
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_feedforward,
            dropout=dropout,
            batch_first=True,
            activation='gelu'  # Using GELU activation
        )
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        # Enhanced output layers with residual connections
        self.fc1 = nn.Linear(d_model, 128)
        self.batch_norm2 = nn.BatchNorm1d(128)
        self.dropout2 = nn.Dropout(dropout)
        self.fc2 = nn.Linear(128, 64)
        self.batch_norm3 = nn.BatchNorm1d(64)
        self.dropout3 = nn.Dropout(dropout)
        self.fc3 = nn.Linear(64, 1)
        
        # Activation functions
        self.activation = nn.GELU()
        
    def forward(self, x):
        # Input projection with residual connection
        identity = x
        x = self.input_proj(x)
        x = self.batch_norm1(x)
        x = self.dropout1(x)
        
        # Add sequence dimension
        x = x.unsqueeze(1)  # [batch_size, 1, d_model]
        
        # Add positional encoding
        x = self.pos_encoder(x)
        
        # Transformer encoder
        x = self.transformer_encoder(x)
        
        # Global pooling
        x = x.squeeze(1)  # Remove sequence dimension
        
        # Enhanced classification layers with residual connections
        x = self.fc1(x)
        x = self.batch_norm2(x)
        x = self.activation(x)
        x = self.dropout2(x)
        
        x = self.fc2(x)
        x = self.batch_norm3(x)
        x = self.activation(x)
        x = self.dropout3(x)
        
        x = self.fc3(x)
        return torch.sigmoid(x)

# Generate synthetic data
def generate_synthetic_data(n_samples=1000):
    # Base features with more noise and complex distributions
    pct_TRMstem = np.random.beta(1.2, 6, n_samples) * 0.4  # Even more skewed
    pct_Exhausted = np.random.beta(1.5, 2.5, n_samples) * 0.3  # More uniform
    TIGIT_score = np.random.normal(0.5, 0.4, n_samples)     # More variance
    PDCD1_expr = np.random.gamma(1.2, 1.5, n_samples)       # More spread out
    
    # Add substantial noise
    noise1 = np.random.normal(0, 0.5, n_samples)    # Increased noise
    noise2 = np.random.normal(0, 0.45, n_samples)   # Increased noise
    noise3 = np.random.normal(0, 0.4, n_samples)    # Increased noise
    
    # Create very complex, highly non-linear relationships
    immune_score = (
        0.2 * pct_TRMstem +                         # Further reduced main effect
        -0.15 * pct_Exhausted + 
        0.1 * TIGIT_score + 
        -0.08 * PDCD1_expr +
        0.3 * np.sin(10 * pct_TRMstem) +           # More oscillations
        0.25 * np.cos(8 * pct_Exhausted) +
        0.2 * np.tan(4 * TIGIT_score) +            # More non-linear terms
        0.15 * np.exp(-3 * PDCD1_expr) +
        0.35 * noise1 +                             # Increased noise influence
        0.3 * noise2 +
        0.25 * noise3
    )
    
    # Add more random flipping of labels (15% of labels)
    immune_label = (immune_score > np.median(immune_score)).astype(int)
    flip_mask = np.random.random(n_samples) < 0.15  # Increased label noise
    immune_label[flip_mask] = 1 - immune_label[flip_mask]
    
    # Create more complex feature interactions with more noise
    interaction_1 = np.sin(pct_TRMstem * TIGIT_score) + noise1 * 0.4
    interaction_2 = np.exp(-pct_Exhausted * PDCD1_expr) + noise2 * 0.35
    ratio_feature = np.clip(
        (pct_TRMstem + noise3 * 0.3) / (pct_Exhausted + 0.01), 
        0, 5
    )
    
    # Add more noisy and irrelevant features
    random_feature1 = np.random.normal(0, 2.0, n_samples)
    random_feature2 = np.random.uniform(-3, 3, n_samples)
    random_feature3 = np.random.exponential(1.5, n_samples)
    random_feature4 = np.random.gamma(1.5, 2.5, n_samples)
    random_feature5 = np.random.lognormal(0, 1, n_samples)
    random_feature6 = np.random.weibull(1.5, n_samples)
    
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
        'Random_Feature5': random_feature5,
        'Random_Feature6': random_feature6,
        'immune_label': immune_label
    })
    
    return data

# Generate synthetic data
data = generate_synthetic_data(4000)  # Use same sample size as test_ML.py

# Features and labels
X = data.drop(columns=['sample_id', 'immune_label'])
y = data['immune_label']

# Scale the features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.4, random_state=42)

# Create datasets and dataloaders
train_dataset = TabularDataset(X_train, y_train)
test_dataset = TabularDataset(X_test, y_test)

train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=64, shuffle=False)

# Initialize model and training components
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = TabularTransformer(
    input_dim=X.shape[1],
    d_model=256,           # Increased model dimension
    nhead=16,              # More attention heads
    num_layers=4,          # More transformer layers
    dim_feedforward=512,   # Larger feed-forward network
    dropout=0.1
).to(device)

# Enhanced training setup
criterion = nn.BCELoss()
optimizer = torch.optim.AdamW(
    model.parameters(),
    lr=0.001,
    weight_decay=0.01,
    betas=(0.9, 0.999),
    eps=1e-8
)

# Training parameters
n_epochs = 300  # Define n_epochs before using it
warmup_steps = 100
total_steps = n_epochs * len(train_loader)

# Learning rate scheduler with warmup
scheduler = torch.optim.lr_scheduler.OneCycleLR(
    optimizer,
    max_lr=0.01,
    total_steps=total_steps,
    pct_start=0.1,  # 10% of steps for warmup
    anneal_strategy='cos',
    div_factor=25.0,
    final_div_factor=10000.0
)

# Training loop with gradient accumulation
gradient_accumulation_steps = 4
best_loss = float('inf')
patience_counter = 0
patience_limit = 30  # Increased patience

# Add learning rate warmup
def warmup_lr_scheduler(optimizer, warmup_steps, warmup_lr):
    def lr_lambda(step):
        if step < warmup_steps:
            return float(step) / float(max(1, warmup_steps))
        return 1.0
    return torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)

warmup_scheduler = warmup_lr_scheduler(optimizer, warmup_steps, 0.001)

for epoch in range(n_epochs):
    model.train()
    total_loss = 0
    optimizer.zero_grad()
    
    for i, (batch_X, batch_y) in enumerate(train_loader):
        batch_X, batch_y = batch_X.to(device), batch_y.to(device)
        
        outputs = model(batch_X).squeeze()
        loss = criterion(outputs, batch_y)
        loss = loss / gradient_accumulation_steps
        loss.backward()
        
        if (i + 1) % gradient_accumulation_steps == 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            
            # Update learning rate
            if epoch < warmup_steps:
                warmup_scheduler.step()
            else:
                scheduler.step()
                
            optimizer.zero_grad()
        
        total_loss += loss.item() * gradient_accumulation_steps
    
    avg_loss = total_loss / len(train_loader)
    
    if avg_loss < best_loss:
        best_loss = avg_loss
        patience_counter = 0
        # Save best model
        torch.save(model.state_dict(), os.path.join(OUTPUT_DIR, 'best_transformer_model.pth'))
    else:
        patience_counter += 1
    
    if patience_counter >= patience_limit:
        print(f"Early stopping at epoch {epoch+1}")
        break
    
    if (epoch + 1) % 10 == 0:
        current_lr = optimizer.param_groups[0]['lr']
        print(f'Epoch [{epoch+1}/{n_epochs}], Loss: {avg_loss:.4f}, LR: {current_lr:.6f}')

# Load best model for evaluation
model.load_state_dict(torch.load(os.path.join(OUTPUT_DIR, 'best_transformer_model.pth')))
model.eval()

# Evaluation
all_preds = []
all_labels = []

with torch.no_grad():
    for batch_X, batch_y in test_loader:
        batch_X, batch_y = batch_X.to(device), batch_y.to(device)
        outputs = model(batch_X).squeeze()
        all_preds.extend(outputs.cpu().numpy())
        all_labels.extend(batch_y.cpu().numpy())

# Calculate metrics
fpr, tpr, _ = roc_curve(all_labels, all_preds)
roc_auc = auc(fpr, tpr)
accuracy = ((np.array(all_preds) > 0.5) == np.array(all_labels)).mean()

print(f"Transformer Model Performance:")
print(f"AUC-ROC: {roc_auc:.4f}")
print(f"Accuracy: {accuracy:.4f}")

# Plot ROC curve
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Transformer Model ROC Curve')
plt.legend(loc="lower right")
plt.savefig(os.path.join(OUTPUT_DIR, 'transformer_roc.png'))
plt.close()

# Save performance metrics
performance_df = pd.DataFrame({
    'Model': ['Transformer'],
    'AUC-ROC': [roc_auc],
    'Accuracy': [accuracy]
})
performance_df.to_csv(os.path.join(OUTPUT_DIR, 'transformer_performance.csv'), index=False)

print("\nTransformer Model Performance:")
print(performance_df)
print(f"\nCheck the output directory: {OUTPUT_DIR} for detailed results.") 