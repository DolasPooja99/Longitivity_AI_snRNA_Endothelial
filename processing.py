import pandas as pd
import anndata as ad
import os

# --- Configuration ---
# IMPORTANT: Replace with the actual path to your .h5ad file
H5AD_FILE_PATH = "/path/to/your/file.h5ad"

# Name of the directory to store the output Parquet files
OUTPUT_DIR_NAME = "snRNA-seq_of_human_optic_nerve_and_optic_nerve_head_endothelial_cells"

# --- 1. Create Output Directory ---
os.makedirs(OUTPUT_DIR_NAME, exist_ok=True)
print(f"Created output directory: {OUTPUT_DIR_NAME}")

# --- 2. Load the H5AD file into an AnnData object ---
try:
    adata = ad.read_h5ad(H5AD_FILE_PATH)
    print(f"Successfully loaded AnnData object from: {H5AD_FILE_PATH}")
    print(f"AnnData object shape: {adata.shape}")
except FileNotFoundError:
    print(f"Error: H5AD file not found at {H5AD_FILE_PATH}. Please check the path.")
    exit()
except Exception as e:
    print(f"An error occurred while loading the H5AD file: {e}")
    exit()

# --- 3. Save adata.X as expression.parquet ---
expression_parquet_path = os.path.join(OUTPUT_DIR_NAME, "expression.parquet")

# Convert adata.X (which can be a NumPy array or sparse matrix) to a pandas DataFrame
# Ensure it has observation names (index) and variable names (columns)
if isinstance(adata.X, (pd.DataFrame, pd.Series)):
    df_expression = adata.X
elif hasattr(adata.X, 'toarray'): # If it's a sparse matrix (e.g., scipy.sparse.csr_matrix)
    df_expression = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
else: # If it's a dense NumPy array
    df_expression = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

df_expression.to_parquet(expression_parquet_path, index=False)
print(f"Saved expression data to: {expression_parquet_path}")

# --- 4. Save adata.var as feature_metadata.parquet ---
feature_metadata_parquet_path = os.path.join(OUTPUT_DIR_NAME, "feature_metadata.parquet")
adata.var.to_parquet(feature_metadata_parquet_path, index=False)
print(f"Saved feature metadata to: {feature_metadata_parquet_path}")

print(f"\nAll required Parquet files have been created in the '{OUTPUT_DIR_NAME}' directory.")
print("You can now use these files for your submission.")