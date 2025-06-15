
#  snRNA-seq Data Preprocessing for Human Optic Nerve and Optic Nerve Head Endothelial Cells

This repository contains a Python script `convert_h5ad_to_parquet.py` that processes a `.h5ad` single-nucleus RNA sequencing (snRNA-seq) dataset and converts it into optimized formats (`.parquet`) suitable for sharing and downstream analysis.

The dataset is hosted on Hugging Face:  
[View on Hugging Face](https://huggingface.co/datasets/longevity-db/snRNAseq_of_human_optic_nerve_and_optic_nerve_head_endothelial_cells)

---

## Dataset Summary

This dataset consists of single-nucleus RNA-sequencing (snRNA-seq) data of **human optic nerve** and **optic nerve head** endothelial cells. It is structured in the `.h5ad` format, commonly used in the AnnData ecosystem for storing large-scale omics data, particularly from `scanpy`.

---

## What the Script Does (`processing.py`)

The script converts the raw `.h5ad` file into more accessible and interoperable `.parquet` files for easier analysis and integration. It performs the following steps:

### 1. **Setup Configuration**
- Defines the input path to the `.h5ad` file.
- Sets the output directory name based on the dataset title.

### 2. **Creates Output Directory**
- A directory is created to hold all the output files, if it doesn't already exist.

### 3. **Loads the `.h5ad` File**
- Reads the input `.h5ad` file using the `anndata` library.
- Extracts the AnnData object which contains:
  - `adata.X`: Expression matrix (cells × genes)
  - `adata.var`: Gene/feature metadata
  - `adata.obs`: Cell-level metadata (not exported here)

### 4. **Saves Expression Matrix**
- Converts `adata.X` (expression data) to a Pandas DataFrame.
- Handles dense and sparse matrix formats.
- Saves the expression data as `expression.parquet`.

### 5. **Saves Feature Metadata**
- Converts `adata.var` (gene annotations) to a DataFrame.
- Saves as `feature_metadata.parquet`.

### 6. **Output**
- Two files are saved:
  - `expression.parquet`: Gene expression matrix
  - `feature_metadata.parquet`: Metadata about each gene

---

## Requirements

Install the required packages before running:

```bash
pip install pandas anndata pyarrow
```

---

## Usage

1. Replace the `H5AD_FILE_PATH` in the script with the path to your `.h5ad` file:

```python
H5AD_FILE_PATH = "/path/to/your/file.h5ad"
```

2. Run the script:

```bash
python convert_h5ad_to_parquet.py
```

3. Find the output `.parquet` files in the directory:

```
snRNA-seq_of_human_optic_nerve_and_optic_nerve_head_endothelial_cells/
│
├── expression.parquet
└── feature_metadata.parquet
```

---

## Notes

- This script **does not export** the `obs` (cell metadata) from the `.h5ad` file. You can modify it if needed.
- Output files are in `.parquet` format, which is efficient for big data workflows and compatible with tools like Apache Spark, Pandas, and cloud platforms.

---

## References

- Source Dataset:https: 1.//datasets.cellxgene.cziscience.com/f5b09167-e4b5-4f32-b7b4-f0b7c402a4c4.h5ad 2. https://cellxgene.cziscience.com/collections/05e3d0fc-c9dd-4f14-9163-2b242b3bb5c2
- Final Output Dataset: [longevity-db/snRNAseq_of_human_optic_nerve_and_optic_nerve_head_endothelial_cells](https://huggingface.co/datasets/longevity-db/snRNAseq_of_human_optic_nerve_and_optic_nerve_head_endothelial_cells)
- File Format: [Anndata Documentation](https://anndata.readthedocs.io/en/latest/)
- Hugging Face Datasets: [https://huggingface.co/datasets](https://huggingface.co/datasets)

  
Contributed by CellVPA Team
Venkatachalam, Pooja, Albert
