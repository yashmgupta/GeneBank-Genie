# GeneBank Genie

**Version 1.0**  
© Dr Yash Munnalal Gupta

GeneBank Genie is a versatile desktop tool for the analysis of GenBank records. This application provides modules for:

- **General Analysis:** Calculate nucleotide composition, GC content, and gene count. Perform Principal Component Analysis (PCA) with outlier detection based on taxonomic groups.
- **Gene Analysis:** Filter GenBank records for a selected gene, compute gene-specific metrics, and perform PCA to visualize deviations.
- **Taxonomic Visualization (Sankey Diagram):** Visualize the flow of taxonomic annotations across hierarchical levels.
- **Additional Visualizations:** Generate correlation heatmaps, run K-Means clustering, and create pair plots to explore relationships between features.
- **Dendrogram Analysis:** Use hierarchical clustering to generate dendrograms that reveal natural groupings within the data.
- **Sequence Extraction:** Extract full or gene-specific sequences (in FASTA format) from GenBank files.

The graphical user interface is built using Tkinter, and the app leverages several key Python libraries for bioinformatics and data visualization.

---

## Features

- **User-Friendly Interface:** An easy-to-use GUI with multiple analysis modules organized in tabs.
- **Versatile Analysis Options:** Perform statistical analysis, PCA, clustering, and visualizations on GenBank record data.
- **Sequence Extraction:** Extract gene-specific or full sequences from GenBank files and save them in FASTA format.
- **Customizable Plots:** Built-in options to choose colors, markers, and taxonomic groups for all visualizations.

---

## Installation

### Prerequisites

- **Python 3.x** (tested with Python 3.8+)
- Tkinter (usually included with Python installations)
- Git (to clone the repository)

### Steps

1. **Clone the repository:**

    ```bash
    git clone https://github.com/yashmgupta/GeneBank-Genie
    cd GeneBank-Genie
    ```

2. **Create and activate a virtual environment (optional but recommended):**

    ```bash
    python -m venv venv
    # On Windows:
    venv\Scripts\activate
    # On macOS/Linux:
    source venv/bin/activate
    ```

3. **Install the required dependencies:**

    ```bash
    pip install -r requirements.txt
    ```

---

## Running the Application

To launch the GeneBank Genie app, run the following command from your project directory:

```bash
python mt9.2.py

##📋 Usage Instructions

Once the app is launched:

### 1. Load a GenBank File
- Click the **"Browse"** button to select a `.gb` or `.gbk` file.
- The selected file will be used across all analysis modules.

### 2. Explore Functional Tabs:
| Tab | Functionality |
|:---|:---|
| **General Analysis** | Perform feature extraction, PCA, and detect outliers based on taxonomy. |
| **Gene Analysis** | Select a gene, analyze its metrics, and identify deviations using PCA. |
| **Sankey Diagram** | Generate a taxonomy flow visualization across hierarchical levels. |
| **Additional Visualizations** | Create correlation heatmaps, KMeans clustering, and pairplots for feature exploration. |
| **Dendrogram Analysis** | Perform hierarchical clustering and generate dendrograms. |
| **Sequence Extraction** | Extract and download full sequences or gene-specific sequences in FASTA format. |
| **Summary** | View the detailed log of all actions and results during the session. |

---

# 📊 Key Visual Outputs

- **PCA Scatter Plots**: Visualize nucleotide composition, sequence length, GC content, and gene count distributions.
- **Outlier Detection**: Identify records that deviate from taxonomic groups based on PCA + Mahalanobis distance.
- **Sankey Diagrams**: Explore taxonomy transitions interactively.
- **Correlation Heatmaps**: Study relationships between multiple sequence features.
- **KMeans Clustering**: Group sequences based on feature similarity.
- **Hierarchical Dendrograms**: Understand phylogenetic or feature-based relationships visually.
- **Sequence FASTA Files**: Export full or gene-specific sequences easily.

---

# 🧹 Troubleshooting

| Problem | Solution |
|:-------|:---------|
| **File not found / file error** | Ensure you load a valid GenBank `.gb` or `.gbk` file first. |
| **No genes detected** | Your GenBank file may not have gene annotations. Try another dataset. |
| **Application freeze** | Large datasets can cause temporary lag during visualization (especially Sankey/Dendrogram). Please wait patiently. |
| **KMeans error** | Ensure enough data points exist for the number of clusters specified. |

---
