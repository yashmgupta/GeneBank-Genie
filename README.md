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
    git clone https://github.com/yourusername/your-repository-name.git
    cd your-repository-name
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
python mt9.1.py
