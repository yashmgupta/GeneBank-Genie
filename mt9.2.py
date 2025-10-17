import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import subprocess
import logging
from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import UndefinedSequenceError
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.stats import chi2
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import plotly.graph_objects as go
from matplotlib.lines import Line2D

# Setup basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# ---------------------------
# Core Classes and Functions
# ---------------------------
class GenBankParser:
    """
    Handles loading and parsing of a GenBank file.
    """
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.records = []

    def load_records(self):
        """
        Loads records from the GenBank file.
        Returns:
            A list of SeqRecord objects.
        """
        if not os.path.isfile(self.filepath):
            raise FileNotFoundError(f"File not found: {self.filepath}")
        self.records = list(SeqIO.parse(self.filepath, "genbank"))
        logging.info(f"Loaded {len(self.records)} records from {self.filepath}")
        return self.records


class AnalysisEngine:
    """
    Provides static analysis methods such as PCA, feature duplication and
    Mahalanobis distance computation.
    """
    @staticmethod
    def perform_pca(features, n_components=2):
        scaler = StandardScaler()
        features_scaled = scaler.fit_transform(features)
        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(features_scaled)
        return pca_result, pca.explained_variance_ratio_

    @staticmethod
    def duplicate_feature(X):
        return np.hstack([X, X])

    @staticmethod
    def compute_mahalanobis(points):
        mean_vec = np.mean(points, axis=0)
        cov_matrix = np.cov(points, rowvar=False)
        inv_cov_matrix = np.linalg.inv(cov_matrix)
        diff = points - mean_vec
        md_squared = np.sum(diff.dot(inv_cov_matrix) * diff, axis=1)
        return np.sqrt(md_squared)

# ---------------------------
# Sequence Extraction Helper Functions
# ---------------------------
def extract_gene_sequences(genbank_file, selected_gene):
    """
    Extracts gene sequences corresponding to the selected gene from the GenBank file.
    Returns a tuple: (list of SeqRecord objects, list of (record_id, reason) tuples for skipped records).
    """
    records = list(SeqIO.parse(genbank_file, "genbank"))
    gene_seqs = []
    skipped_records = []
    for rec in records:
        try:
            gene_found = False
            for feature in rec.features:
                if feature.type.lower() in ["gene", "cds"]:
                    if "gene" in feature.qualifiers:
                        gene_name = feature.qualifiers["gene"][0]
                        if gene_name.upper() == selected_gene.upper():
                            seq = feature.extract(rec.seq)
                            new_rec = SeqRecord(seq, id=rec.id, description=rec.annotations.get("organism", ""))
                            gene_seqs.append(new_rec)
                            gene_found = True
                            break  # extract one gene per record
            if not gene_found:
                skipped_records.append((rec.id, "gene not found in record"))
        except UndefinedSequenceError:
            logging.warning(f"Skipping record {rec.id}: undefined sequence in extract_gene_sequences")
            skipped_records.append((rec.id, "undefined sequence"))
    return gene_seqs, skipped_records

def extract_full_sequences(genbank_file):
    """
    Extracts the complete sequence of each record in the GenBank file and returns them as SeqRecord objects.
    Returns a tuple: (list of SeqRecord objects, list of (record_id, reason) tuples for skipped records).
    """
    records = list(SeqIO.parse(genbank_file, "genbank"))
    full_seqs = []
    skipped_records = []
    for rec in records:
        try:
            new_rec = SeqRecord(rec.seq, id=rec.id, description=rec.description)
            full_seqs.append(new_rec)
        except UndefinedSequenceError:
            logging.warning(f"Skipping record {rec.id}: undefined sequence in extract_full_sequences")
            skipped_records.append((rec.id, "undefined sequence"))
    return full_seqs, skipped_records

def save_fasta(sequences, output_file):
    """
    Saves a list of SeqRecord objects to a FASTA file.
    """
    SeqIO.write(sequences, output_file, "fasta")
    logging.info(f"Saved FASTA to {output_file}")

# ---------------------------
# Main Application Class
# ---------------------------
class GenBankAnalysisApp:
    """
    GeneBank Genie:
    A versatile tool for the analysis of GenBank records.
    Version 1.0 – © Dr. Yash Munnalal Gupta

    This application provides modules for general analysis, gene analysis,
    taxonomic visualization, additional visualizations, dendrogram analysis,
    and sequence extraction.
    """
    def __init__(self, root):
        self.root = root

        # Apply modern styling with ttk and a modern theme.
        self.style = ttk.Style()
        self.style.theme_use('clam')
        self.style.configure('TButton', padding=6, font=("Segoe UI", 10))
        self.style.configure('TLabel', padding=4, font=("Segoe UI", 10))
        self.style.configure('TEntry', padding=4, font=("Segoe UI", 10))
        self.style.configure('TNotebook.Tab', padding=[12, 8])

        # Set up main window properties.
        self.root.title("GeneBank Genie")
        self.root.geometry("1000x700")
        try:
            self.icon_image = tk.PhotoImage(file="symbol.png")
            self.root.iconphoto(False, self.icon_image)
        except Exception as e:
            logging.warning(f"Could not load icon image: {e}")

        # Global data variables
        self.global_gb_file = None  # File path string for the GenBank file
        self.df_general = None      # DataFrame for general analysis results
        self.df_gene = None         # DataFrame for gene analysis results
        self.summary_text = ""      # Summary log text for the Analysis Results tab

        self.setup_ui()
        self.create_menubar()

    # ---------------------------
    # Utility Methods
    # ---------------------------
    def update_summary(self, text):
        """Append new information to the summary text box."""
        self.summary_text += text + "\n"
        self.text_summary.config(state='normal')
        self.text_summary.delete(1.0, tk.END)
        self.text_summary.insert(tk.END, self.summary_text)
        self.text_summary.config(state='disabled')

    def browse_global_file(self):
        """Allow the user to browse for a GenBank file and update the global file."""
        file_path = filedialog.askopenfilename(
            title="Select GenBank File",
            filetypes=[("GenBank Files", "*.gb *.gbk"), ("All Files", "*.*")]
        )
        if file_path:
            self.global_gb_file = file_path
            self.entry_global_file.config(state='normal')
            self.entry_global_file.delete(0, tk.END)
            self.entry_global_file.insert(0, file_path)
            self.entry_global_file.config(state='readonly')
            self.update_summary(f"Global file loaded: {file_path}")

    def populate_tax_groups(self, source):
        """
        Populate the taxonomic groups for either the General Analysis or Gene Analysis view.
        """
        if not self.global_gb_file or not os.path.isfile(self.global_gb_file):
            messagebox.showerror("File Error", "Load a GenBank file first.")
            return
        try:
            level = int(self.entry_tax_level_general.get() if source == "general" 
                        else self.entry_tax_level_gene.get())
        except Exception:
            messagebox.showerror("Parameter Error", "Taxonomy level index must be an integer.")
            return
        records = list(SeqIO.parse(self.global_gb_file, "genbank"))
        groups = set()
        for rec in records:
            taxonomy = rec.annotations.get("taxonomy", [])
            if len(taxonomy) > level:
                groups.add(taxonomy[level])
            elif taxonomy:
                groups.add(taxonomy[-1])
        groups = sorted(list(groups))
        if source == "general":
            self.combo_sel_tax_general['values'] = groups
            if groups:
                self.combo_sel_tax_general.current(0)
        else:
            self.entry_sel_tax_group_gene.delete(0, tk.END)
            if groups:
                self.entry_sel_tax_group_gene.insert(0, groups[0])
        messagebox.showinfo("Tax Groups", f"Found {len(groups)} groups.")

    def load_gene_list(self):
        """
        Load and populate the gene list from the GenBank file.
        """
        if not self.global_gb_file or not os.path.isfile(self.global_gb_file):
            messagebox.showerror("File Error", "Load a GenBank file first.")
            return
        records = list(SeqIO.parse(self.global_gb_file, "genbank"))
        gene_set = set()
        for rec in records:
            for feature in rec.features:
                if feature.type.lower() in ["gene", "cds"]:
                    if "gene" in feature.qualifiers:
                        gene_name = feature.qualifiers["gene"][0]
                        gene_set.add(gene_name.upper())
        gene_list = sorted(list(gene_set))
        self.combo_gene['values'] = gene_list
        if gene_list:
            self.combo_gene.current(0)
        messagebox.showinfo("Gene List Loaded", f"Found {len(gene_list)} genes.")

    def refresh_common_genes(self):
        """
        Populate all gene names (the union of genes from all records) into the common gene dropdown.
        """
        if not self.global_gb_file or not os.path.isfile(self.global_gb_file):
            messagebox.showerror("File Error", "Load a GenBank file first.")
            return
        records = list(SeqIO.parse(self.global_gb_file, "genbank"))
        all_genes = set()
        for rec in records:
            for feature in rec.features:
                if feature.type.lower() in ["gene", "cds"] and "gene" in feature.qualifiers:
                    all_genes.add(feature.qualifiers["gene"][0].upper())
        if not all_genes:
            messagebox.showerror("Gene Extraction", "No genes found in the GenBank file.")
            self.combo_common_genes['values'] = []
            self.combo_common_genes.set("")
        else:
            all_genes = sorted(list(all_genes))
            self.combo_common_genes['values'] = all_genes
            self.combo_common_genes.current(0)
            messagebox.showinfo("Genes Populated", f"Found {len(all_genes)} genes.")

    def download_selected_gene(self):
        """
        Extract sequences for the gene selected from the dropdown and save as a FASTA file.
        """
        selected_gene = self.combo_common_genes.get()
        if not selected_gene:
            messagebox.showerror("Selection Error", "No gene selected.")
            return
        self.update_summary(f"Extracting sequences for gene: {selected_gene}")
        gene_seqs, skipped_records = extract_gene_sequences(self.global_gb_file, selected_gene)
        if not gene_seqs:
            messagebox.showerror("Gene Extraction Error", f"No sequences found for gene {selected_gene}")
            return
        out_dir = filedialog.askdirectory(title="Select Output Directory for Gene Extraction")
        if not out_dir:
            return
        output_file = os.path.join(out_dir, f"{selected_gene}_extracted.fasta")
        save_fasta(gene_seqs, output_file)
        self.update_summary(f"Selected gene sequences saved as FASTA: {output_file}")
        if skipped_records:
            self.update_summary(f"\n⚠ SKIPPED RECORDS ({len(skipped_records)}):")
            for rec_id, reason in skipped_records:
                self.update_summary(f"  - {rec_id}: {reason}")
        messagebox.showinfo("Download Complete", f"FASTA file saved: {output_file}")

    def download_full_sequences(self):
        """
        Extract full sequences from the GenBank file and save as a FASTA file.
        """
        self.update_summary("Extracting full sequences from GenBank records...")
        full_seqs, skipped_records = extract_full_sequences(self.global_gb_file)
        if not full_seqs:
            messagebox.showerror("Extraction Error", "No full sequences could be extracted.")
            return
        out_dir = filedialog.askdirectory(title="Select Output Directory for Full Sequences")
        if not out_dir:
            return
        output_file = os.path.join(out_dir, "full_sequences.fasta")
        save_fasta(full_seqs, output_file)
        self.update_summary(f"Full sequences saved as FASTA: {output_file}")
        if skipped_records:
            self.update_summary(f"\n⚠ SKIPPED RECORDS ({len(skipped_records)}):")
            for rec_id, reason in skipped_records:
                self.update_summary(f"  - {rec_id}: {reason}")
        messagebox.showinfo("Download Complete", f"FASTA file saved: {output_file}")

    def run_general_analysis(self):
        """
        Run general analysis on the GenBank file:
         - Parse records, compute features, perform PCA,
         - Detect outliers within the selected taxonomic group,
         - And create subplots for visual comparison.
        """
        if not self.global_gb_file or not os.path.isfile(self.global_gb_file):
            messagebox.showerror("File Error", "Load a GenBank file first.")
            return
        try:
            tax_level_index = int(self.entry_tax_level_general.get())
        except Exception:
            messagebox.showerror("Parameter Error", "Taxonomy level index must be an integer.")
            return

        selected_tax_group = self.combo_sel_tax_general.get()
        color_palette = self.combo_color_palette_general.get()
        default_marker = self.entry_default_marker_general.get()
        outlier_marker = self.entry_out_marker_general.get()

        records = list(SeqIO.parse(self.global_gb_file, "genbank"))
        logging.info(f"Running general analysis with {len(records)} records.")
        data = []
        all_tax_levels = defaultdict(list)
        skipped_records = []
        for rec in records:
            try:
                seq = str(rec.seq).upper()
            except UndefinedSequenceError:
                logging.warning(f"Skipping record {rec.id}: undefined sequence")
                skipped_records.append((rec.id, "undefined sequence"))
                continue
            total_length = len(seq)
            if total_length == 0:
                skipped_records.append((rec.id, "empty sequence"))
                continue
            countA = seq.count("A")
            countC = seq.count("C")
            countG = seq.count("G")
            countT = seq.count("T")
            propA = countA / total_length
            propC = countC / total_length
            propG = countG / total_length
            propT = countT / total_length
            gc_content = (countG + countC) / total_length
            gene_count = sum(1 for f in rec.features if f.type.lower() == 'gene')
            taxonomy = rec.annotations.get("taxonomy", [])
            for i, tax in enumerate(taxonomy):
                all_tax_levels[i].append(tax)
            if len(taxonomy) > tax_level_index:
                chosen_tax = taxonomy[tax_level_index]
            elif taxonomy:
                chosen_tax = taxonomy[-1]
            else:
                chosen_tax = "Unknown"
            organism = rec.annotations.get("organism", "Unknown")
            data.append([propA, propC, propG, propT, total_length, gc_content, gene_count,
                         organism, " | ".join(taxonomy), chosen_tax])
        df = pd.DataFrame(data, columns=['A','C','G','T','SeqLength','GC','GeneCount',
                                         'Organism','Full_Taxonomy','TaxLevel'])
        self.df_general = df

        # Prepare features and perform various PCAs
        features_nuc = df[['A','C','G','T']].values
        features_add = df[['SeqLength','GC','GeneCount']].values
        features_combined = df[['A','C','G','T','SeqLength','GC','GeneCount']].values
        features_seq = df[['SeqLength']].values
        features_gc = df[['GC']].values
        features_gene = df[['GeneCount']].values

        pca_nuc, var_nuc = AnalysisEngine.perform_pca(features_nuc)
        pca_add, var_add = AnalysisEngine.perform_pca(features_add)
        pca_comb, var_comb = AnalysisEngine.perform_pca(features_combined)
        pca_seq, var_seq = AnalysisEngine.perform_pca(AnalysisEngine.duplicate_feature(features_seq))
        pca_gc, var_gc = AnalysisEngine.perform_pca(AnalysisEngine.duplicate_feature(features_gc))
        pca_gene, var_gene = AnalysisEngine.perform_pca(AnalysisEngine.duplicate_feature(features_gene))

        df['PC1_nuc'] = pca_nuc[:, 0]; df['PC2_nuc'] = pca_nuc[:, 1]
        df['PC1_add'] = pca_add[:, 0]; df['PC2_add'] = pca_add[:, 1]
        df['PC1_comb'] = pca_comb[:, 0]; df['PC2_comb'] = pca_comb[:, 1]
        df['PC1_seq'] = pca_seq[:, 0]; df['PC2_seq'] = pca_seq[:, 1]
        df['PC1_gc'] = pca_gc[:, 0]; df['PC2_gc'] = pca_gc[:, 1]
        df['PC1_gene'] = pca_gene[:, 0]; df['PC2_gene'] = pca_gene[:, 1]

        # Outlier detection within the selected taxonomic group.
        group_mask = df['TaxLevel'] == selected_tax_group
        df.loc[:, 'Outlier'] = False
        group_df = df[group_mask]
        if not group_df.empty:
            group_points = group_df[['PC1_comb', 'PC2_comb']].values
            distances = AnalysisEngine.compute_mahalanobis(group_points)
            threshold = np.sqrt(chi2.ppf(0.95, df=2))
            outlier_flags = distances > threshold
            df.loc[group_mask, 'Outlier'] = outlier_flags
            self.update_summary(f"General Analysis: Outlier threshold for {selected_tax_group}: {threshold:.2f}")
            self.update_summary("General Analysis: Outlier distances: " +
                                ", ".join([f"{d:.2f}" for d in distances[outlier_flags]]))
        else:
            self.update_summary(f"General Analysis: No records found for {selected_tax_group}.")

        outlier_species = df[(df['TaxLevel'] == selected_tax_group) & (df['Outlier'] == True)]['Organism']
        if not outlier_species.empty:
            unique_species = outlier_species.unique()
            self.update_summary("General Analysis: Outlier species in " + selected_tax_group + ": " +
                                ", ".join(unique_species))
        else:
            self.update_summary("General Analysis: No outlier species detected in " + selected_tax_group + ".")

        unique_tax = df['TaxLevel'].unique()
        cmap = plt.cm.get_cmap(color_palette, len(unique_tax))
        color_dict = {tax: cmap(i) for i, tax in enumerate(unique_tax)}

        # Create subplots for comparisons.
        fig, axes = plt.subplots(2, 3, figsize=(24, 12))
        axes = axes.flatten()
        for tax in unique_tax:
            subset = df[df['TaxLevel'] == tax]
            idx = subset.index
            axes[0].scatter(pca_nuc[idx, 0], pca_nuc[idx, 1],
                            label=tax, color=color_dict[tax], edgecolors='k', alpha=0.7)
        axes[0].set_title(f"Nucleotide PCA (var: {var_nuc[0]:.2f}, {var_nuc[1]:.2f})")
        axes[0].set_xlabel("PC1"); axes[0].set_ylabel("PC2"); axes[0].grid(True)
        for tax in unique_tax:
            subset = df[df['TaxLevel'] == tax]
            idx = subset.index
            axes[1].scatter(pca_seq[idx, 0], pca_seq[idx, 1],
                            label=tax, color=color_dict[tax], edgecolors='k', alpha=0.7)
        axes[1].set_title(f"SeqLength PCA (var: {var_seq[0]:.2f}, {var_seq[1]:.2f})")
        axes[1].set_xlabel("PC1"); axes[1].set_ylabel("PC2"); axes[1].grid(True)
        for tax in unique_tax:
            subset = df[df['TaxLevel'] == tax]
            idx = subset.index
            axes[2].scatter(pca_gc[idx, 0], pca_gc[idx, 1],
                            label=tax, color=color_dict[tax], edgecolors='k', alpha=0.7)
        axes[2].set_title(f"GC PCA (var: {var_gc[0]:.2f}, {var_gc[1]:.2f})")
        axes[2].set_xlabel("PC1"); axes[2].set_ylabel("PC2"); axes[2].grid(True)
        for tax in unique_tax:
            subset = df[df['TaxLevel'] == tax]
            idx = subset.index
            axes[3].scatter(pca_gene[idx, 0], pca_gene[idx, 1],
                            label=tax, color=color_dict[tax], edgecolors='k', alpha=0.7)
        axes[3].set_title(f"Gene Count PCA (var: {var_gene[0]:.2f}, {var_gene[1]:.2f})")
        axes[3].set_xlabel("PC1"); axes[3].set_ylabel("PC2"); axes[3].grid(True)
        for tax in unique_tax:
            subset = df[df['TaxLevel'] == tax]
            axes[4].scatter(pca_add[subset.index, 0], pca_add[subset.index, 1],
                            label=tax, color=color_dict[tax], edgecolors='k', alpha=0.7)
        axes[4].set_title(f"Additional PCA (var: {var_add[0]:.2f}, {var_add[1]:.2f})")
        axes[4].set_xlabel("PC1"); axes[4].set_ylabel("PC2"); axes[4].grid(True)
        for tax in unique_tax:
            subset = df[df['TaxLevel'] == tax]
            markers = [default_marker] * len(subset)
            if tax == selected_tax_group:
                for i, rec_idx in enumerate(subset.index):
                    if df.loc[rec_idx, 'Outlier']:
                        markers[i] = outlier_marker
            for i, rec_idx in enumerate(subset.index):
                axes[5].scatter(df.loc[rec_idx, 'PC1_comb'], df.loc[rec_idx, 'PC2_comb'],
                                marker=markers[i],
                                color=color_dict[tax],
                                edgecolors='k', alpha=0.7)
        axes[5].set_title(f"Combined PCA (var: {var_comb[0]:.2f}, {var_comb[1]:.2f})")
        axes[5].set_xlabel("PC1"); axes[5].set_ylabel("PC2"); axes[5].grid(True)

        legend_elements = [Line2D([0], [0], marker=default_marker, color='w', label=tax,
                                    markerfacecolor=color_dict[tax], markersize=8, markeredgecolor='k')
                           for tax in unique_tax]
        fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.92, 0.92))
        fig.suptitle("General Analysis: PCA Comparisons", fontsize=18)
        plt.tight_layout(rect=[0, 0, 0.9, 0.93])
        plt.show()

        df.to_csv("pca_features.csv", index=False)
        self.update_summary("General analysis complete. Results saved to pca_features.csv")
        if skipped_records:
            self.update_summary(f"\n⚠ SKIPPED RECORDS ({len(skipped_records)}):")
            for rec_id, reason in skipped_records:
                self.update_summary(f"  - {rec_id}: {reason}")
        else:
            self.update_summary("✓ All records were processed successfully.")

    def run_gene_analysis(self):
        """
        Run gene analysis:
         - Filter records by the selected gene.
         - Compute features and run PCA.
         - Detect outliers for the chosen taxonomic group.
         - Generate an interactive scatter plot with a custom legend.
        """
        if not self.global_gb_file or not os.path.isfile(self.global_gb_file):
            messagebox.showerror("File Error", "Load a GenBank file first.")
            return
        selected_gene = self.combo_gene.get()
        try:
            tax_level_index = int(self.entry_tax_level_gene.get())
        except Exception:
            messagebox.showerror("Parameter Error", "Taxonomy level index must be an integer.")
            return
        selected_tax_group = self.entry_sel_tax_group_gene.get()
        color_palette = self.entry_color_palette_gene.get()
        default_marker = self.entry_default_marker_gene.get()
        outlier_marker = self.entry_out_marker_gene.get()

        records = list(SeqIO.parse(self.global_gb_file, "genbank"))
        gene_data = []
        all_tax_levels = defaultdict(list)
        skipped_records = []
        for rec in records:
            gene_seq = None
            try:
                for feature in rec.features:
                    if feature.type.lower() in ["gene", "cds"]:
                        if "gene" in feature.qualifiers:
                            gene_name = feature.qualifiers["gene"][0]
                            if gene_name.upper() == selected_gene.upper():
                                gene_seq = str(feature.extract(rec.seq)).upper()
                                break
            except UndefinedSequenceError:
                logging.warning(f"Skipping record {rec.id}: undefined sequence during gene extraction")
                skipped_records.append((rec.id, "undefined sequence"))
                continue
            if gene_seq is None or len(gene_seq) == 0:
                skipped_records.append((rec.id, "gene not found or empty"))
                continue
            gene_length = len(gene_seq)
            countA = gene_seq.count("A")
            countC = gene_seq.count("C")
            countG = gene_seq.count("G")
            countT = gene_seq.count("T")
            propA = countA / gene_length
            propC = countC / gene_length
            propG = countG / gene_length
            propT = countT / gene_length
            gc_content = (countG + countC) / gene_length
            taxonomy = rec.annotations.get("taxonomy", [])
            for i, tax in enumerate(taxonomy):
                all_tax_levels[i].append(tax)
            if len(taxonomy) > tax_level_index:
                chosen_tax = taxonomy[tax_level_index]
            elif taxonomy:
                chosen_tax = taxonomy[-1]
            else:
                chosen_tax = "Unknown"
            organism = rec.annotations.get("organism", "Unknown")
            gene_data.append([propA, propC, propG, propT, gene_length, gc_content,
                              organism, " | ".join(taxonomy), chosen_tax])
        df = pd.DataFrame(gene_data, columns=["A", "C", "G", "T", "GeneLength", "GC",
                                               "Organism", "Full_Taxonomy", "TaxLevel"])
        self.df_gene = df
        features = df[["A", "C", "G", "T", "GeneLength", "GC"]].values
        pca_result, var_gene = AnalysisEngine.perform_pca(features)
        df["PC1"] = pca_result[:, 0]
        df["PC2"] = pca_result[:, 1]
        group_mask = df["TaxLevel"] == selected_tax_group
        df.loc[:, "Outlier"] = False
        group_df = df[group_mask]
        if not group_df.empty:
            group_points = group_df[["PC1", "PC2"]].values
            distances = AnalysisEngine.compute_mahalanobis(group_points)
            threshold = np.sqrt(chi2.ppf(0.95, df=2))
            outlier_flags = distances > threshold
            df.loc[group_mask, "Outlier"] = outlier_flags
            self.update_summary(f"Gene Analysis ({selected_gene}): Outlier threshold for {selected_tax_group} = {threshold:.2f}")
            self.update_summary("Gene Analysis: Outlier distances: " +
                                ", ".join([f"{d:.2f}" for d in distances[outlier_flags]]))
        else:
            self.update_summary(f"Gene Analysis: No records found for {selected_tax_group}.")
        outlier_species = df[(df["TaxLevel"] == selected_tax_group) & (df["Outlier"] == True)]["Organism"]
        if not outlier_species.empty:
            unique_species = outlier_species.unique()
            self.update_summary(f"Gene Analysis: Outlier species in {selected_tax_group}: " +
                                ", ".join(unique_species))
        else:
            self.update_summary(f"Gene Analysis: No outlier species detected in {selected_tax_group}.")

        unique_tax = df["TaxLevel"].unique()
        cmap = plt.cm.get_cmap(color_palette, len(unique_tax))
        color_dict = {tax: cmap(i) for i, tax in enumerate(unique_tax)}
        plt.figure(figsize=(8, 6))
        for tax in unique_tax:
            subset = df[df["TaxLevel"] == tax]
            marker_list = [default_marker] * len(subset)
            if tax == selected_tax_group:
                for i, idx in enumerate(subset.index):
                    if df.loc[idx, "Outlier"]:
                        marker_list[i] = outlier_marker
            for i, idx in enumerate(subset.index):
                plt.scatter(df.loc[idx, "PC1"], df.loc[idx, "PC2"],
                            marker=marker_list[i],
                            color=color_dict[tax],
                            edgecolors="k", alpha=0.7)
        plt.title(f"PCA of Selected Gene {selected_gene}\nColored by Taxonomic Level")
        plt.xlabel("PC1"); plt.ylabel("PC2")
        legend_elements = [Line2D([0], [0], marker=default_marker, color='w', label=tax,
                                    markerfacecolor=color_dict[tax], markersize=8, markeredgecolor='k')
                           for tax in unique_tax]
        plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.tight_layout()
        plt.show()
        df.to_csv("selected_gene_pca_features.csv", index=False)
        self.update_summary(f"Gene analysis complete for {selected_gene}. Results saved to selected_gene_pca_features.csv")
        if skipped_records:
            self.update_summary(f"\n⚠ SKIPPED RECORDS ({len(skipped_records)}):")
            for rec_id, reason in skipped_records:
                self.update_summary(f"  - {rec_id}: {reason}")
        else:
            self.update_summary("✓ All records were processed successfully.")

    def run_sankey_diagram(self):
        """
        Generate and display a Sankey diagram of taxonomy flow starting at the given level.
        """
        if not self.global_gb_file or not os.path.isfile(self.global_gb_file):
            messagebox.showerror("File Error", "Load a GenBank file first.")
            return
        try:
            start_level = int(self.entry_start_level.get())
        except Exception:
            messagebox.showerror("Parameter Error", "Start level must be an integer.")
            return
        records = list(SeqIO.parse(self.global_gb_file, "genbank"))
        import collections
        link_counts = collections.defaultdict(int)
        nodes_set = set()
        for rec in records:
            taxonomy = rec.annotations.get("taxonomy", [])
            if len(taxonomy) <= start_level:
                continue
            for i in range(start_level, len(taxonomy) - 1):
                source_node = f"L{i}-{taxonomy[i]}"
                target_node = f"L{i+1}-{taxonomy[i+1]}"
                nodes_set.add(source_node)
                nodes_set.add(target_node)
                link_counts[(source_node, target_node)] += 1
        def sort_key(label):
            try:
                level = int(label.split("-")[0][1:])
            except:
                level = 9999
            return (level, label)
        nodes_list = sorted(list(nodes_set), key=sort_key)
        node_to_index = {node: i for i, node in enumerate(nodes_list)}
        sources = []
        targets = []
        values = []
        for (src, tgt), count in link_counts.items():
            sources.append(node_to_index[src])
            targets.append(node_to_index[tgt])
            values.append(count)
        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=nodes_list,
                color="blue"
            ),
            link=dict(
                source=sources,
                target=targets,
                value=values
            ))])
        fig.update_layout(title_text=f"Sankey Diagram of Taxonomy Flow (Starting at Level {start_level})", font_size=10)
        fig.show()

    def show_heatmap(self):
        """
        Display a correlation heatmap for selected features from the general analysis.
        """
        if self.df_general is None:
            messagebox.showerror("Data Error", "Run General Analysis first.")
            return
        selections = list(self.listbox_heat.curselection())
        if not selections:
            messagebox.showerror("Selection Error", "Select at least one feature for the heatmap.")
            return
        feature_cols = [self.listbox_heat.get(i) for i in selections]
        plt.figure(figsize=(10, 8))
        corr_matrix = self.df_general[feature_cols].corr()
        sns.heatmap(corr_matrix, annot=True, cmap='coolwarm')
        plt.title("Correlation Heatmap")
        plt.show()

    def run_kmeans(self):
        """
        Run K-Means clustering globally over the general analysis data and visualize the result.
        """
        if self.df_general is None:
            messagebox.showerror("Data Error", "Run General Analysis first.")
            return
        try:
            n_clusters = int(self.entry_n_clusters.get())
        except:
            messagebox.showerror("Parameter Error", "Number of clusters must be an integer.")
            return
        features_combined = self.df_general[['A','C','G','T','SeqLength','GC','GeneCount']].values
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        clusters = kmeans.fit_predict(features_combined)
        self.df_general['KMeans_Cluster'] = clusters
        sil_score = silhouette_score(features_combined, clusters)
        self.update_summary(f"K-Means Clustering: Silhouette score = {sil_score:.2f}")
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x='PC1_comb', y='PC2_comb', hue='KMeans_Cluster', data=self.df_general,
                        palette='tab10', s=80)
        plt.title("K-Means Clustering (Combined PCA Projection)")
        plt.xlabel("PC1_comb"); plt.ylabel("PC2_comb")
        plt.legend(bbox_to_anchor=(1.05, 1))
        plt.show()

    def show_pairplot(self):
        """
        Display a pair plot for selected features (A, C, G, T).
        """
        if self.df_general is None:
            messagebox.showerror("Data Error", "Run General Analysis first.")
            return
        feature_cols = ['A', 'C', 'G', 'T']
        sns.pairplot(self.df_general[feature_cols], diag_kind='kde')
        plt.suptitle("Pair Plot of Selected Features", y=1.02)
        plt.show()

    def run_dendrogram(self):
        """
        Generate and display a dendrogram for the selected set of features.
        """
        if self.df_general is None:
            messagebox.showerror("Data Error", "Run General Analysis first.")
            return
        feature_option = self.combo_dendro.get()
        if feature_option == "Nucleotide Composition":
            features = self.df_general[['A', 'C', 'G', 'T']].values
        elif feature_option == "Additional Features":
            features = self.df_general[['SeqLength', 'GC', 'GeneCount']].values
        elif feature_option == "Combined Features":
            features = self.df_general[['A', 'C', 'G', 'T', 'SeqLength', 'GC', 'GeneCount']].values
        else:
            messagebox.showerror("Selection Error", "Select a valid feature option for dendrogram.")
            return
        scaler = StandardScaler()
        features_scaled = scaler.fit_transform(features)
        linkage_matrix = sch.linkage(features_scaled, method='ward')
        plt.figure(figsize=(12, 8))
        sch.dendrogram(linkage_matrix, labels=self.df_general['TaxLevel'].values, leaf_rotation=90)
        plt.title("Dendrogram Analysis")
        plt.xlabel("Taxonomic Level")
        plt.ylabel("Euclidean Distance")
        plt.show()

    # ---------------------------
    # Help Menu Functions
    # ---------------------------
    def show_help_general(self):
        info = ("General Analysis:\n"
                "This view parses GenBank records and calculates features such as nucleotide composition, GC content, gene count, "
                "and performs PCA along with outlier detection based on a selected taxonomic group. The resulting plots help you "
                "visualize data patterns and detect anomalies.")
        messagebox.showinfo("General Analysis Info", info)

    def show_help_gene(self):
        info = ("Gene Analysis:\n"
                "This view filters records for a selected gene, computes gene-specific metrics, performs PCA, and uses outlier detection "
                "to identify records that deviate from the norm for that gene.")
        messagebox.showinfo("Gene Analysis Info", info)

    def show_help_sankey(self):
        info = ("Sankey Diagram:\n"
                "This view generates a Sankey diagram to visualize the flow of taxonomic annotations across different hierarchical levels. "
                "It helps in understanding the structure and distribution of taxonomic groups in your dataset.")
        messagebox.showinfo("Sankey Diagram Info", info)

    def show_help_visual(self):
        info = ("Additional Visualizations:\n"
                "This view provides tools such as correlation heatmaps, K-Means clustering, and pair plots to explore relationships among "
                "various features extracted from the GenBank records.")
        messagebox.showinfo("Additional Visualizations Info", info)

    def show_help_dendrogram(self):
        info = ("Dendrogram Analysis:\n"
                "This view performs hierarchical clustering on selected feature sets and displays a dendrogram, highlighting natural groupings "
                "within your data.")
        messagebox.showinfo("Dendrogram Analysis Info", info)

    def show_help_extraction(self):
        info = ("Sequence Extraction:\n"
                "This view lets you extract and download sequences from the GenBank file. You can download full sequences for every record "
                "or populate a list of all detected genes and download sequences for a gene of your choice in FASTA format.")
        messagebox.showinfo("Sequence Extraction Info", info)

    def show_help_about(self):
        info = ("GeneBank Genie:\n"
                "Version 1.0 – © Dr Yash Munnalal Gupta\n"
                "Faculty of Science, Department of Biology\n"
                "Naresuan University, Thailand\n\n"
                "GeneBank Genie is a comprehensive tool for analyzing GenBank files. It includes modules for general analysis, gene analysis, "
                "taxonomic visualization, additional visualizations, dendrogram analysis, and sequence extraction. Use the Help menu for detailed info.")
        messagebox.showinfo("About GeneBank Genie", info)

    def create_menubar(self):
        menubar = tk.Menu(self.root)
        helpmenu = tk.Menu(menubar, tearoff=0)
        helpmenu.add_command(label="General Analysis Info", command=self.show_help_general)
        helpmenu.add_command(label="Gene Analysis Info", command=self.show_help_gene)
        helpmenu.add_command(label="Sankey Diagram Info", command=self.show_help_sankey)
        helpmenu.add_command(label="Additional Visualizations Info", command=self.show_help_visual)
        helpmenu.add_command(label="Dendrogram Analysis Info", command=self.show_help_dendrogram)
        helpmenu.add_command(label="Sequence Extraction Info", command=self.show_help_extraction)
        helpmenu.add_separator()
        helpmenu.add_command(label="About GeneBank Genie", command=self.show_help_about)
        menubar.add_cascade(label="Help", menu=helpmenu)
        self.root.config(menu=menubar)

    # ---------------------------
    # Setup UI
    # ---------------------------
    def setup_ui(self):
        """
        Build the main UI using Tkinter and its themed widgets.
        """
        # Global File Selection Frame
        frame_file = ttk.Frame(self.root, padding=(10, 10))
        frame_file.pack(fill='x', padx=10, pady=10)
        ttk.Label(frame_file, text="GenBank File:").pack(side='left')
        self.entry_global_file = ttk.Entry(frame_file, width=70, state='readonly')
        self.entry_global_file.pack(side='left', padx=5)
        ttk.Button(frame_file, text="Browse...", command=self.browse_global_file).pack(side='left', padx=5)

        # Notebook for Tabs
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill='both', expand=True, padx=10, pady=10)

        # ----- Tab 1: General Analysis -----
        frame_general = ttk.Frame(notebook, padding=10)
        notebook.add(frame_general, text="General Analysis")
        ttk.Label(frame_general, text="Taxonomy Level Index:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.entry_tax_level_general = ttk.Entry(frame_general, width=10)
        self.entry_tax_level_general.grid(row=0, column=1, sticky="w", padx=5, pady=5)
        self.entry_tax_level_general.insert(0, "11")
        ttk.Label(frame_general, text="Selected Taxonomic Group:").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.combo_sel_tax_general = ttk.Combobox(frame_general, width=20)
        self.combo_sel_tax_general.grid(row=1, column=1, sticky="w", padx=5, pady=5)
        self.combo_sel_tax_general.set("Acrididea")
        ttk.Label(frame_general, text="Color Palette:").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.combo_color_palette_general = ttk.Combobox(frame_general, width=15,
                                                        values=["tab20", "viridis", "plasma", "inferno", "magma", "cividis"])
        self.combo_color_palette_general.grid(row=2, column=1, sticky="w", padx=5, pady=5)
        self.combo_color_palette_general.current(0)
        ttk.Label(frame_general, text="Default Marker:").grid(row=3, column=0, sticky="w", padx=5, pady=5)
        self.entry_default_marker_general = ttk.Entry(frame_general, width=10)
        self.entry_default_marker_general.grid(row=3, column=1, sticky="w", padx=5, pady=5)
        self.entry_default_marker_general.insert(0, "o")
        ttk.Label(frame_general, text="Outlier Marker:").grid(row=4, column=0, sticky="w", padx=5, pady=5)
        self.entry_out_marker_general = ttk.Entry(frame_general, width=10)
        self.entry_out_marker_general.grid(row=4, column=1, sticky="w", padx=5, pady=5)
        self.entry_out_marker_general.insert(0, "D")
        ttk.Button(frame_general, text="Populate Tax Groups", command=lambda: self.populate_tax_groups("general")).grid(row=0, column=2, padx=5, pady=5)
        ttk.Button(frame_general, text="Run General Analysis", command=self.run_general_analysis).grid(row=5, column=1, padx=5, pady=15)
        info_gen = ("General Analysis: Parses GenBank records to calculate nucleotide composition, GC content, gene count, "
                    "and performs PCA with outlier detection based on taxonomic groups.")
        lbl_info_gen = ttk.Label(frame_general, text=info_gen, wraplength=800, foreground="blue")
        lbl_info_gen.grid(row=6, column=0, columnspan=3, padx=5, pady=10)

        # ----- Tab 2: Gene Analysis -----
        frame_gene = ttk.Frame(notebook, padding=10)
        notebook.add(frame_gene, text="Gene Analysis")
        ttk.Label(frame_gene, text="Select Gene:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.combo_gene = ttk.Combobox(frame_gene, width=20)
        self.combo_gene.grid(row=0, column=1, padx=5, pady=5)
        self.combo_gene.set("COX1")
        ttk.Button(frame_gene, text="Load Gene List", command=self.load_gene_list).grid(row=1, column=0, padx=5, pady=5)
        ttk.Button(frame_gene, text="Populate Tax Groups", command=lambda: self.populate_tax_groups("gene")).grid(row=1, column=1, padx=5, pady=5)
        ttk.Label(frame_gene, text="Taxonomy Level Index:").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.entry_tax_level_gene = ttk.Entry(frame_gene, width=10)
        self.entry_tax_level_gene.grid(row=2, column=1, sticky="w", padx=5, pady=5)
        self.entry_tax_level_gene.insert(0, "11")
        ttk.Label(frame_gene, text="Selected Taxonomic Group:").grid(row=3, column=0, sticky="w", padx=5, pady=5)
        self.entry_sel_tax_group_gene = ttk.Entry(frame_gene, width=20)
        self.entry_sel_tax_group_gene.grid(row=3, column=1, padx=5, pady=5, sticky="w")
        self.entry_sel_tax_group_gene.insert(0, "Acrididea")
        ttk.Label(frame_gene, text="Color Palette:").grid(row=4, column=0, sticky="w", padx=5, pady=5)
        self.entry_color_palette_gene = ttk.Entry(frame_gene, width=15)
        self.entry_color_palette_gene.grid(row=4, column=1, sticky="w", padx=5, pady=5)
        self.entry_color_palette_gene.insert(0, "tab20")
        ttk.Label(frame_gene, text="Default Marker:").grid(row=5, column=0, sticky="w", padx=5, pady=5)
        self.entry_default_marker_gene = ttk.Entry(frame_gene, width=10)
        self.entry_default_marker_gene.grid(row=5, column=1, padx=5, pady=5, sticky="w")
        self.entry_default_marker_gene.insert(0, "o")
        ttk.Label(frame_gene, text="Outlier Marker:").grid(row=6, column=0, sticky="w", padx=5, pady=5)
        self.entry_out_marker_gene = ttk.Entry(frame_gene, width=10)
        self.entry_out_marker_gene.grid(row=6, column=1, padx=5, pady=5, sticky="w")
        self.entry_out_marker_gene.insert(0, "D")
        ttk.Button(frame_gene, text="Run Gene Analysis", command=self.run_gene_analysis).grid(row=7, column=1, padx=5, pady=15)
        info_gene = ("Gene Analysis: Filters records for a selected gene, computes gene-specific metrics, performs PCA, "
                     "and visualizes outlier detection.")
        lbl_info_gene = ttk.Label(frame_gene, text=info_gene, wraplength=800, foreground="blue")
        lbl_info_gene.grid(row=8, column=0, columnspan=3, padx=5, pady=10)

        # ----- Tab 3: Sankey Diagram -----
        frame_sankey = ttk.Frame(notebook, padding=10)
        notebook.add(frame_sankey, text="Sankey Diagram")
        ttk.Label(frame_sankey, text="Start Level (for taxonomy):").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.entry_start_level = ttk.Entry(frame_sankey, width=10)
        self.entry_start_level.grid(row=0, column=1, sticky="w", padx=5, pady=5)
        self.entry_start_level.insert(0, "12")
        ttk.Button(frame_sankey, text="Run Sankey Diagram", command=self.run_sankey_diagram).grid(row=1, column=1, padx=5, pady=15)
        info_sankey = ("Sankey Diagram: Visualizes the flow of taxonomic annotations across hierarchical levels, "
                       "helping you understand group relationships.")
        lbl_info_sankey = ttk.Label(frame_sankey, text=info_sankey, wraplength=800, foreground="blue")
        lbl_info_sankey.grid(row=2, column=0, columnspan=3, padx=5, pady=10)

        # ----- Tab 4: Additional Visualizations -----
        frame_visual = ttk.Frame(notebook, padding=10)
        notebook.add(frame_visual, text="Additional Visualizations")
        ttk.Label(frame_visual, text="Select Features for Heatmap (CTRL+Click):").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.listbox_heat = tk.Listbox(frame_visual, selectmode=tk.MULTIPLE, height=6)
        self.listbox_heat.grid(row=1, column=0, padx=5, pady=5, sticky="w")
        for feat in ['A','C','G','T','SeqLength','GC','GeneCount']:
            self.listbox_heat.insert(tk.END, feat)
        ttk.Button(frame_visual, text="Show Heatmap", command=self.show_heatmap).grid(row=2, column=0, padx=5, pady=5)
        ttk.Label(frame_visual, text="Number of Clusters:").grid(row=0, column=1, sticky="w", padx=5, pady=5)
        self.entry_n_clusters = ttk.Entry(frame_visual, width=10)
        self.entry_n_clusters.grid(row=0, column=2, padx=5, pady=5)
        self.entry_n_clusters.insert(0, "2")
        ttk.Button(frame_visual, text="Run K-Means Clustering", command=self.run_kmeans).grid(row=2, column=2, padx=5, pady=5)
        ttk.Button(frame_visual, text="Show Pair Plot", command=self.show_pairplot).grid(row=3, column=0, padx=5, pady=15)
        info_visual = ("Additional Visualizations: Includes correlation heatmaps, K-Means clustering, and pair plots to explore inter-feature relationships.")
        lbl_info_visual = ttk.Label(frame_visual, text=info_visual, wraplength=800, foreground="blue")
        lbl_info_visual.grid(row=4, column=0, columnspan=3, padx=5, pady=10)

        # ----- Tab 5: Dendrogram Analysis -----
        frame_dendro = ttk.Frame(notebook, padding=10)
        notebook.add(frame_dendro, text="Dendrogram Analysis")
        ttk.Label(frame_dendro, text="Select Feature Set for Dendrogram:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.combo_dendro = ttk.Combobox(frame_dendro, width=25, values=["Nucleotide Composition", "Additional Features", "Combined Features"])
        self.combo_dendro.grid(row=0, column=1, padx=5, pady=5)
        self.combo_dendro.current(0)
        ttk.Button(frame_dendro, text="Run Dendrogram", command=self.run_dendrogram).grid(row=1, column=1, padx=5, pady=15)
        info_dendro = ("Dendrogram Analysis: Uses hierarchical clustering to build dendrograms that reveal natural groupings in the data.")
        lbl_info_dendro = ttk.Label(frame_dendro, text=info_dendro, wraplength=800, foreground="blue")
        lbl_info_dendro.grid(row=2, column=0, columnspan=3, padx=5, pady=10)

        # ----- Tab 6: Sequence Extraction -----
        frame_extract = ttk.Frame(notebook, padding=10)
        notebook.add(frame_extract, text="Sequence Extraction")
        ttk.Label(frame_extract, text="Populate Genes:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        ttk.Button(frame_extract, text="Populate Genes", command=self.refresh_common_genes).grid(row=0, column=1, padx=5, pady=5)
        ttk.Label(frame_extract, text="Select Gene:").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.combo_common_genes = ttk.Combobox(frame_extract, width=20)
        self.combo_common_genes.grid(row=1, column=1, padx=5, pady=5)
        ttk.Button(frame_extract, text="Download Selected Gene Sequences", command=self.download_selected_gene).grid(row=2, column=0, columnspan=2, padx=5, pady=10)
        ttk.Separator(frame_extract, orient='horizontal').grid(row=3, column=0, columnspan=2, sticky="ew", padx=5, pady=5)
        ttk.Button(frame_extract, text="Download Full Sequences", command=self.download_full_sequences).grid(row=4, column=0, columnspan=2, padx=5, pady=10)
        info_extract = ("Sequence Extraction: Allows you to download either full sequences or gene-specific sequences (from all detected genes) in FASTA format.")
        lbl_info_extract = ttk.Label(frame_extract, text=info_extract, wraplength=800, foreground="blue")
        lbl_info_extract.grid(row=5, column=0, columnspan=2, padx=5, pady=10)

        # ----- Tab 7: Summary -----
        frame_summary = ttk.Frame(notebook, padding=10)
        notebook.add(frame_summary, text="Summary")
        ttk.Label(frame_summary, text="Summary of Analysis Results:").pack(anchor="nw", padx=5, pady=5)
        self.text_summary = tk.Text(frame_summary, height=16, wrap='word', state='disabled')
        self.text_summary.pack(fill='both', expand=True, padx=5, pady=5)
        lbl_copyright = ttk.Label(frame_summary,
                                  text="GeneBank Genie v1.0 – © Dr. Yash Munnalal Gupta",
                                  foreground="gray")
        lbl_copyright.pack(side='bottom', pady=5)

if __name__ == "__main__":
    root = tk.Tk()
    app = GenBankAnalysisApp(root)
    root.mainloop()
