# ILC Progenitors in Fetal Liver (FL) and Bone Marrow (BM)

## Project Overview

This project focuses on understanding the transcriptional and functional diversity of innate lymphoid cell (ILC) progenitors across fetal liver (FL) and bone marrow (BM). It employs state-of-the-art single-cell RNA sequencing (scRNA-seq) and computational methods to:

1. Identify ILC subsets and progenitors.
2. Analyze transcriptional trajectories.
3. Integrate datasets across tissues.
4. Perform signature-based scoring to compare populations.
5. Explore regulatory networks driving ILC differentiation.
6. Map the populations defined in mouse data in human data that are publicly accessible

---

## Key Analyses

### 1. **Data Preprocessing and Integration**
- **Preprocessing:**
  - Normalized and scaled datasets for FL, BM, and small intestine (SI).
  - Identified variable features using `Seurat`
  - Integrate the data using `Harmony` for batch correction.
- **Integration:**
  - Merged datasets while preserving unique tissue-specific signatures.
  - Performed dimensionality reduction (PCA, UMAP) and clustering.

### 2. **ILC Subset Identification**
- Defined marker genes for distinct ILC subsets:
  - **ILC1:** Tbx21, Eomes, Ifng, etc.
  - **ILC2:** Gata3, Bcl11b, Il4, etc.
  - **ILC3:** Rorc, Il22, Batf, etc.
  - **ILCP:** Zbtb16, Tcf7, Runx3, etc.
  - **CLP:** Cd34, Bcl11a, Rag1, etc.
- Visualized subset-specific signatures using DotPlots and FeaturePlots.

### 3. **Transcriptional Trajectory Analysis**
- Conducted pseudotime analysis using:
  - **Diffusion maps** (via `destiny`).
  - **Monocle 3** to reconstruct developmental trajectories.
- Identified key branch points and pseudotime progression in ILC development.

### 4. **Regulatory Network Analysis**
- **SCENIC Workflow:**
  - Constructed gene regulatory networks.
  - Identified key transcription factors (TFs) driving differentiation.
  - Visualized regulon activity using heatmaps and t-SNE plots.
- Key regulators:
  - Eomes, Tbx21, Gata3, Runx3, Ikzf2, etc.

### 5. **Cluster and Signature-Based Scoring**
- Performed module scoring to:
  - Compare ILC1 vs. ILC3 populations.
  - Assess known signatures across clusters and tissues.
- Violin plots and heatmaps for visualizing module scores.

### 6. **Visualization**
- PCA and UMAP for dimensionality reduction.
- DotPlots and FeaturePlots for gene expression.
- Heatmaps for regulon activity and pseudotime-ordered gene expression.
- Barplots for proportional comparisons across tissues and clusters.

---

## Key Results

- Identified transcriptionally distinct ILC subsets across FL and BM.
- Reconstructed differentiation trajectories highlighting transitions between progenitor states and mature ILC subsets.
- Discovered regulatory networks and transcription factors critical for ILC fate decisions.
- Demonstrated tissue-specific differences in ILC progenitor populations.

---

## Repository Structure

```
├── Data/
│   ├── Raw_Data/                # Raw scRNA-seq data files
│   ├── Processed_Data/          # Processed Seurat objects
├── Scripts/
│   ├── Preprocessing.R          # Preprocessing and normalization
│   ├── Integration.R            # Data integration and Harmony
│   ├── Trajectory_Analysis.R    # Diffusion maps and Monocle 3
│   ├── SCENIC_Analysis.R        # Regulatory network analysis
│   ├── Signature_Scoring.R      # Module scoring and signature analysis
│   ├── Visualization.R          # PCA, UMAP, heatmaps, and plots
├── Figures/
│   ├── PCA_Analysis/            # PCA-related figures
│   ├── Trajectories/            # Pseudotime and trajectory plots
│   ├── Regulon_Activity/        # SCENIC regulon heatmaps
├── README.md                    # Project description and instructions
└── LICENSE                      # License information
```

---

## How to Reproduce the Analysis

1. **Install Dependencies:**
   - R packages: `Seurat`, `Harmony`, `SCENIC`, `ggplot2`, `ComplexHeatmap`, etc.
   - Ensure all required scripts and raw data are available in the repository.

2. **Run Scripts:**
   - Follow the order in the `Scripts/` directory to reproduce each step of the analysis.

3. **Generate Figures:**
   - Visualization scripts automatically save outputs to the `Figures/` directory.

---

## Acknowledgments

For questions or issues, contact [rebuffet@ciml.univ-mrs.fr]


