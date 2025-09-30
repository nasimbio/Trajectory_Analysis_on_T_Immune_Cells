# Trajectory Analysis on T Immune Cells

---

## 📘 Project Overview

This project performs trajectory inference on tumor-infiltrating T cells derived from 8 experimental groups of mice (2 timepoints × 4 treatments: 1 control and 3 treatments at each timepoint).

Data includes:

scRNA-seq gene expression

TCR repertoire information

A Seurat object containing preprocessed and annotated gene expression data was converted into a Monocle3 trajectory object to reconstruct developmental lineages and explore pseudotime dynamics. The results provide insights into how T cells transition between functional states under different treatments and timepoints.

---

## 📊 Analysis Tasks

### **Task 1: Seurat → Monocle3 cell_data_set Conversion and Trajectory Learning**
- Converted Seurat object into a Monocle3 cell_data_set object
- Added gene metadata, cell annotation, and UMAP embeddings from Seurat
- Learned trajectory graph with learn_graph()
- Ordered cells in pseudotime with a defined root population (C1_Sox5_Stem-like)
- Visualized trajectories colored by cluster annotation and pseudotime
- Applied graph_test() to identify genes significantly associated with trajectory progression


### **Task 2: Pseudotime Analysis and Visualization**
- Integrated pseudotime back into Seurat metadata
- Generated UMAP pseudotime plots within Seurat
- Compared pseudotime distributions across treatments and timepoints with faceted boxplots

---

## 💡 Notes

- Raw data is not publicly available due to client ownership and confidentiality.
- Some example outputs plots are organized by task in the `output/` folder.
- This project is designed for both reproducibility and clarity.

---

## 📬 Contact

*Author:* Nasim Rahmatpour 
*Email:* nasimrahmatpour1@gmail.com 
*GitHub:* (https://github.com/nasimbio)
