#!/bin/bash

# Step 1: Run the Gene Regulatory Network (GRN) inference using pySCENIC
singularity run -B /mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology:/project /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/rebuffet/docker_files/SCENIC/aertslab-pyscenic-0.12.0.sif \
	pyscenic grn \
       --num_workers 20 \
        -o /project/05_Output/02_RegulatoryNetworkAnalysis/BM/expr_mat.adjacencies.tsv \
        /project/05_Output/02_RegulatoryNetworkAnalysis/data/BM/BM.loom \
        /project/01_Reference/TFList/allTFs_mm.txt
# - This step calculates gene regulatory networks (GRNs) using the provided expression matrix (`BM.loom`) and transcription factor list (`allTFs_mm.txt`).
# - Output is saved as `expr_mat.adjacencies.tsv` in the specified directory.
# - Utilizes 20 workers for parallel processing to speed up the computation.



# Step 1 Alternative: Run GRN inference using Arboreto (if needed)
# Replaces pySCENIC's dask-based GRN inference with Arboreto to prevent potential dask failures.
# Uncomment and modify paths if necessary.
#singularity run -B /mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5:/project /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/rebuffet/docker_files/SCENIC/aertslab-pyscenic-0.12.0.sif \
#	arboreto_with_multiprocessing.py \
#    	/project/05_Output/02_RegulatoryNetworkAnalysis/data/AllNK/AllNK.loom \
#    	/project/01_Reference/TFList/allTFs_hg38.txt \
#    	--method grnboost2 \
#    	--output /project/05_Output/02_RegulatoryNetworkAnalysis/AllNK/expr_mat.adjacencies.tsv \
#    	--num_workers 20 \
#    	--seed 777
# - Uses Arboreto with GRNBoost2 for GRN inference, potentially more robust than dask.
# - Designed for human gene data (`allTFs_hg38.txt`) and saves results in the specified path.



# Step 2: Contextualize GRNs using the motif rankings and annotations        
singularity run -B /mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology:/project /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/rebuffet/docker_files/SCENIC/aertslab-pyscenic-0.12.0.sif \
    	pyscenic ctx \
        /project/05_Output/02_RegulatoryNetworkAnalysis/BM/expr_mat.adjacencies.tsv \
        /project/01_Reference/feather_database/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        /project/01_Reference/feather_database/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        --annotations_fname /project/01_Reference/motif/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
        --expression_mtx_fname /project/05_Output/02_RegulatoryNetworkAnalysis/data/BM/BM.loom \
        --mode "custom_multiprocessing" \
        --output /project/05_Output/02_RegulatoryNetworkAnalysis/BM/regulons.csv \
        --num_workers 20
# - Refines the GRN by identifying motifs associated with target genes using the motif rankings (`feather` files).
# - Requires additional annotations (`motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl`) for accurate motif identification.
# - Outputs the resulting regulons (gene modules) in `regulons.csv`.



# Step 3: Calculate AUC (Area Under the Curve) scores for regulons
singularity run -B /mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology:/project /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/rebuffet/docker_files/SCENIC/aertslab-pyscenic-0.12.0.sif \
    pyscenic aucell \
        /project/05_Output/02_RegulatoryNetworkAnalysis/data/BM/BM.loom \
        /project/05_Output/02_RegulatoryNetworkAnalysis/BM/regulons.csv \
        -o /project/05_Output/02_RegulatoryNetworkAnalysis/BM/auc_mtx.loom \
        --num_workers 10
# - Computes AUC scores for each regulon to evaluate their activity in each cell using the input regulons (`regulons.csv`) and expression matrix (`BM.loom`).
# - Output is saved in a `.loom` file (`auc_mtx.loom`), suitable for downstream visualization or analysis.
# - Utilizes 10 workers for parallel processing.