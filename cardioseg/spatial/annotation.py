import numpy as np
from numpy.linalg import norm

MARKER_DICT = {
    "Cardiomyocyte": [
        "Ryr2", "Fgf12", "Mlip", "Ttn", "Fhl2", "Rbm20", "Ankrd1", 
        "Tecrl", "Mybpc3", "Sgcd", "Myh6", "Myl2", "Myl3", "Tnnt2", 
        "Atp2a2", "Actc1", "Pln" 
    ],
    "Fibroblast": [
        "Dcn", "Apod", "Acsm3", "Cfd", "Abca9", "Mgp", "Lama2", 
        "Kazn", "Abca6", "Gsn", "Pdgfra", "Apoe", "Col1a1", "Col3a1", "Vim" 
    ],
    "Endothelial": [
        "Vwf", "Ldb2", "Ano2", "Ifi27", "Emcn", "Ptprb", "Flt1", 
        "Mecom", "Pecam1", "Egfl7", "Igfbp7", "Cldn5", "Tmsb4x"
    ],
    "Lymphatic endothelial": [
        "Ccl21", "Mmrn1", "Pkhd1l1", "Nrg3", "Reln", "Sema3a", 
        "Sntg2", "Ppfibp1", "Stox2", "St6galnac3"
    ],
    "Mast": [
        "Il18r1", "Kit", "Slc24a3", "Ntm", "Cpa3", "Slc8a3", 
        "Cdk15", "Hpgds", "Slc38a11", "Rab27b"
    ],
    "Mesothelial": [
        "Itln1", "Pdzrn4", "Slc39a8", "Prg4", "Gfpt2", "C3", 
        "Wwc1", "Kcnt2", "Has1", "Wt1"
    ],
    "Myeloid": [
        "F13a1", "Rbpj", "Cd163", "Rbm47", "Mrc1", "Fmn1", 
        "Ms4a6a", "Msr1", "Frmd4b", "Mertk"
    ],
    "Neuronal": [
        "Nrxn1", "Cadm2", "Nrxn3", "Xkr4", "Cdh19", "Chl1", 
        "Kirrel3", "Sorcs1", "Ncam2", "Gpm6b"
    ],
    "NK_T": [
        "Parp8", "Il7r", "Themis", "Aoah", "Skap1", "Cd247", 
        "Itk", "Ptprc", "Camk4", "Gnly", "Cd3e"
    ],
    "Pericyte": [
        "Rgs5", "Abcc9", "Gucy1a2", "Egflam", "Frmd3", "Dlc1", 
        "Agt", "Pdgfrb", "Eps8", "Pla2g5"
    ],
    "Smooth muscle": [
        "Myh11", "Itga8", "Acta2", "Tagln", "Carmn", "Kcnab1", 
        "Zfhx3", "Ntrk3", "Prkg1", "Rcan2"
    ],
}

def build_marker_index(marker_dict, gene_names):
    """Standardize names to Title Case and map to indices."""
    clean_gene_names = [g.strip().capitalize() for g in gene_names]
    gene_to_idx = {g: i for i, g in enumerate(clean_gene_names)}

    marker_indices = {}
    for celltype, genes in marker_dict.items():
        valid_indices = [
            gene_to_idx[g.strip().capitalize()]
            for g in genes
            if g.strip().capitalize() in gene_to_idx
        ]
        marker_indices[celltype] = valid_indices
    return marker_indices

def build_reference_profiles(marker_indices, n_genes):
    """Create reference vectors for scoring."""
    profiles = {}
    for celltype, indices in marker_indices.items():
        vec = np.zeros(n_genes)
        if indices:
            vec[indices] = 1.0
        profiles[celltype] = vec
    return profiles

def score_markers(expression_vector, reference_profiles):
    """
    Scores cell types using an optimized similarity check.
    """
    scores = {}
    
    # 1. Pre-calculate expression norm once
    expr_norm = norm(expression_vector)
    if expr_norm == 0:
        return "Unknown", 0.0, 0.0

    # 2. Calculate scores for each cell type
    for celltype, ref_vec in reference_profiles.items():
        # Dot product only focuses on the marker genes
        dot_val = np.dot(expression_vector, ref_vec)
        ref_norm = norm(ref_vec)
        
        # This is the Cosine Similarity calculation
        score = dot_val / (expr_norm * ref_norm + 1e-9)
        scores[celltype] = score

    # 3. Determine the winner
    best_type = max(scores, key=scores.get)
    sorted_scores = sorted(scores.values(), reverse=True)
    
    winner_score = scores[best_type]
    runner_up = sorted_scores[1] if len(sorted_scores) > 1 else 0.0
   
    # 4. Confidence as a ratio of the gap
    confidence = (winner_score - runner_up) / (winner_score + 1e-9)

    # 5. Threshold check (0.02 is usually the sweet spot for spatial)
    if winner_score < 0.02:
        return "Unknown", float(confidence), float(winner_score)

    return best_type, float(confidence), float(winner_score)