# 🧪 Methodological Background

NetPharmPy implements a systematic **Network Pharmacology** workflow. This approach moves beyond the traditional "one drug, one target" paradigm to explore the complex, multi-target mechanisms typical of natural products and polypharmacology.

## 1. The Core Philosophy: "Network as a Target"
Biological systems function as complex networks. Diseases often arise from the perturbation of these networks rather than the malfunction of a single gene. Similarly, effective treatments—especially from natural products like Curcumin—often work by weakly modulating multiple nodes in a network simultaneously, leading to a synergistic therapeutic effect.

## 2. Pipeline Components & Rationale

### Step A: Target Prediction (Ligand-Based)
Since experimental target validation is expensive, we rely on computational prediction based on chemical similarity.

*   **SwissTargetPrediction:** Uses the principle of "similarity determines function." It compares your molecule's chemical structure (2D/3D fingerprints) against a library of ~370,000 known actives.
    *   *Reference:* Daina, A., et al. (2019). *Nucleic Acids Research*.
*   **SuperPred:** Utilizes machine learning models trained on drug-target interaction data (ChEMBL, BindingDB) and Anatomical Therapeutic Chemical (ATC) codes.
    *   *Reference:* Nickel, J., et al. (2014). *Nucleic Acids Research*.

### Step B: Pathway Mapping (Reactome)
To understand biological context, we map gene targets to pathways.
*   **Reactome:** An open-source, curated database of biological pathways. Unlike GO (Gene Ontology), Reactome provides detailed molecular steps (reactions), offering a more mechanistic view.
    *   *Reference:* Gillespie, M., et al. (2022). *Nucleic Acids Research*.

### Step C: Network Construction (STRING)
Proteins rarely act alone. We build a **Protein-Protein Interaction (PPI)** network to visualize the functional neighborhood of the targets.
*   **STRING DB:** Aggregates interactions from multiple sources:
    *   Experimental data (physical binding).
    *   Co-expression (genes active together).
    *   Text mining (co-mention in literature).
*   **Metric:** We use a "Confidence Score" (default > 0.700) to filter out low-reliability predictions.

### Step D: Network Topology Analysis
We use Graph Theory to identify the most important nodes (Mechanism of Action candidates).
*   **Degree Centrality:** Number of direct connections. High degree = **Hub Protein**. Hubs are critical for network integrity.
*   **Betweenness Centrality:** Measures how often a node acts as a bridge along the shortest path between two other nodes. High betweenness = **Bottleneck/Regulator**.

## 3. Workflow Validity
This *in silico* approach generates **hypotheses**. The results (Hub Proteins) are prioritized candidates for experimental validation (e.g., Western Blot, PCR, or Molecular Docking).

## 4. How to Cite
If you use NetPharmPy for your research, please cite the underlying databases:
1.  **NetPharmPy:** (Link to this repo)
2.  **SwissTargetPrediction:** Daina et al. (2019)
3.  **STRING:** Szklarczyk et al. (2023)
4.  **Reactome:** Gillespie et al. (2022)
