# External Workflow: AlphaFold Structure Retrieval & HDOCK Analysis

**Purpose:** This document outlines the standard operating procedure (SOP) for generating the molecular docking scores used in **Figure 4** of the analysis pipeline.

**Scope:** This workflow is **external** to the `netpharmpy` software. The steps described below are performed manually or via external web servers. This repository is responsible only for **consuming** the resulting scores to perform biophysical gating and visualization.

---

## 1. Overview

The "Confirmatory Phase" (Figure 4) relies on biophysical affinity filtering. This requires 3D structures of target proteins and the ligand (Curcumin) to estimate binding energies.

**Workflow Summary:**
1.  **AlphaFold**: Retrieve predicted 3D structures for target proteins.
2.  **Ligand Prep**: Prepare the Curcumin 3D structure.
3.  **HDOCK**: Perform protein-ligand docking to obtain binding energy scores.
4.  **Integration**: Format scores for `netpharmpy`.

---

## 2. Step 1: Protein Structure Retrieval (AlphaFold)

**Goal:** Obtain the most accurate 3D model for each target protein.
**Tool:** [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/)

### Procedure:
1.  Identify the **UniProt ID** for your target gene (e.g., *DYRK2* -> `Q92630`).
2.  Search the UniProt ID in the AlphaFold Database.
3.  **Selection Criteria:**
    *   Select the **primary isoform** (unless specific biology dictates otherwise).
    *   Ensure the model covers the **catalytic domain** or relevant binding site.
4.  **Download:**
    *   Download the PDB file format.
    *   *Naming Convention:* `AF-<UniProtID>-F1-model_v4.pdb`

### Scientific Justification
AlphaFold is utilized here strictly to provide **structural topology** (3D coordinates) for proteins lacking experimental crystal structures. It does **not** perform docking or predict binding affinity itself. The potential impact of local structural uncertainty (e.g., side-chain placement) is mitigated by our **relative scoring approach**: docking scores are normalized against a reference target (DYRK2) rather than interpreted as absolute physical constants. This comparative method focuses on ranking targets by relative interaction potential rather than precise thermodynamic prediction.

---

## 3. Step 2: Ligand Preparation (Curcumin)

**Goal:** Generate a valid PDB file for the small molecule.

### Procedure:
1.  **Source:** Retrieve the 3D structure from [PubChem](https://pubchem.ncbi.nlm.nih.gov/).
    *   Curcumin CID: `969516`
    *   Download format: `SDF` (3D Conformer).
2.  **Conversion:** Convert SDF to PDB format.
    *   *Tools:* OpenBabel, PyMOL, or online converters.
    *   *Command (OpenBabel example):* `obabel curcumin.sdf -O curcumin.pdb --gen3d`
3.  **Verification:** Ensure the PDB contains all hydrogens and correct bond orders.

---

## 4. Step 3: Molecular Docking (HDOCK)

**Goal:** Estimate the binding energy (Docking Score) between the protein and ligand.
**Tool:** [HDOCK Server](http://hdock.phys.hust.edu.cn/)

### Procedure:
1.  **Input Receptor:** Upload the AlphaFold PDB file (from Step 1).
2.  **Input Ligand:** Upload the Curcumin PDB file (from Step 2).
3.  **Parameters:**
    *   **Docking Type:** Default (Blind docking) is standard unless the binding pocket is experimentally known.
    *   **Model:** Default.
4.  **Execution:** Submit the job and wait for completion.
5.  **Output Extraction:**
    *   Navigate to the "Docking Result" table.
    *   Locate **Model 1** (Top ranked model).
    *   Record the **Docking Score** (Dimensionless energy score, typically negative).

### Model Selection Rationale
We select **Model 1** exclusively because HDOCK ranks its output poses internally based on energy scoring functions. Model 1 represents the most energetically favorable predicted binding mode. To maintain methodological consistency and avoid bias, we do **not** average scores across multiple models or cherry-pick lower-ranked poses based on visual inspection. This ensures the ranking process remains algorithmic and reproducible.

### What the Docking Score Is NOT
To ensure appropriate interpretation, please note:
*   It is **NOT** an experimental binding affinity ($K_d$, $K_i$).
*   It is **NOT** a thermodynamic free energy value ($\Delta G$).
*   It is **NOT** directly comparable to scores from other docking engines (e.g., AutoDock Vina scores cannot be compared to HDOCK scores).
*   It is used **ONLY** as a heuristic for ranking and biophysical gating (filtering) within this specific computational pipeline.

---

## 5. Step 4: Data Integration (Handoff to Repository)

**Goal:** Extract docking results from the HDOCK web server and format them into a structured CSV for the `netpharmpy` pipeline.

### Procedure:
1.  **Locate Results Table:** On the HDOCK results page, locate the table titled **"Summary of the Top 10 Models"**.
2.  **Model Selection:**
    *   **Ignore Model 0:** In some cases, Model 0 may be template-based.
    *   **Select Model 1:** Use Model 1 as the source for the docking score. It represents the top-ranked free-docking pose as determined by the HDOCK scoring function.
3.  **Data Extraction:** Record the value from the **"Docking Score"** column for **Model 1** only. Exactly one docking score per protein target must be recorded.
4.  **CSV Formatting:** Create a CSV file (e.g., `docking_results.csv`) following the schema below.

### CSV Schema Requirements
| Column | Description | Example |
| :--- | :--- | :--- |
| `gene_name` | Official Gene Symbol (Must match pipeline targets) | `DYRK2` |
| `docking_score` | Dimensionless score extracted from **HDOCK Model 1** | `-210.53` |
| `structure_id` | Identifier for the receptor (e.g., AlphaFold ID or PDB ID) | `AF-Q92630` |

### Interpretation Constraints
To ensure scientific hygiene, adhere to the following constraints:
*   **Dimensionless Metrics:** Docking scores are internal HDOCK energy scores; they are not measured in kcal/mol or other physical units.
*   **Heuristic Only:** These scores are used strictly for **relative ranking** of targets and are **NOT** experimental binding affinities ($K_d$, $K_i$) or thermodynamic free energies ($\Delta G$).
*   **Normalization:** All scores in this repository will be normalized relative to the reference standard using the formula:
    $$\Delta Score = Score_{DYRK2} - Score_{Target}$$

### Example CSV Content
```csv
gene_name,docking_score,structure_id
DYRK2,-210.53,AF-Q92630
CD19,-245.21,AF-P15391
MAPK1,-198.45,AF-P27361
TNF,-180.12,AF-P01375
```

### Recommended Reproducibility Metadata
To ensure long-term traceability, we recommend recording the following metadata (e.g., in a separate log or notes column):
*   **AlphaFold Version** (e.g., v4)
*   **Download Date** (Structure retrieval)
*   **HDOCK Submission Date**
*   **Ligand Source** (e.g., PubChem CID 969516)

---

## 6. Methodological Context

It is critical to clarify the role of this data in the pipeline:

1.  **Exploratory vs. Confirmatory:**
    *   **NetPharm Pipeline (Figure 3)**: Exploratory. Uses graph topology to find *potential* targets.
    *   **Docking (Figure 4)**: Confirmatory. Uses structural biophysics to *validate* a subset of targets.
2.  **The Reference Standard:**
    *   All scores are normalized against **DYRK2**.
    *   $\Delta Score = Score_{DYRK2} - Score_{Target}$
    *   Only targets with $\Delta Score > 0$ (predicted stronger binding than DYRK2) are advanced to Figure 4B enrichment.
