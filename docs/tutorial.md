# NetPharmPy User Tutorial: A Case Study with Curcumin

This tutorial guides you through a complete Network Pharmacology analysis using **NetPharmPy**. We will use **Curcumin**, a well-studied bioactive compound, to demonstrate how to predict targets, map them to immune pathways, and analyze the protein interaction network.

## 🎯 Objective
Identify the molecular mechanism by which Curcumin modulates the **Adaptive Immune System**.

---

## 🚀 Step 1: Configuration

First, ensure your configuration file (`config_curcumin.yaml`) is set up to focus on the immune system.

```yaml
target_prediction:
  swiss_threshold: 0.0  # Capture broad potential targets
  superpred_threshold: 0.5

pathways:
  search_terms:
    - "Adaptive Immune System"  # Specific Reactome pathway name
    - "Cytokine Signaling in Immune system"

string:
  confidence: 0.700  # High confidence interactions only

network:
  scope: "all"  # Analyze all targets, not just those in the pathway
```

## 🏃 Step 2: Running the Pipeline

Execute the main script. The tool uses the provided Compound ID (CID) for Curcumin (`969516`).

```bash
python examples/curcumin_example.py
```

*Alternatively, run via CLI:*
```bash
python main.py --config config_curcumin.yaml
```

---

## 🔍 Step 3: Following the Workflow

The pipeline runs in 5 automated steps. Here is what happens at each stage:

### 1. Compound Retrieval
*   **Action:** Connects to PubChem.
*   **Output:** `outputs/compound_.../step1_compound_info/compound_info.csv`
*   **Check:** Verify the SMILES string matches your expectation.

### 2. Target Prediction (Manual Interaction)
*   **Action:** The script will pause and ask you to visit **SwissTargetPrediction** and **SuperPred**.
*   **Task:**
    1.  Copy the SMILES printed in the terminal.
    2.  Paste it into the websites.
    3.  Download the CSV results.
    4.  Save them to the `outputs/.../data/` folder as instructed.
*   **Resume:** Press ENTER in the terminal to continue.

### 3. Pathway Mapping
*   **Action:** Queries Reactome for "Adaptive Immune System".
*   **Result:** Finds proteins belonging to this pathway.
*   **Overlap:** Identifies which of your predicted targets (from Step 2) are *also* immune-related proteins.

### 4. Network Construction
*   **Action:** Queries the STRING database for interactions between your target proteins.
*   **Output:** `outputs/.../step4_network/network_metrics.csv`
*   **Key Insight:** Look for the **"Hub Proteins"** log output. These are proteins with the highest "Degree" (connections).
    > *Example Log:*
    > `TNF: 45 connections`
    > `MAPK1: 38 connections`

### 5. Visualization
*   **Action:** Generates interactive HTML files.
*   **Output:** `outputs/.../network_visualization.html`

---

## 📊 Step 4: Interpreting Results

Navigate to your output directory (e.g., `outputs/compound_969516_20260202...`).

### 1. The Network Graph (`network_visualization.html`)
Open this file in your browser.
*   **Nodes:** Proteins.
    *   **Red Nodes:** Hub proteins (highly connected).
    *   **Blue Nodes:** Peripheral proteins.
*   **Edges:** Known interactions (thicker line = higher confidence).
*   **Interaction:** You can drag nodes to rearrange the network.

### 2. The Hub Analysis (`network_metrics.csv`)
Open this CSV in Excel. Sort by `Degree` (descending).
*   **Top 5 Proteins:** These are your primary candidates for the drug's Mechanism of Action (MoA).
*   **Betweenness Centrality:** High values here indicate proteins that act as "bridges" or bottlenecks in the signaling network.

---

## 🔬 Next Steps: Validation

Now that you have identified potential targets (e.g., TNF, IL6), you should validate them:
1.  **Literature Search:** Has Curcumin been shown to inhibit these targets?
2.  **Docking Simulation:** Use the **[AlphaFold & HDOCK Workflow](alphafold_hdock_workflow.md)** to verify if Curcumin can physically bind to these proteins.
