# External Workflow: DAVID Functional Enrichment Analysis

**Purpose:** This document provides the standard operating procedure (SOP) for performing the functional enrichment analysis used in **Figure 4** (Confirmatory Phase).

**Scope:** While `netpharmpy` automates exploratory analysis (Figure 3), the confirmatory phase uses the [DAVID Knowledgebase](https://david.ncifcrf.gov/) to ensure methodological alignment with the reference study. This step requires manual data entry and retrieval from the DAVID web server.

---

## 1. Prerequisites

Before starting, you must have:
1.  Completed the **Biophysical Gating** step (Filtering targets by docking score).
2.  Generated the target gene list file (usually `gene_list_for_david.txt`).

---

## 2. Step-by-Step Procedure

### A. Data Submission
1.  Open the [DAVID Functional Annotation Tool](https://david.ncifcrf.gov/summary.jsp).
2.  **Upload Tab:**
    *   **Step 1:** Paste your gene list (e.g., `CD19`, `MAPK1`, `AKT1`) into the "Box A" area.
    *   **Step 2:** Select Identifier: **OFFICIAL_GENE_SYMBOL**.
    *   **Step 3:** Select List Type: **Gene List**.
    *   **Step 4:** Click **Submit List**.

### B. Species Selection & Background
1.  After submission, click the **Background** tab (if available) or ensure **Homo sapiens** is selected in the Species tool.
2.  **Critical methodology note:** Ensure the background is set to **Homo sapiens (Whole Genome)** unless a custom background is scientifically justified. This affects the calculation of enrichment p-values.

### C. Quality Control (QC) Check
**Before proceeding, check the "Gene List Manager" or summary header:**
*   Verify the number of **Mapped IDs**.
*   *Action:* If a significant number of genes are "Unmapped" or "Unknown", verify your gene symbols in the input file. Poor mapping will skew statistical results.

### D. Extracting Results
Navigate to the **Functional Annotation Chart** (or "Annotation Summary Results") and extract the following categories. Note that DAVID interfaces may vary slightly; look for these specific database terms:

1.  **Gene Ontology (GO):** 
    *   Expand the **Gene_Ontology** section.
    *   Select **GOTERM_BP_DIRECT** (Biological Process).
    *   Click the **Chart** button to view the table.
    *   **Download:** Click "Download File" (top right of popup). Save/convert the content as `david_go_bp.csv`.
    
2.  **KEGG Pathways:**
    *   Expand the **Pathways** section.
    *   Select **KEGG_PATHWAY**.
    *   Click **Chart**.
    *   **Download:** Save as `david_kegg.csv`.

3.  **Reactome Pathways:**
    *   Expand the **Pathways** section.
    *   Select **REACTOME_PATHWAY**.
    *   *Note:* If Reactome is not listed under Pathways, check the "Functional Categories" or "Protein Domains" sections, or ensure Reactome is enabled in your DAVID preferences.
    *   **Download:** Save as `david_reactome.csv`.

---

## 3. Data Integration (Handoff)

To visualize these results in `netpharmpy` (Figure 4B), place the three CSV files into your project's data directory (e.g., `data/reproduction/david_results/`).

**Required Filenames:**
- `david_go_bp.csv`
- `david_kegg.csv`
- `david_reactome.csv`

**Required CSV Columns:**
Ensure your files are comma-separated and include the following headers (standard DAVID export):
- `Term` (e.g., "GO:0006954~inflammatory response")
- `Count`
- `%` (Percentage of list)
- `Benjamini` (Adjusted p-value)

*Tip: If DAVID provides a tab-delimited text file, open it in a spreadsheet editor and "Save As CSV".*

---

## 4. Methodological Note

### Why use DAVID for Figure 4?
*   **Confirmatory Rigor:** Figure 4 focuses on a small, high-confidence subset of targets (those with $\Delta Score > 0$).
*   **Study Alignment:** The original curcumin/CAR-T study utilized DAVID for this specific subset.
*   **Manual Curation:** DAVID's "Direct" GO terms often provide more specific, manually curated functional descriptions suitable for small gene lists compared to automated high-throughput tools.

### What the Results Represent
These results statistically describe the **biological modules** most likely influenced by the high-affinity targets of your compound. They do not prove clinical efficacy but provide a biophysically-grounded mechanism of action.
