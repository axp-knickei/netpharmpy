# 🆘 Troubleshooting & FAQ

This guide addresses common issues encountered while setting up or running NetPharmPy.

---

## 🔌 API & Connectivity Issues

### 1. "ConnectionError" or "HTTP Error 503/504"
*   **Cause:** One of the external services (STRING, Reactome, or PubChem) is temporarily down or under heavy load.
*   **Solution:** 
    *   The pipeline includes automatic retries with exponential backoff for most calls. 
    *   If it fails after all retries, wait 5–10 minutes and try again.
    *   Check your internet connection and verify if the service is reachable in your browser (e.g., [reactome.org](https://reactome.org)).

### 2. "No compound found for given SMILES"
*   **Cause:** The SMILES string provided is invalid or not indexed in PubChem.
*   **Solution:**
    *   Verify the SMILES on [PubChem](https://pubchem.ncbi.nlm.nih.gov/).
    *   Try searching by **CID** (Compound ID) instead of SMILES in your configuration.

### 3. STRING API returns "No interactions found"
*   **Cause:** The interaction confidence threshold is too high, or the gene list is too small/unrecognized.
*   **Solution:**
    *   Lower the `confidence` in `config.yaml` (e.g., from `0.700` to `0.400`).
    *   Ensure the predicted targets use standard **Gene Symbols** (e.g., `TNF` instead of `Tumor Necrosis Factor`).

---

## 🛠️ Step 2: Target Prediction (Manual Workflow)

### 1. "File not found: swiss_results.csv"
*   **Cause:** The file was not saved in the correct directory or with the correct name.
*   **Solution:**
    *   Save the CSV exactly as `swiss_results.csv` inside the `outputs/compound_.../data/` folder.
    *   Do not rename the column headers in the CSV.

### 2. "Could not find probability column"
*   **Cause:** SwissTargetPrediction changed their export format.
*   **Solution:**
    *   Ensure you are downloading the **CSV** format, not the Excel or PDF version.
    *   The parser looks for "Probability", "Probability*", or "Probability (%)". If your file has something else, please report a bug.

---

## 🧬 Scientific & Logic Questions

### 1. "Why are no overlapping targets found?"
*   **Cause:** Your predicted targets don't match any proteins in the Reactome pathways you searched for.
*   **Solution:**
    *   **Broaden Search:** Use more general search terms in `config.yaml` (e.g., use "Immune System" instead of "IL-17 signaling").
    *   **Lower Thresholds:** Lower the prediction thresholds in Step 2 to include more potential targets.
    *   **Check Organism:** Ensure you selected "Homo sapiens" during the manual target prediction steps.

### 2. "What makes a protein a 'Hub'?"
*   **Cause:** Definition of network metrics.
*   **Answer:** Hubs are identified by **Degree Centrality** (total number of connections). In `network_metrics.csv`, the `Degree` column indicates how many other proteins in your network interact with that specific protein. High-degree proteins are often central to the biological mechanism.

### 3. "Can I analyze more than one compound at once?"
*   **Answer:** Currently, the pipeline is designed for single-compound analysis. To compare multiple compounds, run the pipeline separately for each and compare the `network_metrics.csv` files.

---

## 🐛 Still having trouble?

If your issue isn't listed here:
1.  Check the `netpharm.log` file in your output directory for detailed error traces.
2.  Open an issue on the GitHub repository.
3.  Include your `config.yaml` and the SMILES/CID you are using.
