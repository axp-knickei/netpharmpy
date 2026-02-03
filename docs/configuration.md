# Configuration Guide

NetPharmPy uses YAML configuration files (e.g., `config_curcumin.yaml`) to control the behavior of the analysis pipeline. This guide explains each parameter and how to tune it for your research.

## 📋 Configuration Structure

A typical configuration file looks like this:

```yaml
target_prediction:
  swiss_threshold: 0.0
  superpred_threshold: 0.5

pathways:
  search_terms:
    - "immune system"
    - "cytokine"
    - "R-HSA-1280218"

string:
  confidence: 0.700

network:
  scope: "all"  # Options: 'all' or 'overlap'

enrichment:
  method: "gprofiler"
  run_david: false
  david_email: "your@email.com"
```

---

## 🔬 Parameter Details

### 1. Target Prediction
*   `swiss_threshold` (Float, 0.0 to 1.0): 
    *   Filters targets from **SwissTargetPrediction** based on the "Probability" score.
    *   *Tip:* Set to `0.0` to capture all possible binders, or `0.1+` for higher confidence.
*   `superpred_threshold` (Float, 0.0 to 1.0):
    *   Filters predicted targets from **SuperPred**.
    *   *Tip:* SuperPred tends to return many results; a threshold of `0.5` is a common starting point.

### 2. Pathways
*   `search_terms` (List of Strings):
    *   A list of keywords or specific **Reactome IDs** (starting with `R-HSA-`).
    *   If a keyword is used, the pipeline retrieves the top 5 matching pathways.
    *   *Tip:* Use specific IDs like `R-HSA-1280215` (Cytokine Signaling) for more reproducible results.

### 3. STRING PPI Network
*   `confidence` (Float, 0.0 to 1.0):
    *   The minimum interaction score required for an edge to be included in the network.
    *   `0.400`: Medium confidence (includes more discovery-oriented data).
    *   `0.700`: High confidence (recommended for conservative analysis).
    *   `0.900`: Highest confidence (only very well-validated interactions).

### 4. Network Scope
*   `scope` (String):
    *   `all`: Build the network using **all** predicted targets from Step 2.
    *   `overlap`: Build the network using **only** targets that also appeared in the Reactome pathways from Step 3.
    *   *Tip:* Use `overlap` if you want a cleaner network focused strictly on the biological mechanism of interest.

### 5. Enrichment
*   `method` (String):
    *   `gprofiler`: Fast, automated enrichment analysis via g:Profiler API.
    *   `david`: Guides you through manual analysis on the DAVID Bioinformatics website.
*   `run_david` (Boolean): Set to `true` if you wish to use automated DAVID services (requires email).
*   `david_email` (String): Your registered DAVID email address.

---

## 🛠️ Usage

To run the pipeline with a specific configuration:

```bash
python main.py --config config_your_setup.yaml
```
