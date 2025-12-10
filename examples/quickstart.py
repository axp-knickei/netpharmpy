"""
Network Pharmacology Tool - Quick Start Example
Simple usage without configuration file.
"""

from netpharm import NetworkPharmacology

# Initialize with curcumin (CID: 969516)
npharma = NetworkPharmacology(cid=969516)

print("Starting analysis for curcumin...")
print("This will require manual steps for target prediction.")
print("Please have these websites ready:")
print("  - http://www.swisstargetprediction.ch/")
print("  - https://prediction.charite.de/")
print("\nPress Ctrl+C to cancel\n")

# Step 1: Get compound info
compound = npharma.get_compound_info()

# Step 2: Predict targets (manual)
npharma.predict_targets(swiss_threshold=0.0, superpred_threshold=0.5)

# Step 3: Analyze pathways
npharma.analyze_pathways([
    "cytokine signaling immune system",
    "adaptive immune system",
    "T cell activation"
])

# Step 4: Build network
npharma.build_network(confidence=0.700)

# Step 5: Enrichment (using g:Profiler API)
npharma.enrichment_analysis(method='gprofiler')

print("\n✅ Analysis complete!")
print(f"Results in: {npharma.output_dir}")
