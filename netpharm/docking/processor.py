"""
docking/processor.py

Handles loading, normalization, and filtering of molecular docking scores.
"""

import pandas as pd
from pathlib import Path
from typing import List, Optional, Tuple

class DockingProcessor:
    """
    Processes docking scores to identify high-affinity targets 
    relative to a reference protein (Biophysical Gating).
    """

    def __init__(self, filepath: Path):
        """
        Initialize with path to docking scores CSV.
        
        Expected columns:
        - gene_name: str
        - docking_score: float (binding energy, more negative = stronger)
        """
        self.filepath = Path(filepath)
        if not self.filepath.exists():
            raise FileNotFoundError(f"Docking file not found: {self.filepath}")
        
        self.df = pd.read_csv(self.filepath)
        self._validate_schema()
        self.processed_df = None

    def _validate_schema(self):
        required = {'gene_name', 'docking_score'}
        if not required.issubset(self.df.columns):
            raise ValueError(f"CSV missing required columns: {required - set(self.df.columns)}")

    def normalize_scores(self, reference_gene: str = "DYRK2") -> pd.DataFrame:
        """
        Calculate Delta Docking Score relative to reference.
        
        Delta = Score(Reference) - Score(Target)
        
        Since scores are typically negative energies (e.g. -200),
        if Ref = -200 and Target = -250 (stronger):
        Delta = -200 - (-250) = +50 (Positive means stronger than ref)
        """
        ref_row = self.df[self.df['gene_name'] == reference_gene]
        if ref_row.empty:
            raise ValueError(f"Reference gene '{reference_gene}' not found in data.")
        
        ref_score = ref_row['docking_score'].values[0]
        
        # Create a copy to avoid SettingWithCopy warnings
        self.processed_df = self.df.copy()
        
        # Calculate Delta
        self.processed_df['delta_score'] = ref_score - self.processed_df['docking_score']
        
        # Rank by Delta (highest affinity first)
        self.processed_df = self.processed_df.sort_values('delta_score', ascending=False)
        
        return self.processed_df

    def get_high_affinity_targets(self, threshold: float = 0.0) -> List[str]:
        """
        Get list of genes with Delta Score > threshold.
        Default threshold 0.0 means 'stronger than reference'.
        """
        if self.processed_df is None:
            raise RuntimeError("Must call normalize_scores() before getting targets.")
            
        high_affinity = self.processed_df[self.processed_df['delta_score'] > threshold]
        return high_affinity['gene_name'].tolist()

    def get_processed_data(self) -> pd.DataFrame:
        if self.processed_df is None:
            raise RuntimeError("Must call normalize_scores() first.")
        return self.processed_df
