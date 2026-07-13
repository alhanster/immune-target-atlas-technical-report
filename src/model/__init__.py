"""PU target-nomination pipeline for immune/autoimmune disease.

Ranks genes as candidate drug targets using approved targets as positives and
everything else as unlabeled (positive-unlabeled learning). See README.md.
"""

__all__ = ["config", "data", "labels", "features", "cv", "model", "evaluate", "shap_explain"]
