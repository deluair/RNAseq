"""
Alignment module for RNA-seq analysis.
Handles RNA-seq read alignment to reference genomes.
"""

from .star_aligner import STARAligner

__all__ = ['STARAligner'] 