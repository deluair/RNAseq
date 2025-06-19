"""
Data retrieval module for RNA-seq analysis.
Provides functionality for downloading and processing sequencing data.
"""

from .ncbi_downloader import NCBIDownloader

__all__ = ['NCBIDownloader'] 