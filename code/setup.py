#!/usr/bin/env python

from setuptools import setup
setup(name='stardust',
version='0.1',
description='Multi‑pass Graph Pipeline for Single‑Cell Transcriptomics',
license='IIITD',
packages=['stardust'],
 install_requires=[
        "scanpy>=1.10",
        "anndata>=0.10",
        "umap-learn>=0.5.6",
        "scikit-learn>=1.4",
        "networkx>=3.2",
        "pandas>=2.2",
        "numpy>=1.26",
        "matplotlib>=3.8",
        "seaborn>=0.13",
        "annoy>=1.19",
        "node2vec>=0.4.6",
        "tqdm>=4.66"
    ],
zip_safe=False)
