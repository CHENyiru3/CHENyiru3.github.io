---
layout: page
title: FEAST
description: Simulation and interpolation of spatial transcriptomics from parameter cloud
importance: 1
category: work
related_publications: false
---

## FEAST: Simulation and interpolation of spatial transcriptomics from parameter cloud

**Status:** Independent first-author project; manuscript in final preparation

**Role:** First author, Independent Research Project

**Institution:** ZJU-UoE Institute, Zhejiang University; collaboration with Prof. Maizie Zhou Lab, Vanderbilt University

**Duration:** August 2024 – Present (Remote and on-site collaboration)

<a href="{{ '/assets/pdf/preprint/FEAST_recomb2026.pdf' | relative_url }}" class="btn btn-sm btn-outline-primary" target="_blank" rel="noopener noreferrer">
  <i class="fas fa-file-pdf"></i> Preprint PDF
</a>

### Project Overview

FEAST (FEAture-space based modeling for Spatial Transcriptomics) is an independent first-author project on parameter-cloud modeling for spatial transcriptomics. The framework represents gene-level mean, variance, and sparsity in a latent parameter space, then uses that structure to simulate realistic slices and interpolate tissue states across space.

Beyond two dimensions, FEAST extends to 3D interpolation through optimal-transport-guided transitions between parameter clouds, supporting volumetric reconstruction and systematic benchmarking for spatial omics methods.

**Code:** [GitHub](https://github.com/maiziezhoulab/FEAST) · [PyPI](https://pypi.org/project/FEAST-py/)

### Technical Innovation

**Statistical Modeling:**

- **Parameter-cloud representation:** models gene-level mean, variance, and sparsity in a unified latent space
- **C-vine copula modeling:** captures nonlinear and asymmetric dependence among gene parameters
- **Flexible count models:** supports Poisson, Negative Binomial, ZIP, and ZINB distributions

**Simulation Framework:**

- **Single-slice simulation:** generates realistic 2D spatial transcriptomics slices with controllable perturbations
- **Expression alteration:** tunes mean, variance, and sparsity for robustness testing
- **Spatial transformation:** enables geometry-aware benchmarking for alignment methods

**3D Interpolation:**

- **Optimal transport:** Wasserstein-style interpolation between parameter clouds
- **Alignment-guided coordinates:** transport-informed coordinate interpolation
- **Volumetric reconstruction:** generates intermediate slices for continuous 3D tissue architectures

### Key Results

- **High-fidelity simulation** with near-perfect correlation for gene means and variances across multiple ST platforms (10x Visium, MERFISH, OpenST, Slide-seqV2, Xenium, Stereo-seq)
- **Clustering benchmarking** revealed algorithmic sensitivities under controlled expression and sparsity perturbations
- **Alignment evaluation** demonstrated Spateo's robustness over SPACEL under geometric transformations
- **3D interpolation** achieved >0.9 correlation with ground-truth experimental slices in leave-one-out validation

### Research Scope

This project has involved the full cycle of independent method development:

- literature review and problem formulation
- model construction and implementation
- benchmarking and application studies
- manuscript preparation and figure development

### Current Positioning

FEAST is the clearest expression of my independent research direction so far: building mathematically grounded computational frameworks for simulation, benchmarking, and representation of spatial biological data.
