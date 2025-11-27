---
layout: page
title: FEAST
description: Simulation and interpolation of spatial transcriptomics from parameter cloud
importance: 1
category: work
related_publications: true
---

## FEAST: Simulation and interpolation of spatial transcriptomics from parameter cloud

**Status:** Active project, ready for RECOMB 2026 conference submission

**Role:** First author, Independent Research Project

**Institution:** Prof. Maizie Zhou Lab, BME & CS Department, Vanderbilt University

**Duration:** August 2024 – Present (Remote and on-site collaboration)

<a href="/assets/pdf/preprint/FEAST_recomb2026.pdf" class="btn btn-sm btn-outline-primary" target="_blank">
  <i class="fas fa-file-pdf"></i> Preprint PDF
</a>


### Project Overview

FEAST (FEAture-space based modeling for Spatial Transcriptomics) is a computational infrastructure that models ST data within a parameter cloud — a latent manifold encoding gene-level mean, variance, and sparsity. By sampling and perturbing this manifold, FEAST generates high-fidelity synthetic slices with tunable spatial and transcriptional variation, enabling systematic evaluation of clustering, deconvolution, and spatial alignment algorithms.

Beyond two dimensions, FEAST performs 3D parameter-cloud interpolation guided by optimal transport and benchmark alignment, reconstructing continuous tissue architectures while preserving molecular coherence. Together, these capabilities establish FEAST as a foundational platform for standardized benchmarking, data augmentation, and 3D reconstruction in spatial transcriptomics.

**Code:** [GitHub](https://github.com/maiziezhoulab/FEAST) · [PyPI](https://pypi.org/project/FEAST-py/)

### Technical Innovation

**Statistical Modeling:**

- **Parameter Cloud Representation:** Models gene-level mean, variance, and sparsity in a unified latent manifold
- **C-vine Copula:** Captures nonlinear and asymmetric dependencies among gene parameters
- **Flexible Count Models:** Supports Poisson, Negative Binomial, ZIP, and ZINB distributions

**Simulation Framework:**

- **Single Slice Simulation:** Generates high-fidelity 2D ST slices with controllable alterations
- **Expression Alteration:** Tunable perturbations in mean, variance, and sparsity for robustness testing
- **Spatial Transformation:** Controlled geometric transformations for alignment benchmarking

**3D Interpolation:**

- **Optimal Transport:** Wasserstein barycenter-based interpolation between parameter clouds
- **Alignment-guided Coordinates:** Transport-guided spatial coordinate interpolation
- **Volumetric Reconstruction:** Generates intermediate slices for continuous 3D tissue architectures

### Key Results

- **High-fidelity simulation** with near-perfect correlation for gene means and variances across multiple ST platforms (10x Visium, MERFISH, OpenST, Slide-seqV2, Xenium, Stereo-seq)
- **Clustering benchmarking** revealed algorithmic sensitivities under controlled expression and sparsity perturbations
- **Alignment evaluation** demonstrated Spateo's robustness over SPACEL under geometric transformations
- **3D interpolation** achieved >0.9 correlation with ground-truth experimental slices in leave-one-out validation

### Complete Research Cycle

This project represents a comprehensive research experience including:

- **Literature Review:** Extensive survey of spatial transcriptomics simulation methods
- **Idea Exploration:** Creative hypothesis generation and initial concept development
- **Failure Analysis:** Learning from initial approaches that didn't meet expectations
- **Method Validation:** Rigorous testing and validation of novel approaches
- **Model Construction:** Implementation of robust, scalable computational framework
- **Benchmark Development:** Creating comprehensive evaluation metrics and comparison studies
- **Application Studies:** Demonstrating utility across diverse biological scenarios
- **Manuscript Preparation:** Scientific writing and presentation of results

### Professional Development

The project has provided extensive training in:

- **Public Presentation:** Regular presentation of progress at different project milestones
- **Scientific Communication:** Development of clear, compelling research narratives
- **Independent Research:** Self-directed project management and problem-solving
- **International Collaboration:** Remote and on-site work with US-based research team

### Expected Impact

FEAST provides the research community with:

- **Standardized Benchmarking:** Controllable ground-truth datasets for evaluating clustering, deconvolution, and alignment algorithms
- **Data Augmentation:** Generate realistic synthetic ST data for method development and validation
- **3D Reconstruction:** Interpolate missing slices to reconstruct continuous tissue volumes from sparse experimental sections
- **Platform-agnostic Design:** Works across diverse ST technologies with unified parameter-cloud representation

### Conference Submission

The project is currently being prepared for submission to **RECOMB 2026**, one of the premier conferences in computational biology, demonstrating the high quality and significance of this research.

_This project showcases independent research capability and innovative thinking in computational biology._
