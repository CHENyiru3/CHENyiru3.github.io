---
layout: page
title: Projects
permalink: /projects/
description: Research projects, posters, and coursework.
nav: false
nav_order: 2
---

## Highlighted Projects

Selected research efforts in algorithm development, computational frameworks, and spatial omics methodology.

<div class="projects">
{% assign sorted_projects = site.projects | sort: "importance" %}
<div class="row row-cols-1 row-cols-md-3">
{% for project in sorted_projects %}
  {% include projects.liquid %}
{% endfor %}
</div>
</div>

---

## Additional Materials

<details>
<summary><strong>Posters & Presentations</strong></summary>

<div class="materials-grid">

<div class="material-item">
<h4>SpatialTCR — Best Poster Award of GPB conference</h4>
<p>An integrated computational framework for high-resolution spatial sequencing of T-cell receptor repertoires. Awarded Best Poster at the GPB Omics and Bioinformatics Frontiers Symposium in 2025.</p>
<img src="{{ '/assets/pdf/poster/GPB_2025_SpatialTCR.jpg' | relative_url }}" alt="SpatialTCR Poster" style="width:100%; border:1px solid #ddd; border-radius:4px;">
</div>

<div class="material-item">
<h4>Spatial Transcriptomics Database Poster</h4>
<p>A database-oriented project for organizing spatial transcriptomics datasets and analysis tools. Presented at the Undergraduate Research Poster Day at Zhejiang University in 2024.</p>
<object data="{{ '/assets/pdf/poster/STdatabase_poster.pdf' | relative_url }}" type="application/pdf" width="100%" height="400">
<a href="{{ '/assets/pdf/poster/STdatabase_poster.pdf' | relative_url }}" class="btn btn-sm btn-outline-primary">Download PDF</a>
</object>
</div>

<div class="material-item">
<h4>Bulk TCR Repertoire Analysis Poster</h4>
<p>Computational analysis of bulk T-cell receptor repertoire sequencing data, presented at the Undergraduate Research Poster Day at Zhejiang University.</p>
<img src="{{ '/assets/pdf/poster/bulkTCR_poster.jpeg' | relative_url }}" alt="Bulk TCR Poster" style="width:100%; border:1px solid #ddd; border-radius:4px;">
</div>

</div>

</details>

<details>
<summary><strong>Selected Coursework</strong></summary>

<div class="materials-grid">

<div class="material-item">
<h4>Computational Modeling and Machine Learning — ICA 1 (95/100), 2025</h4>
<p>Individual coursework on protein interaction analysis using molecular dynamics simulation and AlphaFold3.</p>
<object data="{{ '/assets/pdf/key_coursework/CMML3_ICA1.pdf' | relative_url }}" type="application/pdf" width="100%" height="400">
<a href="{{ '/assets/pdf/key_coursework/CMML3_ICA1.pdf' | relative_url }}" class="btn btn-sm btn-outline-primary">Download PDF</a>
</object>
</div>

<div class="material-item">
<h4>Computational Modeling and Machine Learning — ICA 2 (97/100), 2025</h4>
<p>A benchmarking study of deconvolution methods for spatial transcriptomics.</p>
<object data="{{ '/assets/pdf/key_coursework/CMML3_ICA2.pdf' | relative_url }}" type="application/pdf" width="100%" height="400">
<a href="{{ '/assets/pdf/key_coursework/CMML3_ICA2.pdf' | relative_url }}" class="btn btn-sm btn-outline-primary">Download PDF</a>
</object>
</div>

<div class="material-item">
<h4>Genomics & Proteomics — ICA (97/100), 2024</h4>
<p>ATAC-seq analysis of <em>Drosophila</em>.</p>
<object data="{{ '/assets/pdf/key_coursework/GP_2137_97.pdf' | relative_url }}" type="application/pdf" width="100%" height="400">
<a href="{{ '/assets/pdf/key_coursework/GP_2137_97.pdf' | relative_url }}" class="btn btn-sm btn-outline-primary">Download PDF</a>
</object>
</div>

</div>

</details>

<style>
.materials-grid {
  display: grid;
  gap: 2rem;
  margin: 1.5rem 0;
}
.material-item {
  padding: 1rem;
  border: 1px solid var(--global-divider-color);
  border-radius: 8px;
}
.material-item h4 {
  margin-top: 0;
  color: var(--global-theme-color);
}
.material-item p {
  margin: 0.5rem 0;
  font-size: 0.95rem;
}
</style>
