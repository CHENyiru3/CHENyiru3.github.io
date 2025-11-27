---
layout: page
title: Projects
permalink: /projects/
description: Research projects, posters, and coursework.
nav: true
nav_order: 2
---

## Highlighted Projects

<div class="projects">
{% assign sorted_projects = site.projects | sort: "importance" %}
<div class="row row-cols-1 row-cols-md-3">
{% for project in sorted_projects %}
  {% include projects.liquid %}
{% endfor %}
</div>
</div>

---

## All Projects

<details>
<summary><strong>Posters & Presentations</strong></summary>

<div class="materials-grid">

<div class="material-item">
<h4>SpatialTCR — Best Poster Award</h4>
<p>An integrated platform for high-resolution spatial sequencing of T cell receptor repertoires. Awarded Best Poster at GPB Omics and Bioinformatics Frontiers Symposium, 2025.</p>
<img src="/assets/pdf/poster/GPB_2025_SpatialTCR.jpg" alt="SpatialTCR Poster" style="width:100%; border:1px solid #ddd; border-radius:4px;">
</div>

<div class="material-item">
<h4>Spatial Transcriptomics Database Poster</h4>
<p>A comprehensive database for spatial transcriptomics datasets and analysis tools.</p>
<object data="/assets/pdf/poster/STdatabase_poster.pdf" type="application/pdf" width="100%" height="400">
<a href="/assets/pdf/poster/STdatabase_poster.pdf" class="btn btn-sm btn-outline-primary">Download PDF</a>
</object>
</div>

<div class="material-item">
<h4>Bulk TCR Repertoire Analysis Poster</h4>
<p>Computational methods for analyzing bulk T cell receptor repertoire sequencing data.</p>
<img src="/assets/pdf/poster/bulkTCR_poster.jpeg" alt="Bulk TCR Poster" style="width:100%; border:1px solid #ddd; border-radius:4px;">
</div>

</div>

</details>

<details>
<summary><strong>Selected Coursework</strong></summary>

<div class="materials-grid">

<div class="material-item">
<h4>Computational Molecular Medicine Lab — ICA 1</h4>
<p>Individual coursework focusing on computational approaches in molecular medicine.</p>
<object data="{{ '/assets/pdf/key_coursework/CMML3_ICA1.pdf' | relative_url }}" type="application/pdf" width="100%" height="400">
<a href="{{ '/assets/pdf/key_coursework/CMML3_ICA1.pdf' | relative_url }}" class="btn btn-sm btn-outline-primary">Download PDF</a>
</object>
</div>

<div class="material-item">
<h4>Computational Molecular Medicine Lab — ICA 2</h4>
<p>Advanced coursework on computational methods for disease analysis.</p>
<object data="{{ '/assets/pdf/key_coursework/CMML3_ICA2.pdf' | relative_url }}" type="application/pdf" width="100%" height="400">
<a href="{{ '/assets/pdf/key_coursework/CMML3_ICA2.pdf' | relative_url }}" class="btn btn-sm btn-outline-primary">Download PDF</a>
</object>
</div>

<div class="material-item">
<h4>Genetics & Population Report</h4>
<p>Population genetics analysis and report on genetic variation patterns.</p>
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
