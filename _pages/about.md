---
layout: about
title: About
permalink: /
subtitle:

profile:
  align: right
  image: cyr.jpg
  image_circular: false
  more_info: >
    <p>Website version: 2026/04/05</p>

selected_papers: false # includes a list of papers marked as "selected={true}"
social: true # includes social icons at the bottom of the page

announcements:
  enabled: false # includes a list of news items
  scrollable: true # adds a vertical scroll bar if there are more than 3 news items
  limit: 5 # leave blank to include all the news in the `_news` folder

latest_posts:
  enabled: false
  scrollable: true # adds a vertical scroll bar if there are more than 3 new posts items
  limit: 3 # leave blank to include all the blog posts
---

Hi, I'm Yiru Chen.

**Decoding biomedical complexity through the lens of AI and algorithms.**

<a class="btn btn-primary btn-sm" href="{{ '/assets/pdf/CHEN_Yiru_CV.pdf' | relative_url }}" target="_blank" rel="noopener noreferrer" role="button">Download CV (PDF)</a>
<span class="ml-3 small">
<a href="https://orcid.org/0009-0002-5114-4947" target="_blank" rel="noopener noreferrer">ORCID</a>
&nbsp;·&nbsp;
<a href="https://scholar.google.com/citations?user=Rfv54HwAAAAJ" target="_blank" rel="noopener noreferrer">Google Scholar</a>
&nbsp;·&nbsp;
<a href="https://github.com/CHENyiru3" target="_blank" rel="noopener noreferrer">GitHub</a>
</span>

## Quick Background

I am currently completing my undergraduate training in Bioinformatics at the Zhejiang University-University of Edinburgh Institute (ZJU-UoE), Zhejiang University. In August 2026, I will join the Quantitative Biology and Medicine (QBM) Program at Duke-NUS Medical School as an incoming PhD student with a full scholarship.

I will work in Associate Professor Jinmiao Chen's laboratory, where I plan to focus on spatial omics and representation learning. You can reach me at yiru.22@intl.zju.edu.cn or yiru2chen@gmail.com.

<div class="journey-strip" aria-label="Academic journey across China, the United Kingdom, the United States, and Singapore">
  <div class="journey-heading">Academic Journey</div>
  <div class="journey-grid">
    <div class="journey-item">
      <span class="journey-flag" role="img" aria-label="China">🇨🇳</span>
      <span class="journey-country">China</span>
      <span class="journey-school">ZJU</span>
    </div>
    <div class="journey-item">
      <span class="journey-flag" role="img" aria-label="United Kingdom">🇬🇧</span>
      <span class="journey-country">United Kingdom</span>
      <span class="journey-school">UoE</span>
    </div>
    <div class="journey-item">
      <span class="journey-flags" aria-label="United States and Singapore">
        <span class="journey-flag" role="img" aria-label="United States">🇺🇸</span>
        <span class="journey-flag" role="img" aria-label="Singapore">🇸🇬</span>
      </span>
      <span class="journey-country">United States + Singapore</span>
      <span class="journey-school">Duke-NUS</span>
    </div>
  </div>
</div>

---

### Short research statement

I develop algorithmic and statistical frameworks for **spatial omics**, **computational immunology**, and **biomedical AI**. My current work centers on multimodal integration, spatially resolved T-cell receptor analysis, simulation of spatial transcriptomics, and rigorous evaluation of computational methods. Across projects, I aim to turn mathematical modeling and machine learning into practical tools for decoding complex biological systems.

---

## Publications

<style>
.journey-strip {
  margin: 1.5rem 0 0.5rem;
  padding: 1.1rem 1.2rem;
  border: 1px solid var(--global-divider-color);
  border-radius: 14px;
  background:
    linear-gradient(135deg, rgba(38, 152, 186, 0.08), transparent 42%),
    linear-gradient(315deg, rgba(255, 165, 0, 0.08), transparent 38%),
    var(--global-card-bg-color);
}
.journey-heading {
  margin-bottom: 0.75rem;
  font-size: 0.82rem;
  font-weight: 700;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--global-text-color-light);
}
.journey-grid {
  display: grid;
  grid-template-columns: repeat(3, minmax(0, 1fr));
  gap: 0.9rem;
}
.journey-item {
  display: flex;
  flex-direction: column;
  gap: 0.15rem;
  padding-left: 0.85rem;
  border-left: 2px solid var(--global-theme-color);
}
.journey-flags {
  display: flex;
  gap: 0.25rem;
  align-items: center;
}
.journey-flag {
  font-size: 1.3rem;
  line-height: 1;
}
.journey-country {
  font-size: 0.82rem;
  font-weight: 700;
}
.journey-school {
  font-size: 0.78rem;
  color: var(--global-text-color-light);
}
.note { margin-bottom: 1em; font-size: 0.95em; color: #666; }
.pub-section { margin-top: 1.2em; margin-bottom: 0.5em; font-size: 1.1em; border-left: 4px solid; padding-left: 10px; }
.pub-section.published { border-color: #2698BA; }
.pub-section.preprint { border-color: #FFA500; }
.pub-section.presentation { border-color: #28A745; }
@media (max-width: 768px) {
  .journey-grid {
    grid-template-columns: repeat(2, minmax(0, 1fr));
  }
}
@media (max-width: 480px) {
  .journey-grid {
    grid-template-columns: 1fr;
  }
}
</style>

{% include publication_sections.liquid %}
