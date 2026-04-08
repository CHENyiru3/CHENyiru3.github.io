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
    <p>Yiru's Version: 2026/04/05</p>
    <p class="profile-credential">
      <span class="profile-credential-degree">BS</span>
      <span class="profile-credential-text">@ Bioinformatics, ZJU-UoE</span>
      <span class="profile-credential-flags" aria-label="China and United Kingdom">🇨🇳 🇬🇧</span>
    </p>
    <p class="profile-credential">
      <span class="profile-credential-degree">PhD</span>
      <span class="profile-credential-text">@ Quantitative Biology and Medicine, Duke-NUS</span>
      <span class="profile-credential-flags" aria-label="Singapore and United States">🇸🇬 🇺🇸</span>
    </p>

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

{% assign cv_pdf_filename = site.cv_pdf_filename %}
{% assign cv_pdf_date = cv_pdf_filename | regex_replace: '^CHEN_Yiru_CV_(\d{4}\.\d{2}\.\d{2})\.pdf$', '\1' %}
{% assign cv_pdf_path = '/assets/pdf/' | append: cv_pdf_filename %}

Hi, I'm Yiru Chen.

**Decoding biomedical complexity through the lens of AI and algorithms.**

<a class="btn btn-primary btn-sm" href="{{ cv_pdf_path | relative_url }}" target="_blank" rel="noopener noreferrer" role="button">Download CV (PDF)</a>
<span class="ml-2 small text-muted">PDF version: {{ cv_pdf_date }}</span>
<span class="ml-3 small">
<a href="https://orcid.org/0009-0002-5114-4947" target="_blank" rel="noopener noreferrer">ORCID</a>
&nbsp;·&nbsp;
<a href="https://scholar.google.com/citations?user=Rfv54HwAAAAJ" target="_blank" rel="noopener noreferrer">Google Scholar</a>
&nbsp;·&nbsp;
<a href="https://github.com/CHENyiru3" target="_blank" rel="noopener noreferrer">GitHub</a>
</span>

## Quick Background

I am currently completing my undergraduate training in Bioinformatics at the Zhejiang University-University of Edinburgh Institute (ZJU-UoE), Zhejiang University. In August 2026, I will join the Quantitative Biology and Medicine (QBM) Program at Duke-NUS Medical School as an incoming PhD student with a full scholarship.

My research interests center on spatial omics and representation learning. You can reach me at yiru.22@intl.zju.edu.cn or yiru2chen@gmail.com.

---

### Short research statement

I develop algorithmic and statistical frameworks for **spatial omics**, **computational immunology**, and **biomedical AI**. My current work centers on multimodal integration, spatially resolved T-cell receptor analysis, simulation of spatial transcriptomics, and rigorous evaluation of computational methods. Across projects, I aim to turn mathematical modeling and machine learning into practical tools for decoding complex biological systems.

---

## Publications

<style>
.post .profile {
  width: 100%;
}

@media (min-width: 576px) {
  .post .profile {
    width: 34%;
  }
}

.profile-credential {
  margin: 0.35rem 0 0;
  font-size: 0.82rem;
  line-height: 1.45;
}
.profile-credential-degree {
  display: inline-block;
  min-width: 2.1rem;
  font-weight: 700;
  color: var(--global-theme-color);
}
.profile-credential-text {
  color: var(--global-text-color);
}
.profile-credential-flags {
  margin-left: 0.25rem;
  letter-spacing: 0.08em;
}
.note { margin-bottom: 1em; font-size: 0.95em; color: #666; }
.pub-year {
  margin-top: 1.1rem;
  margin-bottom: 0.55rem;
  font-size: 1rem;
  font-weight: 700;
  color: var(--global-text-color);
}
.pub-subgroup {
  margin-top: 0.95rem;
  margin-bottom: 0.45rem;
  font-size: 0.92rem;
  font-weight: 700;
  color: var(--global-text-color-light);
  letter-spacing: 0.02em;
}
.pub-section { margin-top: 1.2em; margin-bottom: 0.5em; font-size: 1.1em; border-left: 4px solid; padding-left: 10px; }
.pub-section.published { border-color: #2698BA; }
.pub-section.preprint { border-color: #FFA500; }
.pub-section.manuscript { border-color: #C46BDB; }
.pub-section.presentation { border-color: #28A745; }
</style>

{% include publication_sections.liquid %}
