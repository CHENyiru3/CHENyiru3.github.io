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
    <p>This is Yiru, the version on 2025/11/26</p>

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

Hi there — I’m Yiru CHEN. Glad that you find me here!

<a class="btn btn-primary btn-sm" href="/assets/pdf/CHEN_Yiru_CV.pdf" target="_blank" role="button">Download CV (PDF)</a>
<span class="ml-3 small">
<a href="https://orcid.org/0009-0002-5114-4947" target="_blank">ORCID</a>
&nbsp;·&nbsp;
<a href="https://scholar.google.com/citations?user=Rfv54HwAAAAJ" target="_blank">Google Scholar</a>
&nbsp;·&nbsp;
<a href="https://github.com/CHENyiru3" target="_blank">GitHub</a>
</span>

## Quick Background

I grew up in Hunan Province, China. Currently, I am an undergraduate student in the Bioinformatics program at the Zhejiang University–University of Edinburgh (ZJE) Institute, Zhejiang University, expected to graduate in 2026. You are very welcomed to contact me with yiru.22@intl.zju.edu.cn or yiru2chen@gmail.com

---

**Actively applying for PhD positions (2025–2026).** Open to rotations and funded opportunities in computational biology / bioinformatics!

### Short research statement

I develop computational methods for **spatial multi-omics** and **immunology**. My work focuses on statistical modeling, integrative multimodal analysis, and algorithms for spatial transcriptomics and spatially resolved T-cell receptor profiling. Looking ahead, I am eager to explore new research areas while staying grounded in computational approaches — particularly statistical learning and deep learning for biomedcine.

---

## Publications

<style>
.pub-section { margin-top: 1.2em; margin-bottom: 0.5em; font-size: 1.1em; border-left: 4px solid; padding-left: 10px; }
.pub-section.published { border-color: #2698BA; }
.pub-section.preprint { border-color: #FFA500; }
.pub-section.presentation { border-color: #28A745; }
</style>

<p class="note" style="font-size:0.9em;color:#666;margin-bottom:1em;">* Equal contribution · † Corresponding author</p>

<h3 class="pub-section published">Peer-Reviewed Publications</h3>

<div class="publications">
{% bibliography --query @*[category=published] %}
</div>

<h3 class="pub-section preprint">Manuscripts Under Review / Submitted</h3>

<div class="publications">
{% bibliography --query @*[category=preprint] %}
</div>

<h3 class="pub-section presentation">Presentations & Posters</h3>

<div class="publications">
{% bibliography --query @*[category=presentation] %}
</div>
