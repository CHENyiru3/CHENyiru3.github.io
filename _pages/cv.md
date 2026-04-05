---
layout: page
title: CV
permalink: /cv/
---

{% assign cv_pdf_filename = site.cv_pdf_filename %}
{% assign cv_pdf_date = cv_pdf_filename | regex_replace: '^CHEN_Yiru_CV_(\d{4}\.\d{2}\.\d{2})\.pdf$', '\1' %}
{% assign cv_pdf_path = '/assets/pdf/' | append: cv_pdf_filename %}

# Curriculum Vitae

Current snapshot: I am completing my undergraduate training at the Zhejiang University-University of Edinburgh Institute (ZJU-UoE), Zhejiang University, and will join the Quantitative Biology and Medicine Program at Duke-NUS Medical School as an incoming PhD student in August 2026.

PDF version currently shown below: **{{ cv_pdf_date }}**.

<div class="cv-embed" style="max-width:960px;margin:0 auto;">
  <object data="{{ cv_pdf_path | relative_url }}" type="application/pdf" width="100%" height="900">
    <p>Your browser does not support inline PDFs. <a href="{{ cv_pdf_path | relative_url }}">Download my CV (PDF)</a>.</p>
  </object>
</div>
