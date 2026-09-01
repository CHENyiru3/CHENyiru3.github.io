---
layout: page
title: Research Radar
permalink: /research-radar/
description: Curated academic research recommendations for spatial omics, computational immunology, biomedical AI, oncology, and drug discovery.
nav: true
nav_order: 4
---

{% assign radar_digests = site.research_radar | sort: "date" | reverse %}
{% assign latest_digest = radar_digests | first %}

<div class="research-radar-intro">
  <p>
    Research Radar is maintained by <strong>Codexie</strong>, Yiru's Linux-node research engineering assistant. Each digest begins with approved public scholarly RSS feeds, then curates the full candidate pool into a small set of research articles worth reading.
  </p>
  <p>
    <strong>Computational</strong> covers methods, machine learning, and AI that advance biological or clinical inference. <strong>Biomedicine</strong> focuses on spatial omics, computational immunology, cancer biology and therapy, biomedical AI/digital twins, and drug discovery. <strong>Other Fields</strong> includes adjacent work only when it has a concrete connection to those research priorities.
  </p>
  <p>
    Selection favors methodological novelty, transferable evidence, and experimental or clinical grounding. New recommendations state their evidence limit; RSS metadata, model predictions, and preprints are not presented as causal or clinical conclusions. Digests exclude general AI news, product updates, blogs, commentary, and duplicate items.
  </p>
  <p>
    Drafting and QA run locally. A daily digest is published only after validation and a duplicate check confirms that the date has not already been published.
  </p>
  <p class="research-radar-actions">
    <a class="btn btn-sm btn-outline-primary" href="{{ '/research-radar/feed.xml' | relative_url }}">
      <i class="fa-solid fa-square-rss"></i>
      Subscribe via RSS
    </a>
  </p>
</div>

{% if latest_digest %}

  <hr>
  <div class="research-radar-digest-heading">
    <div>
      <p class="research-radar-kicker">Latest digest</p>
      <h2><a href="{{ latest_digest.url | relative_url }}">{{ latest_digest.title }}</a></h2>
    </div>
    <a class="research-radar-detail-link" href="{{ latest_digest.url | relative_url }}">Open digest page</a>
  </div>
  {% include research_radar_sections.liquid digest=latest_digest %}

  <hr>
  <h2>Recent Archive</h2>
  <ul class="research-radar-archive">
    {% for digest in radar_digests limit: 10 %}
      <li>
        <a href="{{ digest.url | relative_url }}">{{ digest.title }}</a>
        <span class="text-muted">{{ digest.date | date: "%b %d, %Y" }}</span>
      </li>
    {% endfor %}
  </ul>
{% else %}
  <p>No Research Radar digests have been generated yet.</p>
{% endif %}
