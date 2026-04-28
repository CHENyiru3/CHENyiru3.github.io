---
layout: page
title: Research Radar
permalink: /research-radar/
description: Daily academic article recommendations from curated scholarly feeds.
nav: true
nav_order: 4
---

{% assign radar_digests = site.research_radar | sort: "date" | reverse %}
{% assign latest_digest = radar_digests | first %}

<div class="research-radar-intro">
  <p>
    Research Radar is a daily AI-assisted academic reading digest from curated scholarly feeds, maintained by <strong>DeepSeek-V4-Pro</strong>. It is separate from the milestone-oriented News page and focuses only on academic articles, not AI industry news, product updates, blogs, or general commentary.
  </p>
  <p>
    The digest reflects my current academic reading interests through three views: computational methods and AI, biomedical discoveries, and cross-disciplinary breakthroughs from other fields.
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
