---
layout: page
permalink: /publications/
title: Publications
description: ""
nav: false
nav_order: 2
---

<style>
.note { font-size: 0.9em; color: #666; margin-bottom: 1.5em; }
</style>

<p class="note">
* Equal contribution &nbsp;·&nbsp; † Corresponding author
</p>

## Peer-Reviewed Publications

<div class="publications">
{% bibliography --query @*[category=published] %}
</div>

## Manuscripts Under Review

<div class="publications">
{% bibliography --query @*[category=preprint] %}
</div>

## Presentations & Posters

<div class="publications">
{% bibliography --query @*[category=presentation] %}
</div>
