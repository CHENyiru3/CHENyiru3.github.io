---
layout: page
permalink: /publications/
title: Publications
description: ""
nav: false
nav_order: 2
---

<style>
.note { font-size: 0.9em; color: #666; margin-bottom: 1em; }
.pub-section { margin-top: 1.2em; margin-bottom: 0.5em; font-size: 1.1em; border-left: 4px solid; padding-left: 10px; }
.pub-section.published { border-color: #2698BA; }
.pub-section.preprint { border-color: #FFA500; }
.pub-section.presentation { border-color: #28A745; }
</style>

<p class="note">* Equal contribution · † Corresponding author</p>

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
