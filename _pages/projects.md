---
layout: page
title: projects
permalink: /projects/
description: A growing collection of your cool projects.
nav: true
nav_order: 3
display_categories: [work, fun]
horizontal: false
---

<!-- pages/projects.md -->
<!-- Project file index (published / preprints / posters / coursework) -->
{% comment %} Documents section: groups PDFs under assets/pdf/<category>/ and lists files in reverse name order (use date prefixes like 2025-07-01-title.pdf to order by time). {% endcomment %}
## Project documents & materials

Below are PDFs related to projects (published papers, preprints, posters, and course projects). Upload PDFs into the corresponding folder under `assets/pdf/` and they will appear below.

{% assign categories = "published,preprint,poster,key_coursework" | split: ',' %}
{% for cat in categories %}
### {{ cat | replace:'_',' ' | capitalize }}

{% assign base = '/assets/pdf/' | append: cat | append: '/' %}
{% assign files = site.static_files | where_exp: "f", "f.path contains base" %}
{% if files == empty %}
_No items yet. Upload PDFs to `assets/pdf/{{ cat }}/`._
{% else %}
{% assign sorted = files | sort: 'name' | reverse %}
{% for f in sorted %}
{% assign base_name = f.name | remove: '.pdf' %}
{% assign parts = base_name | split: '-' %}
{% assign display_date = '' %}
{% if parts[0] %}
  {% if parts[0] | size == 4 %}
    {% assign display_date = parts[0] %}
    {% if parts.size > 1 and parts[1] | size == 2 %}
      {% assign display_date = display_date | append: '-' | append: parts[1] %}
      {% if parts.size > 2 and parts[2] | size == 2 %}
        {% assign display_date = display_date | append: '-' | append: parts[2] %}
      {% endif %}
    {% endif %}
  {% endif %}
{% endif %}
- {% if display_date != '' %}{{ display_date }} — {% endif %}[{{ f.name }}]({{ f.path | relative_url }}){:target="_blank"}
{% endfor %}
{% endif %}

{% endfor %}

<div class="projects">
{% if site.enable_project_categories and page.display_categories %}
  <!-- Display categorized projects -->
  {% for category in page.display_categories %}
  <a id="{{ category }}" href=".#{{ category }}">
    <h2 class="category">{{ category }}</h2>
  </a>
  {% assign categorized_projects = site.projects | where: "category", category %}
  {% assign sorted_projects = categorized_projects | sort: "importance" %}
  <!-- Generate cards for each project -->
  {% if page.horizontal %}
  <div class="container">
    <div class="row row-cols-1 row-cols-md-2">
    {% for project in sorted_projects %}
      {% include projects_horizontal.liquid %}
    {% endfor %}
    </div>
  </div>
  {% else %}
  <div class="row row-cols-1 row-cols-md-3">
    {% for project in sorted_projects %}
      {% include projects.liquid %}
    {% endfor %}
  </div>
  {% endif %}
  {% endfor %}

{% else %}

<!-- Display projects without categories -->

{% assign sorted_projects = site.projects | sort: "importance" %}

  <!-- Generate cards for each project -->

{% if page.horizontal %}

  <div class="container">
    <div class="row row-cols-1 row-cols-md-2">
    {% for project in sorted_projects %}
      {% include projects_horizontal.liquid %}
    {% endfor %}
    </div>
  </div>
  {% else %}
  <div class="row row-cols-1 row-cols-md-3">
    {% for project in sorted_projects %}
      {% include projects.liquid %}
    {% endfor %}
  </div>
  {% endif %}
{% endif %}
</div>
