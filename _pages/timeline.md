---
layout: page
title: Growth & timeline
permalink: /timeline/
---

# Growth & timeline

This timeline lists notable milestones (education, projects, publications, presentations, posters, course highlights). Add entries in `_data/timeline.yml` as YAML objects with fields: `date` (YYYY or YYYY-MM), `title`, `type`, `description`, and optionally `link`.

{% assign entries = site.data.timeline | sort: 'date' | reverse %}

{% for e in entries %}
## {{ e.date }} — {{ e.title }}

{% if e.type %}_Type: {{ e.type | capitalize }}_{% endif %}

{% if e.description %}{{ e.description }}{% endif %}

{% if e.link %}[Open associated file or repo]({{ e.link }}){% endif %}

---
{% endfor %}

<!-- Try to auto-find uploaded PDFs that match this timeline entry by date prefix or title slug -->
{% for e in entries %}
{% assign title_slug = e.title | downcase | replace: ' ', '-' | replace: ':', '' | replace: ',', '' | replace: '.', '' %}
{% assign match_expr_date = e.date %}
{% assign files = site.static_files | where_exp: "f", "f.path contains '/assets/pdf/' and (f.name downcase contains '#{e.date}' or f.name downcase contains '#{title_slug}')" %}
{% if files and files != empty %}
### Files related to: {{ e.title }}
{% assign sorted = files | sort: 'name' | reverse %}
{% for f in sorted %}
- [{{ f.name }}]({{ f.path | relative_url }}){:target="_blank"}
{% endfor %}
{% endif %}
{% endfor %}
