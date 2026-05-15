---
layout: page
title: Tech Blog
permalink: /tech-blog/
description: Notes on tooling, engineering, and computational workflows.
nav: true
nav_order: 5
---

<div class="tech-blog-shell">
  <section class="tech-blog-lede">
    <p class="tech-blog-kicker">Technical notes</p>
    <p>
      A focused space for practical write-ups on code, research tooling, computational workflows, and the engineering details behind reproducible biomedical analysis.
    </p>
  </section>

  <div class="tech-blog-topic-row" aria-label="Tech blog topics">
    <span><i class="fa-solid fa-code"></i> Code</span>
    <span><i class="fa-solid fa-diagram-project"></i> Workflows</span>
    <span><i class="fa-solid fa-terminal"></i> Tooling</span>
    <span><i class="fa-solid fa-chart-line"></i> Methods</span>
  </div>

  {% assign has_tech_posts = false %}

  <section class="tech-blog-list" aria-label="Tech blog posts">
    {% for post in site.posts %}
      {% assign is_tech_post = false %}
      {% if post.categories contains "tech-blog" or post.categories contains "tech" or post.tags contains "tech" %}
        {% assign is_tech_post = true %}
      {% endif %}

      {% if is_tech_post %}
        {% assign has_tech_posts = true %}
        <article class="tech-blog-post-card">
          <div class="tech-blog-post-meta">
            <time datetime="{{ post.date | date_to_xmlschema }}">{{ post.date | date: "%b %d, %Y" }}</time>
            {% if post.tags %}
              <span>{{ post.tags | join: " / " }}</span>
            {% endif %}
          </div>
          <h2><a href="{{ post.url | relative_url }}">{{ post.title }}</a></h2>
          {% if post.description %}
            <p>{{ post.description }}</p>
          {% endif %}
          <a class="tech-blog-read-link" href="{{ post.url | relative_url }}">Read note</a>
        </article>
      {% endif %}
    {% endfor %}

    {% unless has_tech_posts %}
      <div class="tech-blog-empty">
        <i class="fa-solid fa-pen-nib"></i>
        <p>No tech notes published yet.</p>
      </div>
    {% endunless %}
  </section>
</div>

<style>
.tech-blog-shell {
  display: grid;
  gap: 1.4rem;
}
.tech-blog-lede {
  border-bottom: 1px solid var(--global-divider-color);
  padding-bottom: 1rem;
}
.tech-blog-kicker {
  margin-bottom: 0.35rem;
  color: var(--global-theme-color);
  font-size: 0.78rem;
  font-weight: 700;
  letter-spacing: 0.08em;
  text-transform: uppercase;
}
.tech-blog-lede p:last-child {
  max-width: 46rem;
  margin-bottom: 0;
  color: var(--global-text-color-light);
}
.tech-blog-topic-row {
  display: grid;
  grid-template-columns: repeat(4, minmax(0, 1fr));
  gap: 0.75rem;
}
.tech-blog-topic-row span {
  min-height: 3rem;
  display: inline-flex;
  align-items: center;
  justify-content: center;
  gap: 0.45rem;
  border: 1px solid var(--global-divider-color);
  border-radius: 8px;
  color: var(--global-text-color);
  background: var(--global-bg-color);
  font-size: 0.9rem;
  font-weight: 600;
}
.tech-blog-topic-row i {
  color: var(--global-theme-color);
}
.tech-blog-list {
  display: grid;
  gap: 1rem;
}
.tech-blog-post-card,
.tech-blog-empty {
  border: 1px solid var(--global-divider-color);
  border-radius: 8px;
  padding: 1rem;
  background: var(--global-bg-color);
}
.tech-blog-post-meta {
  display: flex;
  flex-wrap: wrap;
  gap: 0.45rem 0.8rem;
  margin-bottom: 0.4rem;
  color: var(--global-text-color-light);
  font-size: 0.82rem;
}
.tech-blog-post-card h2 {
  margin: 0 0 0.45rem;
  font-size: 1.15rem;
}
.tech-blog-post-card p {
  margin-bottom: 0.65rem;
  color: var(--global-text-color-light);
}
.tech-blog-read-link {
  font-size: 0.88rem;
  font-weight: 700;
}
.tech-blog-empty {
  min-height: 9rem;
  display: grid;
  place-items: center;
  text-align: center;
  color: var(--global-text-color-light);
}
.tech-blog-empty i {
  margin-bottom: 0.5rem;
  color: var(--global-theme-color);
  font-size: 1.35rem;
}
.tech-blog-empty p {
  margin: 0;
}
@media (max-width: 767px) {
  .tech-blog-topic-row {
    grid-template-columns: repeat(2, minmax(0, 1fr));
  }
}
@media (max-width: 420px) {
  .tech-blog-topic-row {
    grid-template-columns: 1fr;
  }
}
</style>
