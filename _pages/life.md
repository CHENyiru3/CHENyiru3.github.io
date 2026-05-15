---
layout: page
permalink: /life/
title: Life
description: Beyond the lab — interests, hobbies, and what drives me
nav: true
nav_order: 6
---

<div class="life-page">
  <section class="life-intro">
    <p class="life-kicker">Beyond Research</p>
    <p>
      Research is central to my life, but not the whole of it. The routines and interests below help me stay grounded, curious, and creatively sharp.
    </p>
  </section>

  <section class="life-section-card life-section-text-only">
    <div class="life-section-copy">
      <p class="life-section-index">01</p>
      <h2>Sports &amp; Wellness</h2>
      <p>
        I stay active through <strong>running</strong>, <strong>swimming</strong>, and ball games. Table tennis is the sport I return to most often; I like the combination of rhythm, reflex, and tactical thinking. Regular exercise is also my best reset during long periods of coding, reading, and analysis.
      </p>
    </div>
  </section>

  <section class="life-section-card life-section-text-only">
    <div class="life-section-copy">
      <p class="life-section-index">02</p>
      <h2>Coffee Enthusiast</h2>
      <p>
        I enjoy the craft of coffee, especially the small decisions that shape a cup: bean choice, grind size, water temperature, and extraction time. Lately I have been spending more time on <strong>pour-over coffee</strong>, which appeals to the same part of me that likes careful experimentation and repeatable process design.
      </p>
    </div>
  </section>

  <section class="life-section-card life-section-text-only">
    <div class="life-section-copy">
      <p class="life-section-index">03</p>
      <h2>History &amp; Culture</h2>
      <p>
        I am drawn to the long arc of how ideas, institutions, and cultures evolve. Current interests include <strong>Medieval Western Europe</strong> and <strong>Chinese history</strong>.
      </p>
      <p>
        Studying history gives me perspective on how knowledge accumulates, how paradigms shift, and how large systems change over time. Those lessons translate naturally to scientific work.
      </p>
    </div>
  </section>

  <section class="life-section-card life-section-text-only">
    <div class="life-section-copy">
      <p class="life-section-index">04</p>
      <h2>Gaming &amp; Technology</h2>
      <p>
        I used to spend more time playing games than I do now, but my interest has shifted toward the <strong>design philosophy</strong> and technical craft behind them. As a longtime Nintendo fan, I admire the way their work combines clarity, creativity, and careful user experience design.
      </p>
    </div>
  </section>

  <section class="life-section-card">
    <div class="life-section-copy">
      <p class="life-section-index">05</p>
      <h2>Take Photos</h2>
      <p>
        I like recording daily life with my iPhone or Nikon Z50 II. I photograph landscapes and buildings that catch my eye, and I also enjoy photographing animals such as birds, cats, dogs, and insects.
      </p>
    </div>
    <div class="life-photo-grid" aria-label="Photography gallery">
      <figure class="life-photo-slot life-photo-image life-photo-large">
        <img src="{{ '/assets/img/life/train-station-01.jpg' | relative_url }}" alt="Green train beside a modern station building">
      </figure>
      <figure class="life-photo-slot life-photo-image">
        <img src="{{ '/assets/img/life/photography-02.jpeg' | relative_url }}" alt="Modern tower rising above yellow flowers">
      </figure>
      <figure class="life-photo-slot life-photo-image life-photo-cat">
        <img src="{{ '/assets/img/life/photography-01.jpeg' | relative_url }}" alt="Cat resting in sunlit grass">
      </figure>
      <figure class="life-photo-slot life-photo-image">
        <img src="{{ '/assets/img/life/photography-03.jpeg' | relative_url }}" alt="Yellow flowers with a visiting insect">
      </figure>
      <figure class="life-photo-slot life-photo-image">
        <img src="{{ '/assets/img/life/photography-04.jpeg' | relative_url }}" alt="Tree-lined road under green leaves">
      </figure>
    </div>
  </section>
</div>

<style>
.life-page {
  display: grid;
  gap: 1.6rem;
}
.life-intro {
  border-bottom: 1px solid var(--global-divider-color);
  padding-bottom: 1rem;
}
.life-kicker,
.life-section-index {
  margin-bottom: 0.35rem;
  color: var(--global-theme-color);
  font-size: 0.78rem;
  font-weight: 700;
  letter-spacing: 0.08em;
  text-transform: uppercase;
}
.life-intro p:last-child {
  max-width: 46rem;
  margin-bottom: 0;
  color: var(--global-text-color-light);
}
.life-section-card {
  display: grid;
  grid-template-columns: minmax(0, 0.92fr) minmax(260px, 1.08fr);
  gap: 1.2rem;
  align-items: stretch;
  border: 1px solid var(--global-divider-color);
  border-radius: 8px;
  padding: 1rem;
  background: var(--global-bg-color);
}
.life-section-text-only {
  grid-template-columns: 1fr;
}
.life-section-copy h2 {
  margin: 0 0 0.65rem;
  font-size: 1.25rem;
}
.life-section-copy p {
  color: var(--global-text-color-light);
}
.life-section-copy p:last-child {
  margin-bottom: 0;
}
.life-photo-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  grid-auto-rows: minmax(7rem, 1fr);
  gap: 0.65rem;
}
.life-photo-slot {
  min-height: 7rem;
  display: grid;
  place-items: end start;
  margin: 0;
  overflow: hidden;
  border: 1px dashed var(--global-divider-color);
  border-radius: 8px;
  padding: 0.7rem;
  background:
    linear-gradient(135deg, rgba(38, 152, 186, 0.12), rgba(40, 167, 69, 0.08)),
    var(--global-code-bg-color);
}
.life-photo-large {
  grid-row: span 2;
}
.life-photo-image {
  display: block;
  padding: 0;
  border-style: solid;
  background: var(--global-code-bg-color);
}
.life-photo-image img {
  width: 100%;
  height: 100%;
  display: block;
  object-fit: cover;
}
.life-photo-cat img {
  object-position: 78% center;
}
@media (max-width: 767px) {
  .life-section-card {
    grid-template-columns: 1fr;
  }
}
@media (max-width: 420px) {
  .life-photo-grid {
    grid-template-columns: 1fr;
  }
  .life-photo-large {
    grid-row: auto;
  }
}
</style>
