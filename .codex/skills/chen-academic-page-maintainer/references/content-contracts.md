# Content Contracts

## Publications

- Data source: `_bibliography/papers.bib`
- Shared renderer: `_includes/publication_sections.liquid`
- Current top-level grouping:
  - Preprints
  - Publications
  - Manuscripts in Preparation, Under Review, or Forthcoming
  - Presentations and Posters at Conferences
- Manuscripts keep accepted/forthcoming above under-review or preparing items.
- Year grouping may be added if the list grows.
- Never shorten or hide long author lists unless the user explicitly requests it.

## CV

- `cv_pdf_filename` in `_config.yml` is the source of truth for the current downloadable CV.
- The visible CV date should be derived from the PDF filename, not manually duplicated.
- Keep the filename convention `CHEN_Yiru_CV_YYYY.MM.DD.pdf` unless the user changes it.
- Treat `_data/cv.yml` as the source of truth for structured CV entries.

## Homepage

- `_pages/about.md` mixes narrative text, profile metadata, and publication rendering.
- In `_pages/about.md` `profile.more_info`, the homepage version line should use the label `Yiru's Version`.
- Treat identity, degree, affiliation, advisor, scholarship, and future-plan statements as protected.
- Check for duplicated facts in `_pages/cv.md`, `_data/cv.yml`, and `_news/*.md` before finalizing changes.

## Projects

- Source: `_projects/*.md`
- Keep each project as its own long-form page unless the user requests structural change.
- Be careful when editing status lines or publication references because they may mirror publication metadata elsewhere.

## News

- Source: `_news/*.md`
- Items can be inline or full posts.
- Milestone announcements are fact-sensitive and may affect homepage or CV consistency.

## Research Radar

- Sources:
  - Daily digest files: `_research_radar/*.md`
  - Index page: `_pages/research-radar.md`
  - RSS feed: `_pages/research-radar-feed.xml`
  - Renderers: `_layouts/research_radar.liquid`, `_includes/research_radar_sections.liquid`, `_includes/research_radar_article.liquid`
- Research Radar is an academic-article recommendation channel and must remain separate from milestone News.
- Keep the route `/research-radar/` and feed `/research-radar/feed.xml`.
- Digest files should use `layout: research_radar`.
- Preserve the two-section structure:
  - `relevant_articles`: Top 5 Relevant Reads
  - `breakthrough_articles`: Top 3 Field Breakthroughs
- Digest files may include an optional `hot_topic` object with `title`, `summary`, and `signals` to highlight the day's strongest academic theme.
- Each article item should include rank, title, authors, source, published date, URL, DOI when available, article type, topics, recommendation, summary, relevance rationale, and Zotero-candidate flag.
- The website renders final Markdown produced by an external agent. It does not perform RSS reading, Zotero inspection, ranking, or automation.
- Exclude AI industry news, blogs, product announcements, newsletters, policy commentary, and Research Radar's own RSS feed.
- Do not manually alter article rankings, summaries, or metadata unless the user requests a correction or regeneration.

## Assets

- Verify referenced PDF, image, and poster paths before finishing.
- Prefer config- or content-driven asset references over duplicated hardcoded paths when the repo already supports derivation.
