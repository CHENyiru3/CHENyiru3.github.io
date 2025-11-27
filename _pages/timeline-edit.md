---
layout: page
title: Edit timeline
permalink: /timeline/edit/
---

# How to add a timeline entry

This site stores timeline entries in `_data/timeline.yml` as a YAML list. To add a new milestone, edit that file and append an object with these fields:

- `date`: `YYYY` or `YYYY-MM` or `YYYY-MM-DD`
- `title`: short title of the milestone
- `type`: one of `education`, `project`, `publication`, `presentation`, `poster`, `coursework`, `other`
- `description`: a short sentence or two
- `link`: optional URL or relative path to a PDF or repo (e.g., `/assets/pdf/key_coursework/2025-07-homework.pdf`)

Example entry:

```
- date: 2025-07
  title: Reproducibility note: optimized preprocessing
  type: project
  description: Improved pipeline I/O and reduced runtime by 30%.
  link: /assets/pdf/key_coursework/2025-07-repro-note.pdf
```

After saving the file, the timeline page (`/timeline/`) will show the new entry automatically.
