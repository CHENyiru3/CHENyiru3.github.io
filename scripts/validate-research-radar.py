#!/usr/bin/env python3
"""Validate Research Radar digest YAML structure and required fields."""
import sys
import yaml

def validate_digest(filepath):
    with open(filepath) as f:
        content = f.read()
    parts = content.split('---')
    if len(parts) < 3:
        print(f"✗ No YAML frontmatter found in {filepath}")
        return 1

    yaml_str = parts[1]
    try:
        data = yaml.safe_load(yaml_str)
    except Exception as e:
        print(f"✗ YAML parse error: {e}")
        return 1

    errors = []

    # Frontmatter required fields
    fm_required = ['layout', 'title', 'date', 'page_last_updated', 'generated_at', 'provider', 'scope']
    for field in fm_required:
        if field not in data:
            errors.append(f"Missing frontmatter field: {field}")

    # hot_topic
    ht = data.get('hot_topic', {})
    if not ht:
        errors.append("Missing hot_topic")
    else:
        for field in ['title', 'summary', 'signals']:
            if field not in ht:
                errors.append(f"hot_topic missing: {field}")
        if 'signals' in ht and not isinstance(ht['signals'], list):
            errors.append("hot_topic.signals is not a list")

    # Article sections
    sections = ['computational_articles', 'biomedicine_articles', 'field_articles']
    required_fields = ['rank', 'title', 'authors', 'source', 'published', 'url',
                       'doi', 'article_type', 'topics', 'recommendation',
                       'summary', 'why_it_matters', 'why_for_me']

    for section in sections:
        articles = data.get(section, [])
        if not articles:
            errors.append(f"{section}: empty or missing (requires >= 1 article)")
            continue
        for i, article in enumerate(articles):
            prefix = f"{section}[{i}] (rank {article.get('rank', '?')})"
            for field in required_fields:
                if field not in article:
                    errors.append(f"{prefix}: missing field '{field}'")
            if 'topics' in article and not isinstance(article['topics'], list):
                errors.append(f"{prefix}: topics is not a list")
            # Check author format - no literal "Last Author"
            if 'authors' in article:
                au = article['authors']
                if 'Last Author' in au or 'Last, F.' in au:
                    errors.append(f"{prefix}: placeholder author detected: '{au[:80]}'")

    # Friday check
    from datetime import date
    today = date.today()
    biotech = data.get('biotech_articles', [])
    if today.weekday() == 4 and not biotech:  # Friday
        errors.append("Friday edition: biotech_articles section missing")
    if today.weekday() != 4 and biotech:
        errors.append(f"Non-Friday edition has biotech_articles ({len(biotech)} items)")

    if errors:
        print(f"✗ {len(errors)} issue(s) found:")
        for e in errors:
            print(f"  • {e}")
        return 1
    else:
        total = sum(len(data.get(s, [])) for s in sections)
        print(f"✓ All checks passed ({total} articles across {len(sections)} sections)")
        return 0


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: validate-research-radar.py <digest.md>")
        sys.exit(1)
    sys.exit(validate_digest(sys.argv[1]))
