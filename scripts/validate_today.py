#!/usr/bin/env python3
"""Validate today's Research Radar digest."""
import yaml, sys, os

digest_file = '_research_radar/2026-05-23.md'

with open(digest_file) as f:
    content = f.read()

all_passed = True
lines = content.split('\n')

# 1. Front matter delimiters
print("=== Structural Checks ===")
if not content.startswith('---'):
    print('FAIL: Missing opening ---')
    all_passed = False
else:
    print('PASS: Opening --- present')

if not lines[1].startswith('layout:'):
    print('FAIL: Line 2 does not start with layout:')
    print(f'  Got: {lines[1][:50]}')
    all_passed = False
else:
    print('PASS: layout: immediately after ---')

# Find closing ---
fm_end = content.find('---', 3)
if fm_end == -1:
    print('FAIL: Missing closing ---')
    all_passed = False
else:
    print('PASS: Closing --- found')

# 2. YAML parse
print("\n=== YAML Parse ===")
fm = content[3:fm_end].strip()
try:
    data = yaml.safe_load(fm)
    print('PASS: YAML parse OK')
except Exception as e:
    print(f'FAIL: YAML parse error — {e}')
    sys.exit(1)

# 3. Content counts
print("\n=== Content Counts ===")
sections = {
    'computational_articles': (3, 6),
    'biomedicine_articles': (3, 6),
    'field_articles': (3, 6),
    'biotech_articles': (0, 6),
}
for sec, (min_c, max_c) in sections.items():
    articles = data.get(sec, [])
    count = len(articles) if articles else 0
    if sec == 'biotech_articles' and count == 0:
        print(f'PASS: {sec} count = {count} (optional, OK)')
    elif count < min_c:
        print(f'FAIL: {sec} has only {count} articles (min {min_c})')
        all_passed = False
    elif count > max_c:
        print(f'FAIL: {sec} has {count} articles (max {max_c})')
        all_passed = False
    else:
        print(f'PASS: {sec} count = {count}')

# 4. Null DOI check
print("\n=== DOI Checks ===")
fail_doi = False
for section in ['computational_articles', 'biomedicine_articles', 'field_articles', 'biotech_articles']:
    for idx, article in enumerate(data.get(section, [])):
        if article.get('doi') is None:
            print(f'FAIL: Null doi in {section}[{idx}]')
            fail_doi = True
            all_passed = False
if not fail_doi:
    print('PASS: No null DOI entries')

# 5. Topics array check
print("\n=== Topics Array Check ===")
fail_topics = False
for section in ['computational_articles', 'biomedicine_articles', 'field_articles', 'biotech_articles']:
    for idx, article in enumerate(data.get(section, [])):
        topics = article.get('topics')
        if not isinstance(topics, list):
            print(f'FAIL: Non-list topics in {section}[{idx}]')
            fail_topics = True
            all_passed = False
        elif topics is None or len(topics) == 0:
            print(f'WARN: Empty/null topics in {section}[{idx}]')
if not fail_topics:
    print('PASS: All topics are lists')

# 6. Required fields check
print("\n=== Required Fields Check ===")
required = ['rank', 'title', 'authors', 'source', 'published', 'url', 'doi', 'article_type', 'topics', 'recommendation', 'summary', 'why_it_matters', 'why_for_me']
for section in ['computational_articles', 'biomedicine_articles', 'field_articles']:
    for idx, article in enumerate(data.get(section, [])):
        missing = [f for f in required if not article.get(f)]
        if missing:
            print(f'FAIL: {section}[{idx}] missing: {missing}')
            all_passed = False

# 7. hot_topic check
print("\n=== Hot Topic Check ===")
ht = data.get('hot_topic')
if ht:
    if ht.get('title'):
        print('PASS: hot_topic.title present')
    else:
        print('FAIL: hot_topic.title missing')
        all_passed = False
    if isinstance(ht.get('signals'), list):
        print(f'PASS: hot_topic.signals is list ({len(ht["signals"])} items)')
    else:
        print('FAIL: hot_topic.signals not a list')
        all_passed = False
else:
    print('FAIL: hot_topic missing')
    all_passed = False

# 8. Section key blank line check
print("\n=== Blank Line After Section Keys ===")
section_keys = ['computational_articles:', 'biomedicine_articles:', 'field_articles:', 'biotech_articles:']
for i, line in enumerate(lines):
    if line.strip() in section_keys:
        if i+1 < len(lines) and lines[i+1].strip() == '':
            print(f'FAIL: Blank line after {line.strip()} at line {i+1}')
            all_passed = False
print('PASS: No blank lines after section keys' if all_passed else '')

# Summary
print(f"\n{'='*40}")
if all_passed:
    print("ALL CHECKS PASSED")
else:
    print("SOME CHECKS FAILED - fix before committing")
    sys.exit(1)
