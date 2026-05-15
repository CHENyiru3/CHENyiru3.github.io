import yaml

with open('_research_radar/2026-05-15.md') as f:
    content = f.read()

parts = content.split('---\n', 2)
fm = yaml.safe_load(parts[1])

# Counts
comp = len(fm.get('computational_articles', []))
biomed = len(fm.get('biomedicine_articles', []))
fields = len(fm.get('field_articles', []))
bio = len(fm.get('biotech_articles', []))

print(f'Counts: comp={comp} biomed={biomed} fields={fields} biotech={bio}')

# Null DOI check
bad_doi = []
for k in ['computational_articles', 'biomedicine_articles', 'field_articles', 'biotech_articles']:
    for a in fm.get(k, []):
        if a.get('doi') is None:
            bad_doi.append(f"{k} rank {a.get('rank')}: doi is None")
print(f'Null DOI check: {"PASS" if not bad_doi else bad_doi}')

# Array topics check
bad_topics = []
for k in ['computational_articles', 'biomedicine_articles', 'field_articles', 'biotech_articles']:
    for a in fm.get(k, []):
        if not isinstance(a.get('topics'), list):
            bad_topics.append(f"{k} rank {a.get('rank')}: topics not a list ({type(a.get('topics')).__name__})")
print(f'Array topics check: {"PASS" if not bad_topics else bad_topics}')

# Structural check
lines = content.split('\n')
assert lines[0] == '---', f"Line 1 not ---: {lines[0]}"
assert lines[1] != '', f"Line 2 is blank after opening ---"
last_content = lines[-1] if lines[-1] != '' else lines[-2]
assert last_content == '---', f"Last content line not ---: {last_content}"

# Section key blank line check
for i, line in enumerate(lines):
    if line.rstrip().endswith(':') and any(line.strip().startswith(k) for k in ['computational_articles', 'biomedicine_articles', 'field_articles', 'biotech_articles', 'hot_topic']):
        if i+1 < len(lines) and lines[i+1] == '':
            print(f"FAIL: Blank line after {line.strip()} at line {i+1}")
            raise SystemExit(1)

print("All structural checks PASSED")
