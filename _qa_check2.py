import yaml
with open('_research_radar/2026-05-05.md') as f:
    content = f.read()
lines = content.split('\n')
first = None
closing = None
for i, l in enumerate(lines):
    if l.strip() == '---':
        if first is None:
            first = i
        else:
            closing = i
            break
fm = '\n'.join(lines[first+1:closing])
data = yaml.safe_load(fm)

# Verify field_articles count
fa = data.get('field_articles', [])
print(f'field_articles count: {len(fa)}')

# Check DOIs that are empty string
for section in ['computational_articles', 'biomedicine_articles', 'field_articles']:
    sec = data.get(section, [])
    for i, a in enumerate(sec):
        doi = a.get('doi')
        if doi == '' or doi is None:
            title = a.get('title', '?')
            print(f'{section}[{i}]: doi={repr(doi)} title={title[:60]}')
