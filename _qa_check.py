import yaml, sys

path = '_research_radar/2026-05-07.md'
with open(path) as f:
    content = f.read()

lines = content.split('\n')
all_passed = True

# 2a: Front matter delimiters
if not content.startswith('---'):
    print('FAIL 2a: Missing opening ---')
    all_passed = False
else:
    second_dash = content.find('---', 3)
    if second_dash == -1:
        print('FAIL 2a: Missing closing ---')
        all_passed = False
    else:
        print('PASS 2a: Front matter delimiters OK')

# 2b: No blank line after opening ---
if lines[1].strip() == '':
    print('FAIL 2b: Blank line after opening ---')
    all_passed = False
else:
    print('PASS 2b: No blank after opening ---')

# 2c: Section key blank lines
section_keys = ['computational_articles:', 'biomedicine_articles:', 'field_articles:', 'biotech_articles:']
for i, line in enumerate(lines):
    if line.strip() in section_keys:
        if i+1 < len(lines) and lines[i+1].strip() == '':
            print(f'FAIL 2c: Blank line after {line.strip()} at line {i+1}')
            all_passed = False
if all_passed:
    print('PASS 2c: No blank lines after section keys (prelim)')
# Re-check after possible failures
fail_2c = False
for i, line in enumerate(lines):
    if line.strip() in section_keys:
        if i+1 < len(lines) and lines[i+1].strip() == '':
            if not fail_2c:
                print(f'FAIL 2c: Blank line after {line.strip()} at line {i+1}')
            fail_2c = True
if not fail_2c:
    print('PASS 2c: No blank lines after section keys')

# 2e: YAML parse test
fm_end = content.find('---', 3)
fm = content[3:fm_end].strip()
try:
    data = yaml.safe_load(fm)
    print('PASS 2e: YAML parse OK')
except Exception as e:
    print(f'FAIL 2e: YAML parse error — {e}')
    sys.exit(1)

# 2d: hot_topic indentation
if 'hot_topic' in data and data['hot_topic'] is not None:
    if 'title' in data['hot_topic']:
        print('PASS 2d: hot_topic.title exists')
    else:
        print('FAIL 2d: hot_topic.title missing')
        all_passed = False
else:
    print('FAIL 2d: hot_topic missing or null')
    all_passed = False

# 2f: Null doi check
fail_2f = False
for section in ['computational_articles', 'biomedicine_articles', 'field_articles', 'biotech_articles']:
    if section in data and data[section]:
        for idx, article in enumerate(data[section]):
            if 'doi' in article and article['doi'] is None:
                print(f'FAIL 2f: Null doi in {section}[{idx}]')
                fail_2f = True
                all_passed = False
if not fail_2f:
    print('PASS 2f: No null doi entries')

# 2g: Topics array check
fail_2g = False
for section in ['computational_articles', 'biomedicine_articles', 'field_articles', 'biotech_articles']:
    if section in data and data[section]:
        for idx, article in enumerate(data[section]):
            if 'topics' in article and not isinstance(article['topics'], list):
                print(f'FAIL 2g: Non-list topics in {section}[{idx}]')
                fail_2g = True
                all_passed = False
if not fail_2g:
    print('PASS 2g: All topics are lists')

# Step 3: Content validation
print('--- Step 3: Content Validation ---')
for section in ['computational_articles', 'biomedicine_articles', 'field_articles']:
    articles = data.get(section, [])
    count = len(articles) if articles else 0
    
    if section in ['computational_articles', 'biomedicine_articles', 'field_articles']:
        min_c, max_c = 3, 6
        if count < min_c:
            print(f'WARN: {section} has only {count} articles (min {min_c})')
            all_passed = False
        elif count > max_c:
            print(f'WARN: {section} has {count} articles (max {max_c})')
            all_passed = False
        else:
            print(f'PASS: {section} count = {count} (3-6 OK)')
    
    for idx, article in enumerate(articles):
        missing = []
        if not article.get('title'):
            missing.append('title')
        if not article.get('url'):
            missing.append('url')
        if not article.get('summary'):
            missing.append('summary')
        if missing:
            print(f'FAIL: {section}[{idx}] missing: {missing}')
            all_passed = False

if all_passed:
    print('ALL CHECKS PASSED')
else:
    print('SOME CHECKS FAILED')
