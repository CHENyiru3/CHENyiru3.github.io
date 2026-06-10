#!/usr/bin/env python3
"""Phase 2: Fetch full abstracts and metadata for selected candidates via Crossref API.

Usage: python3 scripts/fetch_abstracts_today.py [date_str]
Input: _research_radar/candidates_{date_str}.json
Output: _research_radar/enriched_{date_str}.json
"""
import json
import sys
import os
import time
import re
import requests

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def fetch_crossref(doi):
    """Fetch full metadata from Crossref API."""
    if not doi:
        return None
    
    url = f"https://api.crossref.org/works/{doi}"
    try:
        resp = requests.get(url, headers={'User-Agent': 'ResearchRadar/1.0'}, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            msg = data.get('message', {})
            
            # Extract abstract
            abstract = msg.get('abstract', '')
            abstract = re.sub(r'<[^>]+>', ' ', abstract)
            abstract = re.sub(r'\s+', ' ', abstract).strip()
            
            # Extract authors
            authors = []
            for au in msg.get('author', []):
                given = au.get('given', '')
                family = au.get('family', '')
                if family:
                    authors.append(f"{family}, {given[0]}." if given else family)
            
            # Extract published date
            pub_parts = msg.get('published-print', {}) or msg.get('published-online', {}) or {}
            pub_date = ''
            if pub_parts:
                date_parts = pub_parts.get('date-parts', [[None]])
                if date_parts and date_parts[0]:
                    parts = date_parts[0]
                    pub_date = '-'.join(str(p).zfill(2) for p in parts if p)
            
            # Type
            article_type = msg.get('type', '')
            
            # Title
            titles = msg.get('title', [''])
            title = titles[0] if titles else ''
            
            # Container
            container = msg.get('container-title', [''])
            journal = container[0] if container else ''
            
            return {
                'title': title,
                'authors': authors,
                'journal': journal,
                'published': pub_date,
                'abstract': abstract[:3000],
                'type': article_type,
            }
    except Exception as e:
        print(f"  Crossref error for {doi}: {e}", file=sys.stderr)
    
    return None


def main():
    date_str = sys.argv[1] if len(sys.argv) > 1 else '2026-06-03'
    
    input_file = os.path.join(REPO_ROOT, '_research_radar', f'candidates_{date_str}.json')
    
    with open(input_file) as f:
        data = json.load(f)
    
    # Collect all candidate DOIs
    all_candidates = []
    for cat in ['computational_articles', 'biomedicine_articles', 'field_articles']:
        for p in data[cat]:
            if p.get('doi'):
                all_candidates.append(p)
    
    print(f"Fetching abstracts for {len(all_candidates)} candidates with DOIs...")
    
    enriched = {}
    for i, p in enumerate(all_candidates):
        doi = p['doi']
        if doi in enriched:
            continue
        
        print(f"  [{i+1}/{len(all_candidates)}] {doi} - {p.get('title', '')[:80]}")
        result = fetch_crossref(doi)
        if result:
            enriched[doi] = result
        time.sleep(0.1)  # Rate limiting
    
    # Merge enriched data back
    for cat in ['computational_articles', 'biomedicine_articles', 'field_articles']:
        for p in data[cat]:
            doi = p.get('doi', '')
            if doi in enriched:
                cr = enriched[doi]
                p['abstract_full'] = cr['abstract']
                p['authors_full'] = cr['authors']
                p['journal_cr'] = cr['journal']
                p['pub_date_cr'] = cr['published']
                p['article_type_cr'] = cr['type']
    
    output_file = os.path.join(REPO_ROOT, '_research_radar', f'enriched_{date_str}.json')
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print(f"\n✓ Enriched data written to {output_file}")
    print(f"  Enriched {len(enriched)} papers with full abstracts")


if __name__ == '__main__':
    main()
