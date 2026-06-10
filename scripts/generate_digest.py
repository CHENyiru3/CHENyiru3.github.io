#!/usr/bin/env python3
"""Phase 1: Fetch papers from academic RSS feeds, score, classify, and output ranked candidates.

This script queries journal RSS feeds, bioRxiv API, and (on Fridays) biotech RSS feeds,
scores papers for relevance, and outputs a JSON file with structured article data
ready for Phase 2 digest generation.

Usage: python3 scripts/generate_digest.py [--date YYYY-MM-DD]
Output: _research_radar/candidates_YYYY-MM-DD.json
"""

import sys
import os
import json
import re
import time
from datetime import date, datetime, timedelta
from collections import defaultdict

# Add repo root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML not installed. Run: pip3 install pyyaml")
    sys.exit(1)

try:
    import feedparser
except ImportError:
    print("ERROR: feedparser not installed. Run: pip3 install feedparser")
    sys.exit(1)

try:
    import requests
except ImportError:
    print("ERROR: requests not installed. Run: pip3 install requests")
    sys.exit(1)


# === CONFIGURATION ===
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Journal AOP RSS feeds (Nature-branded journals use Nature AOP feed pattern)
RSS_FEEDS = {
    # Cell Press
    "Cell": "https://www.cell.com/cell/inpress.rss",
    "Cell Systems": "https://www.cell.com/cell-systems/inpress.rss",
    "Cancer Cell": "https://www.cell.com/cancer-cell/inpress.rss",
    "Immunity": "https://www.cell.com/immunity/inpress.rss",
    # Nature family
    "Nature": "https://www.nature.com/nature.rss",
    "Nature Biotechnology": "https://www.nature.com/nbt.rss",
    "Nature Genetics": "https://www.nature.com/ng.rss",
    "Nature Immunology": "https://www.nature.com/ni.rss",
    "Nature Medicine": "https://www.nature.com/nm.rss",
    "Nature Methods": "https://www.nature.com/nmeth.rss",
    "Nature Communications": "https://www.nature.com/ncomms.rss",
    "Nature Machine Intelligence": "https://www.nature.com/natmachintell.rss",
    "Nature Computational Science": "https://www.nature.com/natcomputsci.rss",
    # Science
    "Science": "https://www.science.org/rss/express.xml",
    # PLOS Comp Biol
    "PLOS Computational Biology": "https://journals.plos.org/ploscompbiol/feed/atom",
}

# bioRxiv subject areas we care about
BIORXIV_SUBJECTS = ["bioinformatics", "cancer_biology", "immunology"]

# Biotech news RSS (Fridays only)
BIOTECH_RSS = {
    "FierceBiotech": "https://www.fiercebiotech.com/rss/xml",
    "Endpoints News": "https://endpts.com/feed/",
    "STAT": "https://www.statnews.com/feed/",
    "GenomeWeb": "https://www.genomeweb.com/rss.xml",
}

# Keywords for scoring and classification
COMPUTATIONAL_KEYWORDS = [
    "machine learning", "deep learning", "neural network", "computational model",
    "algorithm", "bioinformatics", "single-cell", "spatial transcriptom",
    "multi-omic", "network analysis", "language model", "protein structure prediction",
    "drug design", "molecular dynamics", "docking", "generative model",
    "graph neural", "transformer", "attention mechanism", "representation learning",
    "variational", "bayesian", "statistical method", "computational method",
    "in silico", "simulation", "systems biology", "mathematical model",
    "compartment", "ordinary differential equation", "stochastic model",
    "image analysis", "computer vision", "natural language process",
    "foundation model", "flow matching", "diffusion model", "equivariant",
    "TWAS", "GWAS", "polygenic", "genomic prediction", "variant calling",
    "metagenomic", "phylogenetic", "sequence analysis", "structural bioinformatic",
    "drug repurposing", "virtual screening", "quantitative",
    "CRISPR screen", "perturbation", "synthetic lethality",
    "spatial omics", "cell-cell communication", "ligand-receptor",
    "trajectory inference", "lineage tracing", "RNA velocity",
    "transfer learning", "self-supervised", "contrastive learning",
    "knowledge graph", "ontology", "database", "resource",
]

BIOMEDICINE_KEYWORDS = [
    "cancer", "tumour", "tumor", "immunotherapy", "checkpoint", "CAR-T",
    "T cell", "B cell", "NK cell", "macrophage", "dendritic cell",
    "myeloid", "neutrophil", "Treg", "exhaustion", "microenvironment",
    "TME", "metastasis", "angiogenesis", "chemotherapy", "radiotherapy",
    "targeted therapy", "antibody", "ADC", "bispecific",
    "clinical trial", "phase I", "phase II", "phase III",
    "biomarker", "prognostic", "diagnostic", "companion diagnostic",
    "autoimmune", "inflammation", "infectious disease", "vaccine",
    "cytokine", "chemokine", "interferon", "IL-", "TGF",
    "PD-1", "PD-L1", "CTLA-4", "LAG-3", "TIM-3",
    "neoantigen", "tumour antigen", "cancer vaccine",
    "adoptive cell", "TCR", "engineered T cell",
    "organoid", "patient-derived", "xenograft", "PDX",
    "single-cell atlas", "pan-cancer", "tumour heterogeneity",
    "senescence", "apoptosis", "ferroptosis", "pyroptosis",
    "metabolism", "glycolysis", "OXPHOS", "lipid",
    "microbiome", "gut-brain", "neuroimmune",
    "autoantibody", "immune tolerance", "transplant",
    "regeneration", "wound healing", "fibrosis",
    "stem cell", "iPSC", "differentiation", "organoid",
    "TLS", "tertiary lymphoid", "germinal center",
    "antibody response", "plasma cell", "memory",
    "myeloid-derived suppressor", "MDSC", "CAF",
    "EMT", "epithelial-mesenchymal", "drug resistance",
]

FIELD_KEYWORDS = [
    "cryo-EM", "cryo-electron", "X-ray crystallography", "NMR",
    "optogenetics", "chemogenetics", "CRISPR", "gene editing",
    "base editing", "prime editing", "epigenetic editing",
    "synthetic biology", "genetic circuit", "metabolic engineering",
    "protein engineering", "directed evolution", "yeast display",
    "antibody engineering", "nanobody", "proteolysis",
    "autophagy", "ubiquitin", "proteasome", "degradation",
    "phase separation", "condensate", "intrinsically disordered",
    "chromatin", "epigenetic", "histone", "nucleosome",
    "transcription", "enhancer", "promoter", "super-enhancer",
    "splicing", "alternative splicing", "RNA processing",
    "non-coding RNA", "lncRNA", "miRNA", "circRNA",
    "translation", "ribosome", "tRNA", "codon",
    "DNA repair", "recombination", "replication",
    "cell cycle", "mitosis", "meiosis", "checkpoint",
    "cytoskeleton", "actin", "microtubule", "motor protein",
    "membrane", "ion channel", "GPCR", "receptor",
    "signal transduction", "kinase", "phosphatase",
    "metabolism", "mitochondria", "ER", "Golgi",
    "vesicle", "endocytosis", "exocytosis", "secretion",
    "neuroscience", "synapse", "neuron", "glia",
    "development", "morphogen", "organoid", "embryo",
    "evolution", "comparative genomics", "phylogenetic",
    "plant", "Arabidopsis", "photosynthesis",
    "microbiology", "bacteria", "phage", "virus",
    "ecology", "biodiversity", "climate",
    "imaging", "microscopy", "super-resolution",
    "biophysics", "single-molecule", "FRET",
    "proteomics", "metabolomics", "lipidomics",
    "chemical biology", "chemoproteomic", "activity-based",
    "quantum", "spin", "radical pair", "magnetoreception",
    "optically detected", "covalent inhibitor", "PROTAC",
    "molecular glue", "targeted degradation",
    "structure determination", "structural biology",
    "protein design", "de novo design", "Rosetta",
    "AlphaFold", "RoseTTAFold", "structure prediction",
    "molecular simulation", "free energy", "binding affinity",
    "allostery", "conformational", "dynamics",
]


def is_recent(published_str, max_days=7):
    """Check if a publication date is within max_days of today."""
    # Try various date formats
    formats = [
        "%Y-%m-%d",
        "%a, %d %b %Y %H:%M:%S %z",
        "%a, %d %b %Y %H:%M:%S %Z",
        "%Y-%m-%dT%H:%M:%S%z",
        "%Y-%m-%dT%H:%M:%SZ",
    ]
    for fmt in formats:
        try:
            pub_date = datetime.strptime(published_str.strip(), fmt)
            if pub_date.tzinfo:
                pub_date = pub_date.replace(tzinfo=None)
            diff = (datetime.now() - pub_date).days
            return diff <= max_days
        except (ValueError, TypeError):
            continue

    # Fallback: try to extract YYYY-MM-DD
    m = re.search(r'(\d{4})-(\d{2})-(\d{2})', str(published_str))
    if m:
        try:
            pub_date = datetime(int(m.group(1)), int(m.group(2)), int(m.group(3)))
            diff = (datetime.now() - pub_date).days
            return diff <= max_days
        except ValueError:
            pass
    return False


def score_paper(title, summary, source):
    """Score a paper's relevance for Research Radar inclusion."""
    text = f"{title} {summary}".lower()
    score = 0

    # Bonus for high-impact journals
    high_impact = ["nature", "science", "cell"]
    source_lower = source.lower()
    for j in high_impact:
        if j in source_lower:
            score += 3
            break

    # Nature sub-journals are still high-impact
    if "nature" in source_lower and source_lower != "nature":
        score += 2

    # Bonus for computational/methods journals
    comp_journals = ["nature computational", "nature machine intelligence",
                     "nature methods", "plos computational", "cell systems"]
    for j in comp_journals:
        if j in source_lower:
            score += 2

    return score


def classify_paper(title, summary):
    """Classify paper into computational, biomedicine, or field."""
    text = f"{title} {summary}".lower()

    comp_score = sum(1 for kw in COMPUTATIONAL_KEYWORDS if kw.lower() in text)
    biomed_score = sum(1 for kw in BIOMEDICINE_KEYWORDS if kw.lower() in text)
    field_score = sum(1 for kw in FIELD_KEYWORDS if kw.lower() in text)

    # Normalize by category size
    comp_norm = comp_score / len(COMPUTATIONAL_KEYWORDS)
    biomed_norm = biomed_score / len(BIOMEDICINE_KEYWORDS)
    field_norm = field_score / len(FIELD_KEYWORDS)

    if comp_norm > biomed_norm and comp_norm > field_norm:
        return "computational"
    elif biomed_norm > field_norm:
        return "biomedicine"
    else:
        return "field"


def extract_doi(entry, source_name):
    """Extract DOI from a feed entry."""
    # Try prism:doi
    doi = None
    if hasattr(entry, 'prism_doi'):
        doi = entry.prism_doi
    if not doi and hasattr(entry, 'dc_identifier'):
        dc_id = entry.dc_identifier
        if dc_id and 'doi:' in str(dc_id).lower():
            doi = str(dc_id).lower().replace('doi:', '').strip()

    # Try id field
    if not doi and hasattr(entry, 'id'):
        id_str = entry.id
        if id_str and 'doi.org/' in str(id_str):
            doi = str(id_str).split('doi.org/')[-1]

    # Try link
    if not doi and hasattr(entry, 'link'):
        link = entry.link
        if 'doi.org/' in str(link):
            doi = str(link).split('doi.org/')[-1]
        elif 'doi=' in str(link):
            doi = str(link).split('doi=')[-1].split('&')[0]

    return doi


def extract_authors(entry):
    """Extract authors from feed entry."""
    authors = []
    if hasattr(entry, 'authors'):
        for author in entry.authors:
            name = author.get('name', '')
            if name:
                authors.append(name)

    if not authors and hasattr(entry, 'dc_creator'):
        authors = [entry.dc_creator]

    if not authors:
        # Try to extract from summary
        if hasattr(entry, 'summary'):
            summary = entry.summary
            # Common patterns: "Author et al." or "Author, A. et al."
            m = re.search(r'([A-Z][a-z]+(?:\s+[A-Z]\.)?)\s+et\s+al\.', summary)
            if m:
                authors = [f"{m.group(1)} et al."]

    return authors


def fetch_rss_feed(name, url, max_days=7):
    """Fetch and parse an RSS feed, returning recent papers."""
    papers = []
    try:
        resp = requests.get(url, headers={'User-Agent': 'ResearchRadar/1.0'}, timeout=30)
        resp.raise_for_status()
        feed = feedparser.parse(resp.content)

        for entry in feed.entries:
            title = entry.get('title', '')
            summary = entry.get('summary', entry.get('description', ''))
            published = entry.get('published', entry.get('updated', ''))

            if not title:
                continue

            # Clean HTML from summary
            summary = re.sub(r'<[^>]+>', '', summary).strip()

            # Check if recent
            if published and not is_recent(published, max_days):
                continue

            doi = extract_doi(entry, name)
            url_link = entry.get('link', '')

            authors = extract_authors(entry)
            authors_str = ', '.join(authors) if authors else ''

            papers.append({
                'title': title.strip(),
                'source': name,
                'published': published.strip() if published else '',
                'url': url_link,
                'doi': doi or '',
                'summary': summary[:2000],
                'authors_raw': authors_str,
            })

        print(f"  {name}: {len(papers)} recent papers")
    except Exception as e:
        print(f"  {name}: ERROR - {e}")

    return papers


def fetch_biorxiv(subject, max_days=7):
    """Fetch recent bioRxiv papers for a subject area."""
    papers = []
    try:
        # bioRxiv content API
        url = f"https://api.biorxiv.org/details/biorxiv/2024-01-01/2099-12-31/0/100"
        # Actually, let's use the RSS feed instead
        rss_url = f"https://connect.biorxiv.org/biorxiv_xml.php?subject={subject}"
        resp = requests.get(rss_url, headers={'User-Agent': 'ResearchRadar/1.0'}, timeout=30)
        resp.raise_for_status()
        feed = feedparser.parse(resp.content)

        for entry in feed.entries:
            title = entry.get('title', '')
            summary = entry.get('summary', entry.get('description', ''))
            published = entry.get('published', entry.get('updated', ''))

            if not title:
                continue

            summary = re.sub(r'<[^>]+>', '', summary).strip()

            if published and not is_recent(published, max_days):
                continue

            doi = extract_doi(entry, f"bioRxiv ({subject})")
            url_link = entry.get('link', '')

            papers.append({
                'title': title.strip(),
                'source': 'bioRxiv',
                'published': published.strip() if published else '',
                'url': url_link,
                'doi': doi or '',
                'summary': summary[:2000],
                'authors_raw': '',
            })

        print(f"  bioRxiv/{subject}: {len(papers)} recent papers")
    except Exception as e:
        print(f"  bioRxiv/{subject}: ERROR - {e}")

    return papers


def fetch_biotech_rss(name, url, max_days=3):
    """Fetch biotech news RSS feed (Fridays only)."""
    items = []
    try:
        resp = requests.get(url, headers={'User-Agent': 'ResearchRadar/1.0'}, timeout=30)
        resp.raise_for_status()
        feed = feedparser.parse(resp.content)

        for entry in feed.entries:
            title = entry.get('title', '')
            summary = entry.get('summary', entry.get('description', ''))
            published = entry.get('published', entry.get('updated', ''))

            if not title:
                continue

            summary = re.sub(r'<[^>]+>', '', summary).strip()

            if published and not is_recent(published, max_days):
                continue

            url_link = entry.get('link', '')

            items.append({
                'title': title.strip(),
                'source': name,
                'published': published.strip() if published else '',
                'url': url_link,
                'summary': summary[:1000],
            })

        print(f"  {name} (biotech): {len(items)} recent items")
    except Exception as e:
        print(f"  {name} (biotech): ERROR - {e}")

    return items


def main():
    target_date = date.today()
    if len(sys.argv) > 2 and sys.argv[1] == '--date':
        target_date = date.fromisoformat(sys.argv[2])

    date_str = target_date.isoformat()
    is_friday = target_date.weekday() == 4

    print(f"=== Research Radar Phase 1: Paper Fetching ===")
    print(f"Date: {date_str} ({target_date.strftime('%A')})")
    print(f"Biotech: {'YES (Friday)' if is_friday else 'No'}")
    print()

    all_papers = []

    # 1. Fetch journal RSS feeds
    print("--- Journal RSS Feeds ---")
    for name, url in RSS_FEEDS.items():
        papers = fetch_rss_feed(name, url, max_days=7)
        for p in papers:
            p['score'] = score_paper(p['title'], p['summary'], p['source'])
            category = classify_paper(p['title'], p['summary'])
            p['category'] = category
        all_papers.extend(papers)

    # 2. Fetch bioRxiv
    print("\n--- bioRxiv Feeds ---")
    for subject in BIORXIV_SUBJECTS:
        papers = fetch_biorxiv(subject, max_days=7)
        for p in papers:
            p['score'] = score_paper(p['title'], p['summary'], p['source'])
            category = classify_paper(p['title'], p['summary'])
            p['category'] = category
        all_papers.extend(papers)

    # 3. Fetch biotech news (Fridays only)
    biotech_items = []
    if is_friday:
        print("\n--- Biotech News Feeds ---")
        for name, url in BIOTECH_RSS.items():
            items = fetch_biotech_rss(name, url, max_days=7)
            biotech_items.extend(items)

    # 4. Deduplicate by DOI
    print(f"\n--- Deduplication ---")
    print(f"Total before dedup: {len(all_papers)}")
    seen_dois = set()
    deduped = []
    for p in all_papers:
        doi = p.get('doi', '')
        title_key = p['title'].lower()[:100]
        if doi and doi in seen_dois:
            continue
        if doi:
            seen_dois.add(doi)
        # Also check by title to catch papers without DOIs
        if not doi:
            if title_key in [d['title'].lower()[:100] for d in deduped]:
                continue
        deduped.append(p)
    print(f"After dedup: {len(deduped)}")

    # 5. Sort by score and categorize
    computational = sorted([p for p in deduped if p['category'] == 'computational'],
                          key=lambda x: x['score'], reverse=True)
    biomedicine = sorted([p for p in deduped if p['category'] == 'biomedicine'],
                        key=lambda x: x['score'], reverse=True)
    field = sorted([p for p in deduped if p['category'] == 'field'],
                   key=lambda x: x['score'], reverse=True)

    print(f"\n--- Category Breakdown ---")
    print(f"  Computational: {len(computational)}")
    print(f"  Biomedicine: {len(biomedicine)}")
    print(f"  Field: {len(field)}")

    # 6. Output candidates JSON
    output_dir = os.path.join(REPO_ROOT, '_research_radar')
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'candidates_{date_str}.json')

    output = {
        'date': date_str,
        'is_friday': is_friday,
        'fetched_at': datetime.now().isoformat(),
        'computational_articles': computational[:20],  # Top 20 candidates
        'biomedicine_articles': biomedicine[:20],
        'field_articles': field[:20],
        'biotech_articles': biotech_items[:15] if is_friday else [],
    }

    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    print(f"\n✓ Candidates written to {output_file}")
    print(f"  Computational top candidates: {len(computational[:20])}")
    print(f"  Biomedicine top candidates: {len(biomedicine[:20])}")
    print(f"  Field top candidates: {len(field[:20])}")

    # Print top 5 per category for verification
    for cat_name, cat_papers in [("Computational", computational[:5]),
                                  ("Biomedicine", biomedicine[:5]),
                                  ("Field", field[:5])]:
        print(f"\n  Top {cat_name}:")
        for i, p in enumerate(cat_papers):
            print(f"    {i+1}. [{p['source']}] {p['title'][:80]}... (score={p['score']})")

    return 0


if __name__ == '__main__':
    sys.exit(main())
