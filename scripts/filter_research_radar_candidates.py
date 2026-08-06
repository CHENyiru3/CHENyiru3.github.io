#!/usr/bin/env python3
"""Prevent Research Radar from republishing read or previously covered articles.

The script has two modes:
1. Filter a Phase-1 candidate JSON file before it is handed to the writing agent.
   Candidates are removed if their corresponding NetNewsWire entry is currently
   read, or if their DOI/canonical URL/normalised title is already present in a
   prior published Research Radar digest.
2. Validate a newly written digest against all earlier digests before commit.

This is deliberately conservative: if the NetNewsWire database cannot be read,
filter mode exits non-zero rather than allowing the cron job to republish items
whose current read status is unknown.
"""

from __future__ import annotations

import argparse
import copy
import json
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit

try:
    import yaml
except ImportError as exc:  # pragma: no cover - environment prerequisite
    raise SystemExit(f"PyYAML is required: {exc}")

ARTICLE_SECTIONS = (
    "computational_articles",
    "biomedicine_articles",
    "field_articles",
    "biotech_articles",
)
TRACKING_QUERY_KEYS = {"rss", "rssfeed", "ref", "source"}
DATE_FILE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}\.md$")
DOI_RE = re.compile(r"(10\.\d{4,9}/[-._;()/:a-z0-9]+)", re.IGNORECASE)

DEFAULT_NNW_DB = (
    Path.home()
    / "Library/Containers/com.ranchero.NetNewsWire-Evergreen/Data/Library/Application Support"
    / "NetNewsWire/Accounts/2_iCloud/DB.sqlite3"
)


def normalise_title(value: Any) -> str:
    """Return a stable title key that tolerates punctuation and whitespace."""
    return re.sub(r"[^\w]+", "", str(value or "").casefold())


def normalise_doi(value: Any) -> str:
    match = DOI_RE.search(str(value or ""))
    return match.group(1).rstrip(".,;)\]").casefold() if match else ""


def canonical_url(value: Any) -> str:
    """Normalise RSS tracking variants while retaining meaningful query values."""
    raw = str(value or "").strip()
    if not raw:
        return ""
    parts = urlsplit(raw)
    query = [
        (key, val)
        for key, val in parse_qsl(parts.query, keep_blank_values=True)
        if key.casefold() not in TRACKING_QUERY_KEYS
        and not key.casefold().startswith("utm_")
    ]
    path = parts.path.rstrip("/") or "/"
    return urlunsplit(
        (
            parts.scheme.casefold(),
            parts.netloc.casefold(),
            path,
            urlencode(query, doseq=True),
            "",
        )
    )


def keys_for_article(article: dict[str, Any]) -> set[tuple[str, str]]:
    """Build DOI, URL and title keys, deriving a DOI from a URL when possible."""
    doi = normalise_doi(article.get("doi")) or normalise_doi(article.get("url"))
    url = canonical_url(article.get("url"))
    title = normalise_title(article.get("title"))
    keys: set[tuple[str, str]] = set()
    if doi:
        keys.add(("doi", doi))
    if url:
        keys.add(("url", url))
    if title:
        keys.add(("title", title))
    return keys


def front_matter(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---"):
        return {}
    parts = text.split("---", 2)
    if len(parts) < 3:
        return {}
    parsed = yaml.safe_load(parts[1]) or {}
    return parsed if isinstance(parsed, dict) else {}


def published_articles(history_dir: Path, exclude: Path | None = None) -> list[dict[str, Any]]:
    """Extract article records from all published daily digest Markdown files."""
    records: list[dict[str, Any]] = []
    excluded = exclude.resolve() if exclude and exclude.exists() else None
    for path in sorted(history_dir.glob("*.md")):
        if not DATE_FILE_RE.match(path.name):
            continue
        if excluded and path.resolve() == excluded:
            continue
        try:
            data = front_matter(path)
        except (OSError, yaml.YAMLError) as exc:
            raise RuntimeError(f"Cannot parse digest history {path}: {exc}") from exc
        for section in ARTICLE_SECTIONS:
            for article in data.get(section) or []:
                if isinstance(article, dict):
                    records.append({"path": path.name, "article": article})
    return records


def history_index(records: Iterable[dict[str, Any]]) -> dict[tuple[str, str], set[str]]:
    index: dict[tuple[str, str], set[str]] = defaultdict(set)
    for record in records:
        for key in keys_for_article(record["article"]):
            index[key].add(record["path"])
    return index


def nnw_read_index(db_path: Path) -> tuple[dict[str, int], dict[str, list[int]]]:
    """Read NetNewsWire status using URL first and exact title as a fallback."""
    if not db_path.is_file():
        raise RuntimeError(f"NetNewsWire database not found: {db_path}")
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        rows = conn.execute(
            """
            SELECT a.title, COALESCE(a.externalURL, a.url), s.read
            FROM articles AS a
            INNER JOIN statuses AS s ON a.articleID = s.articleID
            """
        ).fetchall()
    except sqlite3.Error as exc:
        raise RuntimeError(f"Cannot read NetNewsWire database: {exc}") from exc
    finally:
        try:
            conn.close()
        except UnboundLocalError:
            pass

    by_url: dict[str, int] = {}
    by_title: dict[str, list[int]] = defaultdict(list)
    for title, url, read in rows:
        if url:
            by_url[canonical_url(url)] = int(read)
        title_key = normalise_title(title)
        if title_key:
            by_title[title_key].append(int(read))
    return by_url, by_title


def read_reason(
    article: dict[str, Any], by_url: dict[str, int], by_title: dict[str, list[int]]
) -> str | None:
    url = canonical_url(article.get("url"))
    if url and by_url.get(url) == 1:
        return "nnw_read_url"
    title_key = normalise_title(article.get("title"))
    statuses = by_title.get(title_key, [])
    # A title-only match is accepted only when every matching NNW record is read.
    if statuses and all(status == 1 for status in statuses):
        return "nnw_read_title"
    return None


def filter_candidates(
    payload: dict[str, Any], history: dict[tuple[str, str], set[str]], db_path: Path
) -> tuple[dict[str, Any], dict[str, int]]:
    selected = payload.get("selected")
    if not isinstance(selected, list):
        raise RuntimeError("Candidate JSON must contain a list under 'selected'.")
    by_url, by_title = nnw_read_index(db_path)
    kept: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    counts: dict[str, int] = defaultdict(int)

    for article in selected:
        if not isinstance(article, dict):
            continue
        matched_history = sorted(
            {path for key in keys_for_article(article) for path in history.get(key, set())}
        )
        if matched_history:
            reason = "published_history"
            detail = ",".join(matched_history)
        else:
            reason = read_reason(article, by_url, by_title)
            detail = ""
        if reason:
            counts[reason] += 1
            excluded.append(
                {
                    "title": article.get("title", ""),
                    "reason": reason,
                    "history": detail,
                }
            )
        else:
            kept.append(article)

    result = copy.deepcopy(payload)
    result["selected"] = kept
    result["dedup_audit"] = {
        "input_count": len(selected),
        "kept_count": len(kept),
        "excluded_count": len(excluded),
        "excluded_by_reason": dict(sorted(counts.items())),
        "excluded": excluded,
    }
    return result, dict(counts)


def validate_digest(path: Path, history_dir: Path) -> list[str]:
    current = front_matter(path)
    records = published_articles(history_dir, exclude=path)
    prior_index = history_index(records)
    seen_here: set[tuple[str, str]] = set()
    errors: list[str] = []
    for section in ARTICLE_SECTIONS:
        for article in current.get(section) or []:
            if not isinstance(article, dict):
                continue
            keys = keys_for_article(article)
            duplicate_prior = sorted(
                {prior for key in keys for prior in prior_index.get(key, set())}
            )
            duplicate_here = any(key in seen_here for key in keys)
            if duplicate_prior or duplicate_here:
                location = ", ".join(duplicate_prior) if duplicate_prior else "same digest"
                errors.append(f"{article.get('title', '<untitled>')} -> {location}")
            seen_here.update(keys)
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history-dir", type=Path, default=Path("_research_radar"))
    parser.add_argument("--nnw-db", type=Path, default=DEFAULT_NNW_DB)
    parser.add_argument("--input", type=Path, help="Phase-1 candidate JSON input")
    parser.add_argument("--output", type=Path, help="Filtered candidate JSON output")
    parser.add_argument("--check-digest", type=Path, help="Validate a rendered digest before commit")
    args = parser.parse_args()

    if args.check_digest:
        errors = validate_digest(args.check_digest, args.history_dir)
        if errors:
            print("Research Radar duplicate validation failed:")
            for error in errors:
                print(f"- {error}")
            return 1
        print(f"No prior-digest duplicates found: {args.check_digest}")
        return 0

    if not args.input or not args.output:
        parser.error("filter mode requires both --input and --output")
    try:
        payload = json.loads(args.input.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise RuntimeError("Candidate JSON root must be an object.")
        target_date = str(payload.get("date", ""))
        target_digest = args.history_dir / f"{target_date}.md" if target_date else None
        history = history_index(published_articles(args.history_dir, exclude=target_digest))
        filtered, counts = filter_candidates(payload, history, args.nnw_db)
    except (OSError, json.JSONDecodeError, RuntimeError) as exc:
        print(f"Research Radar candidate filter failed: {exc}", file=sys.stderr)
        return 2

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(filtered, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    audit = filtered["dedup_audit"]
    print(
        "Research Radar candidates: "
        f"input={audit['input_count']} kept={audit['kept_count']} "
        f"excluded={audit['excluded_count']} reasons={counts}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
