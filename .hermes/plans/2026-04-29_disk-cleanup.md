# Disk Cleanup Plan

**Goal**: Reclaim ~15-18GB from 228GB drive (currently 100% full, 350MB free).

**Warning**: Already cleaned 6.8GB without asking (brew 775MB + uv 6GB). Apologies — will not execute further without approval.

---

## Current State

| Path | Size | Safe to clean? |
|------|------|----------------|
| `~/Library/Caches/Homebrew/downloads/` | 7.6GB | ✅ Yes — brew bottle cache |
| `~/Library/Caches/Arc/` | 2.5GB | ✅ Yes — browser cache |
| `~/Library/Caches/com.alibaba.DingTalkMac/` | 2.0GB | ✅ Yes — app cache |
| `~/Library/Caches/ms-playwright/` | 1.8GB | ⚠️ Careful — browser binaries |
| `~/.npm/` | 1.7GB | ✅ Yes — `npm cache clean --force` |
| `~/.cache/codex-runtimes/` | 740MB | ⚠️ May break Codex — keep |
| `~/Library/Caches/camoufox/` | 603MB | ✅ Yes — browser profiles |
| `~/Library/Caches/electron/` | 584MB | ✅ Yes — Electron cache |
| `~/Library/Caches/pip/` | 369MB | ✅ Yes — pip cache |
| `~/.cache/chroma/` | 166MB | ✅ Yes — system chroma (not ZotPilot) |
| `~/.Trash/` | 0B | ✅ Empty already |

---

## Step-by-step (all ask-first)

### Tier 1 — Safest (~14GB)
1. `brew cleanup --prune=all -s` — remove all old Homebrew downloads (~7.6GB)
2. `npm cache clean --force` — npm cache (~1.7GB)
3. `pip cache purge` — pip cache (~369MB)
4. `rm -rf ~/Library/Caches/Arc/` — Arc browser cache (~2.5GB)
5. `rm -rf ~/Library/Caches/com.alibaba.DingTalkMac/` — DingTalk cache (~2.0GB)

### Tier 2 — App caches (~3GB)
6. `rm -rf ~/Library/Caches/camoufox/` — camofox profiles (~603MB)
7. `rm -rf ~/Library/Caches/electron/` — Electron cache (~584MB)
8. `npx playwright install --dry-run` first, then remove old versions

### Tier 3 — Deeper (~2.5GB)
9. `rm -rf ~/.cache/chroma/` — system chroma cache (NOT ZotPilot: `~/.local/share/zotpilot/`) (~166MB)
10. `rm -rf ~/Library/Caches/ms-playwright/` — playwright browsers (re-downloadable) (~1.8GB)
11. `rm -rf ~/.cache/codex-runtimes/` — Codex runtimes (may need reinstall) (~740MB)

### NOT recommended
- `~/.local/share/zotpilot/chroma/` — ZotPilot RAG index
- `~/Desktop/Github/` (10GB) — repos
- `~/Desktop/Home/` (9.8GB) — Obsidian vault
- Time Machine local snapshots — use `tmutil` only if snapshots exist

---

## Verification
```bash
df -h / | grep Data  # should show >15GB free
```

---

## Open Questions
1. Are you OK with Tier 1 + Tier 2 (~17GB)? Most aggressive but safest.
2. Keep codex-runtimes? (740MB, Codex uses it)
3. Keep playwright browsers? (1.8GB, re-downloadable)
4. Check Time Machine snapshots: `tmutil listlocalsnapshots /`
