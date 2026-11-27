"""Enumerate MaveDB experiments by URN and record which carry insertion variants.

The search endpoint caps at 100 results, so URNs are probed directly. Writes a catalogue that the
replication work then draws on.

Usage: python -I 12_mavedb_catalogue.py <out.json> [max_urn]
"""
import json
import sys
import time
import urllib.error
import urllib.request

API = "https://api.mavedb.org/api/v1"
UA = "Mozilla/5.0 (compatible; research-script)"


def get(url, timeout=45):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        return None if e.code == 404 else ""
    except Exception:
        return ""


def main(out_path, max_urn=1400):
    max_urn = int(max_urn)
    cat, misses = {}, 0
    for n in range(1, max_urn + 1):
        found_any = False
        for suf in "abcdefgh":
            txt = get(f"{API}/experiments/urn:mavedb:{n:08d}-{suf}")
            if not txt:
                if suf == "a":
                    break
                continue
            try:
                d = json.loads(txt)
            except Exception:
                continue
            found_any = True
            cat[d["urn"]] = {"title": d.get("title"), "short": d.get("shortDescription"),
                             "score_sets": d.get("scoreSetUrns") or [],
                             "pub": [p.get("identifier") for p in (d.get("primaryPublicationIdentifiers") or [])]}
        if not found_any:
            misses += 1
        if n % 50 == 0:
            print(f"probed {n}, experiments found {len(cat)}", flush=True)
            json.dump(cat, open(out_path, "w"), indent=1)
        time.sleep(0.02)
    json.dump(cat, open(out_path, "w"), indent=1)
    print("total experiments", len(cat), "score sets", sum(len(v["score_sets"]) for v in cat.values()))


if __name__ == "__main__":
    main(sys.argv[1], *(sys.argv[2:3] or []))
