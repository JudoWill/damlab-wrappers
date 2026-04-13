"""Wrapper for detecting deletion blocks from BAM files"""

__author__ = "Will Dampier"
__copyright__ = "Copyright 2026"
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "1.3.0"

import logging
import math
import re
import sys
import yaml
import csv
from collections import Counter, defaultdict
from typing import Dict, Set, Tuple, List, Optional, Any
import cigarmath as cm  # type: ignore

if "snakemake" not in locals():
    import snakemake  # type: ignore


def _setup_query_debug_logger(sample_name: str = "") -> logging.Logger:
    """Log to the rule log file (``snakemake.log``) when present, else stderr."""
    log = logging.getLogger("cigarmath.deletion_block_detection")
    log.handlers.clear()
    log.setLevel(logging.DEBUG)
    log.propagate = False
    tag = f" sample={sample_name!r}" if sample_name else ""
    fmt = logging.Formatter(f"[deletion_block_detection{tag}] %(message)s")
    path = None
    if hasattr(snakemake, "log") and snakemake.log:
        try:
            path = snakemake.log[0]
        except (TypeError, IndexError, AttributeError):
            path = None
    if path:
        h: logging.Handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    else:
        h = logging.StreamHandler(sys.stderr)
    h.setFormatter(fmt)
    log.addHandler(h)
    return log


def _log_query_param_introspection(
    log: logging.Logger, params: Any, query_raw: Any, query_tokens: List[str]
) -> None:
    """Emit how Snakemake params resolve for ``deletion_query`` / ``query`` (debugging)."""
    log.info("wrapper version %s", __version__)
    try:
        keys = list(params.keys()) if hasattr(params, "keys") else []
    except Exception as e:
        keys = ["<keys() failed: %s>" % e]
    log.info("snakemake.params keys (%d): %s", len(keys), keys)
    gdq = getattr(params, "deletion_query", "<no such attribute>")
    gq = getattr(params, "query", "<no such attribute>")
    log.info("getattr(params, 'deletion_query'): %r (type=%s)", gdq, type(gdq).__name__)
    log.info("getattr(params, 'query'): %r (type=%s)", gq, type(gq).__name__)
    if hasattr(params, "get"):
        log.info("params.get('deletion_query'): %r", params.get("deletion_query"))
        log.info("params.get('query'): %r", params.get("query"))
    if hasattr(params, "__dict__"):
        # Namedlist stores named params on the instance; helps when .get misses
        d = getattr(params, "__dict__", {})
        log.info(
            "params.__dict__ keys overlapping query: %s",
            [k for k in d if "query" in k.lower() or "deletion" in k.lower()],
        )
    log.info(
        "resolved query_raw after unwrap: %r (type=%s)",
        query_raw,
        type(query_raw).__name__,
    )
    log.info("query_tokens (%d): %s", len(query_tokens), query_tokens)


def _log_reference_names_for_query_debug(
    log: logging.Logger, read_data: List[dict], query_regions: List[Tuple[str, str, int, int]]
) -> None:
    """If query refs may not match BAM @SQ SN, show what references appear in alignments."""
    if not read_data:
        log.info("reference_name debug: no reads in read_data")
        return
    ref_counts: Counter = Counter(
        (r.get("reference_name") or "") for r in read_data
    )
    top = ref_counts.most_common(8)
    log.info(
        "reference_name counts (top 8, empty string = unmapped/unknown): %s",
        top,
    )
    if query_regions:
        wanted = {qref for _, qref, _, _ in query_regions}
        present = set(ref_counts.keys())
        missing = wanted - present
        if missing:
            log.warning(
                "query reference(s) %s never appear in alignments; "
                "reads_covering will be 0 unless @SQ SN matches exactly (check samtools view -H).",
                sorted(missing),
            )


def parse_region(region_str: str) -> Tuple[str, Optional[int], Optional[int]]:
    """Parse region string in format 'ref:start-end' or 'ref' (same convention as cigarmath/slice).

    Reference names may include alphanumerics, underscores, dots, and hyphens (e.g. NC_045512.2).
    """
    region_str = region_str.strip()
    match = re.match(r"([\w.-]+)(?::(\d+)-(\d+))?", region_str)
    if not match:
        raise ValueError(
            f"Invalid region format: {region_str!r}. Expected 'ref:start-end' or 'ref'"
        )
    ref = match.group(1)
    start = int(match.group(2)) if match.group(2) else None
    end = int(match.group(3)) if match.group(3) else None
    return ref, start, end


def _resolve_query_raw_from_params(params: Any) -> Any:
    """Prefer ``deletion_query`` (Snakemake param name); fall back to ``query`` for older rules."""
    v = getattr(params, "deletion_query", None)
    if v is None and hasattr(params, "get"):
        v = params.get("deletion_query")
    if v is None:
        v = getattr(params, "query", None)
    if v is None and hasattr(params, "get"):
        v = params.get("query")
    return v


def _unwrap_scalar(x: Any) -> Any:
    """Extract Python scalar from numpy/pandas 0-d values so ``isinstance(..., str)`` works."""
    if x is None or isinstance(x, (str, list, tuple)):
        return x
    if isinstance(x, float) and math.isnan(x):
        return None
    if hasattr(x, "item"):
        try:
            return x.item()
        except Exception:
            return str(x)
    return x


def normalize_query_param(raw: Any) -> List[str]:
    """Turn params.query into a list of region tokens (semicolon-separated in CSV cells)."""
    if raw is None:
        return []
    if isinstance(raw, float) and math.isnan(raw):
        return []
    if isinstance(raw, str):
        s = raw.strip()
        if not s:
            return []
        parts = []
        for tok in s.split(";"):
            t = tok.strip()
            if t:
                parts.append(t)
        return parts
    if isinstance(raw, (list, tuple)):
        out: List[str] = []
        for item in raw:
            out.extend(normalize_query_param(item))
        return out
    return normalize_query_param(str(raw))


def parse_query_regions(tokens: List[str]) -> List[Tuple[str, str, int, int]]:
    """Return list of (original_token, reference, start, end). start/end required."""
    regions: List[Tuple[str, str, int, int]] = []
    for tok in tokens:
        ref, qs, qe = parse_region(tok)
        if qs is None or qe is None:
            raise ValueError(
                f"Query region {tok!r} must include coordinates (ref:start-end), not reference-only."
            )
        regions.append((tok, ref, qs, qe))
    return regions


def read_overlaps_query_interval(
    read_start: int, read_end: int, q_start: int, q_end: int
) -> bool:
    """Same overlap rule as cigarmath/slice: read interval overlaps [q_start, q_end)."""
    return read_start < q_end and read_end > q_start


def deletion_overlaps_query_interval(
    del_start: int, del_end: int, q_start: int, q_end: int
) -> bool:
    """Half-open deletion [del_start, del_end) vs query [q_start, q_end)."""
    return del_start < q_end and del_end > q_start


def parse_allowedlist(allowedlist_path: str) -> Set[Tuple[int, int]]:
    """Parse allowedlist file containing deletion coordinates.
    
    Expected format: CSV with columns 'start' and 'end' or 'deletion_start' and 'deletion_end'
    """
    allowed = set()
    with open(allowedlist_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'start' in row and 'end' in row:
                allowed.add((int(row['start']), int(row['end'])))
            elif 'deletion_start' in row and 'deletion_end' in row:
                allowed.add((int(row['deletion_start']), int(row['deletion_end'])))
    return allowed


def get_combined_segments(bam_paths: List[str]):
    """Get combined segments from one or more BAM files."""
    all_segments = []
    for bam_path in bam_paths:
        segments = list(cm.io.segment_stream_pysam(bam_path, mode='rb'))
        all_segments.extend(segments)
    
    all_segments.sort(key=lambda x: x.query_name)
    valid_segments = [s for s in all_segments if s.cigartuples]
    return list(cm.io.combined_segment_stream(iter(valid_segments)))


def merge_blocks(
    blocks: Counter,
    n: int,
    *,
    sum_counts: bool = True,
) -> Tuple[Counter, Dict[Tuple[int, int], Tuple[int, int]]]:
    """
    Merge blocks whose start and stop are within n units of each other.
    
    Args:
        blocks: Counter keyed by (start, stop) tuples.
        n: integer proximity threshold.
        sum_counts: if True, representative's count is sum of component counts.
                    if False, representative's count is the original count of the chosen block.
    
    Returns:
        Tuple of:
          - Counter of merged blocks: representative_block -> count
          - Dict mapping each original block to its representative block
    """
    if not blocks:
        return Counter(), {}

    items: List[Tuple[Tuple[int, int], int]] = list(blocks.items())
    m = len(items)

    # Union-find
    parent = list(range(m))
    rank = [0] * m

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[rb] < rank[ra]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    for i in range(m):
        (si, ei), _ = items[i]
        for j in range(i + 1, m):
            (sj, ej), _ = items[j]
            if abs(si - sj) <= n and abs(ei - ej) <= n:
                union(i, j)

    groups: Dict[int, List[int]] = defaultdict(list)
    for idx in range(m):
        groups[find(idx)].append(idx)

    merged: Counter = Counter()
    remap: Dict[Tuple[int, int], Tuple[int, int]] = {}
    for group_idxs in groups.values():
        rep_idx = max(
            group_idxs,
            key=lambda idx: (items[idx][1], -items[idx][0][0], -items[idx][0][1]),
        )
        rep_key = items[rep_idx][0]
        if sum_counts:
            total = sum(items[idx][1] for idx in group_idxs)
            merged[rep_key] += total
        else:
            merged[rep_key] += items[rep_idx][1]
        for idx in group_idxs:
            remap[items[idx][0]] = rep_key

    return merged, remap


def calculate_shannon_entropy(counter: Counter) -> float:
    """Calculate Shannon entropy of a Counter's value distribution."""
    total = sum(counter.values())
    if total == 0:
        return 0.0
    entropy = 0.0
    for count in counter.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log(p)
    return entropy


# Get input/output files
input_bams = snakemake.input.bams if hasattr(snakemake.input, 'bams') else list(snakemake.input)
output_reads_csv = snakemake.output.reads
output_deletions_csv = snakemake.output.deletions
output_summary_yaml = snakemake.output.summary
output_query_stats = snakemake.output.query_stats

# Get parameters
min_deletion_size = snakemake.params.get("min_deletion_size", 50)
merge_distance = snakemake.params.get("merge_distance", 0)
sample_name = snakemake.params.get("sample_name", "sample")
_debug_deletion_query = snakemake.params.get("debug_deletion_query", True)

_dbg = _setup_query_debug_logger(sample_name)
query_raw = _unwrap_scalar(_resolve_query_raw_from_params(snakemake.params))
query_tokens = normalize_query_param(query_raw)

_dbg.info("input_bams=%r", input_bams)
if _debug_deletion_query:
    _log_query_param_introspection(_dbg, snakemake.params, query_raw, query_tokens)

# Get optional allowedlist
allowedlist_input = snakemake.input.get("allowedlist", None)
allowed_deletions: Optional[Set[Tuple[int, int]]] = None
if allowedlist_input:
    allowed_deletions = parse_allowedlist(allowedlist_input)

# Ensure input_bams is a list
if isinstance(input_bams, str):
    input_bams = [input_bams]

# Process all BAM files together
combined_stream = get_combined_segments(input_bams)

# Data collection
read_data: List[dict] = []
deletion_counter: Counter = Counter()
total_reads = 0
reads_with_deletions = 0

# First pass: collect all deletions and read data
for start, cigartuples, segments in combined_stream:
    if not cigartuples:
        continue
    
    total_reads += 1
    read_name = segments[0].query_name
    reference_name = segments[0].reference_name or ""
    ref_end = start + cm.reference_offset(cigartuples)
    
    deletions = list(cm.reference_deletion_blocks(
        cigartuples, 
        reference_start=start, 
        min_size=min_deletion_size
    ))
    
    # Filter by allowedlist if specified
    if allowed_deletions is not None:
        deletions = [d for d in deletions if d in allowed_deletions]
    
    if deletions:
        reads_with_deletions += 1
    
    read_data.append({
        'read_name': read_name,
        'reference_name': reference_name,
        'reference_start': start,
        'reference_end': ref_end,
        'deletions': ';'.join(f"{d[0]}-{d[1]}" for d in deletions) if deletions else ''
    })
    
    for del_block in deletions:
        deletion_counter[del_block] += 1

# Optionally merge deletion blocks with similar endpoints
if merge_distance > 0:
    deletion_counter, block_remap = merge_blocks(deletion_counter, merge_distance)

    # Update deletion coordinates in read_data to use representative blocks
    for read in read_data:
        if read['deletions']:
            old_dels: List[Tuple[int, int]] = [
                tuple(map(int, d.split('-'))) for d in read['deletions'].split(';')  # type: ignore
            ]
            new_dels = [block_remap.get(d, d) for d in old_dels]
            # Deduplicate while preserving order (two old blocks may collapse to the same rep)
            seen_dels: Set[Tuple[int, int]] = set()
            unique_dels = []
            for d in new_dels:
                if d not in seen_dels:
                    seen_dels.add(d)
                    unique_dels.append(d)
            read['deletions'] = ';'.join(f"{d[0]}-{d[1]}" for d in unique_dels)

# Second pass: calculate coverage for each deletion
coverage_counter: Counter = Counter()
for read in read_data:
    read_start = read['reference_start']
    read_end = read['reference_end']
    for del_block in deletion_counter:
        if read_start <= del_block[0] and read_end >= del_block[1]:
            coverage_counter[del_block] += 1

query_regions_parsed: List[Tuple[str, str, int, int]] = []
if query_tokens:
    try:
        query_regions_parsed = parse_query_regions(query_tokens)
    except ValueError as err:
        _dbg.error("parse_query_regions failed (writing header-only query_stats): %s", err)
        query_regions_parsed = []

if _debug_deletion_query and query_regions_parsed:
    _log_reference_names_for_query_debug(_dbg, read_data, query_regions_parsed)

# Per-query-region stats (same rows written to query_stats.csv; also embedded in summary YAML)
query_stats_rows: List[dict] = []
if query_regions_parsed:
    for token, qref, qs, qe in query_regions_parsed:
        reads_covering = 0
        reads_del_ovl = 0
        for read in read_data:
            if (read.get('reference_name') or '') != qref:
                continue
            rs = read['reference_start']
            re = read['reference_end']
            if not read_overlaps_query_interval(rs, re, qs, qe):
                continue
            reads_covering += 1
            del_str = read.get('deletions') or ''
            dels: List[Tuple[int, int]] = []
            if del_str.strip():
                dels = [
                    tuple(map(int, x.split('-')))  # type: ignore[misc]
                    for x in del_str.split(';')
                ]
            if any(
                deletion_overlaps_query_interval(d0, d1, qs, qe) for d0, d1 in dels
            ):
                reads_del_ovl += 1
        row = {
            'region': token,
            'reference': qref,
            'start': qs,
            'end': qe,
            'reads_covering': reads_covering,
            'reads_with_deletion_overlapping': reads_del_ovl,
        }
        query_stats_rows.append(row)
        if _debug_deletion_query:
            _dbg.info("query_stats %s", row)
elif query_tokens and not query_regions_parsed:
    _dbg.warning(
        "query_tokens non-empty but query_regions_parsed empty (parse error above); CSV is header-only"
    )
elif not query_tokens and _debug_deletion_query:
    _dbg.info("no query tokens; query_stats CSV is header-only")

# Write read-centered CSV (exclude internal reference_name; not part of public schema)
_read_csv_fields = ['read_name', 'reference_start', 'reference_end', 'deletions']
with open(output_reads_csv, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=_read_csv_fields)
    writer.writeheader()
    writer.writerows({k: r[k] for k in _read_csv_fields} for r in read_data)

# Write deletion-centered CSV
deletion_data = []
for del_block, read_count in deletion_counter.most_common():
    deletion_data.append({
        'deletion_start': del_block[0],
        'deletion_end': del_block[1],
        'deletion_size': del_block[1] - del_block[0],
        'read_count': read_count,
        'coverage_count': coverage_counter[del_block]
    })

with open(output_deletions_csv, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['deletion_start', 'deletion_end', 'deletion_size', 'read_count', 'coverage_count'])
    writer.writeheader()
    writer.writerows(deletion_data)

# Calculate summary statistics
total_deletion_count = sum(deletion_counter.values())
richness = len(deletion_counter)
shannon_entropy = calculate_shannon_entropy(deletion_counter)

targeted_regions_yaml: List[dict] = []
for r in query_stats_rows:
    targeted_regions_yaml.append(
        {
            'region': r['region'],
            'reference': r['reference'],
            'start': int(r['start']),
            'end': int(r['end']),
            'reads_covering': int(r['reads_covering']),
            'reads_with_deletion_overlapping': int(r['reads_with_deletion_overlapping']),
        }
    )
target_cov_sum = sum(int(r['reads_covering']) for r in query_stats_rows)
target_ovl_sum = sum(int(r['reads_with_deletion_overlapping']) for r in query_stats_rows)

top_deletions: List[dict] = []
for row in deletion_data[:10]:
    top_deletions.append(
        {
            'deletion_start': int(row['deletion_start']),
            'deletion_end': int(row['deletion_end']),
            'deletion_size': int(row['deletion_size']),
            'read_count': int(row['read_count']),
            'coverage_count': int(row['coverage_count']),
        }
    )

# Write summary YAML for MultiQC
summary = {
    'sample_name': sample_name,
    'total_reads': total_reads,
    'reads_with_deletions': reads_with_deletions,
    'unique_deletion_count': richness,
    'total_deletion_count': total_deletion_count,
    'deletion_frequency': reads_with_deletions / total_reads if total_reads > 0 else 0.0,
    'deletion_richness': richness,
    'deletion_shannon_entropy': round(shannon_entropy, 6),
    'min_deletion_size': min_deletion_size,
    'merge_distance': merge_distance,
    'input_bam_count': len(input_bams),
    'allowedlist_used': allowed_deletions is not None,
    'allowedlist_size': len(allowed_deletions) if allowed_deletions else 0,
    'targeted_regions': targeted_regions_yaml,
    'targeted_region_count': len(targeted_regions_yaml),
    'target_reads_covering_sum': target_cov_sum,
    'target_reads_with_deletion_overlapping_sum': target_ovl_sum,
    'top_deletions': top_deletions,
}

with open(output_summary_yaml, 'w') as f:
    f.write('# Cigarmath Deletion Block Detection\n')
    yaml.dump(summary, f, default_flow_style=False)

# Per-query-region coverage vs deletion overlap (optional params.query)
query_fieldnames = [
    'region',
    'reference',
    'start',
    'end',
    'reads_covering',
    'reads_with_deletion_overlapping',
]
with open(output_query_stats, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=query_fieldnames)
    writer.writeheader()
    for row in query_stats_rows:
        writer.writerow(row)
