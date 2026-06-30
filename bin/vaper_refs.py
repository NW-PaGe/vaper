#!/usr/bin/env python3

# vaper_refs.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import csv
import gzip
import json
import os
import subprocess
from pathlib import Path
from collections import defaultdict
from typing import Any, Dict, List, Tuple, Set

import numpy as np
import sourmash
from sklearn.cluster import DBSCAN

from vaper_utils import logging_config, get_ref_name

LOGGER = logging_config()


def _create_ref_map(data: list[dict[str, Any]]):
    """Create ID and metadata maps for quick lookups."""
    LOGGER.info("Creating reference map for quick look ups")
    id_map: dict[str, Any] = {}
    meta_map: dict[str, dict[str, list[str]]] = {}
    
    name_cache = {}
    for rec in data:
        if not isinstance(rec, dict):
            continue

        seq = rec.get("sequence")
        if seq is None:
            continue

        # format name (& create if needed)
        name = get_ref_name(rec, name_cache)
        rec['name'] = name
        id_map[name] = rec

        metadata = rec.get("metadata")
        if metadata is None:
            continue

        for k, v in metadata.items():
            if v is None:
                continue

            # Normalize v to an iterable of strings
            if isinstance(v, str):
                vals = [v]
            elif isinstance(v, (list, tuple, set)):
                vals = [str(x) for x in v if x is not None]
            else:
                vals = [str(v)]

            key_lower = k.lower()
            bucket = meta_map.setdefault(key_lower, {})

            for val in vals:
                val = val.strip()
                if not val:
                    continue
                bucket.setdefault(val, []).append(name)

    return id_map, meta_map


def _get_exceptions(meta_map, pattern_str: str | None) -> tuple[list[str], list[str]]:
    """Parse include/exclude patterns and return matching reference names."""
    exceptions_by_name: list[str] = []
    exceptions_by_meta: list[str] = []

    if pattern_str is None:
        return exceptions_by_name, exceptions_by_meta

    parts = pattern_str.split(",")

    for p in parts:
        p = p.strip()
        if not p:
            continue
        if p.count("=") != 1:
            LOGGER.warning(f"Malformed include / exclude pattern: {p}")
            continue
        key, value = (x.strip() for x in p.split("=", 1))
        if key == "name":
            exceptions_by_name.append(value)
        else:
            names = meta_map.get(key, {}).get(value, [])
            exceptions_by_meta.extend(names)
    
    return list(set(exceptions_by_name)), list(set(exceptions_by_meta))


def _check_file_ext(filepath):
    """Ensure the input file is JSONL or JSONL.GZ."""
    LOGGER.debug(f"Checking file extension for: {filepath}")
    if not (filepath.endswith("jsonl.gz") or filepath.endswith("jsonl")):
        raise ValueError("Inputs must be JSONL format.")
    LOGGER.info(f"Checked extension: {os.path.basename(filepath)}")


def _load_jsonl(filepath):
    """Load JSON Lines (optionally gzipped) into a list of dicts."""
    LOGGER.debug(f"Loading JSONL: {filepath}")

    def _read_lines(file_obj):
        data = []
        n = 0
        for line in file_obj:
            n += 1
            line = line.strip()
            if line:
                data.append(json.loads(line))
        LOGGER.info(f"Parsed {len(data)} non-empty JSONL records from {os.path.basename(filepath)} (raw lines: {n})")
        return data

    opener = gzip.open if filepath.endswith("gz") else open
    with opener(filepath, "rt" if filepath.endswith("gz") else "r", encoding="utf-8") as f:
        data = _read_lines(f)

    LOGGER.info(f"Loaded: {os.path.basename(filepath)} (records={len(data)})")
    return data


def _write_fasta(records, outpath):
    """Write FASTA records to file (optionally gzipped)."""
    opener = gzip.open if outpath.endswith(".gz") else open
    
    with opener(outpath, "wt") as f:
        # Write in chunks to balance memory vs I/O
        chunk_size = 1000
        chunk = []
        
        for name, seq in records:
            chunk.append(f">{name}\n{seq}\n")
            
            if len(chunk) >= chunk_size:
                f.write("".join(chunk))
                chunk = []
        
        # Write remaining
        if chunk:
            f.write("".join(chunk))


def _dict_to_fasta(id_map, filename, exclusions: list = []):
    """Export reference sequences to gzipped FASTA, excluding specified names."""
    outdir = Path(filename).parent
    if str(outdir):
        outdir.mkdir(parents=True, exist_ok=True)

    exclusion_set = set(exclusions) if exclusions else set()
    subset = [(k, v['sequence']) for k, v in id_map.items() if k not in exclusion_set]

    _write_fasta(subset, filename)

    return filename


def _export_jsonl(data, filename):
    """Export list of dicts as gzipped JSONL."""
    LOGGER.debug(f"Exporting JSONL to: {filename} (records={len(data)})")
    outdir = Path(filename).parent
    if str(outdir):
        outdir.mkdir(parents=True, exist_ok=True)        

    with gzip.open(filename, "wt", encoding="utf-8") as f:
        f.write("\n".join(json.dumps(r) for r in data) + "\n")

    LOGGER.info(f"JSONL exported to {filename} (records={len(data)})")
    return filename


def _run_minimap2(reference, query, output="map.paf", preset="asm5", secondary="yes"):
    """Run minimap2 alignment and save PAF output."""
    outdir = Path(output).parent
    if str(outdir):
        outdir.mkdir(parents=True, exist_ok=True)

    threads = os.cpu_count() or 1
    cmd = [
        "minimap2", "-x", preset, "-t", str(threads),
        "-g", "300", "-r", "100,300", "--secondary", secondary, 
        reference, query
    ]

    stderr_file = str(Path(output).with_suffix(".log"))
    with open(output, "w") as out_f, open(stderr_file, "w") as err_f:
        subprocess.run(cmd, stdout=out_f, stderr=err_f, check=True)

    LOGGER.info(f"minimap2 finished successfully (PAF: {output}, LOG: {stderr_file})")
    return output


def _read_paf(path):
    """
    Read a PAF file and return alignment records.
    PAF columns used (0-based): 0 qname, 5 tname, 6 tlen, 7 tstart, 8 tend, 10 alnlen, 11 qual
    """
    LOGGER.debug(f"Reading PAF: {path}")
    rows = []
    with open(path, "r", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        for cols in r:
            if len(cols) < 12:
                continue
            rows.append({
                "query": cols[0],
                "target": cols[5],
                "qlen": int(cols[1]),
                "qstart": int(cols[2]),
                "qend": int(cols[3]),
                "tlen": int(cols[6]),
                "tstart": int(cols[7]),
                "tend": int(cols[8]),
                "matches": int(cols[9]),
                "align": int(cols[10]),
                "qual": int(cols[11])
            })
    
    if rows:
        uniq_targets = len({x["target"] for x in rows})
        LOGGER.info(f"PAF parsed: records={len(rows)}, unique_targets={uniq_targets}, file={os.path.basename(path)}")
    else:
        LOGGER.warning(f"PAF appears empty: {path}")
    return rows


def _create_subset(id_map, names, prefix="subset"):
    """Write subset of records by name to JSONL.GZ and FA.GZ."""
    LOGGER.debug(f"Creating subset for {len(names)} names with prefix '{prefix}'")
    subset = {n: id_map[n] for n in names if n in id_map}
    missing = set(names) - set(list(subset.keys()))
    
    if missing:
        LOGGER.warning(f"{len(missing)} requested names not found in data: {sorted(missing)[:5]}{' ...' if len(missing) > 5 else ''}")

    LOGGER.info(f'Exporting references: {list(subset.keys())}')

    jsonl_path = f"{prefix}.jsonl.gz"
    fasta_path = f"{prefix}.fa.gz"

    _export_jsonl([{k: v for k, v in rec.items() if k != "sequence"} for rec in subset.values()], jsonl_path)
    _dict_to_fasta(subset, fasta_path)

    LOGGER.info(f"Subset exports complete: JSONL={jsonl_path} (records={len(subset)}), FASTA={fasta_path}")


def _compare_refs(data, names, ksize, scaled, threshold, outdir):
    """
    Build MinHashes for selected names and cluster with DBSCAN (min_samples=1).
    Saves grouping results to <outdir>/reference_groups.json.
    Returns list of sets (groups).
    """
    LOGGER.info(f"Comparing refs (names={len(names)}, ksize={ksize}, scaled={scaled}, threshold={threshold})")

    # Ensure outdir exists
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / "reference_groups.json"

    # Build MinHashes for requested names
    subset = {}
    for rec in data:
        nm = rec.get("name")
        if nm in names:
            mh = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
            mh.add_sequence(rec["sequence"], force=True)
            subset[nm] = mh

    LOGGER.info(f"Built MinHashes: {len(subset)} of {len(names)} requested names present")

    if not subset:
        LOGGER.info("Reference grouping complete: groups=0, pairwise=0")

        # Save empty file
        with outfile.open("w", encoding="utf-8") as f:
            json.dump([], f, indent=2)

        LOGGER.info(f"Wrote reference grouping results to {outfile}")
        return []

    ids = list(subset.keys())
    n = len(ids)

    # Precompute pairwise distances for DBSCAN
    mat = np.zeros((n, n), dtype=float)
    pairwise = 0
    for i in range(n):
        m1 = subset[ids[i]]
        for j in range(i + 1, n):
            d = m1.containment_ani(subset[ids[j]]).dist
            mat[i, j] = mat[j, i] = d
            pairwise += 1

    # Cluster with DBSCAN
    db = DBSCAN(eps=threshold, min_samples=1, metric="precomputed").fit(mat)
    labels = db.labels_

    # Build groups from labels
    clusters = {}
    noise = []
    for idx, lbl in enumerate(labels):
        if lbl == -1:
            noise.append({ids[idx]})
        else:
            clusters.setdefault(int(lbl), set()).add(ids[idx])

    groups = list(clusters.values()) + noise

    # Log summary
    if groups:
        sizes = sorted((len(g) for g in groups), reverse=True)
        LOGGER.info(f"Reference grouping complete: groups={len(groups)}, pairwise={pairwise}, largest_group={sizes[0]}")
    else:
        LOGGER.info(f"Reference grouping complete: groups=0, pairwise={pairwise}")

    # Save groups to JSON
    serializable = [sorted(list(g)) for g in groups]
    with outfile.open("w", encoding="utf-8") as f:
        json.dump(serializable, f, indent=2)

    LOGGER.info(f"Wrote reference grouping results to {outfile}")

    return groups


def _select_refs(groups, paf_con, threshold):
    """
    Select references per group:
      1) If any member has gf >= threshold -> choose max by alignment length * ~ pid.
      2) Else if sum(group gfs) >= threshold -> choose max by alignment length * ~ pid.
      3) Else choose none.
    Ties broken by reference name in alphabetical order.
    """
    LOGGER.debug(f"Selecting refs from {len(groups)} groups with GF threshold={threshold}")
    selected = []

    def gf_of(name):
        return float(paf_con.get(name, {}).get("gf", 0.0))

    def aln_of(name):
        return int(paf_con.get(name, {}).get("aligned", 0))

    def pid_of(name):
        return float(paf_con.get(name, {}).get("pid", 0))

    def score(n):
        return aln_of(n) * pid_of(n)

    def best(names):
        # Highest score wins; ties broken by alphabetically first name.
        return sorted(names, key=lambda n: (-score(n), n))[0]

    for g in groups:
        names = list(g)

        if sum(gf_of(n) for n in names) >= threshold and names:
            selected.append(best(names))
            continue

    LOGGER.info(f"Selection complete: selected={len(selected)} (threshold={threshold})")
    return selected

def _consolidate_paf(data, outdir, out_paf="chosen.paf", target_summary="target-summary.csv"):
    """
    Calculate coverage statistics for reference targets from PAF alignments.
    Also exports the selected (chosen) alignments to outdir/outfile as PAF.
    """
    LOGGER.debug(f"Consolidating PAF records (n={len(data)})")

    # Ensure work dir exists early (so we can export chosen alignments)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ----------------------------
    # 1) Group by query
    # ----------------------------
    by_query = defaultdict(list)
    for rec in data:
        by_query[rec["query"]].append(rec)

    chosen = []

    def _aln_len(rec):
        if "align" in rec and rec["align"] is not None:
            return int(rec["align"])
        s, e = rec["tstart"], rec["tend"]
        return abs(int(e) - int(s))

    # 1.1) Choose the best alignment for each query (one per subject)
    for q, recs in by_query.items():
        recs_sorted = sorted(
            recs,
            key=lambda r: (
                _aln_len(r),
                int(r.get("qual", 0)),
                str(r.get("target"))
            ),
            reverse=True,
        )

        used_subjects = set()
        for r in recs_sorted:
            subject = r["target"]
            if subject not in used_subjects:
                chosen.append(r)
                used_subjects.add(subject)

    LOGGER.info(
        f"Selected best alignments: {len(chosen)} "
        f"(queries={len(by_query)}, subject-region restriction)"
    )

    # ----------------------------
    # 1.2) Export chosen alignments to work directory (PAF)
    # ----------------------------
    def _tags_to_str(rec):
        tags = rec.get("tags")
        if not tags:
            return ""
        # allow tags as dict -> "XX:Z:val" etc if already formatted, else coerce to Z
        if isinstance(tags, dict):
            parts = []
            for k, v in tags.items():
                if isinstance(v, tuple) and len(v) == 2:
                    ttype, tval = v
                    parts.append(f"{k}:{ttype}:{tval}")
                else:
                    parts.append(f"{k}:Z:{v}")
            return "\t" + "\t".join(parts) if parts else ""
        if isinstance(tags, (list, tuple)):
            parts = [str(x) for x in tags if x is not None and str(x) != ""]
            return "\t" + "\t".join(parts) if parts else ""
        return "\t" + str(tags)

    def _paf_line(rec):
        # Required 12 PAF columns
        qname = rec.get("query", "*")
        qlen  = int(rec.get("qlen", 0))
        qst   = int(rec.get("qstart", 0))
        qen   = int(rec.get("qend", 0))
        strand = rec.get("strand", "+")  # if you store it; else '+'
        tname = rec.get("target", "*")
        tlen  = int(rec.get("tlen", 0))
        tst   = int(rec.get("tstart", 0))
        ten   = int(rec.get("tend", 0))
        nmatch = int(rec.get("match", rec.get("nmatch", 0)))
        alnlen = int(rec.get("align", rec.get("alnlen", 0)))
        mapq   = int(rec.get("qual", rec.get("mapq", 0)))

        core = [qname, qlen, qst, qen, strand, tname, tlen, tst, ten, nmatch, alnlen, mapq]
        return "\t".join(map(str, core)) + _tags_to_str(rec) + "\n"

    chosen_path = outdir / out_paf
    with chosen_path.open("w", newline="") as fh:
        for rec in chosen:
            fh.write(_paf_line(rec))
    LOGGER.info(f"Wrote selected alignments to {chosen_path}")

    # ----------------------------
    # 2) Compute coverage using only chosen alignments
    # ----------------------------
    out = {}
    cov = {}
    pid = {}

    for rec in chosen:
        q = rec["query"]
        t = rec["target"]
        L = rec["tlen"]
        s, e = rec["tstart"], rec["tend"]
        if e < s:
            s, e = e, s

        if t not in cov:
            cov[t] = np.zeros(L, dtype=int)
            out[t] = {"hits": [q]}
            LOGGER.debug(f"New target: {t} (length={L})")

        if q not in out[t]["hits"]:
            out[t]["hits"].append(q)

        s = max(0, min(int(s), cov[t].size))
        e = max(s, min(int(e), cov[t].size))
        cov[t][s:e] = True

        aligned = rec['align'  ]  
        matches = rec['matches']

        if t not in pid:
            pid[t] = {'matches': matches, 'aligned': aligned}
        else:
            pid[t]['matches'] += matches
            pid[t]['aligned'] += aligned

    LOGGER.debug(f"Computing coverage and pid for {len(cov)} targets")

    for t, v in cov.items():
        length = v.size
        cov_val = int(v.sum())
        gf = cov_val / length if length else 0.0
        pid_val = pid[t]["matches"] / pid[t]["aligned"] if pid[t]["aligned"] > 0 else 0
        out[t]["length"] = float(length)
        out[t]["aligned"] = float(cov_val)
        out[t]["gf"] = float(gf)
        out[t]["pid"] = float(pid_val)
        LOGGER.debug(f"Target {t}: coverage={gf:.3f}, pid~={pid_val:.3f}, hits={len(out[t]['hits'])}")

    summary_path = outdir / target_summary
    with summary_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["target", "length", "aligned", "gf", "pid", "num_hits"])
        for target in sorted(out.keys(), key=lambda t: out[t]["gf"], reverse=True):
            writer.writerow([
                target,
                out[target]["length"],
                out[target]["aligned"],
                out[target]["gf"],
                out[target]["pid"],
                len(out[target]["hits"]),
            ])

    LOGGER.info(f"Wrote target summary to {summary_path}")
    return out

def _validate(data: List[dict]) -> None:
    """Validate reference records for required fields, bases, and name uniqueness."""
    LOGGER.debug(f"Validating {len(data)} reference records")    

    for i, rec in enumerate(data):
        if not isinstance(rec, dict):
            raise ValueError(f"Input is malformed on line {i + 1}")

        req_fields = ["sequence"]
        if any(k not in rec for k in req_fields):
            raise ValueError(f"Line {i + 1} is missing one or more required field: {req_fields}")

        LEGAL_BASES = {"-", "A", "T", "C", "G", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"}
        for p, b in enumerate(rec["sequence"]):
            if b.upper() not in LEGAL_BASES:
                raise ValueError(
                    f"Illegal base ({b}) at position {p} in reference (line {i})"
                )

    LOGGER.info("Input passed validation")


def main():
    """
    Command-line entry point for VAPER reference formatting.
    Processes reference JSONL files and optionally maps a query assembly.
    """
    version = "1.0.1"

    parser = argparse.ArgumentParser(
        description="VAPER reference processing and selection tool"
    )
    parser.add_argument("--refs", required=True, help="Path to reference JSONL file.")
    parser.add_argument("--query", help="Path to query assembly.")
    parser.add_argument("--genfrac", type=float, default=0.50, help="Genome fraction threshold to select a reference.")
    parser.add_argument("--dist", type=float, default=0.20, help="Distance threshold used to cluster references (1 - ANI/100).")
    parser.add_argument("--scaled", type=int, default=1, help="MinHash scaled factor.")
    parser.add_argument("--ksize", type=int, default=31, help="MinHash k-mer size.")
    parser.add_argument("--include", help="Comma separated list of references to include (name=value or key=value).")
    parser.add_argument("--exclude", help="Comma separated list of references to exclude (name=value or key=value).")
    parser.add_argument("--outdir", default='.', help="Output directory")
    parser.add_argument("--validate", action="store_true", help="Validate the JSONL file.")
    parser.add_argument("--version", action="version", version=version)
    args = parser.parse_args()

    LOGGER.info(f"vaper_refs v{version}")
    LOGGER.info("Author: Jared Johnson")

    _check_file_ext(args.refs)
    data = _load_jsonl(args.refs)

    if args.validate:
        _validate(data)
        return

    if not args.query:
        raise ValueError("No query supplied!")

    id_map, meta_map = _create_ref_map(data)

    include_names, include_meta = _get_exceptions(meta_map, args.include)

    # Handle exclude = "*"
    if args.exclude == "*":
        exclude_names = list(id_map.keys())  # everything
        exclude_meta = []
    else:
        exclude_names, exclude_meta = _get_exceptions(meta_map, args.exclude)

    # Build include/exclude sets
    include = set(include_names) | set(include_meta)
    exclude = set(exclude_names) | set(exclude_meta)

    # Remove anything that is explicitly included
    exclude = list(exclude - include)

    if include or exclude:
        LOGGER.info(f"Manual changes to references: excluding={len(exclude)}, including={len(include)}")

    work_dir = os.path.join(args.outdir, 'work')
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(args.outdir, exist_ok=True)
    LOGGER.info(f"Work directory: {work_dir}")
    LOGGER.info(f"Output directory: {args.outdir}")

    ref = _dict_to_fasta(id_map, os.path.join(work_dir, "ref.fa"), exclude)
    paf_file = _run_minimap2(ref, args.query, os.path.join(work_dir, "map.paf"))
    paf_data = _read_paf(paf_file)
    paf_data_con = _consolidate_paf(paf_data, work_dir)

    targets = sorted(list(paf_data_con.keys()))
    LOGGER.info(f"Targets detected from PAF: {len(targets)}")

    groups = _compare_refs(data, targets, args.ksize, args.scaled, args.dist, work_dir)
    selected = _select_refs(groups, paf_data_con, args.genfrac)

    if include_names:
        before = len(selected)
        selected = list(set(selected + include_names))
        LOGGER.info(f"Adding references specified by name: before={before}, added={len(include_names)}, after={len(selected)}")
        
    if selected:
        _create_subset(id_map, selected)
    else:
        LOGGER.info("No references met selection criteria.")


if __name__ == "__main__":
    main()