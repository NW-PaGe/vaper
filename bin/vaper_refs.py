#!/usr/bin/env python3

# vaper_refs.py
# Author: Jared Johnson, jared.johnson@doh.wa.gov

import argparse
import csv
import gzip
import json
import os
import random
import subprocess
from pathlib import Path
from collections import defaultdict
from typing import Any, List
from itertools import combinations

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

def _match_rate(rec):
    """Matches per aligned base for a single alignment record."""
    aln = rec.get("align", 0)
    return rec.get("matches", 0) / aln if aln else 0.0


def estimate_target_matches(recs):
    """
    Estimate total matches for one reference while removing duplicated
    contributions from overlapping contigs.
    """
    queries = list(recs.values())
    total = sum(r.get("matches", 0) for r in queries)

    rates = {id(r): _match_rate(r) for r in queries}

    bounds = sorted(
        {r["tstart"] for r in queries} |
        {r["tend"] for r in queries}
    )

    for a, b in zip(bounds, bounds[1:]):
        seg = b - a
        if seg <= 0:
            continue

        covering = [
            r for r in queries
            if r["tstart"] <= a and r["tend"] >= b
        ]

        if len(covering) < 2:
            continue

        for r in covering:
            total -= rates[id(r)] * seg

        total += max(rates[id(r)] for r in covering) * seg

    return total


def build_maps(paf, min_qcov=0.50):
    """
    Build contig->references and reference->contigs maps.

    Alignments whose query coverage (matches / query length) is below
    ``min_qcov`` are discarded before the maps are built.
    """

    contig_to_refs = defaultdict(set)
    ref_to_contigs = defaultdict(set)
    ref_records = defaultdict(dict)

    kept = 0
    dropped = 0

    for rec in paf:

        q = rec["query"]
        t = rec["target"]

        qlen = rec.get("qlen", 0)
        qcov = rec["matches"] / qlen if qlen else 0.0

        if qcov < min_qcov:
            dropped += 1
            LOGGER.debug(
                f"Dropped alignment {q} -> {t}: query coverage "
                f"{qcov:.3f} < min_qcov={min_qcov}"
            )
            continue

        kept += 1
        contig_to_refs[q].add(t)
        ref_to_contigs[t].add(q)

        ref_records[t][q] = rec

    LOGGER.info(
        f"Query-coverage filter (min_qcov={min_qcov}): kept {kept} / {len(paf)} "
        f"alignments (dropped {dropped}); "
        f"contigs={len(contig_to_refs)}, references={len(ref_to_contigs)}"
    )

    return contig_to_refs, ref_to_contigs, ref_records


def collapse_reference_groups(ref_to_contigs):
    """
    Collapse references that explain identical sets of contigs.

    Returns

    group_members:
        group_id -> set(reference)

    group_cover:
        group_id -> set(contigs)
    """

    coverage_map = defaultdict(set)

    for ref, contigs in ref_to_contigs.items():
        coverage_map[frozenset(contigs)].add(ref)

    group_members = {}
    group_cover = {}

    for i, (cover, refs) in enumerate(coverage_map.items()):
        gid = f"group{i+1}"
        group_members[gid] = refs
        group_cover[gid] = set(cover)

    return group_members, group_cover


def minimum_cover_groups(group_cover, all_contigs):
    """
    Find every minimum-cardinality group combination that covers all contigs.
    """

    groups = list(group_cover)

    for k in range(1, len(groups) + 1):

        solutions = []

        for combo in combinations(groups, k):

            covered = set()

            for g in combo:
                covered |= group_cover[g]

            if covered >= all_contigs:
                solutions.append(combo)

        if solutions:
            return solutions

    return []


def _reference_coverage(ref, ref_records):
    """
    Estimate coverage statistics for a single reference.

    Returns a dict with est_matches (de-duplicated across overlapping contigs),
    ref_len, explained_fraction (est_matches / ref_len), and num_contigs.
    """
    recs = ref_records[ref]
    est_matches = estimate_target_matches(recs)
    ref_len = max(r["tlen"] for r in recs.values())
    frac = est_matches / ref_len if ref_len else 0.0
    return {
        "est_matches": est_matches,
        "ref_len": ref_len,
        "explained_fraction": frac,
        "num_contigs": len(recs),
    }


def _write_selection_metrics(
    scored,
    outdir,
    ranking_file="reference-selection.csv",
    selected_file="selected-references.csv",
):
    """
    Write reference-selection metrics to the work directory for troubleshooting.

    ``scored`` is the ranked list of group-solution metric dicts (best first;
    the selected solution is placed first and flagged via its ``selected`` key).

    - ranking_file:  one row per group solution, in ranked order, with the total
      coverage used to choose between them plus supporting stats (number of
      references, summed estimated matches, and the minimum / average explained
      fraction).
    - selected_file: the per-reference breakdown of the winning set.

    An empty ``scored`` still writes headers so downstream steps always find the
    files.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ranking_path = outdir / ranking_file
    with ranking_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "rank", "selected", "num_refs", "total_coverage", "summed_est_matches",
            "min_explained_fraction", "avg_explained_fraction", "references",
        ])
        for rank, m in enumerate(scored, start=1):
            writer.writerow([
                rank,
                "yes" if m.get("selected", rank == 1) else "no",
                m["num_refs"],
                f"{m['total_coverage']:.4f}",
                f"{m['summed_est_matches']:.4f}",
                f"{m['min_explained_fraction']:.4f}",
                f"{m['avg_explained_fraction']:.4f}",
                ";".join(sorted(m["references"])),
            ])
    LOGGER.info(
        f"Wrote reference-selection ranking to {ranking_path} "
        f"(solutions={len(scored)})"
    )

    selected_path = outdir / selected_file
    with selected_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "reference", "est_matches", "ref_len",
            "explained_fraction", "num_contigs",
        ])
        if scored:
            best = scored[0]
            for ref in sorted(best["references"]):
                d = best["per_ref"][ref]
                writer.writerow([
                    ref,
                    f"{d['est_matches']:.4f}",
                    d["ref_len"],
                    f"{d['explained_fraction']:.4f}",
                    d["num_contigs"],
                ])
    LOGGER.info(f"Wrote selected-reference metrics to {selected_path}")


def choose_best_reference_set(paf, min_query_cov=0.50, min_ref_cov=0.70, outdir=None):
    """
    Select the best set of references explaining the query contigs.

    The candidate reference sets are the minimum-cardinality group combinations
    that cover every contig. Each reference belonging to a group solution has its
    coverage estimated exactly once; within each group the best-covered member
    (meeting ``min_ref_cov``) is used as the representative. The total coverage of
    a group solution is the summed explained fraction of its representatives, and
    the solution with the highest total coverage is selected. Ties are broken at
    random.

    When ``outdir`` is provided, the ranked table of group solutions and the
    winning set's per-reference breakdown are written there for troubleshooting.

    Returns a list of (total_coverage, references) tuples with the selected
    solution first, or an empty list if no complete, qualifying reference set
    exists.
    """
    LOGGER.info(
        f"Selecting reference set: alignments={len(paf)}, "
        f"min_query_cov={min_query_cov}, min_ref_cov={min_ref_cov}"
    )

    (
        contig_to_refs,
        ref_to_contigs,
        ref_records,
    ) = build_maps(paf, min_qcov=min_query_cov)

    all_contigs = set(contig_to_refs)

    if not all_contigs:
        LOGGER.warning(
            f"No contigs passed the query-coverage filter (min_query_cov="
            f"{min_query_cov}); no reference set can be formed. Consider "
            f"lowering --min-query-cov or checking the query assembly."
        )
        if outdir is not None:
            _write_selection_metrics([], outdir)
        return []

    group_members, group_cover = collapse_reference_groups(
        ref_to_contigs
    )
    LOGGER.info(
        f"Collapsed {len(ref_to_contigs)} reference(s) into "
        f"{len(group_members)} group(s) of identical contig coverage"
    )

    group_solutions = minimum_cover_groups(
        group_cover,
        all_contigs,
    )

    if not group_solutions:
        LOGGER.warning(
            f"No combination of references covers all {len(all_contigs)} "
            f"contig(s); cannot form a complete reference set."
        )
        if outdir is not None:
            _write_selection_metrics([], outdir)
        return []

    LOGGER.info(
        f"Minimum cover uses {len(group_solutions[0])} group(s); "
        f"{len(group_solutions)} equivalent solution(s) found"
    )

    # Estimate coverage once for every reference that belongs to a group
    # solution. References within a group are interchangeable only in which
    # contigs they explain, so each is scored on its own merits.
    used_groups = set()
    for solution in group_solutions:
        used_groups.update(solution)

    ref_coverage = {}
    for g in used_groups:
        for ref in group_members[g]:
            if ref not in ref_coverage:
                ref_coverage[ref] = _reference_coverage(ref, ref_records)

    LOGGER.info(
        f"Estimated coverage for {len(ref_coverage)} reference(s) across "
        f"{len(used_groups)} group(s)"
    )

    # Each group's representative is its best-covered member that meets
    # min_ref_cov; a group with no qualifying member contributes nothing.
    group_rep = {}
    for g in used_groups:
        eligible = [
            r for r in group_members[g]
            if ref_coverage[r]["explained_fraction"] >= min_ref_cov
        ]
        if eligible:
            group_rep[g] = max(
                eligible,
                key=lambda r: ref_coverage[r]["explained_fraction"],
            )
        else:
            group_rep[g] = None
            LOGGER.debug(f"Group {g}: no member meets min_ref_cov={min_ref_cov}")

    # Total coverage of each group solution = summed explained fraction of its
    # representatives.
    scored = []
    for solution in group_solutions:
        reps = [group_rep[g] for g in solution if group_rep[g] is not None]
        if not reps:
            continue
        fractions = [ref_coverage[r]["explained_fraction"] for r in reps]
        total_coverage = sum(fractions)
        scored.append({
            "references": reps,
            "num_refs": len(reps),
            "total_coverage": total_coverage,
            "summed_est_matches": sum(ref_coverage[r]["est_matches"] for r in reps),
            "min_explained_fraction": min(fractions),
            "avg_explained_fraction": total_coverage / len(fractions),
            "per_ref": {r: ref_coverage[r] for r in reps},
        })

    LOGGER.info(
        f"{len(scored)} of {len(group_solutions)} group solution(s) yielded a "
        f"qualifying reference set (min_ref_cov={min_ref_cov})"
    )

    if not scored:
        LOGGER.warning(
            f"No group solution met the reference-coverage threshold "
            f"(min_ref_cov={min_ref_cov}); nothing selected. Consider lowering "
            f"--min-ref-cov."
        )
        if outdir is not None:
            _write_selection_metrics([], outdir)
        return []

    # Select the group solution with the best total coverage; break ties at
    # random.
    best_total = max(m["total_coverage"] for m in scored)
    tied = [m for m in scored if m["total_coverage"] == best_total]
    winner = random.choice(tied)

    if len(tied) > 1:
        LOGGER.info(
            f"{len(tied)} group solution(s) tied at total coverage "
            f"{best_total:.4f}; selected one at random"
        )
    else:
        LOGGER.info(f"Best total coverage {best_total:.4f} (single best solution)")

    # Winner first, remaining solutions by total coverage (for metrics / return).
    rest = sorted(
        (m for m in scored if m is not winner),
        key=lambda m: m["total_coverage"],
        reverse=True,
    )
    ranked = [winner] + rest
    for m in ranked:
        m["selected"] = m is winner

    LOGGER.info(
        f"Selected reference set: references={sorted(winner['references'])}, "
        f"num_refs={winner['num_refs']}, "
        f"total_coverage={winner['total_coverage']:.3f}, "
        f"summed_est_matches={winner['summed_est_matches']:.1f}, "
        f"min_explained_fraction={winner['min_explained_fraction']:.3f}, "
        f"avg_explained_fraction={winner['avg_explained_fraction']:.3f}"
    )

    if outdir is not None:
        _write_selection_metrics(ranked, outdir)

    return [(m["total_coverage"], m["references"]) for m in ranked]


def main():
    """
    Command-line entry point for VAPER reference formatting.
    Processes reference JSONL files and optionally maps a query assembly.
    """
    version = "2.0"

    parser = argparse.ArgumentParser(
        description="VAPER reference processing and selection tool"
    )
    parser.add_argument("--refs", required=True, help="Path to reference JSONL file.")
    parser.add_argument("--query", help="Path to query assembly.")
    parser.add_argument("--min-query-cov", type=float, default=0.50, help="Minimum query (contig) coverage for an alignment to be considered.")
    parser.add_argument("--min-ref-cov", type=float, default=0.70, help="Minimum reference coverage (explained fraction) for a reference to be selected.")
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
    ranked = choose_best_reference_set(
        paf_data,
        min_query_cov=args.min_query_cov,
        min_ref_cov=args.min_ref_cov,
        outdir=work_dir,
    )
    selected = ranked[0][1] if ranked else []

    if include_names:
        before = len(selected)
        selected = list(set(selected + include_names))
        LOGGER.info(f"Adding references specified by name: before={before}, added={len(include_names)}, after={len(selected)}")
        
    if selected:
        LOGGER.info(f"Final selected references ({len(selected)}): {sorted(selected)}")
        _create_subset(id_map, selected)
    else:
        LOGGER.info("No references met selection criteria.")


if __name__ == "__main__":
    main()