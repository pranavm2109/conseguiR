#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import io
import lzma
import math
import os
import zipfile
from collections import defaultdict


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
CCRE_BED = os.path.join(REPO_ROOT, "data", "raw", "ENCODE", "GRCh38-cCREs.bed")
GENE_LINKS_ZIP = os.path.join(REPO_ROOT, "data", "raw", "ENCODE", "Human-Gene-Links.zip")
GENE_LOC = os.path.join(REPO_ROOT, "data", "raw", "NCBI38", "NCBI38.gene.loc")
OUTPUT_DIR = os.path.join(REPO_ROOT, "inst", "extdata", "backend")
OUTPUT_PREFIX = os.path.join(OUTPUT_DIR, "gene_reg_graph_no_scores")

SOURCE_MEMBERS = {
    "3d_chromatin": "V4-hg38.Gene-Links.3D-Chromatin.txt",
    "crispr": "V4-hg38.Gene-Links.CRISPR.txt",
    "eqtl": "V4-hg38.Gene-Links.eQTLs.txt",
}

EVIDENCE_ALPHA = 0.6
WEIGHT_3D = 1.0
WEIGHT_CRISPR = 1.25
WEIGHT_EQTL = 1.0
SUPPORT_BONUS = 0.15


def ensure_parent(path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)


def safe_float(value: str | None) -> float | None:
    if value is None:
        return None
    value = value.strip()
    if value == "" or value.upper() == "NA":
        return None
    try:
        out = float(value)
    except ValueError:
        return None
    if not math.isfinite(out):
        return None
    return out


def significance_from_p(p_value: float | None) -> float | None:
    if p_value is None or p_value <= 0:
        return None
    return -math.log10(max(p_value, 1e-300))


def better_best(new_sig: float | None, new_mag: float | None, old_sig: float | None, old_mag: float | None) -> bool:
    new_key = (
        1 if new_sig is not None else 0,
        float("-inf") if new_sig is None else new_sig,
        float("-inf") if new_mag is None else new_mag,
    )
    old_key = (
        1 if old_sig is not None else 0,
        float("-inf") if old_sig is None else old_sig,
        float("-inf") if old_mag is None else old_mag,
    )
    return new_key > old_key


def percentile_ranks(values: list[float | None]) -> list[float | None]:
    out: list[float | None] = [None] * len(values)
    finite = [(idx, val) for idx, val in enumerate(values) if val is not None and math.isfinite(val)]
    if not finite:
        return out
    if len(finite) == 1:
        out[finite[0][0]] = 1.0
        return out

    finite.sort(key=lambda x: x[1])
    n = len(finite)
    pos = 0
    while pos < n:
        end = pos + 1
        while end < n and finite[end][1] == finite[pos][1]:
            end += 1
        avg_rank = ((pos + 1) + end) / 2.0
        pct = (avg_rank - 1.0) / (n - 1.0)
        for group_pos in range(pos, end):
            out[finite[group_pos][0]] = pct
        pos = end
    return out


def read_ccre_map(path: str) -> dict[str, dict[str, object]]:
    ccre = {}
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 6:
                continue
            reg_id = row[4].strip()
            if not reg_id:
                continue
            ccre[reg_id] = {
                "reg_chr": row[0].strip(),
                "reg_start": int(row[1]) + 1,
                "reg_end": int(row[2]),
                "reg_accession": row[3].strip(),
                "reg_element_type": row[5].strip(),
            }
    return ccre


def read_gene_loc_map(path: str) -> dict[str, dict[str, object]]:
    gene_loc: dict[str, dict[str, object]] = {}
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 6:
                continue
            gene_symbol = row[5].strip()
            if not gene_symbol:
                continue
            entry = gene_loc.get(gene_symbol)
            start = int(row[2])
            end = int(row[3])
            if entry is None:
                gene_loc[gene_symbol] = {
                    "gene_chr": row[1].strip(),
                    "gene_start": start,
                    "gene_end": end,
                }
            else:
                entry["gene_start"] = min(int(entry["gene_start"]), start)
                entry["gene_end"] = max(int(entry["gene_end"]), end)
    return gene_loc


def stream_source_collapsed(zip_path: str, member: str, source_name: str) -> list[dict[str, object]]:
    aggregates: dict[tuple[str, str], dict[str, object]] = {}
    line_count = 0
    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(member, "r") as raw:
            reader = io.TextIOWrapper(raw, encoding="utf-8", newline="")
            for line in reader:
                line_count += 1
                if line_count % 1000000 == 0:
                    print(
                        f"{source_name}: processed {line_count:,} rows; "
                        f"{len(aggregates):,} unique pairs so far",
                        flush=True,
                    )
                fields = line.rstrip("\n").split("\t")
                if source_name == "3d_chromatin":
                    fields += [""] * max(0, 9 - len(fields))
                    reg_id, ensembl_gene_id, common_gene_name, gene_type, assay_type, experiment_id, context_label, metric_raw, p_raw = fields[:9]
                elif source_name == "crispr":
                    fields += [""] * max(0, 10 - len(fields))
                    reg_id, ensembl_gene_id, common_gene_name, gene_type, _guide_id, assay_type, experiment_id, context_label, metric_raw, p_raw = fields[:10]
                else:
                    fields += [""] * max(0, 9 - len(fields))
                    reg_id, ensembl_gene_id, common_gene_name, gene_type, _variant_id, assay_type, context_label, metric_raw, p_raw = fields[:9]
                    experiment_id = assay_type

                reg_id = reg_id.strip()
                common_gene_name = common_gene_name.strip()
                ensembl_gene_id = ensembl_gene_id.strip()
                gene_id = common_gene_name if common_gene_name else ensembl_gene_id
                if not reg_id or not gene_id:
                    continue

                metric_value = safe_float(metric_raw)
                metric_magnitude = None if metric_value is None else abs(metric_value)
                p_value = safe_float(p_raw)
                significance_value = significance_from_p(p_value)

                key = (reg_id, gene_id)
                entry = aggregates.get(key)
                if entry is None:
                    entry = {
                        "reg_id": reg_id,
                        "gene_id": gene_id,
                        "ensembl_gene_id": ensembl_gene_id if ensembl_gene_id else None,
                        "gene_type": gene_type.strip() if gene_type.strip() else None,
                        "row_count": 0,
                        "metric_value": metric_value,
                        "p_value": p_value,
                        "best_sig": significance_value,
                        "best_mag": metric_magnitude,
                        "max_significance": significance_value,
                        "max_magnitude": metric_magnitude,
                        "assay_types": set(),
                        "context_labels": set(),
                    }
                    aggregates[key] = entry

                entry["row_count"] += 1
                if assay_type.strip():
                    entry["assay_types"].add(assay_type.strip())
                if context_label.strip():
                    entry["context_labels"].add(context_label.strip())

                if metric_magnitude is not None and (
                    entry["max_magnitude"] is None or metric_magnitude > entry["max_magnitude"]
                ):
                    entry["max_magnitude"] = metric_magnitude
                if significance_value is not None and (
                    entry["max_significance"] is None or significance_value > entry["max_significance"]
                ):
                    entry["max_significance"] = significance_value

                if better_best(significance_value, metric_magnitude, entry["best_sig"], entry["best_mag"]):
                    entry["metric_value"] = metric_value
                    entry["p_value"] = p_value
                    entry["best_sig"] = significance_value
                    entry["best_mag"] = metric_magnitude

    records = list(aggregates.values())
    magnitude_ranks = percentile_ranks([rec["max_magnitude"] for rec in records])
    significance_ranks = percentile_ranks([rec["max_significance"] for rec in records])
    for rec, mag_rank, sig_rank in zip(records, magnitude_ranks, significance_ranks):
        if sig_rank is not None and mag_rank is not None:
            evidence = EVIDENCE_ALPHA * sig_rank + (1.0 - EVIDENCE_ALPHA) * mag_rank
        elif sig_rank is not None:
            evidence = sig_rank
        else:
            evidence = mag_rank
        rec["source_significance_rank"] = sig_rank
        rec["source_magnitude_rank"] = mag_rank
        rec["source_evidence"] = evidence
        rec["assay_types"] = "|".join(sorted(rec["assay_types"]))
        rec["context_labels"] = "|".join(sorted(rec["context_labels"]))
    return records


def merge_sources(
    per_source: dict[str, list[dict[str, object]]],
    ccre_map: dict[str, dict[str, object]],
    gene_loc_map: dict[str, dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    combined: dict[tuple[str, str], dict[str, object]] = {}

    source_spec = {
        "3d_chromatin": ("score_3d", "p_3d", WEIGHT_3D, "support_3d"),
        "crispr": ("effect_crispr", "p_crispr", WEIGHT_CRISPR, "support_crispr"),
        "eqtl": ("slope_eqtl", "p_eqtl", WEIGHT_EQTL, "support_eqtl"),
    }

    for source_name, rows in per_source.items():
        metric_col, p_col, _, support_col = source_spec[source_name]
        sig_col = f"significance_{source_name}"
        mag_col = f"magnitude_{source_name}"
        evidence_col = f"evidence_{source_name}"
        assay_col = f"assay_types_{source_name}"
        context_col = f"context_labels_{source_name}"
        rows_col = f"rows_{source_name}"
        for row in rows:
            key = (row["reg_id"], row["gene_id"])
            target = combined.get(key)
            if target is None:
                ccre = ccre_map.get(row["reg_id"], {})
                gene_loc = gene_loc_map.get(row["gene_id"], {})
                target = {
                    "reg_id": row["reg_id"],
                    "gene_id": row["gene_id"],
                    "ensembl_gene_id": row.get("ensembl_gene_id"),
                    "gene_type": row.get("gene_type"),
                    "reg_chr": ccre.get("reg_chr"),
                    "reg_start": ccre.get("reg_start"),
                    "reg_end": ccre.get("reg_end"),
                    "reg_accession": ccre.get("reg_accession"),
                    "reg_element_type": ccre.get("reg_element_type"),
                    "gene_chr": gene_loc.get("gene_chr"),
                    "gene_start": gene_loc.get("gene_start"),
                    "gene_end": gene_loc.get("gene_end"),
                    "support_3d": False,
                    "support_crispr": False,
                    "support_eqtl": False,
                }
                combined[key] = target
            target[metric_col] = row.get("metric_value")
            target[p_col] = row.get("p_value")
            target[sig_col] = row.get("source_significance_rank")
            target[mag_col] = row.get("source_magnitude_rank")
            target[evidence_col] = row.get("source_evidence")
            target[assay_col] = row.get("assay_types")
            target[context_col] = row.get("context_labels")
            target[rows_col] = row.get("row_count")
            target[support_col] = True
            if not target.get("ensembl_gene_id"):
                target["ensembl_gene_id"] = row.get("ensembl_gene_id")
            if not target.get("gene_type"):
                target["gene_type"] = row.get("gene_type")

    final_edges: list[dict[str, object]] = []
    reg_target_labels: dict[str, list[str]] = defaultdict(list)
    for entry in combined.values():
        support_count = int(entry.get("support_3d", False)) + int(entry.get("support_crispr", False)) + int(entry.get("support_eqtl", False))
        combined_edge_score = (
            WEIGHT_3D * float(entry.get("evidence_3d_chromatin") or 0.0) +
            WEIGHT_CRISPR * float(entry.get("evidence_crispr") or 0.0) +
            WEIGHT_EQTL * float(entry.get("evidence_eqtl") or 0.0) +
            SUPPORT_BONUS * max(support_count - 1, 0)
        )
        if entry.get("support_3d") and entry.get("support_crispr") and entry.get("support_eqtl"):
            link_method = "3d_chromatin|crispr|eqtl"
        elif entry.get("support_3d") and entry.get("support_crispr"):
            link_method = "3d_chromatin|crispr"
        elif entry.get("support_3d") and entry.get("support_eqtl"):
            link_method = "3d_chromatin|eqtl"
        elif entry.get("support_crispr") and entry.get("support_eqtl"):
            link_method = "crispr|eqtl"
        elif entry.get("support_3d"):
            link_method = "3d_chromatin"
        elif entry.get("support_crispr"):
            link_method = "crispr"
        elif entry.get("support_eqtl"):
            link_method = "eqtl"
        else:
            link_method = None

        edge_row = dict(entry)
        edge_row["support_count"] = support_count
        edge_row["combined_edge_score"] = combined_edge_score
        edge_row["link_value"] = combined_edge_score
        edge_row["link_score"] = combined_edge_score
        edge_row["confidence"] = combined_edge_score
        edge_row["weight"] = 1.0 / (1.0 + combined_edge_score)
        edge_row["evidence_count"] = int(entry.get("rows_3d_chromatin") or 0) + int(entry.get("rows_crispr") or 0) + int(entry.get("rows_eqtl") or 0)
        edge_row["link_method"] = link_method
        edge_row["from"] = entry["reg_id"]
        edge_row["to"] = entry["gene_id"]
        final_edges.append(edge_row)
        reg_target_labels[entry["reg_id"]].append(entry["gene_id"])

    gene_nodes = {}
    reg_nodes = {}
    for edge in final_edges:
        gene_id = edge["gene_id"]
        if gene_id not in gene_nodes:
            gene_nodes[gene_id] = {
                "name": gene_id,
                "node_id": gene_id,
                "node_type": "gene",
                "chr": edge.get("gene_chr"),
                "start": edge.get("gene_start"),
                "end": edge.get("gene_end"),
                "ensembl_gene_id": edge.get("ensembl_gene_id"),
                "gene_type": edge.get("gene_type"),
            }
        reg_id = edge["reg_id"]
        if reg_id not in reg_nodes:
            reg_nodes[reg_id] = {
                "name": reg_id,
                "node_id": reg_id,
                "node_type": "reg",
                "chr": edge.get("reg_chr"),
                "start": edge.get("reg_start"),
                "end": edge.get("reg_end"),
                "reg_chr": edge.get("reg_chr"),
                "reg_start": edge.get("reg_start"),
                "reg_end": edge.get("reg_end"),
                "reg_accession": edge.get("reg_accession"),
                "reg_element_type": edge.get("reg_element_type"),
            }

    all_nodes = list(gene_nodes.values()) + list(reg_nodes.values())
    all_nodes.sort(key=lambda row: (row["node_type"], row["node_id"]))
    for idx, row in enumerate(all_nodes, start=1):
        row["node_index"] = idx
    node_index = {row["node_id"]: row["node_index"] for row in all_nodes}

    compact_edges = []
    for edge in final_edges:
        compact_edges.append({
            "from_idx": node_index[edge["from"]],
            "to_idx": node_index[edge["to"]],
            "confidence": edge["confidence"],
            "link_method": edge["link_method"] or "",
        })

    reg_target_rows = [
        {"reg_elem_id": reg_id, "label": "|".join(sorted(set(targets)))}
        for reg_id, targets in sorted(reg_target_labels.items())
    ]

    return all_nodes, compact_edges, reg_target_rows


def write_tsv_xz(rows: list[dict[str, object]], path: str, fieldnames: list[str]) -> None:
    ensure_parent(path)
    with lzma.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_tsv_gz(rows: list[dict[str, object]], path: str, fieldnames: list[str]) -> None:
    ensure_parent(path)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    print("Reading ENCODE cCRE metadata...", flush=True)
    ccre_map = read_ccre_map(CCRE_BED)
    print("Reading gene locations...", flush=True)
    gene_loc_map = read_gene_loc_map(GENE_LOC)

    per_source = {}
    for source_name, member in SOURCE_MEMBERS.items():
        print(f"Streaming {source_name} from {member} ...", flush=True)
        per_source[source_name] = stream_source_collapsed(GENE_LINKS_ZIP, member, source_name)
        print(f"Collapsed {source_name}: {len(per_source[source_name])} unique reg-gene pairs", flush=True)

    print("Merging sources and building compact backend seed...", flush=True)
    nodes, compact_edges, reg_target_rows = merge_sources(per_source, ccre_map, gene_loc_map)

    node_fields = [
        "name", "node_id", "node_type", "chr", "start", "end",
        "ensembl_gene_id", "gene_type", "reg_chr", "reg_start", "reg_end",
        "reg_accession", "reg_element_type", "node_index"
    ]
    edge_fields = ["from_idx", "to_idx", "confidence", "link_method"]
    label_fields = ["reg_elem_id", "label"]

    write_tsv_xz(nodes, OUTPUT_PREFIX + "_nodes_compact.tsv.xz", node_fields)
    write_tsv_xz(compact_edges, OUTPUT_PREFIX + "_edges_compact.tsv.xz", edge_fields)
    write_tsv_gz(reg_target_rows, os.path.join(OUTPUT_DIR, "reg_target_labels.tsv.gz"), label_fields)

    print(f"Wrote {len(nodes)} nodes and {len(compact_edges)} compact edges", flush=True)


if __name__ == "__main__":
    main()
