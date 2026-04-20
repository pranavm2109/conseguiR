#!/usr/bin/env python3

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd


@dataclass
class DiffusionConfig:
    nodes_path: str = "data/processed/gene_reg_graph_scored_nodes.tsv.gz"
    edges_path: str = "data/processed/gene_reg_graph_scored_edges.tsv.gz"
    output_dir: str = "data/processed"
    output_stem: str = "gene_reg_graph_diffusion"
    top_k: int = 3
    confidence_power: float = 2.0
    beta_germline: float = 0.5
    beta_somatic: float = 0.5
    beta_epigenomic: float = 0.7
    integration_weight_germline: float = 1.0
    integration_weight_somatic: float = 1.0
    integration_weight_epigenomic: float = 1.0
    positive_only: bool = False
    reg_signal_clip: float = 5.0
    top_n_to_save: int = 50


def ensure_parent_dir(path: str) -> str:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def read_scored_gene_reg_nodes(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Scored gene-reg node file does not exist: {path}")

    return pd.read_csv(path, sep="\t", low_memory=False)


def read_scored_gene_reg_edges(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Scored gene-reg edge file does not exist: {path}")

    return pd.read_csv(path, sep="\t", low_memory=False)


def validate_scored_gene_reg_nodes(nodes: pd.DataFrame) -> pd.DataFrame:
    required_cols = {
        "node_id",
        "node_type",
        "somatic_score",
        "germline_score",
        "epigenomic_score",
    }
    missing_cols = required_cols.difference(nodes.columns)
    if missing_cols:
        raise ValueError(f"Scored gene-reg node table is missing columns: {sorted(missing_cols)}")

    node_types = set(nodes["node_type"].dropna().unique())
    bad_types = node_types.difference({"gene", "reg"})
    if bad_types:
        raise ValueError(f"Unsupported node types found in node table: {sorted(bad_types)}")

    if nodes["node_id"].isna().any() or not nodes["node_id"].is_unique:
        raise ValueError("`node_id` must be non-missing and unique in the scored node table.")

    score_cols = ["somatic_score", "germline_score", "epigenomic_score"]
    for col in score_cols:
        if nodes[col].isna().any():
            raise ValueError(f"Score column `{col}` contains NA values.")

    return nodes.copy()


def validate_scored_gene_reg_edges(edges: pd.DataFrame) -> pd.DataFrame:
    required_cols = {"from", "to", "confidence"}
    missing_cols = required_cols.difference(edges.columns)
    if missing_cols:
        raise ValueError(f"Scored gene-reg edge table is missing columns: {sorted(missing_cols)}")

    if edges[["from", "to"]].isna().any().any():
        raise ValueError("Edge table contains missing endpoint identifiers.")

    if edges["confidence"].isna().any():
        raise ValueError("Edge table contains missing confidence values.")

    return edges.copy()


def safe_zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    mu = np.mean(values)
    sd = np.std(values)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(values, dtype=float)
    return (values - mu) / sd


def positive_part(values: np.ndarray) -> np.ndarray:
    return np.maximum(np.asarray(values, dtype=float), 0.0)


def magnitude(values: np.ndarray) -> np.ndarray:
    return np.abs(np.asarray(values, dtype=float))


def signed_stouffer(
    first: np.ndarray,
    second: np.ndarray,
    third: np.ndarray,
    w_first: float,
    w_second: float,
    w_third: float,
) -> np.ndarray:
    weights = np.asarray([w_first, w_second, w_third], dtype=float)
    denom = float(np.sqrt(np.sum(np.square(weights))))
    if not np.isfinite(denom) or denom <= 0:
        raise ValueError("Integration weights must contain at least one finite non-zero value.")
    return (
        weights[0] * np.asarray(first, dtype=float)
        + weights[1] * np.asarray(second, dtype=float)
        + weights[2] * np.asarray(third, dtype=float)
    ) / denom


def split_gene_and_reg_nodes(nodes: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    gene_nodes = (
        nodes.loc[nodes["node_type"] == "gene"]
        .copy()
        .sort_values("node_id")
        .reset_index(drop=True)
    )
    reg_nodes = (
        nodes.loc[nodes["node_type"] == "reg"]
        .copy()
        .sort_values("node_id")
        .reset_index(drop=True)
    )

    gene_nodes["gene_local_idx"] = np.arange(len(gene_nodes), dtype=np.int64)
    reg_nodes["reg_local_idx"] = np.arange(len(reg_nodes), dtype=np.int64)
    return gene_nodes, reg_nodes


def build_node_type_map(nodes: pd.DataFrame) -> Dict[str, str]:
    return dict(zip(nodes["node_id"], nodes["node_type"]))


def orient_gene_reg_edges(
    edges: pd.DataFrame,
    node_type_map: Dict[str, str],
    gene_nodes: pd.DataFrame,
    reg_nodes: pd.DataFrame,
    confidence_power: float,
) -> pd.DataFrame:
    work = edges.copy()
    work["from_type"] = work["from"].map(node_type_map)
    work["to_type"] = work["to"].map(node_type_map)

    valid_mask = (
        ((work["from_type"] == "reg") & (work["to_type"] == "gene"))
        | ((work["from_type"] == "gene") & (work["to_type"] == "reg"))
    )
    work = work.loc[valid_mask].copy()

    work["reg_node_id"] = np.where(work["from_type"] == "reg", work["from"], work["to"])
    work["gene_node_id"] = np.where(work["from_type"] == "gene", work["from"], work["to"])

    gene_lookup = dict(zip(gene_nodes["node_id"], gene_nodes["gene_local_idx"]))
    reg_lookup = dict(zip(reg_nodes["node_id"], reg_nodes["reg_local_idx"]))

    work["gene_idx"] = work["gene_node_id"].map(gene_lookup)
    work["reg_idx"] = work["reg_node_id"].map(reg_lookup)
    work = work.dropna(subset=["gene_idx", "reg_idx"]).copy()

    work["gene_idx"] = work["gene_idx"].astype(np.int64)
    work["reg_idx"] = work["reg_idx"].astype(np.int64)
    work["edge_weight"] = np.asarray(work["confidence"], dtype=float) ** confidence_power
    return work


def extract_gene_and_reg_signals(
    gene_nodes: pd.DataFrame,
    reg_nodes: pd.DataFrame,
) -> Dict[str, np.ndarray]:
    return {
        "gene_germline": gene_nodes["germline_score"].to_numpy(dtype=float),
        "gene_somatic": gene_nodes["somatic_score"].to_numpy(dtype=float),
        "gene_epigenomic": gene_nodes["epigenomic_score"].to_numpy(dtype=float),
        "reg_germline": reg_nodes["germline_score"].to_numpy(dtype=float),
        "reg_somatic": reg_nodes["somatic_score"].to_numpy(dtype=float),
        "reg_epigenomic": reg_nodes["epigenomic_score"].to_numpy(dtype=float),
    }


def topk_gene_aggregate(
    oriented_edges: pd.DataFrame,
    reg_signal: np.ndarray,
    n_gene: int,
    top_k: int,
    positive_only: bool = False,
    reg_signal_clip: float = 5.0,
) -> np.ndarray:
    work = oriented_edges.loc[:, ["gene_idx", "reg_idx", "edge_weight"]].copy()
    work = work.sort_values(["gene_idx", "edge_weight"], ascending=[True, False])
    topk_edges = work.groupby("gene_idx", sort=False).head(top_k).copy()

    reg_signal_used = np.clip(np.asarray(reg_signal, dtype=float), -reg_signal_clip, reg_signal_clip)
    topk_edges["reg_signal"] = reg_signal_used[topk_edges["reg_idx"].to_numpy(dtype=np.int64)]

    if positive_only:
        topk_edges = topk_edges.loc[topk_edges["reg_signal"] > 0].copy()

    if topk_edges.empty:
        return np.zeros(n_gene, dtype=float)

    topk_edges["weighted_signal"] = topk_edges["edge_weight"] * topk_edges["reg_signal"]
    signal_sum = topk_edges.groupby("gene_idx", sort=False)["weighted_signal"].sum()
    weight_mean = topk_edges.groupby("gene_idx", sort=False)["edge_weight"].mean()
    aggregated = signal_sum / weight_mean

    out = np.zeros(n_gene, dtype=float)
    out[aggregated.index.to_numpy(dtype=np.int64)] = aggregated.to_numpy(dtype=float)
    return out


def run_controlled_re_to_gene_diffusion(
    gene_nodes: pd.DataFrame,
    reg_nodes: pd.DataFrame,
    oriented_edges: pd.DataFrame,
    config: DiffusionConfig,
) -> pd.DataFrame:
    signals = extract_gene_and_reg_signals(gene_nodes, reg_nodes)
    n_gene = len(gene_nodes)

    incoming_germline_raw = topk_gene_aggregate(
        oriented_edges=oriented_edges,
        reg_signal=signals["reg_germline"],
        n_gene=n_gene,
        top_k=config.top_k,
        positive_only=config.positive_only,
        reg_signal_clip=config.reg_signal_clip,
    )
    incoming_somatic_raw = topk_gene_aggregate(
        oriented_edges=oriented_edges,
        reg_signal=signals["reg_somatic"],
        n_gene=n_gene,
        top_k=config.top_k,
        positive_only=config.positive_only,
        reg_signal_clip=config.reg_signal_clip,
    )
    incoming_epigenomic_raw = topk_gene_aggregate(
        oriented_edges=oriented_edges,
        reg_signal=signals["reg_epigenomic"],
        n_gene=n_gene,
        top_k=config.top_k,
        positive_only=config.positive_only,
        reg_signal_clip=config.reg_signal_clip,
    )

    incoming_germline = safe_zscore(incoming_germline_raw)
    incoming_somatic = safe_zscore(incoming_somatic_raw)
    incoming_epigenomic = safe_zscore(incoming_epigenomic_raw)

    prediff_germline = signals["gene_germline"]
    prediff_somatic = signals["gene_somatic"]
    prediff_epigenomic = signals["gene_epigenomic"]

    post_germline = prediff_germline + config.beta_germline * incoming_germline
    post_somatic = prediff_somatic + config.beta_somatic * incoming_somatic
    post_epigenomic = prediff_epigenomic + config.beta_epigenomic * incoming_epigenomic

    prediff_norm = np.sqrt(prediff_germline**2 + prediff_somatic**2 + prediff_epigenomic**2)
    post_norm = np.sqrt(post_germline**2 + post_somatic**2 + post_epigenomic**2)

    prediff_germline_enrichment = positive_part(prediff_germline)
    prediff_somatic_enrichment = positive_part(prediff_somatic)
    prediff_epigenomic_strength = magnitude(prediff_epigenomic)
    post_germline_enrichment = positive_part(post_germline)
    post_somatic_enrichment = positive_part(post_somatic)
    post_epigenomic_strength = magnitude(post_epigenomic)

    prediff_vulnerability = np.sqrt(
        prediff_germline_enrichment**2
        + prediff_somatic_enrichment**2
        + prediff_epigenomic_strength**2
    )
    post_vulnerability = np.sqrt(
        post_germline_enrichment**2
        + post_somatic_enrichment**2
        + post_epigenomic_strength**2
    )

    prediff_integrated = signed_stouffer(
        prediff_germline,
        prediff_somatic,
        prediff_epigenomic,
        config.integration_weight_germline,
        config.integration_weight_somatic,
        config.integration_weight_epigenomic,
    )
    post_integrated = signed_stouffer(
        post_germline,
        post_somatic,
        post_epigenomic,
        config.integration_weight_germline,
        config.integration_weight_somatic,
        config.integration_weight_epigenomic,
    )

    gene_label_col = "node_id"
    if "name" in gene_nodes.columns:
        gene_label_col = "name"

    out = gene_nodes.loc[:, ["node_id", gene_label_col]].copy()
    if gene_label_col != "gene_name":
        out = out.rename(columns={gene_label_col: "gene_name"})
    elif "gene_name" not in out.columns:
        out["gene_name"] = out["node_id"]

    out["prediff_germline"] = prediff_germline
    out["prediff_somatic"] = prediff_somatic
    out["prediff_epigenomic"] = prediff_epigenomic
    out["prediff_norm"] = prediff_norm
    out["prediff_germline_enrichment"] = prediff_germline_enrichment
    out["prediff_somatic_enrichment"] = prediff_somatic_enrichment
    out["prediff_epigenomic_strength"] = prediff_epigenomic_strength
    out["prediff_vulnerability"] = prediff_vulnerability
    out["prediff_integrated"] = prediff_integrated

    out["incoming_germline_raw"] = incoming_germline_raw
    out["incoming_somatic_raw"] = incoming_somatic_raw
    out["incoming_epigenomic_raw"] = incoming_epigenomic_raw

    out["incoming_germline"] = incoming_germline
    out["incoming_somatic"] = incoming_somatic
    out["incoming_epigenomic"] = incoming_epigenomic

    out["post_germline"] = post_germline
    out["post_somatic"] = post_somatic
    out["post_epigenomic"] = post_epigenomic
    out["post_norm"] = post_norm
    out["post_germline_enrichment"] = post_germline_enrichment
    out["post_somatic_enrichment"] = post_somatic_enrichment
    out["post_epigenomic_strength"] = post_epigenomic_strength
    out["post_vulnerability"] = post_vulnerability
    out["post_integrated"] = post_integrated

    out["prediff_norm_rank"] = out["prediff_norm"].rank(ascending=False, method="average")
    out["post_norm_rank"] = out["post_norm"].rank(ascending=False, method="average")
    out["prediff_vulnerability_rank"] = out["prediff_vulnerability"].rank(ascending=False, method="average")
    out["post_vulnerability_rank"] = out["post_vulnerability"].rank(ascending=False, method="average")
    out["prediff_rank"] = out["prediff_integrated"].rank(ascending=False, method="average")
    out["post_rank"] = out["post_integrated"].rank(ascending=False, method="average")
    out["rank_shift"] = out["prediff_rank"] - out["post_rank"]

    out = out.sort_values(
        ["post_integrated", "post_vulnerability", "post_norm"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    return out


def save_diffusion_outputs(
    diffusion_df: pd.DataFrame,
    config: DiffusionConfig,
) -> Dict[str, str]:
    all_genes_path = ensure_parent_dir(
        os.path.join(config.output_dir, f"{config.output_stem}_all_genes.tsv")
    )
    top_genes_path = ensure_parent_dir(
        os.path.join(config.output_dir, f"{config.output_stem}_top{config.top_n_to_save}.tsv")
    )

    diffusion_df.to_csv(all_genes_path, sep="\t", index=False)
    diffusion_df.head(config.top_n_to_save).to_csv(top_genes_path, sep="\t", index=False)

    return {
        "all_genes_path": all_genes_path,
        "top_genes_path": top_genes_path,
    }


def run_gene_reg_diffusion(config: DiffusionConfig = DiffusionConfig()) -> Dict[str, object]:
    nodes = validate_scored_gene_reg_nodes(read_scored_gene_reg_nodes(config.nodes_path))
    edges = validate_scored_gene_reg_edges(read_scored_gene_reg_edges(config.edges_path))

    gene_nodes, reg_nodes = split_gene_and_reg_nodes(nodes)
    node_type_map = build_node_type_map(nodes)
    oriented_edges = orient_gene_reg_edges(
        edges=edges,
        node_type_map=node_type_map,
        gene_nodes=gene_nodes,
        reg_nodes=reg_nodes,
        confidence_power=config.confidence_power,
    )

    diffusion_df = run_controlled_re_to_gene_diffusion(
        gene_nodes=gene_nodes,
        reg_nodes=reg_nodes,
        oriented_edges=oriented_edges,
        config=config,
    )
    output_paths = save_diffusion_outputs(diffusion_df, config)

    return {
        "nodes": nodes,
        "edges": edges,
        "gene_nodes": gene_nodes,
        "reg_nodes": reg_nodes,
        "oriented_edges": oriented_edges,
        "diffusion": diffusion_df,
        "output_paths": output_paths,
        "config": config,
    }


def main() -> None:
    result = run_gene_reg_diffusion()
    diffusion_df = result["diffusion"]

    print("\n=== Gene-reg diffusion complete ===")
    print(f"Genes processed: {len(diffusion_df)}")
    print("Top diffused genes:")
    print(
        diffusion_df.loc[
            : min(9, len(diffusion_df) - 1),
            ["node_id", "gene_name", "prediff_integrated", "post_integrated", "rank_shift"],
        ]
    )
    print("\nSaved outputs:")
    for label, path in result["output_paths"].items():
        print(f"  {label}: {path}")


if __name__ == "__main__":
    main()
