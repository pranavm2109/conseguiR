#!/usr/bin/env python3

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

try:
    from ortools.sat.python import cp_model
except ImportError as exc:  # pragma: no cover
    cp_model = None
    _ORTOOLS_IMPORT_ERROR = exc
else:
    _ORTOOLS_IMPORT_ERROR = None


@dataclass
class SubgraphConfig:
    diffusion_path: str = "data/processed/gene_reg_graph_diffusion_all_genes.tsv"
    gg_nodes_path: str = "data/processed/gene_gene_graph_nodes.tsv.gz"
    gg_edges_path: str = "data/processed/gene_gene_graph_edges.tsv.gz"
    output_dir: str = "data/processed"
    output_stem: str = "gene_gene_selected_subgraph"
    target_genes: int = 50
    candidate_pool_size: int = 400
    min_confidence: float = 0.0
    max_edges_in_model: int = 12000
    node_prize_weight: float = 1.0
    edge_conf_weight: float = 1.0
    edge_cost_weight: float = 1.0
    node_scale: int = 1000
    edge_scale: int = 1000
    max_time_seconds: int = 600
    num_workers: int = 8
    random_seed: int = 42
    prize_column: str = "post_integrated"
    confidence_column: str = "confidence"
    edge_cost_column: str = "weight"


def require_ortools() -> None:
    if cp_model is None:
        raise ImportError(
            "OR-Tools is required to call the cardinality-constrained subgraph, "
            "but it could not be imported inside the managed conseguiR Python "
            "environment."
        ) from _ORTOOLS_IMPORT_ERROR


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def read_diffusion_results(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Diffusion results file does not exist: {path}")
    return pd.read_csv(path, sep="\t", low_memory=False)


def read_gene_gene_nodes(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Gene-gene node file does not exist: {path}")
    return pd.read_csv(path, sep="\t", low_memory=False)


def read_gene_gene_edges(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Gene-gene edge file does not exist: {path}")
    return pd.read_csv(path, sep="\t", low_memory=False)


def validate_diffusion_results(
    diffusion: pd.DataFrame,
    prize_column: str,
) -> pd.DataFrame:
    required = {"node_id", "gene_name"}
    missing = required.difference(diffusion.columns)
    if missing:
        raise ValueError(f"Diffusion results are missing columns: {sorted(missing)}")

    if prize_column not in diffusion.columns:
        post_cols = {"post_germline", "post_somatic", "post_epigenomic"}
        if not post_cols.issubset(diffusion.columns):
            raise ValueError(
                f"Diffusion results are missing prize column `{prize_column}` and do not "
                "contain all three post-diffusion component columns."
            )

    if diffusion["node_id"].isna().any():
        raise ValueError("Diffusion results contain missing `node_id` values.")

    if diffusion["gene_name"].isna().any():
        raise ValueError("Diffusion results contain missing `gene_name` values.")

    return diffusion.copy()


def validate_gene_gene_nodes(nodes: pd.DataFrame) -> pd.DataFrame:
    required = {"node_id", "node_type"}
    missing = required.difference(nodes.columns)
    if missing:
        raise ValueError(f"Gene-gene node table is missing columns: {sorted(missing)}")

    if nodes["node_id"].isna().any() or not nodes["node_id"].is_unique:
        raise ValueError("Gene-gene node table must have unique, non-missing `node_id` values.")

    bad_types = set(nodes["node_type"].dropna().unique()).difference({"gene"})
    if bad_types:
        raise ValueError(f"Gene-gene node table contains unsupported node types: {sorted(bad_types)}")

    return nodes.copy()


def validate_gene_gene_edges(
    edges: pd.DataFrame,
    confidence_column: str,
    edge_cost_column: str,
) -> pd.DataFrame:
    required = {"from", "to", confidence_column, edge_cost_column}
    missing = required.difference(edges.columns)
    if missing:
        raise ValueError(f"Gene-gene edge table is missing columns: {sorted(missing)}")

    if edges[["from", "to"]].isna().any().any():
        raise ValueError("Gene-gene edge table contains missing endpoints.")

    if edges[confidence_column].isna().any():
        raise ValueError(f"Gene-gene edge table contains missing `{confidence_column}` values.")

    if edges[edge_cost_column].isna().any():
        raise ValueError(f"Gene-gene edge table contains missing `{edge_cost_column}` values.")

    return edges.copy()


def compute_node_prizes(
    diffusion: pd.DataFrame,
    prize_column: str,
) -> pd.DataFrame:
    work = diffusion.copy()
    fallback_column = None
    if prize_column in work.columns:
        fallback_column = prize_column
    elif "post_integrated" in work.columns:
        fallback_column = "post_integrated"
    elif "post_vulnerability" in work.columns:
        fallback_column = "post_vulnerability"
    elif "post_norm" in work.columns:
        fallback_column = "post_norm"

    if fallback_column is not None:
        work["prize"] = pd.to_numeric(work[fallback_column], errors="coerce")
    else:
        work["prize"] = np.sqrt(
            np.square(pd.to_numeric(work["post_germline"], errors="coerce"))
            + np.square(pd.to_numeric(work["post_somatic"], errors="coerce"))
            + np.square(pd.to_numeric(work["post_epigenomic"], errors="coerce"))
        )

    if work["prize"].isna().any():
        raise ValueError("Unable to compute node prizes because some diffusion scores are non-numeric.")

    work["prize"] = work["prize"].astype(float)
    finite_mask = np.isfinite(work["prize"].to_numpy(dtype=float))
    if not finite_mask.any():
        raise ValueError("Unable to compute node prizes because all prize values are non-finite.")
    if not finite_mask.all():
        finite_max = float(work.loc[finite_mask, "prize"].max())
        work.loc[~finite_mask, "prize"] = finite_max

    work["prize"] = np.maximum(work["prize"], 0.0)
    return work


def attach_diffusion_prizes_to_gene_graph(
    gg_nodes: pd.DataFrame,
    diffusion: pd.DataFrame,
) -> pd.DataFrame:
    diffusion = compute_node_prizes(
        diffusion,
        prize_column="prize" if "prize" in diffusion.columns else "post_integrated",
    )
    diff_by_node = diffusion.drop_duplicates(subset=["node_id"]).copy()

    if "name" not in gg_nodes.columns:
        gg_nodes["name"] = gg_nodes["node_id"]

    attached = gg_nodes.merge(
        diff_by_node[
            [
                col
                for col in [
                    "node_id",
                    "gene_name",
                    "prize",
                    "prediff_integrated",
                    "prediff_vulnerability",
                    "prediff_norm",
                    "post_integrated",
                    "post_vulnerability",
                    "post_norm",
                    "post_germline",
                    "post_somatic",
                    "post_epigenomic",
                    "rank_shift",
                ]
                if col in diff_by_node.columns
            ]
        ],
        on="node_id",
        how="left",
    )

    attached = attached.dropna(subset=["prize"]).copy()
    if attached.empty:
        raise ValueError("No gene-gene graph nodes could be matched to the diffusion output.")

    attached["gene_name"] = attached["gene_name"].fillna(attached["name"])
    return attached


def build_candidate_pool(
    scored_nodes: pd.DataFrame,
    target_genes: int,
    candidate_pool_size: int,
) -> pd.DataFrame:
    if target_genes <= 0:
        raise ValueError("`target_genes` must be a positive integer.")

    if candidate_pool_size < target_genes:
        raise ValueError("`candidate_pool_size` must be at least as large as `target_genes`.")

    ordered = (
        scored_nodes.sort_values(["prize", "gene_name"], ascending=[False, True])
        .reset_index(drop=True)
        .copy()
    )
    if candidate_pool_size > ordered.shape[0]:
        raise ValueError(
            f"`candidate_pool_size` ({candidate_pool_size}) exceeds the number of "
            f"available diffusion-ranked genes ({ordered.shape[0]})."
        )
    candidate_nodes = ordered.head(candidate_pool_size).copy()
    if candidate_nodes.shape[0] < target_genes:
        raise ValueError(
            f"Only {candidate_nodes.shape[0]} candidate genes are available, "
            f"but target cardinality is {target_genes}."
        )
    return candidate_nodes


def restrict_edges_to_candidates(
    edges: pd.DataFrame,
    candidate_node_ids: Iterable[str],
    config: SubgraphConfig,
) -> pd.DataFrame:
    candidate_node_ids = set(candidate_node_ids)
    work = edges.loc[
        edges["from"].isin(candidate_node_ids) & edges["to"].isin(candidate_node_ids)
    ].copy()

    if config.min_confidence > 0:
        work = work.loc[work[config.confidence_column] >= config.min_confidence].copy()

    work["confidence_raw"] = pd.to_numeric(work[config.confidence_column], errors="coerce")
    work["edge_cost_raw"] = pd.to_numeric(work[config.edge_cost_column], errors="coerce")

    if work[["confidence_raw", "edge_cost_raw"]].isna().any().any():
        raise ValueError("Some candidate edges have non-numeric confidence or cost values.")

    confidence_scale = float(work["confidence_raw"].max()) if not work.empty else 1.0
    if not math.isfinite(confidence_scale) or confidence_scale <= 0:
        confidence_scale = 1.0

    work["confidence_unit"] = work["confidence_raw"] / confidence_scale
    work["edge_reward_raw"] = (
        config.edge_conf_weight * work["confidence_unit"]
        - config.edge_cost_weight * work["edge_cost_raw"]
    )
    work = work.loc[work["edge_reward_raw"] > 0].copy()

    work = work.sort_values(
        ["edge_reward_raw", "confidence_raw", "edge_cost_raw"],
        ascending=[False, False, True],
    ).head(config.max_edges_in_model)

    return work.reset_index(drop=True)


def reindex_candidate_problem(
    candidate_nodes: pd.DataFrame,
    candidate_edges: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    candidate_nodes = candidate_nodes.sort_values("node_id").reset_index(drop=True).copy()
    candidate_nodes["opt_idx"] = np.arange(candidate_nodes.shape[0], dtype=np.int64)
    node_to_opt = dict(zip(candidate_nodes["node_id"], candidate_nodes["opt_idx"]))

    candidate_edges = candidate_edges.copy()
    candidate_edges["u_opt"] = candidate_edges["from"].map(node_to_opt)
    candidate_edges["v_opt"] = candidate_edges["to"].map(node_to_opt)
    candidate_edges = candidate_edges.dropna(subset=["u_opt", "v_opt"]).copy()
    candidate_edges["u_opt"] = candidate_edges["u_opt"].astype(np.int64)
    candidate_edges["v_opt"] = candidate_edges["v_opt"].astype(np.int64)
    candidate_edges = candidate_edges.loc[candidate_edges["u_opt"] != candidate_edges["v_opt"]].copy()

    return candidate_nodes, candidate_edges


def integerize_objective_terms(
    candidate_nodes: pd.DataFrame,
    candidate_edges: pd.DataFrame,
    config: SubgraphConfig,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    candidate_nodes = candidate_nodes.copy()
    candidate_edges = candidate_edges.copy()

    candidate_nodes["node_obj"] = np.round(
        config.node_prize_weight * candidate_nodes["prize"].astype(float) * config.node_scale
    ).astype(int)
    candidate_edges["edge_obj"] = np.round(
        candidate_edges["edge_reward_raw"].astype(float) * config.edge_scale
    ).astype(int)
    candidate_edges = candidate_edges.loc[candidate_edges["edge_obj"] > 0].copy()

    return candidate_nodes, candidate_edges


def solve_cardinality_subgraph(
    candidate_nodes: pd.DataFrame,
    candidate_edges: pd.DataFrame,
    config: SubgraphConfig,
) -> Dict[str, object]:
    require_ortools()

    model = cp_model.CpModel()

    x = {}
    for row in candidate_nodes.itertuples(index=False):
        x[row.opt_idx] = model.NewBoolVar(f"x_{row.opt_idx}")

    model.Add(sum(x[i] for i in x) == config.target_genes)

    z = {}
    for edge_id, row in enumerate(candidate_edges.itertuples(index=False)):
        z[edge_id] = model.NewBoolVar(f"z_{edge_id}")
        model.Add(z[edge_id] <= x[row.u_opt])
        model.Add(z[edge_id] <= x[row.v_opt])

    node_term = sum(
        int(row.node_obj) * x[row.opt_idx]
        for row in candidate_nodes.itertuples(index=False)
    )
    edge_term = sum(
        int(row.edge_obj) * z[edge_id]
        for edge_id, row in enumerate(candidate_edges.itertuples(index=False))
    )
    model.Maximize(node_term + edge_term)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = config.max_time_seconds
    solver.parameters.num_search_workers = config.num_workers
    solver.parameters.random_seed = config.random_seed

    status = solver.Solve(model)
    status_name_map = {
        cp_model.OPTIMAL: "OPTIMAL",
        cp_model.FEASIBLE: "FEASIBLE",
        cp_model.INFEASIBLE: "INFEASIBLE",
        cp_model.MODEL_INVALID: "MODEL_INVALID",
        cp_model.UNKNOWN: "UNKNOWN",
    }
    status_name = status_name_map.get(status, str(status))
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        raise RuntimeError(f"Subgraph solver failed with status: {status_name}")

    selected_opt_nodes = [idx for idx in x if solver.Value(x[idx]) == 1]
    selected_edge_ids = [eid for eid in z if solver.Value(z[eid]) == 1]

    selected_nodes = candidate_nodes.loc[candidate_nodes["opt_idx"].isin(selected_opt_nodes)].copy()
    selected_edges = (
        candidate_edges.iloc[selected_edge_ids].copy()
        if selected_edge_ids
        else candidate_edges.iloc[[]].copy()
    )

    return {
        "solver_status": status_name,
        "objective_value": solver.ObjectiveValue(),
        "best_objective_bound": solver.BestObjectiveBound(),
        "selected_nodes": selected_nodes,
        "selected_edges": selected_edges,
    }


def build_summary(
    candidate_nodes: pd.DataFrame,
    candidate_edges: pd.DataFrame,
    selected_nodes: pd.DataFrame,
    selected_edges: pd.DataFrame,
    solve_result: Dict[str, object],
    config: SubgraphConfig,
) -> pd.DataFrame:
    n_selected_edges = selected_edges.shape[0]
    return pd.DataFrame(
        [
            {
                "solver_status": solve_result["solver_status"],
                "candidate_pool_size": config.candidate_pool_size,
                "candidate_edges_in_model": candidate_edges.shape[0],
                "target_genes": config.target_genes,
                "n_selected_nodes": selected_nodes.shape[0],
                "n_selected_edges": n_selected_edges,
                "sum_selected_prize": selected_nodes["prize"].sum(),
                "sum_selected_node_obj": selected_nodes["node_obj"].sum(),
                "sum_selected_edge_reward_raw": selected_edges["edge_reward_raw"].sum() if n_selected_edges else 0.0,
                "sum_selected_edge_obj": selected_edges["edge_obj"].sum() if n_selected_edges else 0,
                "mean_selected_edge_confidence": selected_edges["confidence_raw"].mean() if n_selected_edges else np.nan,
                "mean_selected_edge_cost": selected_edges["edge_cost_raw"].mean() if n_selected_edges else np.nan,
                "objective_value": solve_result["objective_value"],
                "best_objective_bound": solve_result["best_objective_bound"],
                "max_time_seconds": config.max_time_seconds,
                "num_workers": config.num_workers,
                "node_prize_weight": config.node_prize_weight,
                "edge_conf_weight": config.edge_conf_weight,
                "edge_cost_weight": config.edge_cost_weight,
                "node_scale": config.node_scale,
                "edge_scale": config.edge_scale,
                "prize_column": config.prize_column,
            }
        ]
    )


def make_selected_node_output(selected_nodes: pd.DataFrame) -> pd.DataFrame:
    keep = [
        col
        for col in [
            "opt_idx",
            "node_id",
            "gene_name",
            "prize",
            "node_obj",
            "prediff_integrated",
            "prediff_vulnerability",
            "prediff_norm",
            "post_integrated",
            "post_vulnerability",
            "post_norm",
            "post_germline",
            "post_somatic",
            "post_epigenomic",
            "rank_shift",
        ]
        if col in selected_nodes.columns
    ]
    return selected_nodes.loc[:, keep].sort_values(["prize", "gene_name"], ascending=[False, True]).reset_index(drop=True)


def make_selected_edge_output(selected_edges: pd.DataFrame) -> pd.DataFrame:
    work = selected_edges.copy()
    work = work.rename(columns={"from": "gene_u", "to": "gene_v"})
    keep = [
        col
        for col in [
            "u_opt",
            "v_opt",
            "gene_u",
            "gene_v",
            "confidence_raw",
            "edge_cost_raw",
            "edge_reward_raw",
            "edge_obj",
            "n_protein_edges",
        ]
        if col in work.columns
    ]
    return work.loc[:, keep].sort_values(
        [col for col in ["edge_reward_raw", "confidence_raw", "edge_cost_raw"] if col in keep],
        ascending=[False, False, True][: len([col for col in ["edge_reward_raw", "confidence_raw", "edge_cost_raw"] if col in keep])],
    ).reset_index(drop=True)


def _graphml_key(parent: ET.Element, key_id: str, domain: str, name: str, attr_type: str) -> None:
    ET.SubElement(
        parent,
        "key",
        id=key_id,
        **{"for": domain, "attr.name": name, "attr.type": attr_type},
    )


def _graphml_type_for_series(series: pd.Series) -> str:
    if pd.api.types.is_integer_dtype(series):
        return "long"
    if pd.api.types.is_float_dtype(series):
        return "double"
    return "string"


def write_graphml(
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
    path: str,
) -> str:
    graphml = ET.Element(
        "graphml",
        xmlns="http://graphml.graphdrawing.org/xmlns",
        **{"xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance"},
    )

    node_cols = [c for c in nodes_df.columns if c != "node_id"]
    edge_cols = [c for c in edges_df.columns if c not in {"gene_u", "gene_v"}]

    for col in node_cols:
        _graphml_key(graphml, f"node_{col}", "node", col, _graphml_type_for_series(nodes_df[col]))
    for col in edge_cols:
        _graphml_key(graphml, f"edge_{col}", "edge", col, _graphml_type_for_series(edges_df[col]))

    graph = ET.SubElement(graphml, "graph", edgedefault="undirected", id="selected_subgraph")

    for row in nodes_df.itertuples(index=False):
        node = ET.SubElement(graph, "node", id=str(row.node_id))
        row_dict = row._asdict()
        for col in node_cols:
            value = row_dict.get(col)
            if pd.isna(value):
                continue
            data = ET.SubElement(node, "data", key=f"node_{col}")
            data.text = str(value)

    for edge_idx, row in enumerate(edges_df.itertuples(index=False)):
        edge = ET.SubElement(graph, "edge", id=f"e{edge_idx}", source=str(row.gene_u), target=str(row.gene_v))
        row_dict = row._asdict()
        for col in edge_cols:
            value = row_dict.get(col)
            if pd.isna(value):
                continue
            data = ET.SubElement(edge, "data", key=f"edge_{col}")
            data.text = str(value)

    tree = ET.ElementTree(graphml)
    ensure_dir(os.path.dirname(path))
    tree.write(path, encoding="utf-8", xml_declaration=True)
    return path


def save_selected_subgraph_outputs(
    selected_nodes: pd.DataFrame,
    selected_edges: pd.DataFrame,
    summary: pd.DataFrame,
    candidates: pd.DataFrame,
    config: SubgraphConfig,
) -> Dict[str, str]:
    ensure_dir(config.output_dir)

    nodes_path = os.path.join(config.output_dir, f"{config.output_stem}_nodes.tsv")
    edges_path = os.path.join(config.output_dir, f"{config.output_stem}_edges.tsv")
    summary_path = os.path.join(config.output_dir, f"{config.output_stem}_summary.tsv")
    candidates_path = os.path.join(config.output_dir, f"{config.output_stem}_candidate_nodes.tsv")
    graphml_path = os.path.join(config.output_dir, f"{config.output_stem}.graphml")

    selected_nodes.to_csv(nodes_path, sep="\t", index=False)
    selected_edges.to_csv(edges_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    candidates.to_csv(candidates_path, sep="\t", index=False)
    write_graphml(selected_nodes, selected_edges, graphml_path)

    return {
        "nodes_path": nodes_path,
        "edges_path": edges_path,
        "summary_path": summary_path,
        "candidates_path": candidates_path,
        "graphml_path": graphml_path,
    }


def run_cardinality_subgraph_calling(
    config: SubgraphConfig = SubgraphConfig(),
) -> Dict[str, object]:
    diffusion = validate_diffusion_results(
        read_diffusion_results(config.diffusion_path),
        prize_column=config.prize_column,
    )
    gg_nodes = validate_gene_gene_nodes(read_gene_gene_nodes(config.gg_nodes_path))
    gg_edges = validate_gene_gene_edges(
        read_gene_gene_edges(config.gg_edges_path),
        confidence_column=config.confidence_column,
        edge_cost_column=config.edge_cost_column,
    )

    diffusion = compute_node_prizes(diffusion, prize_column=config.prize_column)
    scored_nodes = attach_diffusion_prizes_to_gene_graph(gg_nodes, diffusion)
    candidate_nodes = build_candidate_pool(
        scored_nodes=scored_nodes,
        target_genes=config.target_genes,
        candidate_pool_size=config.candidate_pool_size,
    )
    candidate_edges = restrict_edges_to_candidates(
        edges=gg_edges,
        candidate_node_ids=candidate_nodes["node_id"],
        config=config,
    )
    candidate_nodes, candidate_edges = reindex_candidate_problem(candidate_nodes, candidate_edges)
    candidate_nodes, candidate_edges = integerize_objective_terms(
        candidate_nodes=candidate_nodes,
        candidate_edges=candidate_edges,
        config=config,
    )

    solve_result = solve_cardinality_subgraph(
        candidate_nodes=candidate_nodes,
        candidate_edges=candidate_edges,
        config=config,
    )

    selected_nodes = make_selected_node_output(solve_result["selected_nodes"])
    selected_edges = make_selected_edge_output(solve_result["selected_edges"])
    summary = build_summary(
        candidate_nodes=candidate_nodes,
        candidate_edges=candidate_edges,
        selected_nodes=solve_result["selected_nodes"],
        selected_edges=solve_result["selected_edges"],
        solve_result=solve_result,
        config=config,
    )
    output_paths = save_selected_subgraph_outputs(
        selected_nodes=selected_nodes,
        selected_edges=selected_edges,
        summary=summary,
        candidates=make_selected_node_output(candidate_nodes),
        config=config,
    )

    return {
        "diffusion": diffusion,
        "gene_gene_nodes": gg_nodes,
        "gene_gene_edges": gg_edges,
        "candidate_nodes": candidate_nodes,
        "candidate_edges": candidate_edges,
        "selected_nodes": selected_nodes,
        "selected_edges": selected_edges,
        "summary": summary,
        "output_paths": output_paths,
        "config": config,
    }


def main() -> None:
    result = run_cardinality_subgraph_calling()
    selected_nodes = result["selected_nodes"]
    selected_edges = result["selected_edges"]

    print("\n=== Cardinality-constrained subgraph calling complete ===")
    print(f"Selected genes: {selected_nodes.shape[0]}")
    print(f"Selected edges: {selected_edges.shape[0]}")
    print("\nTop selected genes:")
    print(selected_nodes.head(20).to_string(index=False))
    print("\nSummary:")
    print(result["summary"].to_string(index=False))
    print("\nSaved outputs:")
    for label, path in result["output_paths"].items():
        print(f"  {label}: {path}")


if __name__ == "__main__":
    main()
