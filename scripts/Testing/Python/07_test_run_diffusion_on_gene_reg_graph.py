#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import sys
import uuid
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
DIFFUSION_SCRIPT = REPO_ROOT / "scripts" / "Internals" / "Python" / "07_run_diffusion_on_gene_reg_graph.py"
TEST_OUTPUT_ROOT = REPO_ROOT / "data" / "processed" / "test_outputs" / "python_step7"


def default_backend_dir() -> Path:
    return REPO_ROOT / "data" / "processed"


def resolve_no_score_graph_paths() -> tuple[Path, Path]:
    candidates = [
        (
            default_backend_dir() / "gene_reg_graph_no_scores_nodes.tsv.gz",
            default_backend_dir() / "gene_reg_graph_no_scores_edges.tsv.gz",
        ),
        (
            REPO_ROOT / "inst" / "extdata" / "backend" / "gene_reg_graph_no_scores_nodes.tsv.gz",
            REPO_ROOT / "inst" / "extdata" / "backend" / "gene_reg_graph_no_scores_edges.tsv.gz",
        ),
    ]

    for nodes_path, edges_path in candidates:
        if nodes_path.exists() and edges_path.exists():
            return nodes_path, edges_path

    raise FileNotFoundError(
        "Could not locate a materialized no-score gene-reg graph. "
        "Run initialize_backend_graphs(build_gene_reg=TRUE) first."
    )


def load_diffusion_module():
    spec = importlib.util.spec_from_file_location("conseguiR_step7", DIFFUSION_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def make_repo_test_dir(prefix: str) -> Path:
    TEST_OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    path = TEST_OUTPUT_ROOT / f"{prefix}_{uuid.uuid4().hex}"
    path.mkdir(parents=True, exist_ok=False)
    return path


def make_scored_graph_fixture(tmpdir: Path) -> tuple[Path, Path]:
    no_score_nodes, no_score_edges = resolve_no_score_graph_paths()
    nodes = pd.read_csv(no_score_nodes, sep="\t", low_memory=False)
    edges = pd.read_csv(no_score_edges, sep="\t", low_memory=False)

    nodes["somatic_score"] = 0.0
    nodes["germline_score"] = 0.0
    nodes["epigenomic_score"] = 0.0

    gene_mask = nodes["node_type"] == "gene"
    reg_mask = nodes["node_type"] == "reg"

    gene_ids = nodes.loc[gene_mask, "node_id"].head(3).tolist()
    reg_ids = nodes.loc[reg_mask, "node_id"].head(3).tolist()

    nodes.loc[nodes["node_id"] == gene_ids[0], ["somatic_score", "germline_score"]] = [2.1, 1.4]
    nodes.loc[nodes["node_id"] == gene_ids[1], ["somatic_score", "germline_score"]] = [-1.2, 0.8]
    nodes.loc[nodes["node_id"] == gene_ids[2], ["somatic_score", "germline_score"]] = [0.6, -0.5]

    nodes.loc[nodes["node_id"] == reg_ids[0], ["somatic_score", "germline_score", "epigenomic_score"]] = [1.1, -0.4, 3.0]
    nodes.loc[nodes["node_id"] == reg_ids[1], ["somatic_score", "germline_score", "epigenomic_score"]] = [-0.7, 0.9, -2.0]
    nodes.loc[nodes["node_id"] == reg_ids[2], ["somatic_score", "germline_score", "epigenomic_score"]] = [0.4, 1.2, 0.5]

    scored_nodes_path = tmpdir / "gene_reg_graph_scored_nodes.tsv.gz"
    scored_edges_path = tmpdir / "gene_reg_graph_scored_edges.tsv.gz"

    nodes.to_csv(scored_nodes_path, sep="\t", index=False, compression="gzip")
    edges.to_csv(scored_edges_path, sep="\t", index=False, compression="gzip")
    return scored_nodes_path, scored_edges_path


def test_run_gene_reg_diffusion():
    module = load_diffusion_module()

    tmpdir = make_repo_test_dir("conseguiR_step7_test")
    nodes_path, edges_path = make_scored_graph_fixture(tmpdir)

    config = module.DiffusionConfig(
        nodes_path=str(nodes_path),
        edges_path=str(edges_path),
        output_dir=str(tmpdir),
        output_stem="gene_reg_graph_diffusion_test",
        top_k=3,
        top_n_to_save=10,
    )

    result = module.run_gene_reg_diffusion(config=config)
    diffusion = result["diffusion"]

    required_cols = {
        "node_id",
        "gene_name",
        "prediff_germline",
        "prediff_somatic",
        "prediff_epigenomic",
        "incoming_germline",
        "incoming_somatic",
        "incoming_epigenomic",
        "post_germline",
        "post_somatic",
        "post_epigenomic",
        "post_norm",
        "rank_shift",
    }

    assert required_cols.issubset(diffusion.columns), "Diffusion output is missing expected columns."
    assert len(diffusion) > 0, "Diffusion output is empty."
    assert diffusion["post_norm"].notna().all(), "post_norm contains missing values."

    all_genes_path = Path(result["output_paths"]["all_genes_path"])
    top_genes_path = Path(result["output_paths"]["top_genes_path"])
    assert all_genes_path.exists(), "All-genes diffusion output was not written."
    assert top_genes_path.exists(), "Top-genes diffusion output was not written."

    print("Top diffused genes from scored graph fixture:")
    print(
        diffusion.loc[:, ["node_id", "gene_name", "prediff_norm", "post_norm", "rank_shift"]]
        .head(10)
        .to_string(index=False)
    )


def main():
    test_run_gene_reg_diffusion()
    print("Step 7 diffusion test passed.")


if __name__ == "__main__":
    main()
