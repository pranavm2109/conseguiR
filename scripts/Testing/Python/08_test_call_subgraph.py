#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
SUBGRAPH_SCRIPT = REPO_ROOT / "scripts" / "Internals" / "Python" / "08_call_subgraph.py"
GG_NODES = REPO_ROOT / "data" / "processed" / "gene_gene_graph_nodes.tsv.gz"
GG_EDGES = REPO_ROOT / "data" / "processed" / "gene_gene_graph_edges.tsv.gz"


def load_subgraph_module():
    spec = importlib.util.spec_from_file_location("conseguiR_step8", SUBGRAPH_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def ortools_available(module) -> bool:
    return getattr(module, "cp_model", None) is not None


def make_diffusion_fixture(tmpdir: Path, n_genes: int = 80) -> Path:
    nodes = pd.read_csv(GG_NODES, sep="\t", low_memory=False)
    edges = pd.read_csv(GG_EDGES, sep="\t", low_memory=False)

    edge_node_ids = pd.unique(pd.concat([edges["from"], edges["to"]], ignore_index=True))
    candidate_ids = list(edge_node_ids[:n_genes])

    node_lookup = nodes.set_index("node_id")
    candidate_ids = [node_id for node_id in candidate_ids if node_id in node_lookup.index]

    diffusion = pd.DataFrame(
        {
            "node_id": candidate_ids,
            "gene_name": [node_lookup.loc[node_id, "name"] for node_id in candidate_ids],
            "post_germline": [5.0 - (idx * 0.03) for idx in range(len(candidate_ids))],
            "post_somatic": [3.5 - (idx * 0.02) for idx in range(len(candidate_ids))],
            "post_epigenomic": [2.5 - (idx * 0.015) for idx in range(len(candidate_ids))],
            "prediff_norm": [1.0 + (idx * 0.01) for idx in range(len(candidate_ids))],
            "rank_shift": [0.0 for _ in candidate_ids],
        }
    )
    diffusion["post_norm"] = (
        diffusion["post_germline"] ** 2
        + diffusion["post_somatic"] ** 2
        + diffusion["post_epigenomic"] ** 2
    ) ** 0.5

    diffusion_path = tmpdir / "gene_reg_graph_diffusion_all_genes.tsv"
    diffusion.to_csv(diffusion_path, sep="\t", index=False)
    return diffusion_path


def test_run_cardinality_subgraph_calling_live():
    module = load_subgraph_module()

    if not ortools_available(module):
        print(
            "Skipping live step 8 optimization test because OR-Tools is not installed. "
            "Run this test from `lymphoma_graph_env` to exercise the solver path."
        )
        return

    with tempfile.TemporaryDirectory(prefix="conseguiR_step8_test_") as tmp:
        tmpdir = Path(tmp)
        diffusion_path = make_diffusion_fixture(tmpdir, n_genes=80)

        config = module.SubgraphConfig(
            diffusion_path=str(diffusion_path),
            gg_nodes_path=str(GG_NODES),
            gg_edges_path=str(GG_EDGES),
            output_dir=str(tmpdir),
            output_stem="gene_gene_selected_subgraph_test",
            target_genes=12,
            candidate_pool_size=40,
            max_edges_in_model=500,
            max_time_seconds=30,
            num_workers=4,
        )

        result = module.run_cardinality_subgraph_calling(config=config)
        selected_nodes = result["selected_nodes"]
        selected_edges = result["selected_edges"]

        assert selected_nodes.shape[0] == 12, "Selected node cardinality does not match the request."
        assert {"node_id", "gene_name", "prize"}.issubset(selected_nodes.columns)
        assert {"gene_u", "gene_v"}.issubset(selected_edges.columns)
        assert Path(result["output_paths"]["nodes_path"]).exists()
        assert Path(result["output_paths"]["edges_path"]).exists()
        assert Path(result["output_paths"]["summary_path"]).exists()
        assert Path(result["output_paths"]["graphml_path"]).exists()

        print("Selected genes from step 8 test fixture:")
        print(selected_nodes.head(20).to_string(index=False))


def test_validate_diffusion_results_negative_missing_columns():
    module = load_subgraph_module()
    bad = pd.DataFrame({"node_id": ["TP53"], "post_norm": [1.0]})

    try:
        module.validate_diffusion_results(bad, prize_column="post_norm")
    except ValueError as exc:
        assert "missing columns" in str(exc)
    else:
        raise AssertionError("Expected diffusion validation to fail for missing gene_name.")


def test_run_cardinality_subgraph_calling_negative_missing_ortools():
    module = load_subgraph_module()
    if ortools_available(module):
        return

    with tempfile.TemporaryDirectory(prefix="conseguiR_step8_no_ortools_") as tmp:
        tmpdir = Path(tmp)
        diffusion_path = make_diffusion_fixture(tmpdir, n_genes=40)

        config = module.SubgraphConfig(
            diffusion_path=str(diffusion_path),
            gg_nodes_path=str(GG_NODES),
            gg_edges_path=str(GG_EDGES),
            output_dir=str(tmpdir),
            target_genes=10,
            candidate_pool_size=20,
        )

        try:
            module.run_cardinality_subgraph_calling(config=config)
        except ImportError as exc:
            assert "OR-Tools is required" in str(exc)
        else:
            raise AssertionError("Expected step 8 to fail clearly when OR-Tools is unavailable.")


def test_run_cardinality_subgraph_calling_negative_no_overlap():
    module = load_subgraph_module()

    with tempfile.TemporaryDirectory(prefix="conseguiR_step8_bad_overlap_") as tmp:
        tmpdir = Path(tmp)
        diffusion = pd.DataFrame(
            {
                "node_id": [f"NOT_A_REAL_GENE_{i}" for i in range(20)],
                "gene_name": [f"NOT_A_REAL_GENE_{i}" for i in range(20)],
                "post_norm": [1.0 + i for i in range(20)],
            }
        )
        diffusion_path = tmpdir / "bad_diffusion.tsv"
        diffusion.to_csv(diffusion_path, sep="\t", index=False)

        config = module.SubgraphConfig(
            diffusion_path=str(diffusion_path),
            gg_nodes_path=str(GG_NODES),
            gg_edges_path=str(GG_EDGES),
            output_dir=str(tmpdir),
            target_genes=10,
            candidate_pool_size=20,
        )

        try:
            module.run_cardinality_subgraph_calling(config=config)
        except ValueError as exc:
            assert "could be matched" in str(exc)
        else:
            raise AssertionError("Expected step 8 to fail when diffusion genes do not overlap the graph.")


def test_run_cardinality_subgraph_calling_negative_bad_edges():
    module = load_subgraph_module()

    with tempfile.TemporaryDirectory(prefix="conseguiR_step8_bad_edges_") as tmp:
        tmpdir = Path(tmp)
        diffusion_path = make_diffusion_fixture(tmpdir, n_genes=40)

        bad_edges = pd.read_csv(GG_EDGES, sep="\t", low_memory=False).drop(columns=["confidence"])
        bad_edges_path = tmpdir / "bad_edges.tsv"
        bad_edges.to_csv(bad_edges_path, sep="\t", index=False)

        config = module.SubgraphConfig(
            diffusion_path=str(diffusion_path),
            gg_nodes_path=str(GG_NODES),
            gg_edges_path=str(bad_edges_path),
            output_dir=str(tmpdir),
            target_genes=10,
            candidate_pool_size=20,
        )

        try:
            module.run_cardinality_subgraph_calling(config=config)
        except ValueError as exc:
            assert "missing columns" in str(exc)
        else:
            raise AssertionError("Expected step 8 to fail for malformed gene-gene edges.")


def test_run_cardinality_subgraph_calling_negative_small_candidate_pool():
    module = load_subgraph_module()

    with tempfile.TemporaryDirectory(prefix="conseguiR_step8_bad_pool_") as tmp:
        tmpdir = Path(tmp)
        diffusion_path = make_diffusion_fixture(tmpdir, n_genes=30)

        config = module.SubgraphConfig(
            diffusion_path=str(diffusion_path),
            gg_nodes_path=str(GG_NODES),
            gg_edges_path=str(GG_EDGES),
            output_dir=str(tmpdir),
            target_genes=15,
            candidate_pool_size=10,
        )

        try:
            module.run_cardinality_subgraph_calling(config=config)
        except ValueError as exc:
            assert "at least as large" in str(exc)
        else:
            raise AssertionError("Expected step 8 to fail when candidate_pool_size < target_genes.")


def main():
    test_run_cardinality_subgraph_calling_live()
    test_validate_diffusion_results_negative_missing_columns()
    test_run_cardinality_subgraph_calling_negative_missing_ortools()
    test_run_cardinality_subgraph_calling_negative_no_overlap()
    test_run_cardinality_subgraph_calling_negative_bad_edges()
    test_run_cardinality_subgraph_calling_negative_small_candidate_pool()
    print("Step 8 subgraph-calling tests passed.")


if __name__ == "__main__":
    main()
