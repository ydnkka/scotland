#!/usr/bin/env python3
"""
Superspreading Signature Detection
==================================

Run the SSE detection pipeline: load analysis rows, build transition edges,
assemble cluster-node features, score SSE candidates, validate the regression
schema, and write detector plus transition-summary outputs.

From the repository root:
    ``python -m chapter_analyses.sse_detection.lib.sse.detection``
"""

from __future__ import annotations

import logging
from pathlib import Path
import sys

import pandas as pd


_PROJECT_ROOT = Path(__file__).resolve().parents[4]

if __package__ in {None, ""} and str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from utils import load_analysis_columns  # noqa: E402


from .cluster_features import (  # noqa: E402
    build_cluster_attributes,
    build_cluster_stats,
    build_cluster_table,
)
from .config import (  # noqa: E402
    ANALYSIS_COLUMNS,
    SSE_OUTPUT_DIR,
    TRANSITION_WINDOW_STRIDE,
)
from ..model.prep import validate_regression_cluster_columns  # noqa: E402
from .scoring import add_sse_node_metrics  # noqa: E402
from .transition_graph import (  # noqa: E402
    TransitionSummaryOutputs,
    build_transition_network,
    build_transition_summary_outputs,
)


LOGGER = logging.getLogger(__name__)


def load_sequence_data() -> pd.DataFrame:
    """Load sequence-window rows for the retained transition-window stride."""
    return load_analysis_columns(
        ANALYSIS_COLUMNS,
        add_policy=True,
        window_stride=TRANSITION_WINDOW_STRIDE,
    )


def write_sse_output_tables(
    cluster_table: pd.DataFrame,
    edge_table: pd.DataFrame,
    transition_outputs: TransitionSummaryOutputs,
    output_dir: Path,
) -> None:
    """Write detector outputs plus transition-summary tables."""
    output_dir.mkdir(parents=True, exist_ok=True)
    cluster_table.to_parquet(output_dir / "cluster_table.parquet")
    edge_table.to_parquet(output_dir / "edge_table.parquet")

    csv_tables = {
        "transition_graph_summary",
        "transition_window_summary",
        "transition_component_summary",
    }
    for name, table in transition_outputs.tables().items():
        table.to_parquet(output_dir / f"{name}.parquet")
        if name in csv_tables:
            table.to_csv(output_dir / f"{name}.csv", index=False)


def main() -> None:
    LOGGER.info("...........loading data")
    df = load_sequence_data()

    LOGGER.info("...........building transition graph")
    edge_table, _, component_map = build_transition_network(df)

    LOGGER.info("...........building cluster features")
    cluster_stats = build_cluster_stats(df)
    cluster_att = build_cluster_attributes(df)
    cluster_att["connected_components"] = cluster_att["cluster_id"].map(component_map)
    cluster_table = build_cluster_table(
        cluster_att,
        cluster_stats,
        edge_table,
        sequence_df=df,
    )

    LOGGER.info("...........scoring SSE candidates")
    sse_df = add_sse_node_metrics(cluster_table)
    validate_regression_cluster_columns(sse_df)

    LOGGER.info("...........building transition summaries")
    transition_outputs = build_transition_summary_outputs(
        sse_df,
        edge_table,
    )
    write_sse_output_tables(sse_df, edge_table, transition_outputs, SSE_OUTPUT_DIR)

    LOGGER.info("Done.")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    main()
