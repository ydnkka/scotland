"""Regression tests for selecting individual publication artifacts."""

import unittest
import importlib
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from analyses.sse_detection.lib.figs import fig22, tab09
from analyses.surveillance.lib.figs import tab01
from utils.style import add_panel_labels

from . import registry
from results.migrate_asset_names import migration_plan


class IndividualArtifactTests(unittest.TestCase):
    def test_panel_labels_handle_grids_and_explicit_text_styles(self):
        fig, axes = plt.subplots(2, 2)
        add_panel_labels(axes, kwargs={"fontweight": "normal", "fontsize": 11})
        labels = [ax.texts[-1] for ax in axes.flat]
        self.assertEqual([text.get_text() for text in labels], list("ABCD"))
        self.assertTrue(all(text.get_fontsize() == 11 for text in labels))
        self.assertTrue(all(text.get_fontweight() == "normal" for text in labels))
        with self.assertRaisesRegex(ValueError, "one panel label"):
            add_panel_labels(axes, label=["A"])
        plt.close(fig)

    def test_size_figure_selection_writes_only_the_requested_figure(self):
        cumulative = pd.DataFrame(
            {
                "group": ["CSC", "CSC", "Background", "Background"],
                "cluster_size": [6, 10, 6, 20],
                "cumulative_clusters": [1, 2, 1, 2],
                "cluster_percent": [50, 100, 50, 100],
                "cumulative_records": [6, 16, 6, 26],
                "record_percent": [37.5, 100, 23.08, 100],
            }
        )
        with TemporaryDirectory() as tmp, patch.object(
            fig22, "prepare", return_value={"size_cumulative": cumulative}
        ):
            figures = Path(tmp) / "figures"
            result = registry.build_figures(
                names=["sse_detection:fig22_size_cumulative"],
                figure_dir=figures,
                table_dir=Path(tmp) / "tables",
            )
            self.assertEqual(list(result), ["sse_detection:fig22_size_cumulative"])
            self.assertEqual(
                {p.name for p in figures.iterdir()},
                {"fig22_size_cumulative.png", "fig22_size_cumulative.pdf"},
            )
            self.assertEqual(list((Path(tmp) / "tables").iterdir()), [])
        plt.close("all")

    def test_single_size_table_respects_custom_publication_directory(self):
        sizes = pd.DataFrame(
            {
                "scenario": ["route", "route"],
                "group": ["Background", "Burst only"],
                "n_clusters": [2, 1],
                "n_records": [18, 6],
                "size_median": [9, 6],
                "size_q25": [7, 6],
                "size_q75": [11, 6],
            }
        )
        with TemporaryDirectory() as tmp, patch.object(
            tab09, "prepare", return_value={"size_summary": sizes}
        ):
            output = Path(tmp) / "custom-publication-tables"
            registry.build_tables(
                names=["sse_detection:tab09_route_sizes"],
                figure_dir=Path(tmp) / "figures",
                table_dir=output,
            )
            self.assertEqual(
                [p.name for p in output.iterdir()], ["tab09_route_sizes.tex"]
            )
            text = (output / "tab09_route_sizes.tex").read_text()
            self.assertIn("Background & 2 & 18", text)
            self.assertIn("Burst only & 1 & 6", text)
            self.assertFalse((output.parent / "tables").exists())

    def test_surveillance_table_selection_does_not_build_other_tables(self):
        sequences = pd.DataFrame(
            {
                "sequence_id": ["one", "two", "three"],
                "collection_date": pd.to_datetime(["2021-01-04"] * 3),
                "variant": ["Alpha", "Alpha", "Delta"],
            }
        )
        with TemporaryDirectory() as tmp, patch.object(
            tab01, "load_sequences", return_value=sequences
        ), patch.object(
            registry, "SURVEILLANCE_TABLES_DIR", Path(tmp) / "source-tables"
        ):
            registry.build_tables(
                names=["surveillance:tab01_clade_frequency_by_period"],
                figure_dir=Path(tmp) / "figures",
                table_dir=Path(tmp) / "publication-tables",
            )
            files = list((Path(tmp) / "source-tables").iterdir())
            self.assertEqual({p.suffix for p in files}, {".csv", ".parquet"})
            self.assertEqual({p.stem for p in files}, {"tab01_clade_frequency_by_period"})
            table = pd.read_csv(
                Path(tmp) / "source-tables/tab01_clade_frequency_by_period.csv"
            )
            self.assertAlmostEqual(table.loc[0, "Alpha"], 2 / 3)
            self.assertAlmostEqual(table.loc[0, "Delta"], 1 / 3)
            self.assertEqual(list((Path(tmp) / "figures").iterdir()), [])

    def test_every_export_matches_the_physical_producer_script(self):
        builders = (*registry.figure_builders(), *registry.table_builders())
        for builder in builders:
            module = importlib.import_module(builder.module)
            script = Path(module.__file__).stem
            self.assertTrue(builder.name.startswith(script + "_"), builder.key)
            constant = module.FIGURE_NAME if builder.kind == "figure" else module.TABLE_NAME
            self.assertEqual(builder.name, constant)
        assets = registry.publication_assets()
        self.assertEqual(len(assets), len({asset["stem"] for asset in assets}))
        self.assertFalse(any(asset["domain"] == "surveillance" and asset["kind"] == "table" for asset in assets))

    def test_migration_rejects_conflicts_before_any_export_is_moved(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp) / "figures"
            directory.mkdir()
            assets = [
                {"directory": "figures", "formats": ["pdf"],
                 "legacy_stem": "fig_one", "stem": "fig01_one"},
                {"directory": "figures", "formats": ["pdf"],
                 "legacy_stem": "fig_two", "stem": "fig02_two"},
            ]
            (directory / "fig_one.pdf").write_bytes(b"original one")
            (directory / "fig_two.pdf").write_bytes(b"original two")
            (directory / "fig02_two.pdf").write_bytes(b"different authored export")
            with self.assertRaisesRegex(FileExistsError, "Conflicting"):
                migration_plan(Path(tmp), assets)
            self.assertTrue((directory / "fig_one.pdf").is_file())
            self.assertFalse((directory / "fig01_one.pdf").exists())

    def test_migration_accepts_identical_duplicates_and_already_numbered_files(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp) / "tables"
            directory.mkdir()
            asset = {"directory": "tables", "formats": ["tex"],
                     "legacy_stem": "tab_one", "stem": "tab01_one"}
            old, new = directory / "tab_one.tex", directory / "tab01_one.tex"
            old.write_text("same table")
            new.write_text("same table")
            self.assertEqual(migration_plan(Path(tmp), [asset]), [(old, new)])
            old.unlink()
            self.assertEqual(migration_plan(Path(tmp), [asset]), [])


if __name__ == "__main__":
    unittest.main()
