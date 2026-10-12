"""Check caption coverage, dynamic values and LaTeX rendering."""

import unittest

import yaml

from utils import (
    load_table_metadata,
    render_latex_caption,
    render_latex_grouped_column_table,
    render_latex_longtable,
    render_latex_table,
    table_caption,
)
from utils.table_captions import TABLE_CAPTIONS_PATH

from .registry import publication_assets


class TableCaptionTests(unittest.TestCase):
    def test_registry_covers_every_publication_table(self):
        text = TABLE_CAPTIONS_PATH.read_text(encoding="utf-8")
        registry = yaml.safe_load(text)
        # YAML normally accepts duplicate keys silently.
        self.assertEqual(len(yaml.compose(text).value), len(registry))
        sources = {
            asset["module"].replace(".", "/") + ".py"
            for asset in publication_assets()
            if asset["kind"] == "table"
        }
        self.assertEqual({entry["source"] for entry in registry.values()}, sources)
        self.assertEqual(len(registry), len(sources))
        root = TABLE_CAPTIONS_PATH.parents[1]
        for label in registry:
            with self.subTest(label=label):
                entry = load_table_metadata(label)
                self.assertTrue(entry["columns"])
                self.assertIn(label, (root / entry["source"]).read_text())
                caption = table_caption(label, n_summaries="12", n_pairs="1,234")
                self.assertNotIn("{{", caption["caption"])
                rendered = render_latex_caption(**caption)
                self.assertIn(r"\label{" + label + "}", rendered)
                for note in entry.get("notes", []):
                    self.assertNotIn(note, rendered)

    def test_dynamic_counts_and_missing_values(self):
        label = "tab:cluster_pairwise_distance_summary"
        caption = table_caption(label, n_summaries="12", n_pairs="1,234")
        self.assertIn("12 window-lineage strata", caption["caption"])
        self.assertIn("1,234 recorded pairs", caption["caption"])
        with self.assertRaisesRegex(ValueError, "n_pairs"):
            table_caption(label, n_summaries="12")
        with self.assertRaisesRegex(KeyError, "tab:missing"):
            table_caption("tab:missing")

    def test_all_writers_preserve_latex_and_escape_plain_text(self):
        cases = (
            ("tab:sequence_composition_by_policy", r"n (\%)", r"\textless{}5"),
            (
                "tab:bayesian_fixed_effects_focal",
                r"95\%",
                r"\ref{tab:bayesian_model_specifications}",
            ),
            ("tab:sensitivity_sizes", r"$p\leq0.05$", "six sequences"),
            ("tab:size_composition", r"$\Delta$", r"\ref{tab:size_overview}"),
        )
        for label, expected, extra in cases:
            caption = table_caption(label)
            outputs = (
                render_latex_table(**caption, columns=["A"], rows=[[1]]),
                render_latex_longtable(**caption, columns=["A"], rows=[[1]]),
                render_latex_grouped_column_table(
                    **caption,
                    row_columns=["A"],
                    column_groups=[("Group", ["B"])],
                    rows=[[1, 2]],
                    column_spec="lr",
                ),
            )
            for output in outputs:
                with self.subTest(label=label, output=output.splitlines()[0]):
                    self.assertIn(expected, output)
                    self.assertIn(extra, output)
                    self.assertEqual(output.count(r"\caption["), 1)
                    self.assertEqual(output.count(r"\label{" + label + "}"), 1)
                    self.assertNotIn(r"\textbackslash{}", output)

    def test_short_caption_is_plain_text_even_for_latex_captions(self):
        rendered = render_latex_caption(
            label="tab:example",
            caption=r"Effect $x$ with 95\% interval",
            short_caption="Effects & 95% intervals",
            caption_is_latex=True,
        )
        self.assertIn(r"\caption[Effects \& 95\% intervals]", rendered)
        self.assertIn(r"{Effect $x$ with 95\% interval}", rendered)


if __name__ == "__main__":
    unittest.main()
