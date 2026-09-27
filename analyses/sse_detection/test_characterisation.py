"""Regression tests for descriptive weighting, eligibility and missingness."""
import unittest

import numpy as np
import pandas as pd

from .lib.characterisation import ATTRIBUTES, ENTROPY_ATTRIBUTES, build_summaries, candidate_mask, select_groups, validate_summaries


def fixture():
    nodes, records = [], []
    for i, (size, tier, burst, burden) in enumerate([
        (6, "high_priority_burst", .02, np.nan),
        (10, "high_priority_burden", .7, .03),
        (20, "possible_review", .08, .7),
        (6, "background_or_low_information", .7, np.nan),
        (5, "size_ineligible", .001, np.nan),
    ]):
        row = dict(cluster_id=str(i), cluster_size=size, candidate_tier=tier, sse_tested=size>=6,
            burst_score_upper_p=burst, burden_score_upper_p=burden, burden_eligible=np.isfinite(burden))
        for key, _, _ in ENTROPY_ATTRIBUTES:
            row[key+"_entropy_obs"] = .5
            row[key+"_entropy_z"] = -.5
        nodes.append(row)
        for j in range(size):
            record = dict(cluster_id=str(i), sequence_id=f"s{j}", window_id=f"w{i}")
            for _, col, _ in ATTRIBUTES:
                record[col] = None if i == 3 else ("A" if i == 0 or j < size//2 else "B")
            if i == 1 and j == 0:
                record["sex"] = None
            records.append(record)
    return pd.DataFrame(nodes), pd.DataFrame(records)


class CharacterisationTests(unittest.TestCase):
    def test_weighting_and_missing_denominators(self):
        nodes, records = fixture()
        tables = build_summaries(nodes, records)
        validate_summaries(tables)
        comp = tables["composition"]
        row = comp.loc[comp.scenario.eq("primary") & comp.group.eq("CSC") & comp.size_band.eq("All") & comp.attribute.eq("sex") & comp.category.eq("A")].iloc[0]
        self.assertEqual(row.n_nonmissing_records, 15)
        self.assertEqual(row.n_missing_records, 1)
        self.assertAlmostEqual(row.record_proportion, 10/15)
        self.assertAlmostEqual(row.cluster_mean_proportion, (1+4/9)/2)
        bg = comp.loc[comp.scenario.eq("primary") & comp.group.eq("Background") & comp.size_band.eq("All") & comp.attribute.eq("sex")].iloc[0]
        self.assertEqual(bg.n_missing_clusters, 1)
        self.assertEqual(bg.n_nonmissing_clusters, 1)

    def test_fixed_floor_threshold_and_route_comparator(self):
        nodes, _ = fixture()
        self.assertEqual(int(candidate_mask(nodes, .025).sum()), 1)
        self.assertEqual(int(candidate_mask(nodes, .1).sum()), 3)
        groups = select_groups(nodes, "route")
        self.assertEqual(groups["Burden-eligible background"].cluster_id.tolist(), ["2"])
        restricted = select_groups(nodes, "size_10")
        self.assertEqual(restricted["CSC"].cluster_id.tolist(), ["1"])

    def test_duplicate_record_rejected_but_cross_window_recurrence_allowed(self):
        nodes, records = fixture()
        build_summaries(nodes, records)
        with self.assertRaisesRegex(ValueError, "one record per sequence"):
            build_summaries(nodes, pd.concat([records, records.iloc[:1]]))

    def test_size_or_label_mismatch_rejected(self):
        nodes, records = fixture()
        with self.assertRaisesRegex(ValueError, "cluster sizes"):
            build_summaries(nodes, records.iloc[1:])
        nodes.loc[0, "candidate_tier"] = "background_or_low_information"
        with self.assertRaisesRegex(ValueError, "primary labels"):
            build_summaries(nodes, records)


if __name__ == "__main__":
    unittest.main()
