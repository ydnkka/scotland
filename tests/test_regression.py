import pandas as pd

from sse_detection.lib import regression


def test_term_matching_uses_full_patsy_factor_tokens():
    data = pd.DataFrame(
        {
            "y": [0, 0, 0, 1, 1, 1, 0, 1],
            "age": [20, 30, 40, 50, 60, 70, 35, 55],
            "age_band": [
                "20-29",
                "30-39",
                "40-49",
                "50-59",
                "60-69",
                "70-79",
                "30-39",
                "50-59",
            ],
        }
    )
    result = regression.fit_binomial_glm(data, "y ~ age + C(age_band)")

    assert regression.parameter_names_for_term(result, "age") == ["age"]
    assert all(
        "age_band" in term
        for term in regression.parameter_names_for_term(result, "C(age_band)")
    )

    age_or = regression.tidy_odds_ratios(result, term_filter="age")
    assert age_or["term"].tolist() == ["age"]

    wald = regression.robust_wald_for_prefix(result, "age")
    assert int(wald.loc[0, "df"]) == 1


def test_firth_logit_exposes_statsmodels_like_term_metadata():
    data = pd.DataFrame(
        {
            "y": [0, 0, 0, 1, 1, 1, 0, 1],
            "window_idx": [1, 1, 2, 2, 3, 3, 4, 4],
            "age": [20, 30, 40, 50, 60, 70, 35, 55],
            "age_band": [
                "20-29",
                "30-39",
                "40-49",
                "50-59",
                "60-69",
                "70-79",
                "30-39",
                "50-59",
            ],
        }
    )

    result = regression.fit_firth_logit(
        data,
        "y ~ age + C(window_idx) + C(age_band)",
        maxiter=50,
    )

    assert "age" in result.params.index
    assert regression.parameter_names_for_term(result, "age") == ["age"]
    assert regression.tidy_odds_ratios(result, term_filter="age")["term"].tolist() == ["age"]
