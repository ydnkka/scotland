"""The diagnostics suggest that while there is some zero inflation,
it's not extreme enough to necessarily require a zero-inflated model.
The observed excess zeros are about 10% more than expected under the
fitted NegBin, which could be due to random variation or mild zero inflation.
Given this, we can proceed with a standard negative binomial  model for `n_sequences - 1`,
which will be more straightforward to interpret and implement in Bambi.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import nbinom
from scipy.optimize import minimize

from utils import data, style

def main(y: np.ndarray):
    """Diagnostic checks for zero inflation in a count data set."""

    # --- Fit NegBin to ALL data via MLE ---
    def negbin_nll(params, counts):
        mu, alpha = params
        if mu <= 0 or alpha <= 0:
            return np.inf
        # NegBin parameterised as (r, p) where r=alpha, p=alpha/(alpha+mu)
        r = alpha
        p = alpha / (alpha + mu)
        return -nbinom.logpmf(counts, r, p).sum()

    result = minimize(
        negbin_nll,
        x0=[y.mean(), 1.0],          # starting guess: mean, dispersion
        args=(y,),
        method='Nelder-Mead'
    )
    mu_fit, alpha_fit = result.x
    print(f"Fitted mu:    {mu_fit:.3f}")
    print(f"Fitted alpha: {alpha_fit:.3f}")

    # --- Expected vs observed zeros ---
    r = alpha_fit
    p = alpha_fit / (alpha_fit + mu_fit)

    p_zero_negbin   = nbinom.pmf(0, r, p)          # P(Y=0) under fitted NegBin
    expected_zeros  = p_zero_negbin * len(y)
    observed_zeros  = (y == 0).sum()

    print(f"\nObserved zeros:  {observed_zeros}  ({observed_zeros/len(y):.1%})")
    print(f"Expected zeros:  {expected_zeros:.1f}  ({p_zero_negbin:.1%})")
    print(f"Excess zeros:    {observed_zeros - expected_zeros:.1f}")
    print(f"\nZero inflation genuine: {observed_zeros > expected_zeros * 1.1}")

    # --- Visual check ---
    style.set_theme()
    fig, axes = style.new_figure(
        width="onehalf", nrows=1, ncols=2,
        sharey=True, constrained_layout=True)

    # Left: observed vs NegBin expected counts
    max_val = min(y.max(), 30)
    obs_counts    = np.array([(y == k).sum() for k in range(max_val + 1)])
    exp_counts    = np.array([nbinom.pmf(k, r, p) * len(y) for k in range(max_val + 1)])

    x = np.arange(max_val + 1)
    axes[0].bar(x - 0.2, obs_counts, width=0.4, label='Observed', alpha=0.7, color='steelblue')
    axes[0].bar(x + 0.2, exp_counts, width=0.4, label='NegBin expected', alpha=0.7, color='coral')
    axes[0].axvline(0.5, color='red', linestyle='--', alpha=0.5, label='Zero boundary')
    axes[0].set_xlabel('cluster size - 1')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Observed vs NegBin Expected')
    axes[0].legend()

    # Right: zoom on zero bar only
    categories = ['Observed zeros', 'NegBin expected zeros']
    values     = [observed_zeros, expected_zeros]
    colors     = ['steelblue', 'coral']
    axes[1].bar(categories, values, color=colors, alpha=0.8)
    axes[1].set_ylabel('Count')
    axes[1].set_title(f'Zero excess: {observed_zeros - expected_zeros:.0f} extra zeros')
    for i, v in enumerate(values):
        axes[1].text(i, v + 5, f'{v:.0f}', ha='center', fontweight='bold')

    plt.show()
    _ = style.save_figure(
        fig, data.Paths.from_config().root / "analysis/figures/zero_inflation_diagnostics",
        width="onehalf", save_png=True, save_pdf=True
    )

    plt.close(fig)

if __name__ == "__main__":
    cluster_data = data.load_cluster_features()
    shifted = np.asarray(cluster_data['n_sequences'] - 1)
    main(shifted)
