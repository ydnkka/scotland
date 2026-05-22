import pandas as pd

from utils.data import apply_window_stride, select_window_stride


def test_select_window_stride_uses_sorted_window_positions():
    windows = sorted(list(range(1, 25)))

    assert select_window_stride(windows, stride=1) == windows[::1]
    assert select_window_stride(windows, stride=2) == windows[::2]
    assert select_window_stride(windows, stride=3) == windows[::3]


def test_apply_window_stride_filters_and_renumbers_windows():
    df = pd.DataFrame(
        {
            "window_idx": [1, 1, 2, 3, 3, 4, 5],
            "window_id": ["W001", "W001", "W002", "W003", "W003", "W004", "W005"],
            "value": range(7),
        }
    )

    out = apply_window_stride(df, stride=2)

    assert out["value"].tolist() == [0, 1, 3, 4, 6]
    assert out["window_idx"].tolist() == [1, 1, 2, 2, 3]
    assert out["window_id"].tolist() == ["W001", "W001", "W002", "W002", "W003"]
