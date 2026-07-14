#!/usr/bin/env python3
"""Build the sliding three-week window methods schematic as SVG and PNG.

The SVG generation uses only the standard library; PNG export uses CairoSVG.
Run the script from any directory. By default, both files are written beside it.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from xml.sax.saxutils import escape


DEFAULT_OUTPUT = Path(__file__).with_name("sliding_3week_windows.svg")


def text(x: float, y: float, value: str, *, anchor: str = "start", css: str = "") -> str:
    return (
        f'<text x="{x:g}" y="{y:g}" text-anchor="{anchor}" class="{css}">'
        f"{escape(value)}</text>"
    )


def build_svg() -> str:
    width, height = 1000, 500
    plot_x, cell_w = 150, 72
    weeks = 10
    row_y = [118, 174, 230, 286, 342]

    parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">
  <title id="title">Sliding three-week temporal windows with a one-week step</title>
  <desc id="desc">Five overlapping, half-open three-week windows. Samples are partitioned by Pango lineage and clustered separately within every window.</desc>
  <defs>
    <style>
      text {{ font-family: Arial, Helvetica, sans-serif; fill: #17202a; }}
      .header {{ font-size: 15px; font-weight: 700; }}
      .label {{ font-size: 15px; }}
      .small {{ font-size: 13px; fill: #475569; }}
      .note {{ font-size: 14px; }}
      .window {{ fill: #dbeafe; stroke: #2563eb; stroke-width: 1.5; }}
      .grid {{ stroke: #cbd5e1; stroke-width: 1; }}
    </style>
    <marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto">
      <path d="M 1 1 L 8 5 L 1 9" fill="none" stroke="#334155" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
    </marker>
  </defs>
  <rect width="100%" height="100%" fill="#ffffff"/>
''']

    parts.append(text(42, 38, "Rolling-window construction", css="header"))
    parts.append(text(42, 62, "Three-week duration · one-week step", css="small"))

    for i in range(weeks):
        x = plot_x + (i + 0.5) * cell_w
        parts.append(text(x, 92, f"Week {i + 1}", anchor="middle", css="header"))
        boundary = plot_x + i * cell_w
        parts.append(f'<line x1="{boundary}" y1="101" x2="{boundary}" y2="374" class="grid"/>')
    parts.append(f'<line x1="{plot_x + weeks * cell_w}" y1="101" x2="{plot_x + weeks * cell_w}" y2="374" class="grid"/>')

    for index, y in enumerate(row_y, start=1):
        start_week = index
        x = plot_x + (index - 1) * cell_w
        parts.append(text(plot_x - 16, y + 24, f"Window {index}", anchor="end", css="label"))
        parts.append(f'<rect x="{x}" y="{y}" width="{3 * cell_w}" height="38" rx="5" class="window"/>')
        parts.append(text(x + 1.5 * cell_w, y + 25, f"Weeks {start_week}–{start_week + 2}", anchor="middle", css="label"))

    parts.extend([
        '<line x1="150" y1="389" x2="222" y2="389" stroke="#334155" stroke-width="1.5" marker-end="url(#arrow)"/>',
        text(150, 411, "one-week step", css="small"),
        '<rect x="42" y="432" width="916" height="47" rx="7" fill="#f8fafc" stroke="#cbd5e1"/>',
        text(58, 452, "Membership rule:", css="header"),
        text(184, 452, "start ≤ collection_date < end (half-open interval)", css="note"),
        text(58, 470, "A sequence may occur in up to three overlapping windows; clustering is performed separately by Pango lineage within each window.", css="small"),
        "</svg>\n",
    ])
    return "\n".join(parts)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", type=Path, default=DEFAULT_OUTPUT, help="output SVG path")
    parser.add_argument("--png-output", type=Path, help="output PNG path (default: SVG path with .png suffix)")
    parser.add_argument("--png-scale", type=float, default=2.0, help="PNG resolution multiplier (default: 2)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.png_scale <= 0:
        raise SystemExit("--png-scale must be greater than zero")
    png_output = args.png_output or args.output.with_suffix(".png")
    svg = build_svg()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    png_output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(svg, encoding="utf-8")
    try:
        import cairosvg
    except ImportError as exc:
        raise SystemExit("PNG export requires CairoSVG; run with: conda run -n PhD python <script>") from exc
    cairosvg.svg2png(bytestring=svg.encode("utf-8"), write_to=str(png_output), scale=args.png_scale)
    print(f"Wrote {args.output}")
    print(f"Wrote {png_output}")


if __name__ == "__main__":
    main()
