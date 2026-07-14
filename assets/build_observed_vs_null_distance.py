#!/usr/bin/env python3
"""Build the scenario-specific EpiLink scoring schematic as SVG and PNG.

Despite the legacy filename, the generated figure deliberately uses neither
"null distribution" nor "p-value": the simulated distributions are specific to
each configured epidemiological scenario and shaded areas are percentiles.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from xml.sax.saxutils import escape


DEFAULT_OUTPUT = Path(__file__).with_name("observed_vs_null_distance.svg")


def text(x: float, y: float, value: str, *, anchor: str = "start", css: str = "") -> str:
    return f'<text x="{x:g}" y="{y:g}" text-anchor="{anchor}" class="{css}">{escape(value)}</text>'


def density_path(x0: float, y0: float, width: float, height: float) -> str:
    points = []
    for i in range(81):
        t = i / 80
        density = math.exp(-0.5 * ((t - 0.5) / 0.18) ** 2)
        points.append((x0 + t * width, y0 - density * height))
    return "M " + " L ".join(f"{x:.1f},{y:.1f}" for x, y in points)


def distribution_panel(x: int, y: int, label: str, observed_t: float, q: float, colour: str) -> list[str]:
    w, h = 380, 96
    obs_x = x + observed_t * w
    curve = density_path(x, y + h, w, 68)
    # The shaded area is illustrative and explicitly labelled as an empirical percentile.
    shade_points = []
    steps = max(2, round(observed_t * 80))
    for i in range(steps + 1):
        t = observed_t * i / steps
        density = math.exp(-0.5 * ((t - 0.5) / 0.18) ** 2)
        shade_points.append((x + t * w, y + h - density * 68))
    shade = f"M {x},{y+h} L " + " L ".join(f"{px:.1f},{py:.1f}" for px, py in shade_points) + f" L {obs_x:.1f},{y+h} Z"
    return [
        text(x, y - 18, label, css="panel-title"),
        f'<path d="{shade}" fill="{colour}" opacity="0.28"/>',
        f'<path d="{curve}" fill="none" stroke="#5b6b7c" stroke-width="2"/>',
        f'<line x1="{x}" y1="{y+h}" x2="{x+w}" y2="{y+h}" class="axis"/>',
        f'<line x1="{obs_x:.1f}" y1="{y+11}" x2="{obs_x:.1f}" y2="{y+h}" stroke="{colour}" stroke-width="2.5"/>',
        text(obs_x, y + 7, "observed", anchor="middle", css="observed"),
        text(x, y + h + 20, "smaller distance", css="tiny"),
        text(x + w, y + h + 20, "larger distance", anchor="end", css="tiny"),
        text(x + w + 30, y + 42, f"q = {q:.2f}", css="q"),
        text(x + w + 30, y + 66, "empirical percentile", css="tiny"),
    ]


def build_svg() -> str:
    width, height = 1180, 780
    parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">
  <title id="title">EpiLink scenario-specific distance compatibility scoring</title>
  <desc id="desc">Temporal and genetic observed distances are located in separate scenario-specific simulated distributions. Each percentile is transformed to a two-sided compatibility, multiplied within scenario, then summed over direct and co-primary target scenarios.</desc>
  <defs>
    <style>
      text {{ font-family: Arial, Helvetica, sans-serif; fill: #17202a; }}
      .title {{ font-size: 22px; font-weight: 700; }} .subtitle {{ font-size: 14px; fill: #526273; }}
      .step {{ font-size: 14px; font-weight: 700; fill: #2563eb; letter-spacing: .6px; }}
      .panel-title {{ font-size: 17px; font-weight: 700; }} .body {{ font-size: 15px; }}
      .tiny {{ font-size: 12px; fill: #64748b; }} .observed {{ font-size: 12px; font-weight: 700; fill: #b42318; }}
      .q {{ font-size: 17px; font-weight: 700; }} .formula {{ font-family: 'Courier New', monospace; font-size: 18px; }}
      .formula-big {{ font-family: 'Courier New', monospace; font-size: 21px; font-weight: 700; }}
      .formula-small {{ font-family: 'Courier New', monospace; font-size: 16px; }}
      .axis {{ stroke: #334155; stroke-width: 1.2; }} .card {{ fill: #ffffff; stroke: #cbd5e1; stroke-width: 1.2; }}
    </style>
    <marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto">
      <path d="M 1 1 L 8 5 L 1 9" fill="none" stroke="#64748b" stroke-width="1.5"/>
    </marker>
  </defs>
  <rect width="100%" height="100%" fill="#f8fafc"/>
''']
    parts += [
        text(42, 42, "From observed distances to EpiLink compatibility", css="title"),
        text(42, 67, "10,000 Monte Carlo draws per configured scenario under the natural-history and mutation model", css="subtitle"),
        text(42, 108, "1  SCORE DISTANCES SEPARATELY WITHIN EACH SCENARIO", css="step"),
        '<rect x="42" y="126" width="700" height="392" rx="10" class="card"/>',
        text(66, 160, "Example: one configured target scenario", css="panel-title"),
    ]
    parts += distribution_panel(86, 205, "Temporal distance", 0.35, 0.35, "#d14d41")
    parts += distribution_panel(86, 374, "Genetic distance", 0.72, 0.72, "#7c3aed")
    parts += [
        '<rect x="770" y="126" width="368" height="392" rx="10" class="card"/>',
        text(794, 160, "Percentile → compatibility", css="panel-title"),
        text(954, 213, "c(q) = 1 − 2|q − 0.5|", anchor="middle", css="formula-big"),
        '<path d="M 824 354 L 954 239 L 1084 354" fill="#dbeafe" stroke="#2563eb" stroke-width="2"/>',
        '<line x1="824" y1="354" x2="1084" y2="354" class="axis"/>',
        text(824, 376, "q = 0", anchor="middle", css="tiny"),
        text(954, 376, "0.5", anchor="middle", css="tiny"),
        text(1084, 376, "1", anchor="middle", css="tiny"),
        text(954, 231, "maximum compatibility", anchor="middle", css="tiny"),
        text(954, 417, "Highest near the simulated median;", anchor="middle", css="body"),
        text(954, 440, "decreasing toward either tail.", anchor="middle", css="body"),
        text(954, 482, "q is an empirical percentile.", anchor="middle", css="observed"),
        text(42, 558, "2  COMBINE MARGINAL COMPATIBILITIES", css="step"),
        '<rect x="42" y="578" width="1096" height="160" rx="10" class="card"/>',
        '<line x1="590" y1="598" x2="590" y2="718" stroke="#e2e8f0" stroke-width="1.5"/>',
        text(66, 611, "Within each scenario", css="panel-title"),
        text(66, 648, "scenario score = c(q_temporal) × c(q_genetic)", css="formula"),
        text(66, 693, "Multiply the temporal and genetic marginal compatibilities.", css="body"),
        text(620, 611, "Across configured target scenarios", css="panel-title"),
        text(620, 648, "target score = score_ad(0) + score_ca(0,0)", css="formula-small"),
        text(620, 682, "ad(0): direct ancestor–descendant", css="body"),
        text(620, 710, "ca(0,0): co-primary / common-ancestor", css="body"),
        "</svg>\n",
    ]
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
