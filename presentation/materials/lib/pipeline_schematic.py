"""audience-friendly pipeline schematic for the RIS talk.

Five stages, left to right:
  1. Genomes + metadata
  2. Windowed clusters
  3. Transition graph
  4. Candidate nodes
  5. Association models
"""

from __future__ import annotations

import math
import os
import random

import cairosvg

# ----------------------------------------------------------------- palette
INK = "#1f2937"
INK_SOFT = "#374151"
GRAY = "#64748b"
GRAY_LIGHT = "#94a3b8"
BORDER = "#d9dee4"
PANEL = "#ffffff"

TEAL = "#0f766e"
TEAL_DARK = "#115e54"
TEAL_SOFT = "#5eada3"
TEAL_PALE = "#e6f4f1"

ORANGE = "#dd7a33"
ORANGE_DARK = "#b45309"
ORANGE_PALE = "#fdf0e3"
ORANGE_INK = "#8a5a22"

NODE_GRAY = "#aab9c2"

FONT = (
    "'Helvetica Neue', Helvetica, Arial, 'Liberation Sans', 'DejaVu Sans', sans-serif"
)
FONT_SCALE = 1.25

# ----------------------------------------------------------------- geometry
MARGIN_X = 36
TOP = 30
CARD_W = 240
CARD_H = 372
GAP = 56
N = 5

VBW = 2 * MARGIN_X + N * CARD_W + (N - 1) * GAP  # 1496
VBH = 492

CARD_Y0 = TOP
CARD_Y1 = TOP + CARD_H  # 402
CARD_R = 16
TITLE_H = 64

GLYPH_PAD_X = 16
GLYPH_TOP = CARD_Y0 + TITLE_H + 14  # 108
GLYPH_BOT = CARD_Y1 - 16  # 386
GLYPH_MIDY = (GLYPH_TOP + GLYPH_BOT) / 2  # 247

CAP_Y0 = CARD_Y1 + 34
CAP_DY = 22

STAGES = [
    (
        ("Genomes", "+ metadata"),
        ("Sequences with age, place,", "deprivation and date"),
    ),
    (
        ("Windowed", "clusters"),
        ("group genomes inside sliding", "epidemic-time windows"),
    ),
    (
        ("Cluster", "graph"),
        ("link clusters that share sequence", "continuity across windows"),
    ),
    (
        ("Candidate", "nodes"),
        ("flag nodes with unusual", "amplification or onward spread"),
    ),
    (("Association", "models"), ("ask where candidate signatures", "are enriched")),
]


def card_left(i: int) -> float:
    return MARGIN_X + i * (CARD_W + GAP)


# ------------------------------------------------------------- svg helpers
def _fmt(v) -> str:
    if isinstance(v, float):
        return f"{v:.2f}".rstrip("0").rstrip(".")
    return str(v)


def rect(x, y, w, h, *, rx=0, fill="none", stroke="none", sw=0, opacity=1.0, dash=None):
    a = (
        f'<rect x="{_fmt(x)}" y="{_fmt(y)}" width="{_fmt(w)}" '
        f'height="{_fmt(h)}" rx="{_fmt(rx)}" fill="{fill}" '
    )
    if stroke != "none":
        a += f'stroke="{stroke}" stroke-width="{_fmt(sw)}" '
    if dash:
        a += f'stroke-dasharray="{dash}" '
    if opacity != 1.0:
        a += f'opacity="{_fmt(opacity)}" '
    return a + "/>"


def rounded_top(x, y, w, h, r, fill, opacity=1.0):
    d = (
        f"M{_fmt(x)},{_fmt(y + r)} Q{_fmt(x)},{_fmt(y)} {_fmt(x + r)},{_fmt(y)} "
        f"H{_fmt(x + w - r)} Q{_fmt(x + w)},{_fmt(y)} {_fmt(x + w)},{_fmt(y + r)} "
        f"V{_fmt(y + h)} H{_fmt(x)} Z"
    )
    op = f' opacity="{_fmt(opacity)}"' if opacity != 1.0 else ""
    return f'<path d="{d}" fill="{fill}"{op}/>'


def circle(cx, cy, r, *, fill="none", stroke="none", sw=0, opacity=1.0):
    a = f'<circle cx="{_fmt(cx)}" cy="{_fmt(cy)}" r="{_fmt(r)}" fill="{fill}" '
    if stroke != "none":
        a += f'stroke="{stroke}" stroke-width="{_fmt(sw)}" '
    if opacity != 1.0:
        a += f'opacity="{_fmt(opacity)}" '
    return a + "/>"


def ellipse(cx, cy, rx, ry, *, fill="none", stroke="none", sw=0, opacity=1.0):
    a = (
        f'<ellipse cx="{_fmt(cx)}" cy="{_fmt(cy)}" rx="{_fmt(rx)}" '
        f'ry="{_fmt(ry)}" fill="{fill}" '
    )
    if stroke != "none":
        a += f'stroke="{stroke}" stroke-width="{_fmt(sw)}" '
    if opacity != 1.0:
        a += f'opacity="{_fmt(opacity)}" '
    return a + "/>"


def line(
    x1, y1, x2, y2, *, stroke, sw, dash=None, cap="round", marker=None, opacity=1.0
):
    a = (
        f'<line x1="{_fmt(x1)}" y1="{_fmt(y1)}" x2="{_fmt(x2)}" y2="{_fmt(y2)}" '
        f'stroke="{stroke}" stroke-width="{_fmt(sw)}" stroke-linecap="{cap}" '
    )
    if dash:
        a += f'stroke-dasharray="{dash}" '
    if marker:
        a += f'marker-end="url(#{marker})" '
    if opacity != 1.0:
        a += f'opacity="{_fmt(opacity)}" '
    return a + "/>"


def text(
    x, y, s, *, size, fill, weight="normal", anchor="middle", italic=False, spacing=None
):
    style = ""
    if spacing is not None:
        style = f' letter-spacing="{_fmt(spacing)}"'
    it = ' font-style="italic"' if italic else ""
    return (
        f'<text x="{_fmt(x)}" y="{_fmt(y)}" font-family="{FONT}" '
        f'font-size="{_fmt(size * FONT_SCALE)}" font-weight="{weight}" fill="{fill}" '
        f'text-anchor="{anchor}"{it}{style}>{s}</text>'
    )


def text2(x, y, line1, line2, dy, *, size, fill, weight="normal", anchor="middle"):
    return (
        f'<text x="{_fmt(x)}" y="{_fmt(y)}" font-family="{FONT}" '
        f'font-size="{_fmt(size * FONT_SCALE)}" font-weight="{weight}" fill="{fill}" '
        f'text-anchor="{anchor}">'
        f'<tspan x="{_fmt(x)}">{line1}</tspan>'
        f'<tspan x="{_fmt(x)}" dy="{_fmt(dy * FONT_SCALE)}">{line2}</tspan></text>'
    )


# ------------------------------------------------------------------- card
def draw_card(i, title_lines):
    out = []
    x0 = card_left(i)
    # body
    out.append(rect(x0, CARD_Y0, CARD_W, CARD_H, rx=CARD_R, fill=PANEL))
    # title strip (rounded top only)
    out.append(rounded_top(x0, CARD_Y0, CARD_W, TITLE_H, CARD_R, TEAL_PALE))
    # hairline under strip
    out.append(
        line(
            x0 + 0.5,
            CARD_Y0 + TITLE_H,
            x0 + CARD_W - 0.5,
            CARD_Y0 + TITLE_H,
            stroke=BORDER,
            sw=1,
            cap="butt",
        )
    )
    # crisp outline on top of everything
    out.append(
        rect(x0, CARD_Y0, CARD_W, CARD_H, rx=CARD_R, fill="none", stroke=BORDER, sw=1)
    )
    # number chip
    cy = CARD_Y0 + TITLE_H / 2
    cx = x0 + 30
    out.append(circle(cx, cy, 16, fill=TEAL))
    out.append(text(cx, cy + 6, str(i + 1), size=18, fill="#ffffff", weight="bold"))
    # title (two lines), vertically centred on the strip
    tx = x0 + 56
    out.append(
        text2(
            tx,
            cy - 2,
            title_lines[0],
            title_lines[1],
            20,
            size=19,
            fill=TEAL_DARK,
            weight="bold",
            anchor="start",
        )
    )
    return "\n".join(out)


def draw_caption(i, cap_lines):
    cx = card_left(i) + CARD_W / 2
    return text2(cx, CAP_Y0, cap_lines[0], cap_lines[1], CAP_DY, size=14, fill=GRAY)


def draw_arrow(i):
    x_end = card_left(i) + CARD_W
    x_next = card_left(i + 1)
    y = (CARD_Y0 + CARD_Y1) / 2
    return line(x_end + 12, y, x_next - 12, y, stroke=ORANGE, sw=3.6, marker="arrow")


# ----------------------------------------------------------------- glyphs
def glyph_inputs(i):
    out = []
    x0 = card_left(i)
    cx = x0 + CARD_W / 2

    # "genomes" label + three barcode bars
    out.append(text(cx, GLYPH_TOP + 14, "genomes", size=12.5, fill=GRAY, italic=True))
    bar_w, bar_h = 168, 19
    bx = cx - bar_w / 2
    rng = random.Random(3)
    for k, by in enumerate((GLYPH_TOP + 24, GLYPH_TOP + 52, GLYPH_TOP + 80)):
        out.append(
            rect(bx, by, bar_w, bar_h, rx=3, fill="#f1f5f4", stroke=BORDER, sw=1)
        )
        tick_col = TEAL_DARK if k == 1 else TEAL
        xx = bx + 5
        while xx < bx + bar_w - 6:
            w = rng.uniform(3.0, 9.5)
            if rng.random() < 0.62:
                out.append(rect(xx, by + 3, w, bar_h - 6, fill=tick_col))
            xx += w + rng.uniform(2.0, 6.0)

    # "metadata" label + chips
    out.append(text(cx, GLYPH_BOT - 56, "metadata", size=12.5, fill=GRAY, italic=True))
    chips = ["age", "place", "SIMD", "date"]
    cols = [ORANGE, TEAL, ORANGE, TEAL]
    cw, ch, cg = 46, 28, 7
    total = len(chips) * cw + (len(chips) - 1) * cg
    sx = cx - total / 2
    chip_y = GLYPH_BOT - 44
    for c, col in zip(chips, cols):
        out.append(rect(sx, chip_y, cw, ch, rx=7, fill=col, opacity=0.16))
        txt_col = ORANGE_INK if col == ORANGE else TEAL_DARK
        out.append(
            text(
                sx + cw / 2,
                chip_y + ch / 2 + 4,
                c,
                size=12,
                fill=txt_col,
                weight="bold",
            )
        )
        sx += cw + cg
    return "\n".join(out)


def glyph_window(i):
    out = []
    x0 = card_left(i)
    cx = x0 + CARD_W / 2
    gx0 = x0 + GLYPH_PAD_X
    gw = CARD_W - 2 * GLYPH_PAD_X

    band_top = GLYPH_TOP + 4
    band_h = (GLYPH_BOT - 22) - band_top
    out.append(rect(gx0 + 8, band_top, gw - 16, band_h, rx=6, fill=TEAL, opacity=0.07))
    out.append(
        rect(
            gx0 + 8,
            band_top,
            gw - 16,
            band_h,
            rx=6,
            fill="none",
            stroke=TEAL_SOFT,
            sw=1,
            dash="5 4",
        )
    )

    rng = random.Random(11)
    bmid = band_top + band_h / 2
    clusters = [(cx - 36, bmid - 26, TEAL), (cx + 34, bmid + 30, ORANGE)]
    for mx, my, col in clusters:
        out.append(ellipse(mx, my, 44, 40, fill=col, opacity=0.11, stroke=col, sw=1))
        for _ in range(7):
            px = mx + rng.gauss(0, 16)
            py = my + rng.gauss(0, 15)
            # keep dots inside the ellipse
            while ((px - mx) / 40) ** 2 + ((py - my) / 36) ** 2 > 1:
                px = mx + rng.gauss(0, 16)
                py = my + rng.gauss(0, 15)
            out.append(circle(px, py, 4.2, fill=col, stroke="#ffffff", sw=1))

    out.append(
        text(cx, GLYPH_BOT - 2, "one 3-week window", size=12.5, fill=GRAY, italic=True)
    )
    return "\n".join(out)


def _graph_layout(i):
    x0 = card_left(i)
    gx0 = x0 + GLYPH_PAD_X + 18
    gx1 = x0 + CARD_W - GLYPH_PAD_X - 18
    cols = [gx0, (gx0 + gx1) / 2, gx1]
    m = GLYPH_MIDY
    layout = {
        (0, 0): (cols[0], m - 30),
        (0, 1): (cols[0], m + 30),
        (1, 0): (cols[1], m - 42),
        (1, 1): (cols[1], m + 8),
        (1, 2): (cols[1], m + 50),
        (2, 0): (cols[2], m - 22),
        (2, 1): (cols[2], m + 26),
    }
    edges = [
        ((0, 0), (1, 0)),
        ((0, 0), (1, 1)),
        ((0, 1), (1, 2)),
        ((1, 0), (2, 0)),
        ((1, 1), (2, 1)),
        ((1, 2), (2, 1)),
    ]
    cand = (1, 0)
    return cols, layout, edges, cand


def _trim(x1, y1, x2, y2, r1, r2):
    dx, dy = x2 - x1, y2 - y1
    d = math.hypot(dx, dy) or 1.0
    ux, uy = dx / d, dy / d
    return x1 + ux * r1, y1 + uy * r1, x2 - ux * r2, y2 - uy * r2


def glyph_graph(i, highlight=False):
    out = []
    cols, layout, edges, cand = _graph_layout(i)
    node_col = NODE_GRAY if highlight else TEAL_SOFT
    r = 9

    # edges (drawn under the nodes)
    for a, b in edges:
        ra = 18 if (highlight and a == cand) else r
        rb = 18 if (highlight and b == cand) else r
        xa, ya = layout[a]
        xb, yb = layout[b]
        x1, y1, x2, y2 = _trim(xa, ya, xb, yb, ra + 3, rb + 6)
        out.append(
            line(x1, y1, x2, y2, stroke=GRAY_LIGHT, sw=1, marker="gtip", cap="butt")
        )

    # nodes
    for key, (x, y) in layout.items():
        if highlight and key == cand:
            out.append(circle(x, y, 30, fill=ORANGE, opacity=0.16))
            out.append(circle(x, y, 16, fill=ORANGE, stroke="#ffffff", sw=1))
        else:
            out.append(circle(x, y, r, fill=node_col, stroke="#ffffff", sw=1))

    # candidate label, placed clearly above the glow (no overlap)
    if highlight:
        hx, hy = layout[cand]
        out.append(
            text(
                hx, hy - 30 - 8, "candidate", size=12.5, fill=ORANGE_DARK, weight="bold"
            )
        )

    # time axis ticks
    labels = ["t", "t+1", "t+2"]
    for xx, lab in zip(cols, labels):
        out.append(text(xx, GLYPH_BOT - 2, lab, size=12, fill=GRAY_LIGHT))
    return "\n".join(out)


def glyph_models(i):
    out = []
    x0 = card_left(i)
    gx0 = x0 + GLYPH_PAD_X
    axis_x = x0 + 150
    scale = 24.0

    ests = [
        (1.6, ORANGE, 1.0),
        (-0.6, TEAL_DARK, 0.9),
        (0.5, GRAY, 0.7),
        (1.1, ORANGE, 0.9),
    ]
    labels = ["health board", "Borders", "age", "rural"]
    rows_y = [GLYPH_TOP + 34, GLYPH_TOP + 92, GLYPH_TOP + 150, GLYPH_TOP + 208]

    # OR = 1 reference line
    out.append(
        line(
            axis_x,
            GLYPH_TOP + 12,
            axis_x,
            GLYPH_BOT - 22,
            stroke=GRAY_LIGHT,
            sw=1.2,
            dash="5 4",
            cap="butt",
        )
    )

    for y, (e, col, w), lab in zip(rows_y, ests, labels):
        out.append(
            line(
                axis_x + (e - w) * scale,
                y,
                axis_x + (e + w) * scale,
                y,
                stroke=col,
                sw=5.0,
                opacity=0.45,
            )
        )
        out.append(circle(axis_x + e * scale, y, 6, fill=col, stroke="#ffffff", sw=1))
        out.append(text(gx0, y + 4, lab, size=13, fill=INK_SOFT, anchor="start"))

    out.append(
        text(axis_x, GLYPH_BOT - 2, "odds ratio", size=12.5, fill=GRAY, italic=True)
    )
    return "\n".join(out)


# ------------------------------------------------------------------- build
def build_svg() -> str:
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {VBW} {VBH}" '
        f'width="{VBW}" height="{VBH}" font-family="{FONT}">',
        "<defs>",
        '<marker id="arrow" viewBox="0 0 10 10" refX="8.5" refY="5" '
        'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
        f'<path d="M0,0 L10,5 L0,10 z" fill="{ORANGE}"/></marker>',
        '<marker id="gtip" viewBox="0 0 10 10" refX="8" refY="5" '
        'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
        f'<path d="M0,1.5 L9,5 L0,8.5" fill="none" stroke="{GRAY_LIGHT}" '
        'stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/>'
        "</marker>",
        "</defs>",
        rect(0, 0, VBW, VBH, fill=PANEL),
    ]

    glyphs = [
        glyph_inputs,
        glyph_window,
        lambda i: glyph_graph(i, highlight=False),
        lambda i: glyph_graph(i, highlight=True),
        glyph_models,
    ]

    for i, (title_lines, cap_lines) in enumerate(STAGES):
        parts.append(draw_card(i, title_lines))
        parts.append(glyphs[i](i))
        parts.append(draw_caption(i, cap_lines))
        if i < N - 1:
            parts.append(draw_arrow(i))

    parts.append("</svg>")
    return "\n".join(parts)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    # repo layout: presentation/materials/lib/<this file>
    fig_dir = os.path.abspath(os.path.join(here, "..", "figures"))

    svg = build_svg()

    svg_path = os.path.join(fig_dir, "pipeline_schematic.svg")
    with open(svg_path, "w") as fh:
        fh.write(svg)
    print("wrote", svg_path)
    cairosvg.svg2png(
        bytestring=svg.encode(),
        write_to=os.path.join(fig_dir, "pipeline_schematic.png"),
        output_width=2400,
        background_color="#ffffff",
    )
    cairosvg.svg2pdf(
        bytestring=svg.encode(),
        write_to=os.path.join(fig_dir, "pipeline_schematic.pdf"),
    )
    print("  + png + pdf")


if __name__ == "__main__":
    main()
