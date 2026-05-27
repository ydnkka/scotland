"""fig16 - validation bridge: from graph signature to confirmed event.

Left  : what the graph gives us (a candidate SSE-like signature).
Middle: the validation gap (a graph signature is not a confirmed event).
Right : the external evidence that would close that gap.
"""

from __future__ import annotations

import math
import os

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

CAP_BG = "#f5f7f7"

FONT = (
    "'Helvetica Neue', Helvetica, Arial, 'Liberation Sans', 'DejaVu Sans', sans-serif"
)
FONT_SCALE = 1.25

# ----------------------------------------------------------------- geometry
VBW, VBH = 1240, 370

# columns
L_X, L_W, L_Y, L_H = 24, 322, 34, 300  # left card
R_X, R_W = 790, 426  # right stack
R_Y, R_H = 34, 300
M_X0, M_X1 = L_X + L_W, R_X  # middle gap span
CXM = (M_X0 + M_X1) / 2  # middle centre x
YMID = 178  # bridge centre line

CARD_R = 12
HEADER_H = 38


# ------------------------------------------------------------- svg helpers
def _fmt(v) -> str:
    if isinstance(v, float):
        return f"{v:.2f}".rstrip("0").rstrip(".")
    return str(v)


def _esc(s: str) -> str:
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


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


def path(
    d,
    *,
    fill="none",
    stroke="none",
    sw=0,
    dash=None,
    cap="round",
    join="round",
    opacity=1.0,
):
    a = f'<path d="{d}" fill="{fill}" '
    if stroke != "none":
        a += (
            f'stroke="{stroke}" stroke-width="{_fmt(sw)}" '
            f'stroke-linecap="{cap}" stroke-linejoin="{join}" '
        )
    if dash:
        a += f'stroke-dasharray="{dash}" '
    if opacity != 1.0:
        a += f'opacity="{_fmt(opacity)}" '
    return a + "/>"


def text(
    x, y, s, *, size, fill, weight="normal", anchor="middle", italic=False, spacing=None
):
    extra = ""
    if spacing is not None:
        extra += f' letter-spacing="{_fmt(spacing)}"'
    if italic:
        extra += ' font-style="italic"'
    return (
        f'<text x="{_fmt(x)}" y="{_fmt(y)}" font-family="{FONT}" '
        f'font-size="{_fmt(size * FONT_SCALE)}" font-weight="{weight}" fill="{fill}" '
        f'text-anchor="{anchor}"{extra}>{_esc(s)}</text>'
    )


def text_lines(x, y, lines, dy, *, size, fill, weight="normal", anchor="middle"):
    spans = "".join(
        f'<tspan x="{_fmt(x)}"{"" if k == 0 else f" dy={chr(34)}{_fmt(dy * FONT_SCALE)}{chr(34)}"}>{_esc(ln)}</tspan>'
        for k, ln in enumerate(lines)
    )
    return (
        f'<text x="{_fmt(x)}" y="{_fmt(y)}" font-family="{FONT}" '
        f'font-size="{_fmt(size * FONT_SCALE)}" font-weight="{weight}" fill="{fill}" '
        f'text-anchor="{anchor}">{spans}</text>'
    )


def _trim(x1, y1, x2, y2, r1, r2):
    dx, dy = x2 - x1, y2 - y1
    d = math.hypot(dx, dy) or 1.0
    ux, uy = dx / d, dy / d
    return x1 + ux * r1, y1 + uy * r1, x2 - ux * r2, y2 - uy * r2


# ------------------------------------------------------------------- icons
def icon_chain(cx, cy, col, s):
    r = 0.33 * s
    off = 0.42 * s
    out = [
        circle(cx - off, cy, r, stroke=col, sw=2),
        circle(cx + off, cy, r, stroke=col, sw=2),
        line(cx - 0.10 * s, cy, cx + 0.10 * s, cy, stroke=col, sw=2),
    ]
    return "".join(out)


def icon_building(cx, cy, col, s):
    bw, bh = 0.9 * s, 0.62 * s
    bx, by = cx - bw / 2, cy - 0.10 * s
    roof = (
        f"M{_fmt(cx - 0.58 * s)},{_fmt(by)} L{_fmt(cx)},{_fmt(cy - 0.56 * s)} "
        f"L{_fmt(cx + 0.58 * s)},{_fmt(by)} Z"
    )
    out = [
        rect(bx, by, bw, bh, fill="none", stroke=col, sw=2),
        path(roof, fill=col),
        rect(cx - 0.13 * s, by + bh - 0.30 * s, 0.26 * s, 0.30 * s, fill=col),
    ]
    return "".join(out)


def icon_people(cx, cy, col, s):
    out = []
    for dx in (-0.42 * s, 0.42 * s):
        hx, hy, hr = cx + dx, cy - 0.30 * s, 0.20 * s
        out.append(circle(hx, hy, hr, stroke=col, sw=2))
        # shoulders: half ellipse below the head
        sw_ = 0.34 * s
        sh = 0.34 * s
        d = (
            f"M{_fmt(hx - sw_)},{_fmt(cy + 0.42 * s)} "
            f"A{_fmt(sw_)},{_fmt(sh)} 0 0 1 {_fmt(hx + sw_)},{_fmt(cy + 0.42 * s)}"
        )
        out.append(path(d, stroke=col, sw=2))
    return "".join(out)


def icon_route(cx, cy, col, s):
    d = (
        f"M{_fmt(cx - 0.55 * s)},{_fmt(cy + 0.36 * s)} "
        f"L{_fmt(cx - 0.12 * s)},{_fmt(cy - 0.30 * s)} "
        f"L{_fmt(cx + 0.18 * s)},{_fmt(cy + 0.10 * s)}"
    )
    out = [
        path(d, stroke=col, sw=2, dash="4 3"),
        line(
            cx + 0.18 * s,
            cy + 0.10 * s,
            cx + 0.52 * s,
            cy - 0.36 * s,
            stroke=col,
            sw=2.2,
            marker="rtip",
        ),
    ]
    return "".join(out)


# ------------------------------------------------------------- left card
def draw_left():
    out = []
    out.append(rect(L_X, L_Y, L_W, L_H, rx=CARD_R, fill=PANEL, stroke=TEAL, sw=2))
    out.append(rounded_top(L_X, L_Y, L_W, HEADER_H, CARD_R, TEAL_PALE))
    out.append(
        line(
            L_X + 0.5,
            L_Y + HEADER_H,
            L_X + L_W - 0.5,
            L_Y + HEADER_H,
            stroke=BORDER,
            sw=1,
            cap="butt",
        )
    )
    out.append(rect(L_X, L_Y, L_W, L_H, rx=CARD_R, fill="none", stroke=TEAL, sw=2))
    out.append(
        text(
            L_X + L_W / 2,
            L_Y + HEADER_H / 2 + 4,
            "WHAT GENOMES SHOW",
            size=12.5,
            fill=TEAL_DARK,
            weight="bold",
            spacing=0.6,
        )
    )

    # mini graph with highlighted candidate hub
    gx = L_X + L_W / 2
    gy = L_Y + HEADER_H + 70
    hub_r, glow_r, sat_r = 16, 30, 7
    sats = [
        (gx - 58, gy - 32),
        (gx - 58, gy + 32),
        (gx + 58, gy - 40),
        (gx + 58, gy),
        (gx + 58, gy + 40),
    ]
    lefts, rights = sats[:2], sats[2:]
    for sx, sy in lefts:
        x1, y1, x2, y2 = _trim(sx, sy, gx, gy, sat_r + 2, glow_r - 6)
        out.append(
            line(x1, y1, x2, y2, stroke=GRAY_LIGHT, sw=1.3, marker="gtip", cap="butt")
        )
    for sx, sy in rights:
        x1, y1, x2, y2 = _trim(gx, gy, sx, sy, glow_r - 6, sat_r + 5)
        out.append(
            line(x1, y1, x2, y2, stroke=GRAY_LIGHT, sw=1.3, marker="gtip", cap="butt")
        )
    for sx, sy in sats:
        out.append(circle(sx, sy, sat_r, fill=TEAL_SOFT, stroke="#ffffff", sw=1))
    out.append(circle(gx, gy, glow_r, fill=ORANGE, opacity=0.16))
    out.append(circle(gx, gy, hub_r, fill=ORANGE, stroke="#ffffff", sw=2))
    out.append(
        text(
            gx,
            gy + glow_r + 16,
            "candidate",
            size=12.5,
            fill=ORANGE_DARK,
            weight="bold",
        )
    )

    # bullets
    bullets = [
        "large relative to upstream burden",
        "high novelty / onward dissemination",
        "read within sampling + time context",
    ]
    by0 = L_Y + L_H - 92
    for k, txt in enumerate(bullets):
        yy = by0 + k * 30
        out.append(circle(L_X + 20, yy - 4, 3.2, fill=TEAL))
        out.append(text(L_X + 32, yy, txt, size=12.5, fill=INK_SOFT, anchor="start"))
    return "\n".join(out)


# ------------------------------------------------------------- middle gap
def draw_middle():
    out = []
    # dashed bridge arch over the gap
    ax0, ax1 = M_X0 + 14, M_X1 - 14
    arc = (
        f"M{_fmt(ax0)},{_fmt(YMID - 6)} "
        f"C{_fmt(CXM - 150)},{_fmt(YMID - 150)} "
        f"{_fmt(CXM + 150)},{_fmt(YMID - 150)} "
        f"{_fmt(ax1)},{_fmt(YMID - 6)}"
    )
    out.append(path(arc, stroke=GRAY_LIGHT, sw=2, dash="6 5"))

    # question title
    out.append(
        text_lines(
            CXM,
            YMID - 55,
            ["Is this a real", "superspreading event?"],
            18,
            size=15.5,
            fill=INK,
            weight="bold",
        )
    )

    # "?" gap circle
    qr = 26
    out.append(circle(CXM, YMID, qr, fill="#ffffff", stroke=ORANGE, sw=2))
    out.append(text(CXM, YMID + 11, "?", size=33, fill=ORANGE_DARK, weight="bold"))

    # subtitle
    out.append(
        text_lines(
            CXM,
            YMID + 50,
            [
                "both strands inform the question;",
                "absence of either is not",
                "evidence against an event.",
            ],
            15,
            size=11.5,
            fill=GRAY,
        )
    )

    # both arrows point inward — genomes and field evidence converge on the question
    out.append(
        line(M_X0 + 8, YMID, CXM - qr - 8, YMID, stroke=TEAL, sw=2.6, marker="teal_tip")
    )
    out.append(
        line(
            M_X1 - 8,
            YMID,
            CXM + qr + 8,
            YMID,
            stroke=ORANGE,
            sw=2.6,
            marker="orange_tip",
        )
    )
    return "\n".join(out)


# ------------------------------------------------------------- right card
def draw_right():
    out = []
    out.append(rect(R_X, R_Y, R_W, R_H, rx=CARD_R, fill=PANEL, stroke=ORANGE, sw=2))
    out.append(rounded_top(R_X, R_Y, R_W, HEADER_H, CARD_R, ORANGE_PALE))
    out.append(
        line(
            R_X + 0.5,
            R_Y + HEADER_H,
            R_X + R_W - 0.5,
            R_Y + HEADER_H,
            stroke=BORDER,
            sw=1,
            cap="butt",
        )
    )
    out.append(rect(R_X, R_Y, R_W, R_H, rx=CARD_R, fill="none", stroke=ORANGE, sw=2))
    out.append(
        text(
            R_X + R_W / 2,
            R_Y + HEADER_H / 2 + 4,
            "WHAT FIELD EVIDENCE SHOWS",
            size=12.5,
            fill=ORANGE_DARK,
            weight="bold",
            spacing=0.6,
        )
    )

    # central 2x2 icon grid (parallel to the left card's mini-graph)
    gx = R_X + R_W / 2
    gy = R_Y + HEADER_H + 70
    tile = 52
    off = 32
    icon_specs = [
        (-off, -off, icon_chain),
        (off, -off, icon_building),
        (-off, off, icon_people),
        (off, off, icon_route),
    ]
    for dx, dy, ic in icon_specs:
        tx = gx + dx - tile / 2
        ty = gy + dy - tile / 2
        out.append(rect(tx, ty, tile, tile, rx=10, fill=ORANGE, opacity=0.14))
        out.append(ic(gx + dx, gy + dy, ORANGE, 30))

    # bullets mirror the left card's bullet block
    bullets = [
        "outbreak and exposure links between cases",
        "setting context — venues, contacts, behaviour",
        "movement and travel between regions",
    ]
    by0 = R_Y + R_H - 92
    for k, txt in enumerate(bullets):
        yy = by0 + k * 30
        out.append(circle(R_X + 20, yy - 4, 3.2, fill=ORANGE))
        out.append(text(R_X + 32, yy, txt, size=12.5, fill=INK_SOFT, anchor="start"))
    return "\n".join(out)


# ------------------------------------------------------------------- build
def build_svg() -> str:
    def marker(mid, col):
        return (
            f'<marker id="{mid}" viewBox="0 0 10 10" refX="8.5" refY="5" '
            f'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path d="M0,0 L10,5 L0,10 z" fill="{col}"/></marker>'
        )

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {VBW} {VBH}" '
        f'width="{VBW}" height="{VBH}" font-family="{FONT}">',
        "<defs>",
        marker("teal_tip", TEAL),
        marker("orange_tip", ORANGE),
        marker("rtip", ORANGE),  # placeholder, recoloured per use below
        '<marker id="gtip" viewBox="0 0 10 10" refX="8" refY="5" '
        'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
        f'<path d="M0,1.5 L9,5 L0,8.5" fill="none" stroke="{GRAY_LIGHT}" '
        'stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"/>'
        "</marker>",
        "</defs>",
        rect(0, 0, VBW, VBH, fill=PANEL),
        draw_left(),
        draw_middle(),
        draw_right(),
        "</svg>",
    ]
    return "\n".join(parts)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    # repo layout: presentation/materials/lib/<this file>
    fig_dir = os.path.abspath(os.path.join(here, "..", "figures"))

    svg = build_svg()
    if not os.path.isdir(fig_dir):
        fig_dir = here
    svg_path = os.path.join(fig_dir, "validation_bridge.svg")
    with open(svg_path, "w") as fh:
        fh.write(svg)
    print("wrote", svg_path)

    cairosvg.svg2png(
        bytestring=svg.encode(),
        write_to=os.path.join(fig_dir, "validation_bridge.png"),
        output_width=2400,
        background_color="#ffffff",
    )
    cairosvg.svg2pdf(
        bytestring=svg.encode(),
        write_to=os.path.join(fig_dir, "validation_bridge.pdf"),
    )
    print("  + png + pdf")


if __name__ == "__main__":
    main()
