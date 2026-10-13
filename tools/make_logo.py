#!/usr/bin/env python3
"""Generate the Ephemera logo files in assets/logo/ (development only).

The mark is a lowercase e drawn as a compass construction: the bowl and bar
as one stroke, dashed construction circles and cross-hairs, a ring of 24
inward ticks, and a gold planet with a small orbit at the end of the bar.
Small sizes use a bold, simplified e with the planet.  The wordmark is
Jost Light (SIL Open Font License), letter-spaced and converted to outlines,
so the files need no font.

    python3 tools/make_logo.py --font Jost[wght].ttf

Needs fonttools and cairosvg (for the PNG files).
"""
import argparse
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "assets", "logo")

INK_LIGHT = "#1B1F33"   # ink on light backgrounds
INK_DARK = "#E8E4DA"    # ink on dark backgrounds
GOLD = "#D4A23F"
PAPER = "#F5F2EA"
NIGHT = "#10131F"


def ticks(r1, r2, step):
    d = ""
    for a in range(0, 360, step):
        t = math.radians(a)
        c, s = math.cos(t), math.sin(t)
        d += "M%.2f %.2fL%.2f %.2f" % (100 + r1 * c, 100 + r1 * s, 100 + r2 * c, 100 + r2 * s)
    return d


def mark(ink, accent):
    """The full mark, 200 x 200 units."""
    return f'''<g fill="none" stroke="{ink}" stroke-linecap="round">
<circle cx="100" cy="100" r="90" stroke-width="1"/>
<path d="{ticks(90, 84, 15)}" stroke-width="0.7"/>
<circle cx="100" cy="100" r="60" stroke-width="0.6" stroke-dasharray="2 3"/>
<circle cx="100" cy="100" r="30" stroke-width="0.6" stroke-dasharray="2 3"/>
<line x1="100" y1="22" x2="100" y2="178" stroke-width="0.6" stroke-dasharray="2 3"/>
<line x1="22" y1="100" x2="178" y2="100" stroke-width="0.6" stroke-dasharray="2 3"/>
<path d="M160 100 A60 60 0 1 0 138.57 145.96" stroke-width="2"/>
<line x1="40" y1="100" x2="160" y2="100" stroke-width="2"/>
<circle cx="160" cy="100" r="11" stroke-width="0.8"/>
</g>
<circle cx="138.57" cy="145.96" r="2.2" fill="{ink}"/>
<circle cx="40" cy="100" r="2.2" fill="{ink}"/>
<circle cx="160" cy="100" r="5" fill="{accent}"/>'''


def small_mark(ink, accent, stroke=18):
    """The bold e for favicons and small sizes, 200 x 200 units."""
    return f'''<g fill="none" stroke="{ink}" stroke-linecap="round" stroke-width="{stroke}">
<path d="M164 100 A64 64 0 1 0 141.14 149.03"/>
<line x1="36" y1="100" x2="146" y2="100"/>
</g>
<circle cx="164" cy="100" r="20" fill="{accent}"/>'''


class Wordmark:
    """Text as outlines from Jost Light, with letter spacing."""

    def __init__(self, font_path, weight=300):
        from fontTools.ttLib import TTFont
        from fontTools.varLib import instancer
        font = TTFont(font_path)
        if "fvar" in font:
            font = instancer.instantiateVariableFont(font, {"wght": weight})
        self.font = font
        self.upem = font["head"].unitsPerEm
        self.cap = font["OS/2"].sCapHeight
        self.glyphs = font.getGlyphSet()
        self.cmap = font.getBestCmap()

    def path(self, text, cap_height, tracking=0.0):
        """SVG path data with its baseline at y = 0 and the given cap
        height; tracking in em.  Returns (d, width)."""
        from fontTools.pens.svgPathPen import SVGPathPen
        from fontTools.pens.transformPen import TransformPen
        scale = cap_height / self.cap
        pen = SVGPathPen(self.glyphs)
        x = 0.0
        for i, ch in enumerate(text):
            name = self.cmap[ord(ch)]
            glyph = self.glyphs[name]
            glyph.draw(TransformPen(pen, (scale, 0, 0, -scale, x, 0)))
            x += glyph.width * scale
            if i < len(text) - 1:
                x += tracking * self.upem * scale
        return pen.getCommands(), x


def svg(width, height, body, background=None):
    bg = f'<rect width="{width}" height="{height}" fill="{background}"/>' if background else ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}">{bg}{body}</svg>\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--font", required=True, help="Jost variable font (Jost[wght].ttf)")
    args = ap.parse_args()
    import cairosvg
    os.makedirs(OUT, exist_ok=True)
    words = Wordmark(args.font)

    def write(name, text):
        with open(os.path.join(OUT, name), "w") as f:
            f.write(text)

    def png(name, text, w, h):
        cairosvg.svg2png(bytestring=text.encode(), write_to=os.path.join(OUT, name),
                         output_width=w, output_height=h)

    # The mark.
    for name, ink, accent in (("mark-light.svg", INK_LIGHT, GOLD), ("mark-dark.svg", INK_DARK, GOLD),
                              ("mark-black.svg", "#000000", "#000000"),
                              ("mark-white.svg", "#FFFFFF", "#FFFFFF")):
        write(name, svg(200, 200, mark(ink, accent)))

    # Small mark: follows the browser's light or dark scheme.
    write("mark-small.svg",
          '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="200" viewBox="0 0 200 200">'
          f'<style>.ink{{stroke:{INK_LIGHT}}}@media (prefers-color-scheme:dark){{.ink{{stroke:{INK_DARK}}}}}</style>'
          '<g class="ink" fill="none" stroke-linecap="round" stroke-width="18">'
          '<path d="M164 100 A64 64 0 1 0 141.14 149.03"/><line x1="36" y1="100" x2="146" y2="100"/></g>'
          f'<circle cx="164" cy="100" r="20" fill="{GOLD}"/></svg>\n')

    # Lockups: horizontal (mark beside the word) and stacked.
    for theme, ink in (("light", INK_LIGHT), ("dark", INK_DARK)):
        d, w = words.path("EPHEMERA", 34, tracking=0.32)
        pad = 10
        width = int(200 + 36 + w + pad)
        body = (f'<g transform="translate(0 0)">{mark(ink, GOLD)}</g>'
                f'<path transform="translate(236 117)" d="{d}" fill="{ink}"/>')
        write(f"lockup-horizontal-{theme}.svg", svg(width, 200, body))
        d, w = words.path("EPHEMERA", 22, tracking=0.32)
        width = int(max(200, w) + 2 * pad)
        body = (f'<g transform="translate({(width - 200) / 2:.1f} 0)">{mark(ink, GOLD)}</g>'
                f'<path transform="translate({(width - w) / 2:.1f} 252)" d="{d}" fill="{ink}"/>')
        write(f"lockup-stacked-{theme}.svg", svg(width, 266, body))

    # Favicons: the small mark on a paper tile; app icons: the full mark.
    def tile(content_svg, size, radius_frac=0.22, inset=0.1, fill=PAPER):
        r = size * radius_frac
        s = size * (1 - 2 * inset) / 200
        return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" '
                f'viewBox="0 0 {size} {size}"><rect width="{size}" height="{size}" rx="{r:.1f}" fill="{fill}"/>'
                f'<g transform="translate({size * inset:.2f} {size * inset:.2f}) scale({s:.4f})">'
                f'{content_svg}</g></svg>')

    for px, stroke in ((16, 22), (32, 19), (48, 17)):
        png(f"favicon-{px}.png", tile(small_mark(INK_LIGHT, GOLD, stroke), 256, inset=0.06), px, px)
    png("apple-touch-icon.png", tile(mark(INK_LIGHT, GOLD), 512, radius_frac=0, inset=0.08), 180, 180)
    png("icon-512.png", tile(mark(INK_LIGHT, GOLD), 512, inset=0.08), 512, 512)

    # Social preview (GitHub, link cards): 1280 x 640, dark.
    d, w = words.path("EPHEMERA", 46, tracking=0.32)
    t, tw = words.path("A pocket-sized, MIT-licensed ephemeris", 17, tracking=0.08)
    t2, tw2 = words.path("for astrology and astronomy", 17, tracking=0.08)
    block_w = 280 + 48 + max(w, tw)
    x0 = (1280 - block_w) / 2
    body = (f'<g transform="translate({x0:.1f} 180) scale(1.4)">{mark(INK_DARK, GOLD)}</g>'
            f'<path transform="translate({x0 + 328:.1f} 300)" d="{d}" fill="{INK_DARK}"/>'
            f'<path transform="translate({x0 + 330:.1f} 360)" d="{t}" fill="#A9ADBD"/>'
            f'<path transform="translate({x0 + 330:.1f} 392)" d="{t2}" fill="#A9ADBD"/>')
    social = svg(1280, 640, body, background=NIGHT)
    write("social-preview.svg", social)
    png("social-preview.png", social, 1280, 640)
    print("wrote", sorted(os.listdir(OUT)))


if __name__ == "__main__":
    main()
