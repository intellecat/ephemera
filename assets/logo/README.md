# Ephemera logo

A lowercase *e* drawn as a compass construction: the bowl and bar as one
stroke, dashed construction circles and cross-hairs, a ring of inward
ticks, and a gold planet with a small orbit at the end of the bar.

| File | Use |
|---|---|
| `mark-light.svg`, `mark-dark.svg` | the mark on light or dark backgrounds |
| `mark-black.svg`, `mark-white.svg` | one colour |
| `mark-small.svg` | the bold small-size mark (favicons, below ~64 px); follows the browser's light or dark scheme |
| `lockup-horizontal-*.svg`, `lockup-stacked-*.svg` | mark and wordmark |
| `favicon-16/32/48.png`, `apple-touch-icon.png`, `icon-512.png` | icons |
| `social-preview.png` | link previews, e.g. the GitHub repository's social preview (1280 × 640) |

Colours: ink `#1B1F33` (on light) or `#E8E4DA` (on dark), gold `#D4A23F`,
paper `#F5F2EA`, night `#10131F`. The wordmark is Jost Light with wide
letter spacing, converted to outlines (Jost: SIL Open Font License), so
the files need no font. Regenerate everything with
`python3 tools/make_logo.py --font Jost[wght].ttf`.
