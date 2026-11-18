#!/usr/bin/env python3
"""Build the published documentation figures from their Inkscape sources.

Each docs/figures/<name>.svg is an Inkscape file kept for editing. This script writes two cleaned copies to
docs/src/assets/figures/: <name>.svg for the light scheme and <name>-dark.svg for the dark one. Use them in the
docs as a pair:

    ![alt](../assets/figures/<name>.svg#only-light)
    ![alt](../assets/figures/<name>-dark.svg#only-dark)

The cleanup drops editor metadata, hidden elements and unused definitions, and shortens long decimals.

The dark copy recolours only what is marked: add the token @theme to an object's label in Inkscape (Object >
Object Properties, or double-click the name in the Layers and Objects panel), e.g. "outline @theme". The object
and everything inside it then has its colours swapped according to DARK_COLOURS; unmarked objects keep theirs.
Mark what borders the transparent background (text, outer outlines) and leave inner lines alone. A marked <use>
clone does not recolour the original it points to; mark the original instead.

Usage: python docs/scripts/build_figures.py
"""
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

DOCS = Path(__file__).resolve().parent.parent
SOURCES = DOCS / "figures"
OUTPUT = DOCS / "src" / "assets" / "figures"

SVG = "http://www.w3.org/2000/svg"
XLINK = "http://www.w3.org/1999/xlink"
ET.register_namespace("", SVG)
ET.register_namespace("xlink", XLINK)

# Light-scheme colour -> dark-scheme colour. The dark page background is #0f141b (extra.css).
DARK_COLOURS = {
    "#1a1a1a": "#e6edf3",  # text and operators
    "#000000": "#c9d1d9",  # outlines
    "#dbdbdb": "#6e7681",  # excluded volume, so it does not glare
}
THEME_MARK = re.compile(r"(^|\s)@theme(\s|$)")
LABEL = "{http://www.inkscape.org/namespaces/inkscape}label"
COLOUR_ATTRIBUTES = ("style", "fill", "stroke", "stop-color", "color")


def strip_editor_data(root: ET.Element) -> None:
    """Remove Inkscape/Sodipodi elements and attributes, which only the editor uses."""
    for e in root.iter():
        for key in [k for k in e.attrib if "inkscape" in k or "sodipodi" in k]:
            del e.attrib[key]
        for child in [c for c in e if not c.tag.startswith("{" + SVG)]:
            e.remove(child)


def remove_hidden(root: ET.Element) -> None:
    parents = {c: p for p in root.iter() for c in p}
    for e in list(root.iter()):
        if "display:none" in (e.get("style") or "").replace(" ", "") and e in parents:
            parents[e].remove(e)


def remove_unused_defs(root: ET.Element) -> None:
    """Drop definitions nothing refers to, repeating until none are left (a definition can use another)."""
    for defs in root.iter(f"{{{SVG}}}defs"):
        while True:
            text = ET.tostring(root, encoding="unicode")
            used = set(re.findall(r"url\(#([^)]+)\)", text)) | set(re.findall(r'href="#([^"]+)"', text))
            unused = [d for d in defs if d.get("id") not in used]
            if not unused:
                break
            for d in unused:
                defs.remove(d)


def shorten_numbers(text: str) -> str:
    """Round to 6 significant digits; Inkscape writes up to 17, which the rendering never needs."""
    def fmt(m: re.Match) -> str:
        s = f"{float(m.group(0)):.6g}"
        return "0" if s == "-0" else s
    return re.sub(r"-?\d+\.\d{5,}(?:e-?\d+)?", fmt, text)


def recolour_marked(root: ET.Element) -> int:
    """Swap DARK_COLOURS in every object labelled with @theme, and in everything inside it. Returns the count."""
    pattern = re.compile("|".join(re.escape(k) for k in DARK_COLOURS), re.IGNORECASE)
    marked = [e for e in root.iter() if THEME_MARK.search(e.get(LABEL) or "")]
    for top in marked:
        for e in top.iter():
            for attr in COLOUR_ATTRIBUTES:
                if attr in e.attrib:
                    e.set(attr, pattern.sub(lambda m: DARK_COLOURS[m.group(0).lower()], e.get(attr)))
    return len(marked)


def clean(root: ET.Element) -> str:
    strip_editor_data(root)
    remove_hidden(root)
    remove_unused_defs(root)
    return shorten_numbers(ET.tostring(root, encoding="unicode"))


def build(source: Path) -> None:
    light = clean(ET.parse(source).getroot())
    dark_root = ET.parse(source).getroot()
    if recolour_marked(dark_root) == 0:
        print(f"warning: nothing in {source.name} is labelled @theme, so the dark copy equals the light one",
              file=sys.stderr)
    dark = clean(dark_root)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, text in [(source.stem, light), (source.stem + "-dark", dark)]:
        out = OUTPUT / f"{name}.svg"
        out.write_text(text + "\n", encoding="utf-8")
        print(f"{source.relative_to(DOCS)} -> {out.relative_to(DOCS)} ({len(text) // 1024} KB)")


def main() -> int:
    sources = sorted(SOURCES.glob("*.svg"))
    if not sources:
        print(f"No sources in {SOURCES}", file=sys.stderr)
        return 1
    for source in sources:
        build(source)
    return 0


if __name__ == "__main__":
    sys.exit(main())
