#!/usr/bin/python3
"""Draw a protein as a stylized cartoon for the fit figures, by tracing a PyMOL rendering.

PyMOL renders the cartoon twice in the same view: once as black outlines only, and once as a flat silhouette. Both
are traced into vector paths: the silhouette becomes the yellow fill, with its small interior holes removed, and the
outlines are drawn in black on top. This is the style of the lysozyme drawing in figures/saxs-decomposition.svg.

Writes figures/<name>-molecule.svg. Needs PyMOL, so run it with the system Python that has it:

    /usr/bin/python3 docs/scripts/build_molecule_drawing.py <structure file> <name>
"""
import sys
import tempfile
from pathlib import Path

import contourpy
import numpy as np
from PIL import Image
from pymol import cmd

DOCS = Path(__file__).resolve().parent.parent
FILL = "#e3b341"
SIZE = 1000          # width and height of the rendering, in pixels
PAD = 6              # margin around the drawing, in pixels
MIN_HOLE = 600       # interior holes in the fill smaller than this (pixels²) are filled in
MIN_SPECK = 8        # line fragments smaller than this (pixels²) are dropped
TOLERANCE = 0.7      # largest deviation of the traced paths from the rendering, in pixels


def render(structure: str, folder: Path) -> tuple[np.ndarray, np.ndarray]:
    """The outline and silhouette renderings, as greyscale arrays (0 black, 1 white)."""
    cmd.load(structure, "m")
    cmd.remove("not polymer.protein or hydro")
    cmd.remove(f"not chain {cmd.get_chains('m')[0]}")
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.set("cartoon_oval_width", 0.35)
    cmd.set("cartoon_oval_length", 1.6)
    cmd.set("cartoon_loop_radius", 0.5)
    cmd.orient("m")
    cmd.zoom("m", 0)
    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    cmd.set("antialias", 2)

    cmd.set("ray_trace_mode", 2)
    cmd.set("ray_trace_gain", 0.4)
    cmd.png(str(folder / "lines.png"), width=SIZE, height=SIZE, ray=1)

    cmd.set("ray_trace_mode", 0)
    cmd.color("black", "m")
    for setting, value in (("ambient", 1), ("direct", 0), ("reflect", 0), ("specular", 0), ("ray_shadows", 0), ("depth_cue", 0)):
        cmd.set(setting, value)
    cmd.png(str(folder / "mask.png"), width=SIZE, height=SIZE, ray=1)

    read = lambda name: np.asarray(Image.open(folder / name).convert("L"), dtype=float) / 255
    return read("lines.png"), read("mask.png")


def simplify(ring: np.ndarray, tolerance: float) -> np.ndarray:
    """Ramer-Douglas-Peucker simplification of a closed ring."""
    keep = np.zeros(len(ring), dtype=bool)
    keep[[0, len(ring) // 2, -1]] = True
    stack = [(0, len(ring) // 2), (len(ring) // 2, len(ring) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        a, b = ring[i], ring[j]
        d = b - a
        norm = np.hypot(*d)
        seg = ring[i + 1:j] - a
        dist = np.abs(d[0] * seg[:, 1] - d[1] * seg[:, 0]) / norm if norm > 0 else np.hypot(seg[:, 0], seg[:, 1])
        k = int(np.argmax(dist))
        if dist[k] > tolerance:
            keep[i + 1 + k] = True
            stack += [(i, i + 1 + k), (i + 1 + k, j)]
    return ring[keep]


def area(ring: np.ndarray) -> float:
    x, y = ring[:, 0], ring[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def trace(image: np.ndarray, min_outer: float, min_hole: float) -> str:
    """SVG path data for the dark regions of image, with holes, as even-odd rings."""
    generator = contourpy.contour_generator(z=image, fill_type="OuterOffset")
    points, offsets = generator.filled(-1.0, 0.5)
    rings = []
    for polygon, offset in zip(points, offsets):
        parts = [polygon[a:b] for a, b in zip(offset[:-1], offset[1:])]
        if area(parts[0]) < min_outer:
            continue
        rings.append(parts[0])
        rings += [hole for hole in parts[1:] if area(hole) >= min_hole]
    data = []
    for ring in rings:
        ring = simplify(ring, TOLERANCE)
        if len(ring) < 3:
            continue
        coordinates = " ".join(f"{x - PAD_X:.1f} {y - PAD_Y:.1f}" for x, y in ring[:-1])
        data.append(f"M{coordinates}Z")
    return "".join(data)


def main() -> int:
    global PAD_X, PAD_Y
    structure, name = sys.argv[1], sys.argv[2]
    with tempfile.TemporaryDirectory() as folder:
        lines, mask = render(structure, Path(folder))

    rows, cols = np.nonzero(mask < 0.5)
    y0, y1 = rows.min() - PAD, rows.max() + PAD
    x0, x1 = cols.min() - PAD, cols.max() + PAD
    PAD_X, PAD_Y = x0, y0
    fill = trace(mask, min_outer=MIN_SPECK, min_hole=MIN_HOLE)
    outline = trace(lines, min_outer=MIN_SPECK, min_hole=0)

    width, height = x1 - x0, y1 - y0
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}">\n'
           f'<g id="molecule">\n'
           f'<path fill="{FILL}" fill-rule="evenodd" d="{fill}"/>\n'
           f'<path fill="#000000" fill-rule="evenodd" d="{outline}"/>\n'
           f'</g>\n</svg>\n')
    out = DOCS / "figures" / f"{name}-molecule.svg"
    out.write_text(svg, encoding="utf-8")
    print(f"-> {out.relative_to(DOCS)} ({len(svg) // 1024} kB, {width}x{height})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
