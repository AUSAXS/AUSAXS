#!/usr/bin/env python3
"""Build the animation shown on the rigid-body refinement application page.

Refines the five domains of LAR1-2 against its SAXS data, saving the structure after every accepted move as an .xyz
trajectory, and animates the Cα trace of the trajectory. (A save directly inside optimize_once would also record the
trial conformations of rejected moves.) Writes a pair to docs/src/assets/figures/:

    rigidbody-trajectory.gif, rigidbody-trajectory-dark.gif

Use them as a pair:

    ![alt](../assets/figures/rigidbody-trajectory.gif#only-light)
    ![alt](../assets/figures/rigidbody-trajectory-dark.gif#only-dark)

The refinement is stochastic, so the animation changes with every build. The data are read from data/rigidbody/LAR1-2-full,
which is not under version control.

Needs the ausaxs command-line tool, matplotlib and pillow. Set AUSAXS to the executable if it is not on the path.
Usage: python docs/scripts/build_rigidbody_animation.py
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from build_figures import OUTPUT
from build_fit_figure import THEMES

DOCS = Path(__file__).resolve().parent.parent
DATA = DOCS.parent / "data" / "rigidbody" / "LAR1-2-full"
AUSAXS = os.environ.get("AUSAXS", "ausaxs")
SPLITS = [7, 99, 199, 292]  # first residue of every body but the first
ITERATIONS = 1000

BACKGROUND = {"light": "#ffffff", "dark": "#0f141b"}
BODY_COLOURS = {
    "light": ["#57606a", "#b07c12", "#2f6f9f", "#8a4f7d", "#4e8a5b"],
    "dark":  ["#8b949e", "#e3b341", "#6cb6ff", "#d2a8ff", "#7ee787"],
}

CONFIG = f"""\
output out
load {{
	pdb {DATA / "LAR1-2-full.pdb"}
	saxs {DATA / "LAR1-2-full.dat"}
	split {" ".join(map(str, SPLITS))}
}}
autoconstrain backbone
save trajectory.xyz
parameters {{
	iterations {ITERATIONS}
	translate 10
	rotate 0.5
}}
loop
	optimize_once
		on_improvement
			save trajectory.xyz
		end
	end
end
"""


def ca_atoms(pdb: Path) -> tuple[np.ndarray, np.ndarray]:
    """The indices of the Cα atoms among the ATOM records, and the body each belongs to."""
    indices, bodies, i = [], [], 0
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() == "CA":
            indices.append(i)
            bodies.append(int(np.searchsorted(SPLITS, int(line[22:26]), side="right")))
        i += 1
    return np.array(indices), np.array(bodies)


def read_trajectory(path: Path, indices: np.ndarray) -> np.ndarray:
    """The coordinates of the given atoms in every frame, as an array of shape (frames, atoms, 3)."""
    lines = path.read_text().splitlines()
    n = int(lines[0])
    frames = []
    for start in range(0, len(lines), n + 2):
        rows = [lines[start + 2 + i].split() for i in indices]
        frames.append(np.array([[float(x) for x in r[1:4]] for r in rows]))
    return np.array(frames)


def project(frames: np.ndarray) -> np.ndarray:
    """Project every frame onto the two principal axes of the first frame."""
    centre = frames[0].mean(axis=0)
    _, _, axes = np.linalg.svd(frames[0] - centre, full_matrices=False)
    return (frames - centre) @ axes[:2].T


def animate(xy: np.ndarray, bodies: np.ndarray, steps: list[int], scheme: str, out: Path) -> None:
    # Size the canvas to the extent of the molecule over the whole trajectory.
    lo, hi = xy[steps].reshape(-1, 2).min(axis=0) - 4, xy[steps].reshape(-1, 2).max(axis=0) + 4
    width = 4
    fig, ax = plt.subplots(figsize=(width, width * (hi[1] - lo[1]) / (hi[0] - lo[0])), dpi=100)
    fig.patch.set_facecolor(BACKGROUND[scheme])
    ax.set_facecolor(BACKGROUND[scheme])
    ax.axis("off")
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    fig.subplots_adjust(0, 0, 1, 1)

    lines = [ax.plot([], [], lw=1.6, color=BODY_COLOURS[scheme][b], solid_capstyle="round")[0]
             for b in range(bodies.max() + 1)]

    # Hold the first and last frames a little longer.
    sequence = [steps[0]] * 6 + steps + [steps[-1]] * 15

    def draw(step):
        for b, line in enumerate(lines):
            line.set_data(*xy[step, bodies == b].T)
        return lines

    FuncAnimation(fig, draw, frames=sequence, blit=True).save(out, writer=PillowWriter(fps=8))
    plt.close(fig)


def main() -> int:
    indices, bodies = ca_atoms(DATA / "LAR1-2-full.pdb")
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        (work / "lar.conf").write_text(CONFIG)
        subprocess.run([AUSAXS, "rigidbody", "lar.conf"], cwd=work, check=True, stdout=subprocess.DEVNULL)
        frames = read_trajectory(work / "out" / "trajectory.xyz", indices)

    xy = project(frames)
    steps = list(range(len(frames)))
    print(f"{len(frames) - 1} accepted moves")
    for scheme in THEMES:
        out = OUTPUT / ("rigidbody-trajectory.gif" if scheme == "light" else "rigidbody-trajectory-dark.gif")
        animate(xy, bodies, steps, scheme, out)
        print(f"-> {out.relative_to(DOCS)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
