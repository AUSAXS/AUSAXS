"""
Reduce an assembled GPU runtime package to what is loaded at run time.

    python prune_gpu_package.py <package directory>

acpp --acpp-deploy recreates every development and version alias of a shared library as a symlink beside it.
Neither upload-artifact nor a wheel preserves symlinks, so each alias would ship as a full copy: libLLVM four
times over, which is most of what made the Linux package 400 MB. Each library is therefore collapsed into one
regular file under the name the dynamic loader asks for: its SONAME on Linux, its install name on macOS.
Windows deploys plain DLLs without aliases, so there the collapse finds nothing to do.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

# Link-time and debug artefacts. Nothing in the package links against these at run time.
DROP_SUFFIXES = {".a", ".lib", ".exp", ".pdb"}

# glibc's libmvec. acpp links it into kernels JIT-compiled for its OpenMP host device, which the backend never
# selects (it asks for sycl::gpu_selector_v), and a copy of glibc taken from the build machine must not be paired
# with whatever glibc the user has.
DROP_PREFIXES = ("libmvec.so",)


def is_shared_library(path: Path) -> bool:
    name = path.name
    if sys.platform == "darwin":
        return name.endswith(".dylib")
    return name.endswith(".so") or ".so." in name


def loader_name(path: Path) -> str | None:
    """The file name the dynamic loader looks the library up by, or None if it records none."""
    if sys.platform == "darwin":
        out = subprocess.run(["otool", "-D", str(path)], capture_output=True, text=True).stdout.splitlines()
        return os.path.basename(out[1].strip()) if len(out) > 1 else None
    out = subprocess.run(["readelf", "-d", str(path)], capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "(SONAME)" in line:
            return line.split("[", 1)[1].rstrip("]")
    return None


def size_of(package: Path) -> int:
    return sum(f.stat().st_size for f in package.rglob("*") if f.is_file() and not f.is_symlink())


def collapse_aliases(package: Path) -> list[str]:
    libs = sorted(f for f in package.rglob("*") if (f.is_file() or f.is_symlink()) and is_shared_library(f))
    materialise: dict[Path, Path] = {}  # loader-name path -> the real file to put there
    aliases: list[Path] = []
    for lib in libs:
        real = lib.resolve()
        if not real.is_file():
            continue  # dangling; removed below
        name = loader_name(real)
        if name is None:
            continue
        target = lib.parent / name
        if lib != target or lib.is_symlink():
            aliases.append(lib)
        if not target.is_file() or target.is_symlink():
            materialise.setdefault(target, real)

    # Copied aside before anything is deleted, since the real file is usually one of the aliases.
    staged = {target: target.with_name(target.name + ".staging") for target in materialise}
    for target, staging in staged.items():
        shutil.copy2(materialise[target], staging)
    for alias in aliases:
        alias.unlink()
    for target, staging in staged.items():
        staging.rename(target)
    return [str(a.relative_to(package)) for a in aliases if not (a.exists() or a.is_symlink())]


def drop_unneeded(package: Path) -> list[str]:
    dropped = []
    for f in sorted(package.rglob("*")):
        if not (f.is_file() or f.is_symlink()):
            continue
        dangling = f.is_symlink() and not f.exists()
        if dangling or f.suffix.lower() in DROP_SUFFIXES or f.name.startswith(DROP_PREFIXES):
            f.unlink()
            dropped.append(str(f.relative_to(package)))
    return dropped


def strip_symbols(package: Path) -> list[str]:
    """Linux only: macOS would need every binary re-signed afterwards, and Windows keeps symbols in .pdb files."""
    if not sys.platform.startswith("linux"):
        return []
    stripped = []
    for f in sorted(package.rglob("*")):
        if not f.is_file() or f.is_symlink():
            continue
        with open(f, "rb") as fh:
            if fh.read(4) != b"\x7fELF":
                continue
        sections = subprocess.run(["readelf", "-S", str(f)], capture_output=True, text=True).stdout
        if ".symtab" in sections:
            subprocess.run(["strip", "--strip-unneeded", str(f)], check=True)
            stripped.append(str(f.relative_to(package)))
    return stripped


def main() -> int:
    if len(sys.argv) != 2 or not Path(sys.argv[1]).is_dir():
        print(__doc__)
        return 2
    package = Path(sys.argv[1])
    before = size_of(package)

    report = [
        ("Dropped", drop_unneeded(package)),
        ("Collapsed aliases", collapse_aliases(package)),
        ("Stripped", strip_symbols(package)),
    ]
    remaining = sorted(str(f.relative_to(package)) for f in package.rglob("*") if f.is_symlink())

    after = size_of(package)
    lines = [f"## Pruned {package.name}", "", f"{before / 1e6:.1f} MB -> {after / 1e6:.1f} MB uncompressed", ""]
    for title, files in report:
        lines.append(f"**{title}** ({len(files)}): " + (", ".join(f"`{f}`" for f in files) or "none"))
        lines.append("")
    if remaining:
        # upload-artifact would ship each of these as a full copy of its target.
        lines.append(f"**Symlinks left** ({len(remaining)}): " + ", ".join(f"`{f}`" for f in remaining))
    text = "\n".join(lines)
    print(text)
    if "GITHUB_STEP_SUMMARY" in os.environ:
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as fh:
            fh.write(text + "\n")
    if remaining:
        print(f"::warning::{len(remaining)} symlinks are left in the package and will ship as copies")
    return 0


if __name__ == "__main__":
    sys.exit(main())
