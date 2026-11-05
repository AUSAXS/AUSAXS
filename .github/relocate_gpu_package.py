"""
Make an assembled Linux GPU runtime package loadable from wherever it is unpacked.

    python relocate_gpu_package.py <package directory>

The package is unpacked into a directory of the user's choosing and opened by absolute path, with no LD_LIBRARY_PATH
to help. Each library is then found only through the RUNPATH of the binary that needs it, so every ELF file must
carry one that reaches the top of the package. As built, several do not: the backend's points at the CI build tree,
libLLVM's at a ../lib that does not exist here, llvm-spirv and the bundled system libraries have none at all. Any of
these either fails to load, or worse, silently picks up a mismatched copy from the user's system.

Each RUNPATH is rewritten to its $ORIGIN-relative entries plus the top of the package. Absolute entries are dropped:
they name directories on the build machine, and on a user's machine can only ever find the wrong library.
Must run before UPX, since patchelf cannot edit a compressed library.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def is_elf(path: Path) -> bool:
    with open(path, "rb") as fh:
        return fh.read(4) == b"\x7fELF"


def runpath(path: Path) -> list[str]:
    out = subprocess.run(["patchelf", "--print-rpath", str(path)], capture_output=True, text=True, check=True).stdout
    return [entry for entry in out.strip().split(":") if entry]


def relocate(package: Path) -> list[str]:
    changed = []
    for f in sorted(package.rglob("*")):
        if not f.is_file() or f.is_symlink() or f.suffix == ".bc" or not is_elf(f):
            continue
        top = os.path.relpath(package, f.parent)
        to_top = "$ORIGIN" if top == "." else f"$ORIGIN/{top}"

        old = runpath(f)
        new = [entry for entry in old if entry.startswith("$ORIGIN")]
        if to_top not in (entry.rstrip("/") for entry in new):
            new.append(to_top)
        if new != old:
            subprocess.run(["patchelf", "--set-rpath", ":".join(new), str(f)], check=True)
            changed.append(f"{f.relative_to(package)}: {':'.join(old) or '(none)'} -> {':'.join(new)}")
    return changed


def main() -> int:
    if len(sys.argv) != 2 or not Path(sys.argv[1]).is_dir():
        print(__doc__)
        return 2
    if not sys.platform.startswith("linux"):
        print("relocate_gpu_package.py: only Linux packages need relocating")
        return 0
    package = Path(sys.argv[1])

    changed = relocate(package)
    lines = [f"## Relocated {package.name}", "", f"**RUNPATH rewritten** ({len(changed)}):", ""]
    lines += [f"- `{line}`" for line in changed] or ["none"]
    text = "\n".join(lines)
    print(text)
    if "GITHUB_STEP_SUMMARY" in os.environ:
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as fh:
            fh.write(text + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
