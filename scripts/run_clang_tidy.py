# SPDX-License-Identifier: LGPL-3.0-or-later
# Author: Kristian Lytje

"""
Script to run clang-tidy over AUSAXS sources.

The project is normally configured with GCC, whose compile_commands.json contains
flags clang rejects (-fconstexpr-ops-limit, some -march values). This script writes
a sanitized copy of the compilation database and drives clang-tidy against it, so
no separate clang build directory is needed.

Findings in headers are reported once per translation unit that includes them, so
identical diagnostics are collapsed and each one is shown only the first time.

Requires clang-tidy 19 or newer; see .clang-tidy for why. The binary is taken from
--clang-tidy, else $CLANG_TIDY, else the newest clang-tidy on PATH.

Usage:
    python run_clang_tidy.py                            # whole project
    python run_clang_tidy.py source/rigidbody           # a subtree
    python run_clang_tidy.py source/core/utility/Logging.cpp
    python run_clang_tidy.py --fix source/core/utility  # apply fixes
    python run_clang_tidy.py --checks='-*,performance-*' source/core
    CLANG_TIDY=clang-tidy-21 python run_clang_tidy.py   # pin the toolchain
"""

import argparse
import json
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# Flags that GCC accepts but clang does not, or that only slow the parse down.
UNSUPPORTED_FLAGS = re.compile(
    r"^(-fconstexpr-ops-limit=.*"
    r"|-fno-finite-math-only"
    r"|-flto.*"
    r"|-pipe"
    r"|-w)$"
)

# clang-tidy always parses with clang, so clang-specific workarounds must be applied
# even when the build directory was configured with GCC.
EXTRA_COMPILER_FLAGS = [
    # dlib trips this on clang >= 17; CMakeLists.txt applies the same workaround.
    "-Wno-missing-template-arg-list-after-template-kw",
    # The compile database may name warnings the running clang does not know.
    "-Wno-unknown-warning-option",
]

# misc-include-cleaner.MissingIncludes, which .clang-tidy relies on to keep the
# check usable alongside the *Fwd.h headers, landed in clang-tidy 19. Older
# binaries silently ignore the option and bury the run in false positives.
MINIMUM_VERSION = 19


def find_project_root():
    """Find the project root directory (where CMakeLists.txt is)."""
    current = Path(__file__).resolve().parent
    while current != current.parent:
        if (current / "CMakeLists.txt").exists():
            return current
        current = current.parent
    raise RuntimeError("Could not find project root (CMakeLists.txt)")


def clang_tidy_version(binary):
    """
    Query the major version of a clang-tidy binary.

    Args:
        binary: Path to the binary.

    Returns:
        The major version as an int, or None if it could not be determined.
    """
    try:
        out = subprocess.run(
            [binary, "--version"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            text=True, timeout=30
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return None

    match = re.search(r"version (\d+)\.", out)
    return int(match.group(1)) if match else None


def discover_candidates():
    """
    Find every clang-tidy binary reachable from PATH.

    Returns:
        List of (major version, path) tuples, newest first. Unversioned binaries
        that fail to report a version sort last.
    """
    names = ["clang-tidy"] + [f"clang-tidy-{v}" for v in range(30, MINIMUM_VERSION - 1, -1)]

    seen = set()
    found = []
    for name in names:
        path = shutil.which(name)
        if not path or path in seen:
            continue
        seen.add(path)
        found.append((clang_tidy_version(path) or -1, path))

    return sorted(found, key=lambda x: -x[0])


def find_clang_tidy(override=None):
    """
    Locate a clang-tidy binary new enough for this project's configuration.

    The binary is chosen in this order:
        1. --clang-tidy=<path> on the command line
        2. the CLANG_TIDY environment variable
        3. the newest clang-tidy[-N] on PATH

    Args:
        override: Explicit path supplied on the command line, or None. An explicit
            choice is honoured even if it is older than MINIMUM_VERSION, with a
            warning, so unusual toolchains stay usable.

    Returns:
        Tuple of (path to the binary, major version or None).
    """
    override = override or os.environ.get("CLANG_TIDY")

    if override:
        if not shutil.which(override) and not Path(override).is_file():
            raise RuntimeError(f"clang-tidy not found at: {override}")
        version = clang_tidy_version(override)
        if version is not None and version < MINIMUM_VERSION:
            print(
                f"Warning: {override} is clang-tidy {version}; this project's .clang-tidy "
                f"expects {MINIMUM_VERSION} or newer.\n"
                f"         misc-include-cleaner will report a large number of false "
                f"positives against the *Fwd.h headers.\n"
            )
        return override, version

    candidates = discover_candidates()
    if not candidates:
        raise RuntimeError(
            "No clang-tidy binary found on PATH.\n"
            "Install one with e.g.:\n"
            "  sudo apt install clang-tidy-21        # Debian/Ubuntu\n"
            "  brew install llvm                     # macOS\n"
            "or point at an existing binary with --clang-tidy=<path>."
        )

    version, path = candidates[0]
    if version < MINIMUM_VERSION:
        raise RuntimeError(
            f"Found clang-tidy {version if version > 0 else '(unknown version)'} at {path}, "
            f"but this project requires {MINIMUM_VERSION} or newer.\n"
            "Install a newer one with e.g.:\n"
            "  sudo apt install clang-tidy-21\n"
            "or override the check with --clang-tidy=<path>."
        )

    return path, version


def sanitize_compile_commands(project_root, build_dir):
    """
    Write a clang-compatible copy of the compilation database.

    Entries for fetched dependencies are dropped, and GCC-only flags are stripped.

    Args:
        project_root: Path to the project root.
        build_dir: Path to the CMake build directory.

    Returns:
        Path to the directory holding the sanitized compile_commands.json.
    """
    source_db = build_dir / "compile_commands.json"
    if not source_db.exists():
        raise RuntimeError(
            f"No compilation database at {source_db}\n"
            "Configure the project first:\n"
            "  cmake -B build -S . -DCMAKE_EXPORT_COMPILE_COMMANDS=ON"
        )

    with source_db.open() as f:
        entries = json.load(f)

    sanitized = []
    for entry in entries:
        if "_deps" in entry["file"]:
            continue

        if "command" in entry:
            parts = [p for p in entry["command"].split() if not UNSUPPORTED_FLAGS.match(p)]
            entry["command"] = " ".join(parts)
        else:
            entry["arguments"] = [
                a for a in entry["arguments"] if not UNSUPPORTED_FLAGS.match(a)
            ]
        sanitized.append(entry)

    if not sanitized:
        raise RuntimeError(f"No project translation units found in {source_db}")

    out_dir = build_dir / "clang-tidy"
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "compile_commands.json").open("w") as f:
        json.dump(sanitized, f, indent=1)

    return out_dir, {Path(e["file"]).resolve() for e in sanitized}


def select_files(project_root, targets, known_files):
    """
    Resolve command line targets to translation units present in the database.

    Args:
        project_root: Path to the project root.
        targets: List of file or directory arguments; empty means everything.
        known_files: Set of resolved paths present in the compilation database.

    Returns:
        Sorted list of paths to analyse.
    """
    if not targets:
        return sorted(known_files)

    selected = set()
    for target in targets:
        path = (project_root / target).resolve() if not Path(target).is_absolute() else Path(target)
        if path.is_dir():
            selected |= {f for f in known_files if f.is_relative_to(path)}
        elif path in known_files:
            selected.add(path)
        else:
            print(f"Warning: not in the compilation database, skipping: {target}")

    return sorted(selected)


# clang-tidy's own progress chatter, which carries no information once the
# header filter is in place.
NOISE_LINES = re.compile(
    r"^(\d+ warnings? generated\."
    r"|Suppressed \d+ warnings? .*"
    r"|Use -header-filter=.*"
    r"|Use -system-headers .*)$"
)


# Start of a diagnostic block: path:line:col: <severity>: message
DIAGNOSTIC_HEAD = re.compile(r"^(?P<loc>.+?:\d+:\d+): (?P<severity>warning|error|note): ")


def split_diagnostics(output):
    """
    Split clang-tidy output into individual diagnostics.

    Each diagnostic keeps its trailing source snippet and any 'note:' lines that
    belong to it.

    Args:
        output: Raw clang-tidy stdout.

    Returns:
        List of (key, severity, text) tuples, where key identifies the diagnostic
        by location and message so the same header finding can be recognised
        across translation units.
    """
    diagnostics = []
    for line in output.splitlines():
        if NOISE_LINES.match(line.strip()):
            continue

        match = DIAGNOSTIC_HEAD.match(line)
        if match and match.group("severity") != "note":
            diagnostics.append([line, match.group("severity"), [line]])
        elif diagnostics:
            diagnostics[-1][2].append(line)

    return [(key, severity, "\n".join(body)) for key, severity, body in diagnostics]


def run_one(clang_tidy, db_dir, path, project_root, extra_args):
    """
    Run clang-tidy on a single translation unit.

    Returns:
        Tuple of (path, list of (key, severity, text) diagnostics).
    """
    cmd = [
        clang_tidy, "-p", str(db_dir), *extra_args,
        *(f"--extra-arg={flag}" for flag in EXTRA_COMPILER_FLAGS),
        str(path),
    ]
    result = subprocess.run(
        cmd, cwd=project_root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    return path, split_diagnostics(result.stdout)


def main():
    parser = argparse.ArgumentParser(
        description="Run clang-tidy over AUSAXS sources.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "targets", nargs="*",
        help="Files or directories to analyse. Defaults to the whole project."
    )
    parser.add_argument(
        "--fix", action="store_true",
        help="Apply clang-tidy's suggested fixes in place."
    )
    parser.add_argument(
        "--checks",
        help="Override the check list from .clang-tidy for this run."
    )
    parser.add_argument(
        "--clang-tidy", metavar="PATH",
        help="clang-tidy binary to use. Also settable via the CLANG_TIDY environment "
             "variable; otherwise the newest clang-tidy on PATH is used."
    )
    parser.add_argument(
        "--build-dir", default="build",
        help="CMake build directory containing compile_commands.json (default: build)."
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=multiprocessing.cpu_count(),
        help="Number of files to analyse in parallel."
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Only print files that produced diagnostics."
    )
    args = parser.parse_args()

    project_root = find_project_root()
    build_dir = project_root / args.build_dir

    try:
        clang_tidy, version = find_clang_tidy(args.clang_tidy)
        db_dir, known_files = sanitize_compile_commands(project_root, build_dir)
    except RuntimeError as e:
        print(f"Error: {e}")
        return 1

    files = select_files(project_root, args.targets, known_files)
    if not files:
        print("No matching translation units to analyse.")
        return 1

    extra_args = []
    if args.fix:
        extra_args += ["--fix", "--fix-errors"]
    if args.checks:
        extra_args += [f"--checks={args.checks}"]

    label = f" (version {version})" if version and version > 0 else ""
    print(f"Using: {clang_tidy}{label}")
    print(f"Analysing {len(files)} file(s) with {args.jobs} job(s)\n")

    # A finding in a header is reported once per translation unit that includes it,
    # so the same issue can surface hundreds of times. Report each one once.
    seen = set()
    flagged = []
    duplicates = 0
    per_check = {}

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [
            pool.submit(run_one, clang_tidy, db_dir, f, project_root, extra_args)
            for f in files
        ]
        for i, future in enumerate(futures, start=1):
            path, diagnostics = future.result()
            rel = path.relative_to(project_root)

            fresh = []
            for key, severity, text in diagnostics:
                if key in seen:
                    duplicates += 1
                    continue
                seen.add(key)
                fresh.append(text)

                check = re.search(r"\[([\w.-]+)\]$", text.splitlines()[0])
                if check:
                    per_check[check.group(1)] = per_check.get(check.group(1), 0) + 1

            if not args.quiet or fresh:
                print(f"[{i}/{len(files)}] {rel}")
            if fresh:
                print("\n".join(fresh))
                print()
            if fresh:
                flagged.append((rel, len(fresh)))

    if not flagged:
        print(f"\nNo diagnostics in {len(files)} file(s).")
        return 0

    total = sum(count for _, count in flagged)
    print(f"\n{total} unique diagnostic(s) across {len(flagged)} of {len(files)} file(s)")
    if duplicates:
        print(f"({duplicates} repeat(s) of header findings already shown were suppressed)")

    print("\nBy check:")
    for check, count in sorted(per_check.items(), key=lambda x: -x[1]):
        print(f"  {count:5d}  {check}")

    return 1


if __name__ == "__main__":
    sys.exit(main())
