# SPDX-License-Identifier: LGPL-3.0-or-later
# Author: Kristian Lytje

"""
Run clang-tidy over the AUSAXS sources, exiting non-zero if anything is reported.

The project is normally configured with GCC, whose compile_commands.json contains
flags clang rejects (-fconstexpr-ops-limit, -flto). A sanitized copy is written to
build/clang-tidy, so no separate clang build directory is needed. Only a configure
step is required, not a build:

    cmake -B build -S .
    python scripts/run_clang_tidy.py                            # whole project
    python scripts/run_clang_tidy.py source/rigidbody           # a subtree
    python scripts/run_clang_tidy.py --checks='-*,modernize-*'  # one group at a time
    python scripts/run_clang_tidy.py --fix include/api          # fix one batch
    python scripts/run_clang_tidy.py --summary=report.md        # also write markdown

Findings are printed twice: once as they arrive, so a long run says something while it
works, and once at the end grouped by the file they point at, which is the form that
says what to go and fix. --summary=FILE appends that same grouping as markdown, and is
how CI turns a failure into a report -- see .github/workflows/build-and-test-template.yml,
which passes $GITHUB_STEP_SUMMARY.

Naming an include/ folder also selects the sources mirroring it, so a batch is one
folder: --fix include/api rewrites include/api and source/api and nothing else.

Requires clang-tidy 19 or newer; see .clang-tidy. Set $CLANG_TIDY to pick a binary.
"""

import json
import multiprocessing
import os
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from sys import argv

root = Path(__file__).resolve().parents[1]
build_dir = root / "build"
clang_tidy = os.environ.get("CLANG_TIDY", "clang-tidy")

if not shutil.which(clang_tidy):
    exit(f"{clang_tidy} not found on PATH. Install it, or set $CLANG_TIDY to a binary.")

# misc-include-cleaner.MissingIncludes, which .clang-tidy sets to keep the check usable
# alongside the *Fwd.h headers, landed in clang-tidy 19. Older binaries silently ignore
# the option and bury the run in false positives, so refuse them rather than mislead.
version = re.search(r"version (\d+)\.", subprocess.run(
    [clang_tidy, "--version"], capture_output=True, text=True).stdout)
if version and int(version.group(1)) < 19:
    exit(f"{clang_tidy} is clang-tidy {version.group(1)}; .clang-tidy requires 19 or newer.")

fix = "--fix" in argv[1:]
checks = [a for a in argv[1:] if a.startswith("--checks=")]
targets = [a for a in argv[1:] if not a.startswith("-")]
summary = next((a.split("=", 1)[1] for a in argv[1:] if a.startswith("--summary=")), None)

# flags GCC accepts but clang does not, or that only slow the parse down
unsupported = re.compile(r"^(-fconstexpr-ops-limit=.*|-fno-finite-math-only|-flto.*|-pipe|-w)$")

# clang-tidy always parses with clang, so clang-specific workarounds are needed even
# though the build directory was configured with GCC: dlib trips the first on
# clang >= 17 (CMakeLists.txt applies the same workaround), and the compile database
# may name warnings the running clang does not know
extra_args = [
    "--extra-arg=-Wno-missing-template-arg-list-after-template-kw",
    "--extra-arg=-Wno-unknown-warning-option",
]

# a diagnostic starts here; the lines after it (source snippet, "note:") belong to it.
# the trailing [check] is optional because clang-tidy omits it on the compiler errors it
# forwards, and the message is matched lazily so a message containing its own brackets
# still yields the *last* one as the check name.
diagnostic = re.compile(
    r"^(?P<file>.+?):(?P<line>\d+):\d+: (?:warning|error): "
    r"(?P<message>.*?)(?: \[(?P<check>[\w.,-]+)\])?$"
)

# progress chatter, which carries no information once the header filter is in place
noise = re.compile(
    r"^(\d+ warnings? generated\.|Suppressed \d+ warnings? .*"
    r"|Use -header-filter=.*|Use -system-headers .*)$"
)

source_db = build_dir / "compile_commands.json"
if not source_db.exists():
    exit(f"No compilation database at {source_db}. Configure first: cmake -B build -S .")

entries = [e for e in json.loads(source_db.read_text()) if "_deps" not in e["file"]]
if not entries:
    exit(f"No project translation units found in {source_db}")

for entry in entries:
    if "command" in entry:
        entry["command"] = " ".join(p for p in entry["command"].split() if not unsupported.match(p))
    else:
        entry["arguments"] = [a for a in entry["arguments"] if not unsupported.match(a)]

db_dir = build_dir / "clang-tidy"
db_dir.mkdir(parents=True, exist_ok=True)
(db_dir / "compile_commands.json").write_text(json.dumps(entries, indent=1))

files = sorted(Path(e["file"]).resolve() for e in entries)
prefixes = set()
if targets:
    prefixes = {(root / t).resolve() for t in targets}
    # a header is not a translation unit, so pull in the source that implements it: the
    # mirrored directory for a folder, the mirrored file for a single header. Header-only
    # headers have neither, and are instead covered by any selected unit that includes
    # them, as HeaderFilterRegex in .clang-tidy lets those findings through.
    for t in (Path(t) for t in targets if t.startswith("include/")):
        mirror = root / "source" / t.relative_to("include")
        prefixes |= {mirror, mirror.with_suffix(".cpp")}
    files = [f for f in files if any(f.is_relative_to(p) for p in prefixes)]
if not files:
    # not an error: a change may legitimately touch no translation unit at all
    print(f"Nothing to analyse in: {' '.join(targets)}")
    exit(0)


def run(file):
    """Run clang-tidy on a single translation unit, returning its diagnostics."""
    cmd = [clang_tidy, "-p", str(db_dir), *extra_args, *checks]
    if fix:
        cmd += ["--fix", "--fix-errors"]
        # confine header rewrites to the batch, so fixing one folder leaves the headers
        # of every other folder untouched and each batch stays reviewable on its own
        if prefixes:
            scope = "|".join(sorted(re.escape(str(p)) for p in prefixes))
            cmd.append(f"-header-filter=({scope})/.*")

    # clang-tidy reports diagnostics on stderr, interleaved with its own chatter
    out = subprocess.run([*cmd, str(file)], cwd=root, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True).stdout

    found = []
    for line in out.splitlines():
        if noise.match(line.strip()):
            continue
        if diagnostic.match(line):
            found.append([line])
        elif found:
            found[-1].append(line)
    return found


# fixes to a shared header must not be applied from several translation units at once
jobs = 1 if fix else multiprocessing.cpu_count()
print(f"Using {clang_tidy} on {len(files)} file(s) with {jobs} job(s)\n")

# a finding in a header is reported once per translation unit that includes it, so
# identical diagnostics are collapsed and each is shown only the first time
seen = set()
per_check = Counter()
# the file a finding points at is not the unit that produced it -- a header finding is
# attributed to the header -- so the grouping below is keyed on the former
per_file = defaultdict(lambda: defaultdict(list))
with ThreadPoolExecutor(max_workers=jobs) as pool:
    for i, (file, found) in enumerate(zip(files, pool.map(run, files)), start=1):
        print(f"[{i}/{len(files)}] {file.relative_to(root)}")
        for block in found:
            if block[0] in seen:
                continue
            seen.add(block[0])
            print("\n".join(block))
            head = diagnostic.match(block[0])
            check = head["check"] or "(no check)"
            per_check[check] += 1
            where = Path(head["file"])
            if where.is_absolute() and where.is_relative_to(root):
                where = where.relative_to(root)
            per_file[str(where)][check].append((int(head["line"]), head["message"]))

if not per_check:
    print(f"\nNo diagnostics in {len(files)} file(s).")
    exit(0)


def by_file():
    """The findings grouped per file and check, worst file and worst check first."""
    for name in sorted(per_file, key=lambda f: (-sum(map(len, per_file[f].values())), f)):
        checks = sorted(per_file[name].items(), key=lambda kv: (-len(kv[1]), kv[0]))
        yield name, sum(len(hits) for _, hits in checks), checks


headline = (f"{sum(per_check.values())} unique diagnostic(s) in {len(per_file)} file(s), "
            f"from {len(files)} translation unit(s) analysed")
print(f"\n{headline}:")
for name, total, checks in by_file():
    print(f"\n{name}  ({total})")
    for check, hits in checks:
        print(f"    {check}: {len(hits)}")
        for line, message in hits:
            print(f"        L{line}: {message}")

print("\ntotals per check:")
for check, count in per_check.most_common():
    print(f"  {count:5d}  {check}")

if summary:
    report = ["## clang-tidy", "", headline + ".", "", "| count | check |", "| ----: | :---- |"]
    report += [f"| {count} | `{check}` |" for check, count in per_check.most_common()]
    for name, total, checks in by_file():
        report += ["", f"<details open><summary><code>{name}</code> — {total}</summary>", "", "```"]
        for check, hits in checks:
            report.append(f"{check}: {len(hits)}")
            report += [f"    L{line}: {message}" for line, message in hits]
        report += ["```", "", "</details>"]
    with open(summary, "a", encoding="utf-8") as f:
        f.write("\n".join(report) + "\n")

exit(1)
