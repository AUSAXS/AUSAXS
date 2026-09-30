# SPDX-License-Identifier: LGPL-3.0-or-later
# Author: Kristian Lytje

import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys
import time


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Run each unit-test executable with a per-file timeout."
    )
    parser.add_argument(
        "test_directory",
        type=Path,
        help="Directory containing utest_* executables.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=5.0,
        help="Maximum seconds allowed per executable (default: 5).",
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        required=True,
        help="Directory for timing.csv and per-executable logs.",
    )
    return parser.parse_args()


def test_executables(test_directory):
    return sorted(
        path for path in test_directory.iterdir()
        if path.name.startswith("utest_") and path.is_file() and os.access(path, os.X_OK)
    )


def main():
    arguments = parse_arguments()
    test_directory = arguments.test_directory.resolve()
    output_directory = arguments.output_directory.resolve()
    logs_directory = output_directory / "logs"

    if not test_directory.is_dir():
        print(f"Unit-test directory does not exist: {test_directory}", file=sys.stderr)
        return 2

    executables = test_executables(test_directory)
    if not executables:
        print(f"No unit-test executables found in: {test_directory}", file=sys.stderr)
        return 2

    logs_directory.mkdir(parents=True, exist_ok=True)
    timing_path = output_directory / "timing.csv"
    results = []

    for executable in executables:
        log_path = logs_directory / f"{executable.name}.log"
        start = time.monotonic()
        status = "passed"

        with log_path.open("w") as log_file:
            try:
                completed = subprocess.run(
                    [str(executable), "--reporter", "compact"],
                    cwd=test_directory,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    timeout=arguments.timeout,
                    check=False,
                )
                if completed.returncode != 0:
                    status = f"failed ({completed.returncode})"
            except subprocess.TimeoutExpired:
                status = "timed out"

        elapsed_seconds = time.monotonic() - start
        results.append((executable.name, status, elapsed_seconds, log_path))
        print(f"{executable.name}: {status} ({elapsed_seconds:.3f}s)")

    with timing_path.open("w", newline="") as timing_file:
        writer = csv.writer(timing_file)
        writer.writerow(("executable", "status", "elapsed_seconds", "log"))
        for executable, status, elapsed_seconds, log_path in results:
            writer.writerow((executable, status, f"{elapsed_seconds:.3f}", log_path))

    unsuccessful = [result for result in results if result[1] != "passed"]
    print(f"\nTimed {len(results)} unit-test executables; {len(unsuccessful)} did not pass.")
    return 1 if unsuccessful else 0


if __name__ == "__main__":
    raise SystemExit(main())