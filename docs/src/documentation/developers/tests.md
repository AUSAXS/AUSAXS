---
title: Running the tests - AUSAXS
description: Build and run the AUSAXS unit and feature tests with CMake and the run_test.py script.
---

# Running the tests

The tests are written with [Catch2](https://github.com/catchorg/Catch2), which CMake fetches automatically.

The suite is split in two:

- **Unit tests**, in `tests/unit/`, test a single class or function each. Their folders mirror `include/`.
- **Feature tests**, in `tests/feature/`, test larger parts of the program, often on real structures and data, and take longer.

Run them with the `scripts/run_test.py` script from the repository root, after configuring the build in `build/` as described in [compiling from source](../../getting-started/compiling.md). It builds the tests it needs, and runs every test case in its own process:

```bash
python scripts/run_test.py                   # every unit test
python scripts/run_test.py <test file>       # one test file, unit or feature
python scripts/run_test.py hist/             # every test in a folder
python scripts/run_test.py "<test case>"     # the file containing this test case
```

Prefix the argument with `utest` or `ftest` to choose between the unit and feature tests explicitly. Running the whole feature suite takes a long time, so run the tests relevant to your change locally and leave the rest to the CI.
