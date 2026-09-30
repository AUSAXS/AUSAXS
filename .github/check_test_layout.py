# SPDX-License-Identifier: LGPL-3.0-or-later
# Author: Kristian Lytje

"""
Checks that every unit and feature test file corresponds to a header, so each header's tests can be found from its path alone.

The rule, for a header include/<module>/<dir>/<ClassName>.h (included as <dir>/<ClassName>.h):
    tests/<suite>/<dir>/<class_name>.cpp           the tests of that header, or of its header/source pair
    tests/<suite>/<dir>/<class_name>.<topic>.cpp   an additional file for the same header, split out by topic
where <suite> is unit or feature. The file name is the snake_case form of the header name; a lowercase name is what tells a test apart
from the implementation in a file search. Every name must also be unique within its suite, since each file becomes an executable named
utest_<file name> or ftest_<file name>. Shared test helpers belong in tests/support, never in a suite directory.

Run with --untested to list every header that has no unit test file yet. Feature tests are deliberately left out of that list: they test
aggregate behaviour and the connections between units, so they are not expected to exist for every header.
"""

import os
import re
import sys

include_root = "include"
suites = {"unit": os.path.join("tests", "unit"), "feature": os.path.join("tests", "feature")}
snake_case = re.compile(r"^[a-z0-9]+(_[a-z0-9]+)*$")

# headers that only gather other headers or declare names; they have no behaviour of their own to test
untestable = re.compile(r"(Fwd|All|ExportMacro)$")

def normalize(name):
    return name.lower().replace("_", "")

# (directory as seen from an #include, normalized header name) -> header path
headers = {}
for module in sorted(os.listdir(include_root)):
    module_root = os.path.join(include_root, module)
    if not os.path.isdir(module_root):
        continue
    for root, _, files in os.walk(module_root):
        for file in files:
            if file.endswith(".h"):
                directory = os.path.relpath(root, module_root)
                headers[(directory, normalize(file[:-2]))] = os.path.join(root, file)

flag_fail = False
tested = {suite: set() for suite in suites}
for suite, test_root in suites.items():
    test_names = {}
    for root, _, files in os.walk(test_root):
        for file in sorted(files):
            file_path = os.path.join(root, file)
            if file == "CMakeLists.txt":
                continue
            if not file.endswith(".cpp"):
                flag_fail = True
                print(f"{file_path}: only test sources belong in {test_root}; move shared helpers to tests/support")
                continue

            stem = file[:-4]
            name, _, topic = stem.partition(".")
            if not snake_case.match(name) or (topic and not snake_case.match(topic)):
                flag_fail = True
                print(f"{file_path}: must be named <header_name>.cpp or <header_name>.<topic>.cpp, in lowercase snake_case")
                continue

            if stem in test_names:
                flag_fail = True
                print(f"{file_path}: its executable would clash with {test_names[stem]}")
            test_names[stem] = file_path

            directory = os.path.relpath(root, test_root)
            header = headers.get((directory, normalize(name)))
            if header is None:
                flag_fail = True
                print(f"{file_path}: no matching header {directory}/<{name} in CamelCase>.h under {include_root}/<module>/")
                continue
            tested[suite].add(header)

if "--untested" in sys.argv[1:]:
    testable = sorted(h for h in headers.values() if not untestable.search(os.path.basename(h)[:-2]))
    untested = [h for h in testable if h not in tested["unit"]]
    for header in untested:
        print(f"untested: {header}")
    print(f"{len(untested)} of {len(testable)} headers have no unit test file.")

if flag_fail:
    print("Some test files do not follow the naming convention; see the top of .github/check_test_layout.py.")
    exit(1)

print("All unit and feature test files correspond to a header.")
exit(0)
