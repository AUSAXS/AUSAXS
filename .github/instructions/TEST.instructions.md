---
applyTo: "tests/**"
---

Watch out for potential issues in the production code. It is up to you what you consider issues in the library. Segfaults is definitely an issue that must be fixed, but you can also fix function signature and API issues if you think they are severe enough. Missing includes leading to inconsistent behaviour (such as triggering missing includes for common use-cases) should also be evaluated on a case-to-case basis; on one hand we do not want to inflate the dependency trees, but it should also be relatively simple for users to actually use the main library features. Log any such changes to the GitHub PR. 

For the test file structure, see e.g. "tests/feature/utility/string_utils.cpp" and "include/core/utility/StringUtils.h" for inspiration. Note the use of Section headers and the intentional lack of comments; tests should be more or less self-descriptive. This is not a hard rule; if you believe something is non-trivial, feel free to add comments.

We currently have two test directories: feature tests and unit tests. The first is meant for larger integration tests relying on multiple components working seamlessly together. The latter are smaller tests validating individual methods or simple collective functionality. If you add new tests, consider to which class it belongs and should be added. 

For easily running tests, use the "python scripts/run_test.py" script. See the first few lines for documentation. This will both reinitialize CMake, compile the tests (with updates), and run the provided test. If you do not use this script, CMake must be reinitialized manually every time you add a new test file, and you must compile the tests before running them. Prefer not to run the whole feature test suite, as it is quite large and time-consuming: run only your own feature tests. Unit tests are fine to run in their entirety, as they are quite fast.