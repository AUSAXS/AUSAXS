---
applyTo: "tests/**"
---

Watch out for potential issues in the production code. It is up to you what you consider issues in the library. Segfaults is definitely an issue that must be fixed, but you can also fix function signature and API issues if you think they are severe enough. Missing includes leading to inconsistent behaviour (such as triggering missing includes for common use-cases) should also be evaluated on a case-to-case basis; on one hand we do not want to inflate the dependency trees, but it should also be relatively simple for users to actually use the main library features. Log any such changes to the GitHub PR. 

For the test file structure, see e.g. "tests/feature/utility/string_utils.cpp" and "include/core/utility/StringUtils.h" for inspiration. Note the use of Section headers and the intentional lack of comments; tests should be more or less self-descriptive. This is not a hard rule; if you believe something is non-trivial, feel free to add comments.

We currently have two test directories: feature tests and unit tests. The first is meant for larger integration tests relying on multiple components working seamlessly together. The latter are smaller tests validating individual methods or simple collective functionality. If you add new tests, consider to which class it belongs and should be added. 

CMake uses "glob" to discover new test files. You should therefore make sure to re-initialize CMake every time you add a new test to make it compilable. You should compile all tests to ensure the build works, but prefer only to run unit tests, as the whole feature test suite is too large and time-consuming for you. If you add new feature tests, run only the test file you updated to verify it works.
