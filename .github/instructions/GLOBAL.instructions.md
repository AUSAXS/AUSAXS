---
applyTo: "**"
---

AUSAXS is a modern C++20 project for small-angle scattering analysis. The project has only curl as a dependency which you must install yourself, and may be compiled using standard CMake commands: "cmake -B build -S .", "cmake --build build --target [tests,unit_tests,feature_tests,saxs_fitter,em_fitter,rigidbody_optimizer] -jX". Prefer not to delete build artefacts; if something is not as you would expect, first try to reinitialize CMake. Compiled tests can be run using CTest. If you want to run tests, prefer to only run unit tests unless otherwise instructed, as the feature tests are time-consuming and only meant for the CI. If you run into any compilation issues, you can trace command invocations from the base CMakeFiles.txt.

You may find yourself being invoked from a PR as a reviewer. If the instructions are unclear, it may therefore help to compare the current branch against `master` to better understand the context of the review. 

Do not add documentation files, as it is hosted online on the project wiki. If you have relevant comments, additional considerations or issues discovered during your work, include it in the PR description.
