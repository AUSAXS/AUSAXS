# SPDX-License-Identifier: LGPL-3.0-or-later
# Author: Kristian Lytje

"""
Script to build and run AUSAXS tests.

Usage:
    python run_test.py                  # Run all unit tests
    python run_test.py <test_file>      # Run specific test (auto-detect unit/feature)
    python run_test.py <test_folder>    # Run all tests in a folder, recursively (auto-detect unit/feature)
    python run_test.py <test_name>      # Run file containing given test (auto-detect host test file & unit/feature)
    python run_test.py <path>           # Run a test file or folder by path, e.g. hist/histogram_manager/

Bare names are resolved as a test file first, then a folder, then a test case name.
Since folders often share their name with a file, use a path (anything containing a '/') to select a folder.
All arguments can be prefixed with 'utest' or 'ftest' to specify test type explicitly.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path


def find_project_root():
    """Find the project root directory (where CMakeLists.txt is)."""
    current = Path(__file__).resolve().parent
    while current != current.parent:
        if (current / "CMakeLists.txt").exists():
            return current
        current = current.parent
    raise RuntimeError("Could not find project root (CMakeLists.txt)")


def test_base_dirs(test_type=None):
    """Return the (type, base directory) pairs to search for the given test type."""
    project_root = find_project_root()
    search_dirs = []
    if test_type == 'utest' or test_type is None:
        search_dirs.append(('utest', project_root / "tests" / "unit"))
    if test_type == 'ftest' or test_type is None:
        search_dirs.append(('ftest', project_root / "tests" / "feature"))
    return search_dirs


def find_path(path_str, test_type=None):
    """
    Resolve an explicit path to a test file or folder.
    The path may be relative to the working directory, the project root, or tests/unit|feature.
    
    Returns:
        List of (found_type, found_path) tuples
    """
    project_root = find_project_root()
    rel = Path(path_str)
    found = []
    for ttype, base_dir in test_base_dirs(test_type):
        for candidate in [Path.cwd() / rel, project_root / rel, base_dir / rel]:
            candidate = candidate.resolve()
            if not candidate.exists() or not candidate.is_relative_to(base_dir):
                continue
            if candidate.is_dir() or candidate.suffix == ".cpp":
                if (ttype, candidate) not in found:
                    found.append((ttype, candidate))
    return found


def find_test_file(test_name, test_type=None):
    """
    Find a test file by name.
    
    Args:
        test_name: Name of the test (without prefix/suffix)
        test_type: Either 'utest', 'ftest', or None for auto-detection
    
    Returns:
        Tuple of (found_type, found_path) or (None, None) if not found
    """
    found = []
    for ttype, base_dir in test_base_dirs(test_type):
        # Search recursively for test files
        for cpp_file in base_dir.rglob("*.cpp"):
            if cpp_file.stem == test_name:
                found.append((ttype, cpp_file))
    
    return found


def build_test(test_target, jobs=12):
    """
    Build a test target using CMake.
    
    Args:
        test_target: The CMake target name (e.g., 'utest_histogram_manager')
        jobs: Number of parallel jobs for building
    
    Returns:
        True if build succeeded, False otherwise
    """
    project_root = find_project_root()
    build_dir = project_root / "build"
    
    if not build_dir.exists():
        print(f"Error: Build directory does not exist: {build_dir}")
        print("Please run CMake configuration first:")
        print("  cmake -B build -S .")
        return False
    
    cmd = ["cmake", "--build", str(build_dir), "--target", test_target, f"-j{jobs}"]
    print(f"Building: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, cwd=project_root)
    return result.returncode == 0


def run_test_executable(test_type, test_name, test_case_filter=None):
    """
    Run a specific test executable or test case using CTest.
    
    Args:
        test_type: Either 'utest' or 'ftest'
        test_name: Name of the test file (without prefix)
        test_case_filter: Optional specific test case name to run
    
    Returns:
        The return code from the test execution
    """
    project_root = find_project_root()
    
    # Determine the correct directory based on test type
    if test_type == 'utest':
        test_dir = "unit"
    elif test_type == 'ftest':
        test_dir = "feature"
    else:
        raise ValueError(f"Invalid test type: {test_type}")
    
    test_path = project_root / "build" / "tests" / test_dir
    
    # Use CTest to run the test
    if test_case_filter:
        # Run only the specific test case
        # CTest test names are just the test case name (from TEST_CASE macro)
        cmd = [
            "ctest",
            "--output-on-failure",
            "-R", f"^{test_case_filter}$",
            "--test-dir", str(test_path)
        ]
        print(f"Running specific test case: {test_case_filter}")
    else:
        # Run all test cases in the file
        # For running all tests from a specific file, we can't easily filter by file
        # so we'll fall back to running the executable directly
        test_exe = project_root / "build" / "tests" / test_dir / "bin" / f"{test_type}_{test_name}"
        
        if not test_exe.exists():
            print(f"Error: Test executable not found: {test_exe}")
            return 1
        
        print(f"Running all test cases in: {test_name}")
        print(f"Executable: {test_exe}")
        result = subprocess.run([str(test_exe)], cwd=test_exe.parent)
        return result.returncode
    
    print(f"Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=project_root)
    return result.returncode


def run_all_tests(test_type, jobs=6, repeat=3):
    """
    Run all tests of a given type using CTest.
    
    Args:
        test_type: Either 'utest' or 'ftest'
        jobs: Number of parallel jobs
        repeat: Number of times to repeat until pass
    
    Returns:
        The return code from CTest
    """
    project_root = find_project_root()
    
    # Determine the correct directory based on test type
    if test_type == 'utest':
        test_dir = "unit"
        build_target = "unit_tests"
    elif test_type == 'ftest':
        test_dir = "feature"
        build_target = "feature_tests"
    else:
        raise ValueError(f"Invalid test type: {test_type}")
    
    # Build all tests
    if not build_test(build_target, jobs=8):
        print(f"Failed to build {test_type}")
        return 1
    
    # Run tests with CTest
    jobs = "12" if test_type == 'utest' else "1"
    test_path = project_root / "build" / "tests" / test_dir
    cmd = [
        "ctest",
        "--output-on-failure",
        "--parallel", jobs,
        "--repeat", f"until-pass:{repeat}",
        "--test-dir", str(test_path)
    ]
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=project_root)
    return result.returncode


def find_folder(folder_name, test_type=None):
    """
    Find a test folder by name.
    
    Args:
        folder_name: Name of the folder to search for
        test_type: Either 'utest', 'ftest', or None for auto-detection
    
    Returns:
        List of (found_type, found_path) tuples
    """
    found = []
    for ttype, base_dir in test_base_dirs(test_type):
        for subdir in base_dir.rglob("*"):
            if subdir.is_dir() and subdir.name == folder_name:
                found.append((ttype, subdir))
    return found


def run_tests_in_folder(folder_path, test_type, jobs=8):
    """
    Run all tests in a given folder (recursively) by building and running each test file.

    Args:
        folder_path: Path to the folder containing test .cpp files.
        test_type: Either 'utest' or 'ftest'.
        jobs: Number of parallel jobs for building.

    Returns:
        The return code (0 if all tests pass, 1 otherwise).
    """
    project_root = find_project_root()
    
    # Find all .cpp files in the folder and its subfolders
    test_files = sorted(folder_path.rglob("*.cpp"))
    
    if not test_files:
        print(f"No test files found in {folder_path}")
        return 1
    
    print(f"Found {len(test_files)} test file(s) in {folder_path.relative_to(project_root)}")
    
    # Build and run each test
    failed_tests = []
    for test_file in test_files:
        test_name = test_file.stem
        target_name = f"{test_type}_{test_name}"
        
        print(f"\n{'='*60}")
        print(f"Building and running: {test_name}")
        print(f"{'='*60}")
        
        # Build the test
        if not build_test(target_name, jobs=jobs):
            print(f"Failed to build test: {target_name}")
            failed_tests.append(test_name)
            continue
        
        # Run the test
        result = run_test_executable(test_type, test_name)
        if result != 0:
            failed_tests.append(test_name)
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Summary: {len(test_files) - len(failed_tests)}/{len(test_files)} tests passed")
    if failed_tests:
        print(f"Failed tests: {', '.join(failed_tests)}")
    print(f"{'='*60}")
    
    return 1 if failed_tests else 0


def find_test_case(test_case_name, test_type=None):
    """
    Search for a test case name in test files.

    Args:
        test_case_name: The TEST_CASE name to search for (e.g., "Atom::coordinates")
        test_type: Either 'utest', 'ftest', or None for auto-detection

    Returns:
        Tuple of (found_type, found_path) or (None, None) if not found
    """
    # Search pattern for TEST_CASE("name")
    test_pattern = re.compile(rf'TEST_CASE\s*\(\s*"({re.escape(test_case_name)})"')
    
    for ttype, base_dir in test_base_dirs(test_type):
        for cpp_file in base_dir.rglob("*.cpp"):
            try:
                with open(cpp_file, "r", encoding="utf-8") as file:
                    content = file.read()
                    if test_pattern.search(content):
                        return (ttype, cpp_file)
            except Exception as e:
                # Skip files that can't be read
                continue
    
    return (None, None)


def build_and_run_file(test_type, test_name, jobs=8):
    """Build and run the test executable of a single test file."""
    target_name = f"{test_type}_{test_name}"
    if not build_test(target_name, jobs=jobs):
        print(f"Failed to build test: {target_name}")
        return 1
    return run_test_executable(test_type, test_name)


def report_ambiguous(message, found, test_name):
    """Print the ambiguous matches and how to disambiguate them. Always returns 1."""
    project_root = find_project_root()
    print(f"Error: {message}:")
    for ttype, path in found:
        print(f"  - {ttype}: {path.relative_to(project_root)}")
    print("\nPlease disambiguate with a test type (utest or ftest) and/or a path relative to tests/unit or tests/feature, e.g.:")
    ttype, path = found[0]
    base_dir = dict(test_base_dirs())[ttype]
    rel = path.relative_to(base_dir).as_posix() + ("/" if path.is_dir() else "")
    print(f"  python {sys.argv[0]} {ttype} {rel}")
    return 1


def main():
    parser = argparse.ArgumentParser(
        description="Build and run AUSAXS tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                           # Run all unit tests
  %(prog)s histogram_manager         # Run specific test (auto-detect)
  %(prog)s utest histogram_manager   # Run specific unit test
  %(prog)s ftest histogram_manager   # Run specific feature test
  %(prog)s ftest hist/histogram_manager/  # Run all feature tests in a folder
  %(prog)s utest                     # Run all unit tests
  %(prog)s ftest                     # Run all feature tests
        """
    )
    
    parser.add_argument(
        "args",
        nargs="*",
        help="Test type (utest/ftest) and/or test name, folder, path or test case name"
    )
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=8,
        help="Number of parallel jobs for building (default: 8)"
    )
    
    args = parser.parse_args()
    
    # Parse arguments
    test_type = None
    test_name = None
    
    if len(args.args) == 0:
        # No arguments: run all unit tests
        print("Error: No arguments provided. Please specify test type or name.")
        parser.print_help()
        return 1
    elif len(args.args) == 1:
        arg = args.args[0]
        if arg in ['utest', 'ftest']:
            # Only test type specified: run all tests of that type
            test_type = arg
        else:
            # Only test name specified: auto-detect
            test_name = arg
    elif len(args.args) == 2:
        # Both test type and name specified
        if args.args[0] not in ['utest', 'ftest']:
            print(f"Error: First argument must be 'utest' or 'ftest', got: {args.args[0]}")
            return 1
        test_type = args.args[0]
        test_name = args.args[1]
    else:
        print("Error: Too many arguments")
        parser.print_help()
        return 1
    
    # Case 1: Run all tests of a specific type
    if test_name is None:
        return run_all_tests(test_type, jobs=6, repeat=3)
    
    project_root = find_project_root()

    # Case 2: An explicit path (e.g. 'hist/histogram_manager' or 'tests/feature/grid/grid.cpp')
    if "/" in test_name or test_name.endswith(".cpp"):
        found = find_path(test_name, test_type)
        if len(found) == 0:
            print(f"Error: Path '{test_name}' not found under tests/unit/ or tests/feature/")
            return 1
        if len(found) > 1:
            return report_ambiguous(f"Multiple paths match '{test_name}'", found, test_name)
        found_type, found_path = found[0]
        if found_path.is_dir():
            print(f"Found test folder: {found_path.relative_to(project_root)}")
            return run_tests_in_folder(found_path, found_type, jobs=args.jobs)
        print(f"Found test: {found_path.relative_to(project_root)}")
        return build_and_run_file(found_type, found_path.stem, jobs=args.jobs)

    # Case 3: A test file name. Takes precedence over folders, since many folders share a name with a file.
    found_tests = find_test_file(test_name, test_type)
    if len(found_tests) > 1:
        return report_ambiguous(f"Multiple tests found with name '{test_name}'", found_tests, test_name)
    if len(found_tests) == 1:
        found_type, found_path = found_tests[0]
        print(f"Found test: {found_path.relative_to(project_root)}")
        return build_and_run_file(found_type, test_name, jobs=args.jobs)

    # Case 4: A folder name
    found_folders = find_folder(test_name, test_type)
    if len(found_folders) > 1:
        return report_ambiguous(f"Multiple folders found with name '{test_name}'", found_folders, test_name)
    if len(found_folders) == 1:
        found_type, found_folder = found_folders[0]
        print(f"Found test folder: {found_folder.relative_to(project_root)}")
        return run_tests_in_folder(found_folder, found_type, jobs=args.jobs)

    # Case 5: Search for test case name within files
    found_type, found_file = find_test_case(test_name, test_type)
    if found_file:
        test_file_name = found_file.stem
        print(f"Found test case '{test_name}' in file: {found_file.relative_to(project_root)}")
        
        # Build the test
        target_name = f"{found_type}_{test_file_name}"
        if not build_test(target_name, jobs=args.jobs):
            print(f"Failed to build test: {target_name}")
            return 1
        
        # Run only the specific test case using CTest filter
        return run_test_executable(found_type, test_file_name, test_case_filter=test_name)
    
    # Nothing found
    print(f"Error: Test '{test_name}' not found")
    if test_type:
        print(f"Searched in: tests/{('unit' if test_type == 'utest' else 'feature')}/")
    else:
        print("Searched in: tests/unit/ and tests/feature/")
    print("Searched for: test files, folders, and test case names")
    return 1

if __name__ == "__main__":
    sys.exit(main())