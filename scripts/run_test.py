# SPDX-License-Identifier: LGPL-3.0-or-later
# Author: Kristian Lytje

"""
Script to build and run AUSAXS tests.

Usage:
    python run_test.py                  # Run all unit tests
    python run_test.py <test_file>      # Run specific test (auto-detect unit/feature)
    python run_test.py <test_folder>    # Run all tests in a folder (auto-detect unit/feature)
    python run_test.py <test_name>      # Run file containing given test (auto-detect host test file & unit/feature)

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


def find_test_file(test_name, test_type=None):
    """
    Find a test file by name.
    
    Args:
        test_name: Name of the test (without prefix/suffix)
        test_type: Either 'utest', 'ftest', or None for auto-detection
    
    Returns:
        Tuple of (found_type, found_path) or (None, None) if not found
    """
    project_root = find_project_root()
    
    search_dirs = []
    if test_type == 'utest' or test_type is None:
        search_dirs.append(('utest', project_root / "tests" / "unit"))
    if test_type == 'ftest' or test_type is None:
        search_dirs.append(('ftest', project_root / "tests" / "feature"))
    
    found = []
    for ttype, base_dir in search_dirs:
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


def run_test_executable(test_type, test_name):
    """
    Run a specific test executable.
    
    Args:
        test_type: Either 'utest' or 'ftest'
        test_name: Name of the test (without prefix)
    
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
    
    test_exe = project_root / "build" / "tests" / test_dir / "bin" / f"{test_type}_{test_name}"
    
    if not test_exe.exists():
        print(f"Error: Test executable not found: {test_exe}")
        return 1
    
    print(f"Running: {test_exe}")
    result = subprocess.run([str(test_exe)], cwd=test_exe.parent)
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
        Tuple of (found_type, found_path) or (None, None) if not found
    """
    project_root = find_project_root()
    
    search_dirs = []
    if test_type == 'utest' or test_type is None:
        search_dirs.append(('utest', project_root / "tests" / "unit"))
    if test_type == 'ftest' or test_type is None:
        search_dirs.append(('ftest', project_root / "tests" / "feature"))
    
    for ttype, base_dir in search_dirs:
        # Check direct subfolder
        folder_path = base_dir / folder_name
        if folder_path.is_dir():
            return (ttype, folder_path)
        
        # Check recursively
        for subdir in base_dir.rglob("*"):
            if subdir.is_dir() and subdir.name == folder_name:
                return (ttype, subdir)
    
    return (None, None)


def run_tests_in_folder(folder_path, test_type, jobs=8):
    """
    Run all tests in a given folder by building and running each test file.

    Args:
        folder_path: Path to the folder containing test .cpp files.
        test_type: Either 'utest' or 'ftest'.
        jobs: Number of parallel jobs for building.

    Returns:
        The return code (0 if all tests pass, 1 otherwise).
    """
    project_root = find_project_root()
    
    # Find all .cpp files in the folder
    test_files = list(folder_path.glob("*.cpp"))
    
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
    project_root = find_project_root()
    
    # Search pattern for TEST_CASE("name")
    test_pattern = re.compile(rf'TEST_CASE\s*\(\s*"({re.escape(test_case_name)})"')
    
    search_dirs = []
    if test_type == 'utest' or test_type is None:
        search_dirs.append(('utest', project_root / "tests" / "unit"))
    if test_type == 'ftest' or test_type is None:
        search_dirs.append(('ftest', project_root / "tests" / "feature"))
    
    for ttype, base_dir in search_dirs:
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
  %(prog)s utest                     # Run all unit tests
  %(prog)s ftest                     # Run all feature tests
        """
    )
    
    parser.add_argument(
        "args",
        nargs="*",
        help="Test type (utest/ftest) and/or test name"
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
    
    # Case 2: Check if the provided name matches a folder
    found_type, found_folder = find_folder(test_name, test_type)
    if found_folder:
        print(f"Found test folder: {found_folder.relative_to(find_project_root())}")
        return run_tests_in_folder(found_folder, found_type, jobs=args.jobs)
    
    # Case 3: Check if it's a test file name
    found_tests = find_test_file(test_name, test_type)
    
    if len(found_tests) == 0:
        # Case 4: Search for test case name within files
        found_type, found_file = find_test_case(test_name, test_type)
        if found_file:
            test_file_name = found_file.stem
            print(f"Found test case '{test_name}' in file: {found_file.relative_to(find_project_root())}")
            
            # Build the test
            target_name = f"{found_type}_{test_file_name}"
            if not build_test(target_name, jobs=args.jobs):
                print(f"Failed to build test: {target_name}")
                return 1
            
            # Run the test
            return run_test_executable(found_type, test_file_name)
        
        # Nothing found
        print(f"Error: Test '{test_name}' not found")
        if test_type:
            print(f"Searched in: tests/{('unit' if test_type == 'utest' else 'feature')}/")
        else:
            print("Searched in: tests/unit/ and tests/feature/")
        print("Searched for: test files, folders, and test case names")
        return 1
    
    if len(found_tests) > 1:
        print(f"Error: Multiple tests found with name '{test_name}':")
        for ttype, path in found_tests:
            rel_path = path.relative_to(find_project_root())
            print(f"  - {ttype}: {rel_path}")
        print("\nPlease specify the test type (utest or ftest) to disambiguate:")
        print(f"  python {sys.argv[0]} utest {test_name}")
        print(f"  python {sys.argv[0]} ftest {test_name}")
        return 1
    
    # Exactly one test found
    found_type, found_path = found_tests[0]
    rel_path = found_path.relative_to(find_project_root())
    print(f"Found test: {rel_path}")
    
    # Build the test
    target_name = f"{found_type}_{test_name}"
    if not build_test(target_name, jobs=args.jobs):
        print(f"Failed to build test: {target_name}")
        return 1
    
    # Run the test
    return run_test_executable(found_type, test_name)


if __name__ == "__main__":
    sys.exit(main())