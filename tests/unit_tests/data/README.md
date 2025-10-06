# Unit Tests for AUSAXS Core/Data

This directory contains unit tests for the core data structures and classes in the AUSAXS library.

## Test Files

- **unit_atom.cpp** - Tests for `Atom`, `AtomFF`, and `Water` classes
  - Constructor tests with various parameter combinations
  - Coordinate and weight accessor tests
  - Equality operator tests
  
- **unit_body.cpp** - Tests for the `Body` class
  - Constructor tests (default, vectors, copy, move)
  - Accessor methods (`get_atoms`, `get_atom`, `get_waters`)
  - Center of mass calculations
  - Mass and charge calculations
  - Translation operations
  - Size methods
  
- **unit_state_manager.cpp** - Tests for the `StateManager` class
  - Constructor and initialization
  - Internal and external modification tracking
  - Hydration layer modification tracking
  - Probe management
  - State reset functionality
  
- **unit_signaller.cpp** - Tests for `BoundSignaller` and `UnboundSignaller` classes
  - Constructor tests
  - Modification signaling (internal and external)
  - Integration with StateManager
  
- **unit_symmetry.cpp** - Tests for the `Symmetry` class and related structures
  - Constructor tests for Symmetry, _Relation, and _Repeat
  - Closed symmetry detection
  - Transform generation
  - Equality operators

## Building and Running

To build a specific unit test:
```bash
cmake --build build --target test_unit_atom
```

To run a specific unit test:
```bash
./build/tests/bin/test_unit_atom
```

To build all unit tests:
```bash
cmake --build build --target test_unit_atom test_unit_body test_unit_state_manager test_unit_signaller test_unit_symmetry
```

## Test Structure

These unit tests follow the structure established in the repository:
- Each test file corresponds to one or more related header files
- Tests use Catch2 framework with Section-based organization
- No comments unless necessary for clarity
- Tests are self-descriptive through naming and structure

## Notes

- Unit tests focus on individual class methods and basic functionality
- Feature tests (integration tests) remain in the main `tests/data` directory
- Default constructors are not tested as they may result in uninitialized memory for trivial types
