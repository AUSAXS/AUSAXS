This directory contains all of the unit and feature tests for the entire project. They can be built using `cmake --build build --target tests` and then run through either CMake itself or the accompanying makefile. 

When contributing new tests, typical pitfalls are the following settings:
 * `settings::molecule::center`: This setting automatically centers the Molecule upon construction, thus modifying the atomic coordinates.
 * `settings::molecule::implicit_hydrogens`: This setting modifies the atomic charges to account for attached hydrogens. This will often throw exceptions when used in combination with simple dummy atoms. 
