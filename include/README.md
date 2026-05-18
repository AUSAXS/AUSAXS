AUSAXS is split into several subprojects to reduce interdependencies. These are compiled as individual libraries.

`core`: the basic functionality common to all applications.
`api`: the external API exposing AUSAXS to other languages and tools (e.g. the pyausaxs Python wrapper).
`em`: classes related to electron-microscopy validation.
`gui`: helper-files for the GUI executables.
`math`: a standalone math library used in all areas of AUSAXS.
`rigidbody`: files required for rigidbody optimization.
