file(REMOVE_RECURSE
  "../../lib/libausaxs_static.a"
  "../../lib/libausaxs_static.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/ausaxs_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
