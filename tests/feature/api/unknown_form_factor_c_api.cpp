// Test case to reproduce the segfault with UNKNOWN form factors and Fraser ExV model
// This uses the C API directly, similar to how Python bindings work

#include <catch2/catch_test_macros.hpp>

#include <api/api_pyausaxs.h>
#include <api/pyausaxs/api_settings.h>

TEST_CASE("test_unknown_form_factor_c_api: UNKNOWN form factors with Fraser exv model via C API") {
    // Create atoms without form factor information (like molecule_from_arrays does)
    double x[] = {0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, -1.0};
    double y[] = {0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 1.0, -1.0};
    double z[] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 1.0, -1.0};
    double w[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    int n_atoms = 9;

    int status = 0;

    // Create molecule from arrays
    int mol_id = molecule_from_arrays(x, y, z, w, n_atoms, &status);
    REQUIRE(status == 0);
    REQUIRE(mol_id >= 0);

    // Set problematic settings via C API
    set_exv_settings("Fraser", &status);
    REQUIRE(status == 0);

    set_fit_settings(10, 100, true, true, true, false, false, &status);
    REQUIRE(status == 0);

    // This should fail when trying to create the histogram because Fraser ExV requires form factor info
    double *aa, *aw, *ww, *axis;
    int n_bins;
    molecule_distance_histogram(mol_id, &aa, &aw, &ww, &axis, &n_bins, &status);

    // We expect this to fail because the Fraser model requires form factor information
    REQUIRE(status != 0);

    char* error_msg = nullptr;
    int error_status = 0;
    get_last_error_msg(&error_msg, &error_status);
    REQUIRE(error_msg != nullptr);

    std::string error_str(error_msg);
    CHECK(error_str.find("UNKNOWN form factor") != std::string::npos);
}
