// Test case to reproduce the segfault with UNKNOWN form factors and Fraser ExV model
// This uses the C API directly, similar to how Python bindings work

#include <catch2/catch_test_macros.hpp>

#include <api/api_pyausaxs.h>
#include <api/api_settings.h>

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
    
    // This should trigger the crash if the bug exists
    double *aa, *aw, *ww;
    int n_bins;
    double delta_r;
    int hist_id = molecule_distance_histogram(mol_id, &aa, &aw, &ww, &n_bins, &delta_r, &status);
    
    if (status != 0) {
        char* error_msg = nullptr;
        int error_status = 0;
        get_last_error_msg(&error_msg, &error_status);
        FAIL("Failed to get histogram: " << (error_msg ? error_msg : "Unknown error"));
    }
    
    REQUIRE(hist_id >= 0);
    REQUIRE(n_bins > 0);
    
    // Try debye transform - this should fail with an informative error about UNKNOWN form factors
    double *q, *I;
    int n_points;
    molecule_debye(mol_id, &q, &I, &n_points, &status);
    
    // We expect this to fail because the Fraser model requires form factor information
    REQUIRE(status != 0);
    
    char* error_msg = nullptr;
    int error_status = 0;
    get_last_error_msg(&error_msg, &error_status);
    REQUIRE(error_msg != nullptr);
    
    std::string error_str(error_msg);
    CHECK(error_str.find("UNKNOWN form factor") != std::string::npos);
}
