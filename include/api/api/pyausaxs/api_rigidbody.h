// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <api/api_helper.h>

extern "C" API int rigidbody_load_script(
    const char* script,
    int* status
);

// Explicit structure (symmetries realized) annotated with the per-atom metadata needed to draw
// a preview: which body each atom belongs to, which symmetry copy (0 = original), its residue
// number, and whether it is a Cα.
// Also returns constraints as a flat array of (n_constraints * 3) ints: [idx1, idx2, type, ...],
// where idx1/idx2 are indices into the returned atom arrays, and type is:
//   0 = backbone/bond, 1 = center-of-mass, 2 = attractor, 3 = repulsor.
// constraint_data is null when n_constraints == 0.
extern "C" API int rigidbody_get_preview_structure(
    int rigidbody_id,
    double** x, double** y, double** z,
    int** body_index, int** copy_index, int** residue_seq, int** is_ca,
    int* n_atoms,
    int** constraint_data, int* n_constraints,
    int* status
);

// Latest structure published by an `update structure` element during a run.
extern "C" API int rigidbody_get_live_structure(
    double** x, double** y, double** z,
    int* n_atoms, int* version,
    int* status
);

// Register or unregister as a live consumer. `update` elements are no-ops unless this is set to true. 
extern "C" API void rigidbody_register_live_consumer(bool connected, int* status);

extern "C" API void rigidbody_validate(
    int rigidbody_id,
    int* status
);

extern "C" API int rigidbody_run(
    int rigidbody_id,
    double** q, double** I, double** I_err, double** I_interp, int* n_points,
    int* status
);

extern "C" API void rigidbody_get_valid_elements(
    const char*** elements,
    int* size,
    int* status
);

extern "C" API void rigidbody_get_valid_arguments(
    const char* element_name,
    const char*** arguments,
    int* size,
    int* status
);