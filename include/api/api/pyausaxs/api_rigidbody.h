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
extern "C" API int rigidbody_get_preview_structure(
    int rigidbody_id,
    double** x, double** y, double** z,
    int** body_index, int** copy_index, int** residue_seq, int** is_ca,
    int* n_atoms,
    int* status
);

// Latest structure published by an `update structure` element during a run (coordinates only, in
// the same atom order as rigidbody_get_preview_structure). `version` is bumped on each publish and
// reset to 0 when a new sequence is parsed; poll it to detect new frames. Thread-safe.
extern "C" API int rigidbody_get_live_structure(
    double** x, double** y, double** z,
    int* n_atoms, int* version,
    int* status
);

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