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
//! Deprecated!
extern "C" API int rigidbody_get_live_structure(
    double** x, double** y, double** z,
    int* n_atoms, int* version,
    int* status
);

// Create a live-structure poller. The returned id owns the buffers that rigidbody_poll_live_structure
// writes into, and must be deallocated once when the consumer is done polling.
extern "C" API int rigidbody_create_live_poller(
    int* status
);

// Copy the latest structure published by an `update structure` element into the poller's own buffers, and
// hand out pointers to them. The pointers stay valid until the next poll on the same poller, or until the
// poller is deallocated.
extern "C" API void rigidbody_poll_live_structure(
    int poller_id,
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

// Ask a running refinement to stop. The request is picked up at the start of the next loop iteration, after which
// rigidbody_run returns normally with the best conformation found so far. Safe to call from another thread than the
// one running the refinement, and a no-op if nothing is running.
extern "C" API void rigidbody_stop_run(
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

// Display names of the inline (positional) arguments an element accepts, in script order, e.g. {"old name", "new name"}
// for "rename". A name ending in "..." marks the slot absorbing the remaining arguments. min_count/max_count give the
// accepted range; a max_count of 100 means the element takes an open-ended list.
extern "C" API void rigidbody_get_valid_inline_arguments(
    const char* element_name,
    const char*** arguments,
    int* size,
    int* min_count,
    int* max_count,
    int* status
);

// Display names of the bodies that remain after the setup elements (merge/delete/convert_to_symmetry)
// have been applied, ordered by body index so that name i corresponds to the body with body_index == i
// reported by rigidbody_get_preview_structure. Each name is the body's custom alias if set, else "bN".
extern "C" API void rigidbody_get_body_names(
    int rigidbody_id,
    const char*** names,
    int* size,
    int* status
);

// Symmetry layout of every replica (copy > 0) in the structure, keyed to rigidbody_get_preview_structure's
// (body, copy) pairs so the caller never has to guess the copy -> symmetry mapping.
// Parallel arrays, one entry per replica:
//   body[]     matches preview_structure's body_index.
//   copy[]     matches preview_structure's copy_index (always > 0; copy 0 is the unreplicated original).
//   symmetry[] 0-based index of the symmetry within its body (a body can have several stacked symmetries).
//   replica[]  1-based repetition index of this copy within its symmetry.
//   type[]     predefined-name string of the symmetry (e.g. "c4", "p2", "d3-c2"), matching the names
//              accepted by the "symmetry" sequencer element.
//   name[]     current addressable name of this replica, e.g. "b1s1r1" or its custom alias after rename.
extern "C" API int rigidbody_get_symmetry_layout(
    int rigidbody_id,
    int** body, int** copy, int** symmetry, int** replica,
    const char*** type, const char*** name,
    int* n_replicas,
    int* status
);