# Extending the AUSAXS symmetry system

## Context

AUSAXS rigid-body refinement currently supports **cyclic symmetries** (c2–c6) and a free
dimer (p2). A symmetry generates symmetric duplicates of a body and refines their relative
positions against SAS data. The critical efficiency property is **distance reuse**: in a c3
group A–B–C the inter-copy distances AB, BC, CA are identical, so only one needs to be
histogrammed (with an integer multiplicity). This reuse is currently hard-coded for the
cyclic case in two histogram managers.

We want to add three symmetry classes:

1. **Polyhedral symmetry** — tetrahedral (T, 12 elements) and, generically, octahedral
   (O, 24) and icosahedral (I, 60) chiral rotation groups. No mirror operations.
2. **Reference-symmetry** — one symmetry that replicates a *group* of bodies together
   (a molecule split across several PDB files). It behaves exactly like the existing
   cyclic symmetry, except each "spot" of the group holds multiple bodies; the combined
   centre-of-mass of the participating bodies fixes their relative transforms.
3. **Nested symmetry** — e.g. `p2-3`: a body has a p2 partner, and the pair is itself one
   member of a larger c3 group → 6 copies {A,A',B,B',C,C'}.

The hard part for all three is preserving and generalizing the distance-reuse machinery.

## Current architecture (verified)

- `ISymmetry` (`include/core/data/symmetry/ISymmetry.h`) — abstract. `get_transform(cm, rep)`
  returns the affine map for copy `rep`; `repetitions()`, `is_closed()`, `add()`, `clone()`,
  `span_translation()`/`span_rotation()`.
- `CyclicSymmetry`, `PointSymmetry` — concrete types. `PredefinedSymmetries.h` maps the
  `enum class type` and name strings to factory-built objects.
- Each `Body` owns a `SymmetryStorage` (`vector<unique_ptr<ISymmetry>>`); symmetries in the
  list are applied **independently** to the base atoms (flat, not nested).
- `generate_transformed_data` (`SymmetryHelpers.cpp`) and `explicit_structure`
  (`BodySymmetryFacade.cpp`) build copy `k` *purely* by calling `get_transform(cm, k+1)` in a
  flat loop. **`get_transform` is the single chokepoint** — any group a symmetry can express
  through it becomes a valid set of copies with no change to those call sites.
- Distance reuse is duplicated verbatim in `SymmetryManagerMT.cpp:60-75` and
  `PartialSymmetryManagerMT.cpp:562-577` (`calc_aa`, `isym2==0` branch). The scale formula
  `scale = repetitions - irepeat (+1 if irepeat==0 && closed)`, loop bound
  `repetitions - closed`. Every *other* cross-pairing in both managers uses `scale=1` and
  full enumeration — i.e. no reuse for sym-vs-sym or sym-vs-other-body pairs today.
- `SimpleCalculator` dispatches a runtime `scale` to compile-time templates for scale 1..30.
- `PartialSymmetryManagerMT` is the manager actually installed during refinement
  (`SymmetryElement.cpp:28`); `SymmetryManagerMT` is the whole-histogram reference. Both
  must stay bit-for-bit consistent.

## Core design: move reuse bookkeeping into `ISymmetry`

The central change. Both managers stop computing scales themselves and instead ask each
symmetry for a **pair schedule** — the list of representative `(repA, repB, scale)` jobs;
every pair not listed is a duplicate of a listed one. One mechanism then covers cyclic,
polyhedral, nested, and reference cases.

Add to `ISymmetry` (`include/core/data/symmetry/ISymmetry.h`):

```cpp
namespace ausaxs::symmetry {
    struct CopyPair { int repA; int repB; int scale; }; // rep 0 = original body, 1..N = copies

    class ISymmetry {
        // ... existing ...
        // distinct inter-copy pairs within {original, copy_1, ..., copy_N} of THIS symmetry
        virtual std::vector<CopyPair> internal_pair_schedule() const;          // base = cyclic logic
        // distinct pairs between copies of two bodies sharing the SAME symmetry object
        virtual std::vector<CopyPair> cross_pair_schedule() const;             // base = full enum, scale 1
    };
}
```

Both are non-pure; their base implementations live in a new
`source/core/data/symmetry/ISymmetry.cpp`. The base `internal_pair_schedule()` reproduces
`SymmetryManagerMT.cpp:60-75` exactly, so `CyclicSymmetry`/`PointSymmetry` need no override:

```cpp
std::vector<CopyPair> ISymmetry::internal_pair_schedule() const {
    std::vector<CopyPair> out;
    bool closed = is_closed();
    int reps = (int)repetitions();
    for (int k = 0; k < reps - int(closed); ++k) {
        int scale = reps - k;
        if (k == 0 && closed) scale += 1;
        out.push_back({0, k+1, scale});
    }
    return out;
}
```

**Shared bucketer.** New `include/core/data/symmetry/PairSchedule.h` +
`source/core/data/symmetry/PairScheduleUtil.cpp`:
`compute_pair_schedule(const std::vector<AffineTransform>& placements)` enumerates all
`C(n,2)` unordered placement pairs, computes each pair's relative transform, canonicalizes
it (rotation matrix + translation rounded to ~1e-4, taking the lexicographically smaller of
the transform and its inverse since distance is symmetric), buckets, and emits one
`CopyPair` representative per bucket with `scale = bucket size`. Polyhedral, composite, and
reference symmetries delegate to this. Schedules depend only on fixed group structure (not
on optimisable offsets), so each symmetry computes its schedule **once at construction** and
caches it.

**Manager rewiring.** Replace `SymmetryManagerMT.cpp:60-75` and the `calc_aa` `isym2==0`
branch in `PartialSymmetryManagerMT.cpp:562-577` with a loop over
`sym->internal_pair_schedule()`. Replace the hard-coded sym-vs-sym double loops
(`SymmetryManagerMT.cpp:90-96,117-123`; `PartialSymmetryManagerMT.cpp:619-635`) with
`sym->cross_pair_schedule()`. Add a helper `atomic_at(bodydata, isym, rep)` returning
`atomic[0][0]` for `rep==0` else `atomic[isym][rep-1]`. After this step, both managers share
one code path. **This refactor is behavior-preserving** for c2–c6/p2 and should be the first
commit, gated by the existing histogram tests passing unchanged.

## Feature 1: Polyhedral symmetry (T / O / I)

New `PolyhedralSymmetry : public ISymmetry`
(`include/core/data/symmetry/PolyhedralSymmetry.h`, `source/.../PolyhedralSymmetry.cpp`):

- Constructed from a fixed `std::vector<Matrix<double>>` of the group's rotation matrices
  (12 for T, 24 for O, 60 for I), built once from canonical generators. `repetitions()` =
  `group size - 1` (identity excluded).
- Optimisable parameters identical in count to a cyclic symmetry: `translation` (offset of
  the body from the group centre) + a 3-vector axis-angle **frame orientation**. The 12/24/60
  internal matrices are fixed; only frame placement is free. `span_translation()` /
  `span_rotation()` expose these two 3-vectors.
- `get_transform(cm, rep)` = `v -> F·G_rep·F⁻¹·(v - cm - t) + cm + t`, with `F` the frame
  rotation, `G_rep` the rep-th fixed group matrix, `t` the offset.
- `is_closed()` returns **false** — `repetitions()` already excludes the identity; there is
  no "extra copy that coincides". All reuse is carried by `internal_pair_schedule()`, which
  overrides the base to call `compute_pair_schedule` on the 12/24/60 placements. This avoids
  the cyclic `+1` closed hack. Verify self-scaling `1 + size_symmetry_total()` (= 12/24/60)
  is correct against a brute-force histogram.
- `add()` accumulates `translation` + frame deltas via `dynamic_cast<const PolyhedralSymmetry*>`.
- Registration: add `type::t`, `type::o`, `type::i` to `PredefinedSymmetries`
  (`get(type)`, `get(string_view)` for `"t"`/`"o"`/`"i"`); add cases to
  `OptimizableSymmetryStorage::add` and the `dynamic_cast` ladder in `TransformStrategy.cpp`.

## Feature 2: Reference-symmetry (one symmetry, a group of bodies)

A reference-symmetry behaves like an ordinary cyclic symmetry whose unit is a *group* of
bodies. New `ReferenceSymmetry : public ISymmetry`
(`include/core/data/symmetry/ReferenceSymmetry.h`, `source/.../ReferenceSymmetry.cpp`):

- Holds the same parameters as `CyclicSymmetry` plus a list of participating body indices.
  The rotation reference is the **combined centre-of-mass** of the participating bodies
  (independent of body count — only the combined CM is needed to fix relative transforms).
- The shared object is owned by one designated *primary* body's `SymmetryStorage`. Other
  participating bodies hold a non-owning `ReferenceSymmetryView : public ISymmetry` that
  delegates every call (`get_transform`, `repetitions`, schedules, `span_*`) to the shared
  object via `observer_ptr`. This keeps `body.symmetry().get(i)` returning an `ISymmetry*`
  so all existing iteration in the managers works unchanged.
- `get_transform` ignores the per-body `cm` argument and uses the stored combined CM, so
  every participating body's atoms get transformed by the **same** affine map — the group
  moves rigidly. `generate_transformed_data` needs **no change**: it already forwards
  whatever `get_transform` produces, per body.
- `internal_pair_schedule()` reuses the cyclic base logic (it *is* a cyclic group of spots).
  `cross_pair_schedule()` is overridden: copies of two different participating bodies at the
  same repetition index are rigidly linked, and `dist(A_j, B_k)` depends only on `j-k`, so
  the cross schedule is the cyclic-difference schedule rather than the default full
  enumeration. The managers consume this where they currently double-loop sym-vs-sym.
- Sequencer: extend `SymmetryElement` (`SymmetryElement.cpp` constructor + `_parse`) with a
  multi-body form, e.g. `symmetry { b1 b2 ref-c3 }` — create one `ReferenceSymmetry` on the
  primary body, `ReferenceSymmetryView`s on the others, register names, seed the offset.

**Risk — clone/sharing.** The refinement pipeline clones whole bodies / storages
(`initial_conformation`, `absolute_parameters`). A naive `clone()` duplicates the shared
object per body and breaks sharing. Clone must be **two-phase at Molecule level**: clone all
bodies, then re-link each `ReferenceSymmetryView` to the cloned primary's `ReferenceSymmetry`.
Add a post-clone fix-up in the Molecule copy path / `MoleculeSymmetryFacade`.

## Feature 3: Nested symmetry (`p2-3`)

New `CompositeSymmetry : public ISymmetry`
(`include/core/data/symmetry/CompositeSymmetry.h`, `source/.../CompositeSymmetry.cpp`)
holding an *inner* and an *outer* `unique_ptr<ISymmetry>`. The inner symmetry's
`{original + copies}` set is treated as one unit that the outer symmetry replicates.

- `repetitions()` = `(1+inner.repetitions())·(1+outer.repetitions()) - 1`
  (p2-3 → `2·3 - 1 = 5`, i.e. 6 placements).
- `get_transform(cm, rep)`: decode `rep` into `(outer_k, inner_j)` via
  `placement = outer_k·(1+inner.reps) + inner_j`, then return
  `T_outer(outer_k) ∘ T_inner(inner_j)` (trivial `std::function` composition).
- `generate_transformed_data` and `explicit_structure` need **no change** — their flat loop
  over `repetitions()` calling `get_transform` produces all 6 placements (the payoff of the
  get_transform chokepoint design).
- `internal_pair_schedule()` overrides the base by feeding the 6 placement transforms to
  `compute_pair_schedule` — the product-group reuse (outer-cyclic × inner) emerges
  automatically from the bucketer.
- Parameters: a composite has two parameter sets and `std::span` cannot represent a
  non-contiguous pair. Handle it the way `TransformStrategy` already special-cases by type:
  add a `CompositeSymmetry` branch to the `dynamic_cast` ladder in `TransformStrategy.cpp`
  (`apply_symmetry`/`add_symmetries`) and to `ParameterGenerationStrategies.cpp:39-44` that
  **recurses** into inner+outer. Localizes the special-casing to the two refinement files
  that already special-case by type.
- Registration: `PredefinedSymmetries::get(string_view)` splits names on `-`, recursively
  builds each part, wraps in `CompositeSymmetry`. Add composite cases to
  `OptimizableSymmetryStorage::add`.

## Files

### Create
- `include/core/data/symmetry/PolyhedralSymmetry.h`, `source/core/data/symmetry/PolyhedralSymmetry.cpp`
- `include/core/data/symmetry/CompositeSymmetry.h`, `source/core/data/symmetry/CompositeSymmetry.cpp`
- `include/core/data/symmetry/ReferenceSymmetry.h`, `source/core/data/symmetry/ReferenceSymmetry.cpp` (incl. `ReferenceSymmetryView`)
- `include/core/data/symmetry/PairSchedule.h`, `source/core/data/symmetry/PairScheduleUtil.cpp`
- `source/core/data/symmetry/ISymmetry.cpp` (base `internal_pair_schedule` / `cross_pair_schedule`)
- Test files (see below)

### Modify
- `include/core/data/symmetry/ISymmetry.h` — add `CopyPair`, the two schedule methods
- `include/core/data/symmetry/PredefinedSymmetries.h`, `source/core/data/symmetry/PredefinedSymmetries.cpp` — `type::t/o/i`, `-`-split composite parsing, reference handling
- `source/core/hist/histogram_manager/SymmetryManagerMT.cpp` — consume schedules; `atomic_at` helper
- `source/core/hist/histogram_manager/PartialSymmetryManagerMT.cpp` — consume schedules; verify incremental `symmetry_modified` tracking recomputes every schedule entry touching a modified symmetry
- `include/rigidbody/rigidbody/parameters/OptimizableSymmetryStorage.h` — `add()` switch: new types
- `source/rigidbody/transform/TransformStrategy.cpp` — `dynamic_cast` ladder + composite recursion
- `source/rigidbody/parameters/ParameterGenerationStrategies.cpp` — composite recursion, reference single-delta
- `source/rigidbody/sequencer/elements/setup/SymmetryElement.cpp` — multi-body `symmetry { ... }` form
- Molecule clone path / `MoleculeSymmetryFacade` — post-clone re-link of `ReferenceSymmetryView`s
- `source/core/data/symmetry/CyclicSymmetry.cpp:91` — fix the `span_rotation()` bug
  (`.begin(), .begin()` → `.begin(), .end()`; currently returns an empty span) while in the area

## Riskiest parts

1. **Reference-symmetry clone/sharing** across refinement conformations — needs a
   Molecule-level post-clone fix-up; getting it wrong silently desyncs the optimiser.
2. **Manager consistency** — `PartialSymmetryManagerMT` (refinement) vs `SymmetryManagerMT`
   (reference) must produce identical histograms; the incremental modification tracking in
   the partial manager must recompute every schedule entry involving a changed symmetry.
3. **Bucketer numerical robustness** — canonicalizing relative transforms by rounding can
   misbucket near-degenerate cases; round to a tested tolerance and unit-test bucket counts
   against known group orders.
4. **Composite parameter handling** — non-contiguous params don't fit `std::span`; the
   recurse-in-`TransformStrategy` approach avoids it but means composites are invisible to
   generic span-based code (audited: only the two refinement files use `span_*`).

## Verification

Layered, each layer gating the next:

1. **Refactor-only commit first.** Relocate cyclic scale logic into
   `ISymmetry::internal_pair_schedule()` and rewire both managers. Run
   `tests/feature/hist/symmetry_manager.cpp` and `tests/feature/hist/partial_symmetry_manager.cpp`
   — they must pass **unchanged**, proving the relocation is behavior-neutral.
2. **Unit tests** (extend `tests/unit/data/cyclic_symmetry.cpp`, `point_symmetry.cpp`; new
   `polyhedral_symmetry.cpp`, `composite_symmetry.cpp`, `reference_symmetry.cpp`): assert
   `get_transform` yields the expected placements; assert schedule multiplicities sum to
   `C(reps+1, 2)` (every unordered pair covered exactly once); assert the bucketer agrees
   with brute-force pairwise grouping.
3. **Golden-histogram cross-check** (decisive). For each new symmetry, compute the histogram
   two ways: (a) `SymmetryManagerMT` with reuse, (b) `explicit_structure()` materialized and
   fed to the plain non-symmetry histogram manager (naive O(n²)). Must be bin-for-bin
   identical — this is ground truth for all scale arithmetic. Add to
   `tests/feature/hist/symmetry_manager.cpp`.
4. **Manager-equivalence test** — `PartialSymmetryManagerMT` vs `SymmetryManagerMT` on the
   same molecule for tetrahedral, `p2-3`, and a 2-body reference-c3. Extend
   `tests/feature/hist/partial_symmetry_manager.cpp`.
5. **Incremental-refinement test** — perturb a symmetry parameter, recompute via the partial
   manager, compare against a from-scratch `SymmetryManagerMT` (catches stale modification
   bookkeeping).
6. **Parser tests** (`tests/unit/rigidbody/elements/sequence_parser_symmetry.cpp`):
   `symmetry t`, `symmetry b1 p2-3`, `symmetry { b1 b2 ref-c3 }` produce the expected storage
   and name registrations.
7. **Reference-symmetry clone test** — clone a Molecule containing a reference-symmetry,
   mutate the shared parameter on the clone, assert the original is untouched and both
   participating bodies in the clone see the mutation.

End-to-end: build with the existing CMake setup and run the rigid-body refinement on a
symmetry config (cf. `tests/files/rigidbody/symmetry.conf`) extended with `t` / `p2-3` /
multi-body cases, confirming the optimiser converges and χ² is finite.

