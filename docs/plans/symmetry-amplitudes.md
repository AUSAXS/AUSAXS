# Symmetry amplitudes in the `parameter` element

**Status:** planned, not started.
**Scope:** give symmetry optimisation the same amplitude control that body translations and rotations
already have, and make its translation and rotation components independently selectable.

---

## Background: what prompted this

The `parameter` element accepts `iterations`, `translate`, `rotate`, `mode`/`strategy` and `decay`.
Translations and rotations each take an amplitude; symmetries take nothing. Symmetry optimisation can
only be switched on wholesale via `mode`, and its translation and rotation components cannot be
separated.

A `symmetry` argument was documented on the wiki but never existed in the parser. Because the parser
did not validate argument keys, it was silently discarded. That hole is being closed separately (see
the argument-whitelist work); this document covers only the amplitudes.

## Findings from the investigation

### 1. Symmetry amplitudes are hardcoded

`create_distributions` in `source/rigidbody/parameters/ParameterGenerationStrategy.cpp` builds:

```cpp
translation_symmetry_dist(-2*Rg, 2*Rg);
rotation_symmetry_dist(-3, 3);
angle_symmetry_dist(-pi/2, pi/2);
```

`angle_symmetry_dist` is dead — constructed, tied and copied, never read. `set_max_translation_distance()`
and `set_max_rotation_angle()` only touch the *body* distributions, so even the C++ API cannot tune the
symmetry amplitudes.

### 2. The translation/rotation split already exists — twice

`ISymmetry::span_translation()` and `span_rotation()` are separate spans, and
`LimitedParameterGenerator::next()` already loops them independently, gated on
`OptimizableSymmetryStorage::optimize_translate` / `optimize_rot_axis`. Those two flags are set to
`true` at every site that adds a symmetry (`SymmetryElement`, `SplitElement`,
`ConvertToSymmetryElement`, `OptimizableSymmetryStorage::add`) and to `false` nowhere outside a test.
They are dead switches.

Independently, `selection::ParameterMask` carries `sym_translation` / `sym_axis` with ready-made
`symmetry_only_trans()` / `symmetry_only_axis()` factories, but the `select` element's `parameters`
keyword never exposes them.

So there are three overlapping control points and none is reachable from a script.

### 3. "Symmetry rotation" means different things per symmetry type

This is the part that makes a single amplitude number meaningless today.

| type | `span_rotation()` is | units |
|---|---|---|
| `CyclicSymmetry` (c2–c12) | `_repeat_relation.axis`, an **unnormalised direction vector**; the angle is fixed by the type | dimensionless |
| `PointSymmetry` (p2) | `rotation`, **Euler angles** | radians |
| `IPolyhedralSymmetry` (T/O/I), dihedral | group-frame orientation, **Euler angles** | radians |

`rotation_symmetry_dist(-3, 3)` is therefore ±172° for p2 and polyhedral — essentially a uniform
re-orientation every step — but "±3 axis units" for cyclic.

Worse, for cyclic symmetries `get_transform()` normalises the axis
(`source/core/data/symmetry/CyclicSymmetry.cpp`) while `add()` accumulates raw components. The axis
magnitude therefore random-walks upward over a run, and the *effective* angular step shrinks on its
own, independently of the decay strategy.

`span_translation()` is consistently Å everywhere, except `PlanarDihedralSymmetry`, which exposes only
2 of the 3 components.

---

## Proposed script surface

Mirror the existing translate/rotate pattern:

```
parameter {
    iterations 100
    translate 5
    rotate 0.1
    symmetry_translate 10     # Å
    symmetry_rotate 0.2       # rad
    decay exponential
}
```

- Keep the existing "deduce the mode from which amplitudes are non-zero" logic, extended to four
  amplitudes; `mode` stays an explicit override.
- Add `symmetry_translate_only` and `symmetry_rotate_only` mode names.
- Keep `symmetry_only` meaning both symmetry components, so existing scripts are unaffected.
- Defaults when unspecified stay at today's values (`2*Rg`, and the current rotation behaviour).
- Accept a plain `symmetry <float>` as an alias enabling both symmetry components, since that spelling
  was documented on the wiki and appears in existing user scripts.

Note the deduction quirk that must be preserved or deliberately changed: when both `translate` and
`rotate` are non-zero, none of the current deduction branches fires and the mode falls through to its
`"both"` default, which enables symmetry too.

## Implementation steps

1. **Collapse the template into a runtime flag set.**
   Replace `LimitedParameterGenerator<TRANSLATE, ROTATE, SYMMETRY>` with a single class holding
   `{bool translate, rotate, sym_translate, sym_rotate;}` plus four amplitudes. Four template bools
   would need 16 explicit instantiations, and `next()` runs once per optimisation step — the
   `if constexpr` buys nothing. The call surface is small: the factory plus five test sites
   (`tests/feature/rigidbody/parameter_generation.cpp`, `symmetry_backup.cpp`,
   `tests/unit/rigidbody/elements/sequence_parser_symmetry.cpp`). Keep `AllParameters`,
   `TranslationsOnly`, `RotationsOnly`, `SymmetryOnly` as factory functions so those tests keep
   compiling.

2. **Add symmetry amplitude setters** to `ParameterGenerationStrategy`, alongside the existing
   `set_max_translation_distance` / `set_max_rotation_angle`, and surface them on `ParameterElement`.

3. **Extend `ParameterElement::_parse`** with the new arguments, the extended deduction, and the new
   mode names. Add every new key to `_valid_arguments()` — the parser now rejects anything not
   declared there.

4. **Delete `OptimizableSymmetryStorage::optimize_translate` / `optimize_rot_axis`.** They duplicate
   what `ParameterMask` already does, and this change would give them a third rival. The clean split
   is: the parameter generator owns *what can move and how far* (global), `ParameterMask` owns *what
   moves this step* (per-selection). Removing them touches the four sites listed in finding 2 plus
   `tests/feature/rigidbody/symmetry_backup.cpp`.

5. **Delete or wire up `angle_symmetry_dist`.**

6. **Optionally expose the mask split** in the `select` element's `parameters` keyword
   (`symmetry_translation`, `symmetry_axis`), since `ParameterMask` already implements it.

## Deliberately kept separate

Making `symmetry_rotate` mean an actual **angle for every symmetry type** requires perturbing the
cyclic axis as a rotation of the normalised axis rather than an additive component delta. Without it,
the same number means ±3 rad for p2 and something drifting for c4, so the knob is not portable across
a script's symmetries.

This is a behaviour change to existing runs and should not be smuggled into an API extension — it
deserves its own change, its own tests and its own note in the changelog. Recommended order: land the
API extension first with today's semantics preserved, then normalise the units.

## Testing

- Parser tests: each new argument accepted, rejected when malformed, mode deduction for all
  combinations of the four amplitudes, and the `symmetry` alias.
- Generator tests: with only `symmetry_translate` set, `span_rotation()` deltas must be zero across
  many draws, and vice versa. Extend the existing pattern in
  `tests/feature/rigidbody/parameter_generation.cpp`.
- Regression: a script with no symmetry arguments must produce byte-identical parameter sequences to
  the current implementation for a fixed seed.
