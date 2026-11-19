# Docs site TODO

Review of the `docs-site` branch, checked against `origin/master` at 0accd315 (2026-11-19).
`mkdocs build` is clean: no link or anchor warnings.

## Wrong or broken now
- [x] **GPU troubleshooting is missing.** The old wiki's Troubleshooting page (`~/projects/wiki/Troubleshooting.md`) has a full GPU section: the error messages, NVIDIA/AMD/Intel driver requirements, clusters and containers, "chosen GPU not present", "incomplete installation". Port it to `docs/src/getting-started/troubleshooting.md`. Both `getting-started/gpu.md` and `GPULoader.h::troubleshooting_url` send users there for GPU problems.
- [x] **Broken README link.** `README.md` links to `documentation/cli/saxs-fitter/`, but the page is built at `documentation/cli/fit/`. The README also lists SANS as "(in development)" with no link, though there is now a SANS page.
- [x] **Quote `pip install "pyausaxs[all]"`.** Unquoted, it fails in zsh (the macOS default). Fix it in `documentation/cli/fit.md`, `em-validation.md` and `rigidbody.md`. The rigidbody.md code block is also tagged `python` instead of `bash`.
- [x] **Unsupported claim in `applications/rigid-body-refinement.md`.** It says the refined dimer fits "better than the deposited dimer model", but the figure only shows initial vs. refined. Either add the deposited model's curve or reword. Also note that χ² varies from run to run.
- [x] **`--gpu`/`--no-gpu` coverage.** `documentation/cli/gpu-driver.md` says they exist for `fit`, `em` and `rigidbody`; `sans` and `debye` have them too. The base-group tables in `fit.md` and `em-validation.md` list only `--gpu`, without `--no-gpu` or the "on once a GPU is chosen" default.
- [x] **rigidbody.md `print` argument names.** The reference documents `message`/`color`; the examples use the aliases `msg`/`colour`. Use one form or mention the aliases. Replace the SHA-pinned link to `ConsoleColor.h` with an inline list of colours.
- [x] **Speed numbers disagree.** `why-ausaxs.md` says about 300k atoms in 10 s on a CPU; the README says "hundreds of thousands of atoms in a few seconds". The benchmark figure on `why-ausaxs.md` is still a commented-out placeholder; the README already embeds `benchmark.png`.

## Missing pages
- [x] ~~**Python API reference.**~~ Deferred: it will come from the pyAUSAXS repo once that is set up. The application pages use `create_molecule`, `debye_raw`, `settings.set`, `prepare_rigidbody_refinement`, … but the site has no reference for them; the "Python API" tile leaves for GitHub. The README advertises point scatterers and the gradient, which appear nowhere on the site.
- [x] ~~**Graphical interface page.**~~ Rejected: the GUI should not need documenting. `ausaxs setup` is now in the command table. `ausaxs gui` and the landing-page tile both go to a three-line section in `fit.md`. `ausaxs setup` is only documented under `compiling.md` and is missing from the command table in `documentation/index.md`.
- [x] ~~**Input and output formats.**~~ Rejected: not useful to readers.
  - Structure formats: PDB, mmCIF, CIF and XYZ, as the landing page promises.
  - Data files: columns, headers, the q-unit guess.
  - Accepted EM map formats.
  - The columns of `ausaxs.fit` and the fields of `report.txt`.
- [x] ~~**Settings file reference.**~~ Rejected: not useful to readers. Several pages point to `config --save` and keys like `format`, but no page describes the file or its keys. The wiki's Configuration-files page was empty too.
- [x] ~~**Model/method overview.**~~ Moved: the exv models are documented in the 2025 paper; the rest belongs in the SAXS tutorial (BL-338).
  - The default fit parameters: scaling, background, hydration c, plus optional exv/density.
  - What the exv models mean in practice.
  - How to read χ²ᵣ.

  `fit.md` currently just says "see the main AUSAXS article".

## Gaps in existing pages
- [x] **Troubleshooting (non-GPU).** Add the common runtime errors:
  - unknown atoms or residues (`--allow-unknown-*`, `--force-field-names`, `--offline`);
  - a wrong q-unit guess;
  - the Windows libomp crash fixed in #530;
  - `ausaxs` not found (currently only on the installation page).
- [x] ~~**State defaults.**~~ Rejected: the docs do not state default values unless they are critical.
  - `grid --width`: 1.75 Å since #537.
  - `grid --rvol`, `--threads`.
  - The EM level range, `--frequency` and `--charge-levels`. `em-validation.md` warns against going "below the default limits" without ever stating them.
- [x] **EM upper level.** Since #529, a `--levelmax` that leaves the map empty is refused. Mention it in the level-range section of `em-validation.md`.
- [x] **SANS application page.**
  - It is the only example with placeholder files (`model.pdb`, `data.dat`) and no figure.
  - Mention that deuterium in structure files is now parsed (#532).
  - Clarify `--exchange`: the CLI help says it applies to *implicit* hydrogens, so say what happens with `hydrogens --keep`.
  - Warn that `--d2o` defaults to pure D₂O.
  - The Python example sets a raw `neutron_solvent_density` while the CLI uses `--d2o`; use a D₂O helper if pyAUSAXS has one.
- [x] **`vacuum` exv model.** It is still accepted on master but missing from the `fit.md` list. Put it back under "Additional models" unless dropping it was deliberate because `ausaxs debye` replaces it.
- [x] **Citing page.** Review added; no specific guidance for rigid-body or SANS work (the main paper applies). Add the Lytje, Pedersen & Hub (2026) review that the SAXS and SANS pages cite for background, and guidance for rigid-body and SANS work.
- [x] **Developers section.** C API rejected: not documented. Currently only the macOS debug build. Add:
  - a link to CONTRIBUTING.md;
  - how to build and run the tests;
  - the C API, now with `settings::reset` exported (#535);
  - the vcpkg toolchain flag for Windows in `compiling.md`.
- [x] **Tutorials.** pyAUSAXS examples linked; new tutorials are BL-338. There are only three, all external. There is nothing for SANS, crystals/debye, CLI rigidbody scripts or Python. Link the pyAUSAXS `examples/` folder from `tutorials.md`.
- [x] **Tone.** `rigidbody.md` still has the wiki's "_**NEW**_" banner and italic notice. `troubleshooting.md` still uses raw `<details>` instead of the admonitions used elsewhere.

## Housekeeping
- [x] Remove the unreferenced repo-root `decomposition.svg`, about 5k lines. `docs/figures/saxs-decomposition.svg` is the source the figure scripts use.
- [x] Make page `title:` values consistent: `installation`, `compiling` and `citing` lack the " - AUSAXS" suffix the others have.
