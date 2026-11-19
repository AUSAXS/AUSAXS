# Docs site TODO

Pre-release review of the `docs-site` branch at d939e7b7, checked against `origin/master` at 0accd315 (2026-11-19).
`mkdocs build --strict` is clean, every external link resolves (the IUCr DOIs answer 403 to scripts only), and every
GitHub path linked from the site exists on master. The trivial fixes from the review are applied; this is what is left.

## Before release
- [ ] **Rebase onto master.** `docs-site` is 10 commits behind. The pages already describe master: `ausaxs sans` (#536), `ausaxs debye`, the `output/saxs|sans|em/` folders (#540) and the 0.02 Å debye bin width (#539). Master does not touch `docs/`, so the rebase should be clean.
- [ ] **Benchmark link.** `why-ausaxs.md` says "The benchmark itself is in `tests/benchmark/`", but `benchmark_paper` is still only on PR #531. Merge it first, or drop the sentence.
- [ ] **Regenerate figures after the rebase.**
  - `build_debye_figures.py` now uses 0.02 Å, the new `ausaxs debye` default, to match the pages; the crystal and powder figures were made at 0.05 Å.
  - `build_em_rigidbody_figures.py` now reads `output/em/` (#540) instead of `output/em_fitter/`; it would have failed on master.

## Decisions
- [ ] **`hydrogens --keep` and the exv tables.** `fit.md` says `--keep` "will also switch the excluded volume tables to a version containing non-zero hydrogen volumes". Nothing in the code does this: `SettingsValidation.cpp` only sets `implicit_hydrogens = false`, and the default `MinimumFluctutation_implicit_H` set gives H zero volume. Either remove the sentence or make the code switch to the `_explicit_H` set.
- [ ] **The SANS example.** `applications/sans.md` keeps the explicit hydrogens of SASDWB3's model (4590 H, all protium) for a measurement in 100% D₂O. Exchange only applies to implicit hydrogens, so none of its labile hydrogens become D. A caveat is now on the page, but the example itself contradicts it. Drop `--keep`, or use a model with deuterated exchangeable sites, then regenerate `sans-fit.svg`.
- [ ] **Rigid-body refinement with SANS.** `applications/rigid-body-refinement.md` and `index.md` say "SAXS or SANS data", but no page explains how to refine against SANS, and `load` only takes `saxs`. Document it (probably a `settings` element before `load`; unverified) or say SAXS.
- [ ] **GPU on macOS.** The table in `getting-started/gpu.md` says macOS is "work in progress", but `GPUInstaller.cpp` knows `ausaxs-gpu-macos-arm64.zip` and the release workflow builds a macOS package. Make the page match what is actually released.
- [ ] **"GPU driver" or "GPU backend".** The pages say driver, the CLI messages and the GPU troubleshooting headings say backend, and "driver" collides with the vendor's graphics driver (`troubleshooting.md`, GPU section). Pick one, or say once on `gpu.md` that the CLI calls it the backend.
- [ ] **Defaults in body text.** The convention is to omit them, but these remain: `fit.md` lines 71, 75, 94, 95, 98, 105, 149, 150, 155; `em-validation.md` 61, 62, 88, 103; `debye.md` 36–38, 176; `sans.md` `--exchange`. The `--d2o` default and its warning are the intended exception.
- [ ] **"Clang v15+"** in `compiling.md` is not verified by CI: there is no Linux Clang build, only AppleClang on macOS. Soften it or say AppleClang.
- [ ] **tkinter.** `pyausaxs[all]` cannot provide `_tkinter`. Without `python3-tk` (Linux) or `python-tk` (Homebrew), `ausaxs gui` says to `pip install pyausaxs[gui]`, which does not help. Add a troubleshooting entry, and maybe fix the message in pyAUSAXS.
- [ ] **Unused asset.** `assets/logo-light.svg` is not referenced anywhere. Use it or delete it.

## Code issues found while reviewing
- [ ] **macOS deployment target.** Release and CI builds run on `macos-26` with no `CMAKE_OSX_DEPLOYMENT_TARGET`, so binaries target the build host's macOS. The pyAUSAXS wheel is tagged `macosx_11_0_universal2` but bundles this library, so it may fail to load on older macOS. Set e.g. `-DCMAKE_OSX_DEPLOYMENT_TARGET=11.0`. `troubleshooting.md` now only says the binaries are built on a recent macOS.
- [ ] **`ausaxs rigidbody config --save`** writes `settings.txt` and then fails with "A configuration script must be supplied" (exit 1). `fit`, `em` and `debye` exit cleanly. Throw `CLI::Success` when there is no script.
- [ ] **`load` resolves `saxs` from the launch folder only.** Structure paths fall back to the script's folder, the data path does not (`io::ExistingFile(saxs.value)` in LoadElement), so scripts are not portable. Use `lookup_file` for it too, or document it in `rigidbody.md`.

## SEO
The site itself is in good shape: a unique title and description per page, canonical URLs, Open Graph tags and a sitemap of all 25 pages. The two duplicate titles (rigid-body and EM, application page vs. CLI reference) are fixed.
- [ ] **GitHub repo "About".** The homepage field is empty and the description is "Repository for the AUSAXS software package". Set the homepage to https://ausaxs.github.io/AUSAXS/, use a keyword description (SAXS, SANS, Debye, small-angle scattering) and add topics. This is the strongest link to the site.
- [ ] **pyAUSAXS `[project.urls]`.** Only `Homepage` is set. Add `Documentation = "https://ausaxs.github.io/AUSAXS/"`, so PyPI links to the site.
- [ ] **Google Search Console.** Verify the site (meta tag in `overrides/main.html`) and submit `sitemap.xml`. A project page cannot serve its own `robots.txt`, so crawlers will not find the sitemap by themselves.
- [ ] **GitHub wiki.** Still enabled. Replace its pages with a pointer to the site, so the stale copy does not compete with it.
- [ ] **Long descriptions.** 14 of the 25 page descriptions are over ~160 characters and are cut off in search results. Put the keywords first.
- [ ] **Nice to have.** An `og:image` for link previews, and JSON-LD `SoftwareApplication` structured data.

## Small things outside the site
- [ ] The Rep. Prog. Phys. review in `citing.md` has no volume or article number yet (2026 should be vol. 89). Add them once assigned.
- [ ] `README.md` line 56: the EM citation has no authors, unlike `citing.md`.
