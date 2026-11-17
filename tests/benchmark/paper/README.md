# Paper benchmarks

`benchmark_paper` produces the benchmark figures of the AUSAXS paper:

| Tag | What it measures |
|---|---|
| `[sas]` | Time of a complete fit to SAS data against the number of atoms: the AUSAXS CLI on CPU and GPU, AUSAXS in-process through the library, Pepsi-SAXS, FoXS and CRYSOL. 18 structures, from 602 to 43k atoms (`manifest.tsv`). |
| `[debye]` | Time of the pure Debye sum of fcc silver spheres, from 225 atoms to about 830k: the AUSAXS forward sum on CPU and GPU at bin widths 0.25 and 0.05 Å and its gradient, DebyeCalculator forward and gradient on CPU and GPU, and AES-Debye forward on CPU and GPU. |
| `[accuracy]` | Deviation of every program from the exact Debye sum (DebyeCalculator's direct sum in double precision), and accuracy against cost at fixed size. |

Each measurement is one warm-up run and `AUSAXS_BENCH_REPS` (default 5) timed runs. All runs are written to `results.jsonl`. A restarted run skips what is already there, so an interrupted run can be resumed by restarting it. Once a single run of a program takes longer than `AUSAXS_BENCH_BUDGET` (default 300 s), that program is not run on larger structures.

## Setup on a new machine

Linux only (`posix_spawn`).

1. **AUSAXS.** Build in Release and install the GPU backend:

   ```sh
   cmake -B build-release -S . -DCMAKE_BUILD_TYPE=Release
   cmake --build build-release --target benchmark_paper -j8
   ./build-release/bin/ausaxs gpu install          # pick the GPU; the benchmark checks that it actually loads
   ```

2. **SAS inputs.** The manifest paths are relative to `AUSAXS_BENCH_DATA` (default `data/`, which is not tracked). On a machine without `data/`, unpack the bundle `paper_inputs.tar.gz` (36 files, 4.8 MB) into an empty folder and point `AUSAXS_BENCH_DATA` at it.

3. **SAS programs.** Each is optional; a series whose variable is unset is left out.
   - Pepsi-SAXS: the Linux binary from the Grudinin group's website → `AUSAXS_BENCH_PEPSI`.
   - CRYSOL: from ATSAS (needs a licence) → `AUSAXS_BENCH_CRYSOL`.
   - FoXS: `conda install -c conda-forge imp` provides `foxs`, or build IMP → `AUSAXS_BENCH_FOXS`.
   - If an older binary no longer finds its shared libraries (`ldd <binary> | grep "not found"`), put them in a separate conda env and point `AUSAXS_BENCH_<NAME>_LIBS` at its `lib/`. The harness prepends that folder to `LD_LIBRARY_PATH` for that program only. For example, ATSAS 3/4 needs `libxml2.so.2` (`conda create -n compat-crysol -c conda-forge libxml2=2.13`), and an IMP build against boost 1.83 needs `libboost=1.83`.

4. **Python competitors.** Install them in one env with a CUDA build of torch, and point `AUSAXS_BENCH_PYTHON` at its interpreter:

   ```sh
   conda create -n ausaxs-bench -c conda-forge python=3.12 numpy matplotlib ase cmake ninja pybind11 compilers
   conda activate ausaxs-bench
   pip install torch                                    # the CUDA wheel matching the driver
   pip install debyecalculator==1.1.1
   pip install git+https://gitlab.cs.fau.de/iq23adyz/debye -v   # AES-Debye; builds with CUDA when nvcc is found
   # if nvcc is not found: add -Ccmake.define.CMAKE_CUDA_COMPILER=$(which nvcc)
   python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name())"
   python -c "import aesdebye; aesdebye.DebyeCalculator(useGPU=True, verbose=False)"   # raises if built without CUDA
   ```

## Running

From the repository root:

```sh
export AUSAXS_BENCH_DATA=/path/to/inputs          # or leave unset to use data/
export AUSAXS_BENCH_PEPSI=... AUSAXS_BENCH_FOXS=... AUSAXS_BENCH_CRYSOL=...
export AUSAXS_BENCH_PYTHON=$CONDA_PREFIX/bin/python
B=./build-release/tests/benchmark/bin/benchmark_paper

$B "[sas]"          # roughly 30-60 min; FoXS on the largest structures dominates
$B "[debye]"        # roughly 1-2 h; DebyeCalculator's CPU gradient dominates until it reaches the budget
$B "[accuracy]"     # roughly 10 min; the exact reference is O(N^2 Q)
python3 tests/benchmark/plot_paper.py output/paper_benchmark
```

Before the real run, check that everything starts with a quick run into a separate folder:

```sh
AUSAXS_BENCH_OUT=output/paper_smoke AUSAXS_BENCH_REPS=1 AUSAXS_BENCH_RADII=10,15 AUSAXS_BENCH_ACCURACY_RADII=10 \
    $B "[paper]"
grep -v '"status":"ok"' output/paper_smoke/results.jsonl     # should print nothing
```

Keep the machine otherwise idle while the real run is in progress. GPU timings in particular change with what else uses the GPU.

## Variables

| Variable | Default | Meaning |
|---|---|---|
| `AUSAXS_BENCH_OUT` | `output/paper_benchmark` | results, figures, intermediate files |
| `AUSAXS_BENCH_DATA` | `data` | root of the manifest's paths |
| `AUSAXS_BENCH_MANIFEST` | `tests/benchmark/paper/manifest.tsv` | the SAS structures |
| `AUSAXS_BENCH_AUSAXS` | the `ausaxs` of this build | CLI for the end-to-end series |
| `AUSAXS_BENCH_REPS` | 5 | timed runs after the warm-up |
| `AUSAXS_BENCH_THREADS` | all hardware threads | threads for AUSAXS in-process, DebyeCalculator and AES-Debye. The CLIs run with their defaults. |
| `AUSAXS_BENCH_BUDGET` | 300 | seconds; a single run past this ends that program's series |
| `AUSAXS_BENCH_GPU` | 1 | 0 leaves out every GPU variant |
| `AUSAXS_BENCH_DEVICE` | `cuda` | torch device of the competitors' GPU runs |
| `AUSAXS_BENCH_COMPETITORS` | `dc,aes` | which Python libraries to run |
| `AUSAXS_BENCH_PYTHON` | `python3` | interpreter with the competitors installed |
| `AUSAXS_BENCH_RADII` | 10,15,20,30,40,50,58,75,93,120,150 | Å; particle sizes of `[debye]`. 58 and 93 give about 48k and 198k atoms, the sizes of the paper's table. |
| `AUSAXS_BENCH_ACCURACY_RADII` | 10,20,30 | Å; particle sizes of `[accuracy]` |
| `AUSAXS_BENCH_BIN_WIDTHS` | 0.5,0.25,0.1,0.05,0.02,0.01 | Å; AUSAXS bin widths of `[accuracy]` |
| `AUSAXS_BENCH_<NAME>`, `AUSAXS_BENCH_<NAME>_LIBS` | unset | executable and extra library folder of `PEPSI`, `FOXS`, `CRYSOL` |

## Conventions

- **Structures.** Every SAS program gets the same structure, with only the heavy ATOM records and no waters, in its own working directory. Data files with q in nm⁻¹ (marked in the manifest) are converted to Å⁻¹ first, so every program reads the same q values. The q range of each data set is part of the benchmark.
- **Program options.** AUSAXS: `fit --no-gpu` or `--gpu`, plus `--offline --allow-unknown-atoms --allow-unknown-residues`. CRYSOL: `--constant --implicit-hydrogen 1`, since the structures carry no hydrogens. Pepsi-SAXS and FoXS: defaults. Without an explicit `--no-gpu`, AUSAXS uses an installed GPU, which costs about 0.1 s of start-up per run.
- **Pure Debye.** Every program is timed from coordinates in memory to I(q), or to dL/dr for the gradient, on the q grid q = 0.01 i, i = 1…1000 (up to 10 Å⁻¹). DebyeCalculator runs with `biso=0`, since AUSAXS applies no Debye-Waller factor, and otherwise with its defaults (float32, grid pair sum). The gradient is the vector-Jacobian product for dL/dI = 1. For DebyeCalculator it includes the forward pass that autograd needs; AUSAXS's does not need one. AES-Debye has no gradient.
- **Accuracy.** The reference is DebyeCalculator's `pair_sum='direct'` in float64, which is independent of AUSAXS. The plot script removes each program's normalisation convention: AUSAXS sums over all ordered pairs, twice DebyeCalculator's sum. It checks that the convention is exactly 1/2, 1 or 2 rather than fitting a scale.
- **Statistics.** `summary.csv` has the min, median, mean and standard deviation of every measurement. The figures show the mean, with one standard deviation as the error bar.
