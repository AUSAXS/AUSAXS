#!/usr/bin/env python3
"""
Times one pure Debye calculation of DebyeCalculator or AES-Debye inside its own interpreter, for benchmark_paper.

This is not a benchmark driver: benchmark_paper.cpp decides what to run and calls this once per measurement. It exists
because both programs are Python libraries. Timing them as processes would add the interpreter start and the import of
torch, so they are timed in-process like AUSAXS is. The timed region is from coordinates in memory to I(q) (forward) or
to dL/dr (gradient), as for AUSAXS; reading the file is excluded.

Conventions, set to match AUSAXS's vacuum model: no Debye-Waller factor (DebyeCalculator's biso defaults to 0.3), and
the shared q grid q = q_step*i, i = 1..n_q. DebyeCalculator returns half of AUSAXS's sum over all ordered pairs;
plot_paper.py takes care of the normalisation, so the intensities are written as returned.

One record is appended to --results in the same format as benchmark_paper's own. Exit code 0 on success, 3 if the
warm-up run alone exceeded --budget seconds (the caller then stops the series), 1 on failure.
"""
import argparse
import json
import os
import sys
import time
import traceback
from importlib import metadata

import numpy as np


def read_xyz(path):
    """Plain XYZ: atom count, comment, then element x y z."""
    with open(path) as f:
        n = int(f.readline())
        f.readline()
        elements, coords = [], []
        for _ in range(n):
            parts = f.readline().split()
            elements.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
    return elements, np.array(coords, dtype=np.float64)


def version(dist):
    try:
        return metadata.version(dist)
    except metadata.PackageNotFoundError:
        return "unknown"


def setup_debyecalculator(args, elements, coords, q):
    import torch
    from debyecalculator import DebyeCalculator

    torch.set_num_threads(args.threads)
    dtype = torch.float64 if args.dtype == "float64" else torch.float32
    kwargs = dict(
        qmin=q[0], qmax=q[-1] + args.q_step/2, qstep=args.q_step, biso=0.0, device=args.device, dtype=dtype,
        pair_sum=args.pair_sum, num_threads=args.threads,
    )
    if args.spread_step_factor is not None:
        kwargs["_spread_step_factor"] = args.spread_step_factor
    calc = DebyeCalculator(**kwargs)
    if calc.device != args.device:
        raise RuntimeError(f"DebyeCalculator fell back from {args.device} to {calc.device}")
    q_dc = calc.q.squeeze(-1).detach().cpu().numpy().astype(np.float64)
    if q_dc.size != q.size or np.max(np.abs(q_dc - q)) > 1e-5:
        raise RuntimeError(f"DebyeCalculator's q grid differs: {q_dc.size} points, {q_dc[:3]}...{q_dc[-3:]}")

    pos = torch.tensor(coords, dtype=dtype, device=args.device)
    sync = torch.cuda.synchronize if args.device.startswith("cuda") else (lambda: None)
    info = {
        "version": version("debyecalculator"), "torch": torch.__version__, "dtype": args.dtype, "pair_sum": args.pair_sum,
        "device_name": torch.cuda.get_device_name() if args.device.startswith("cuda") else "cpu",
    }

    if args.mode == "forward":
        def run():
            return calc.iq((elements, pos)).i  # numpy on the host, so the device is synchronised
    else:
        v = torch.ones(q.size, dtype=dtype, device=args.device)  # dL/dI, as for AUSAXS: the gradient of sum_q I(q)

        def run():
            p = pos.clone().requires_grad_(True)
            I = calc.iq((elements, p), keep_on_device=True).i
            (I*v).sum().backward()
            sync()
            return None
    return run, info


def setup_aesdebye(args, elements, coords, q):
    import aesdebye

    if args.mode != "forward":
        raise RuntimeError("AES-Debye has no gradient")
    calc = aesdebye.DebyeCalculator(
        nThreads=args.threads, useGPU=args.device != "cpu", binsResolution=args.bins_resolution, verbose=False,
    )
    pos = aesdebye.Positions(chemicalSymbols=elements, coordinates=coords.tolist())
    # AES-Debye's grid is start + i*(end - start)/steps for i < steps, so this reproduces q exactly
    start, end = q[0], q[0] + args.q_step*q.size

    def run():
        result = calc.calculateProfile(pos, start=start, end=end, steps=q.size)
        _, profile = result["total"] if "total" in result else next(iter(result.values()))
        q_aes = np.asarray(profile.q)
        if q_aes.size != q.size or np.max(np.abs(q_aes - q)) > 1e-9:
            raise RuntimeError(f"AES-Debye's q grid differs: {q_aes.size} points, {q_aes[:3]}...{q_aes[-3:]}")
        return np.asarray(profile.intensity)

    info = {
        "version": version("aesdebye"), "bins_resolution": args.bins_resolution,
        "device_name": "gpu" if args.device != "cpu" else "cpu",
    }
    return run, info


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--suite", required=True)
    p.add_argument("--program", required=True, choices=["debyecalculator", "aes-debye"])
    p.add_argument("--variant", required=True)
    p.add_argument("--structure", required=True)
    p.add_argument("--atoms", type=int, required=True)
    p.add_argument("--xyz", required=True)
    p.add_argument("--results", required=True, help="results.jsonl to append the record to")
    p.add_argument("--iq", help="write q and I(q) of the last run here")
    p.add_argument("--reps", type=int, default=5, help="timed runs after the warm-up; 0 keeps only the single first run")
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--budget", type=float, default=300)
    p.add_argument("--n-q", type=int, default=1000)
    p.add_argument("--q-step", type=float, default=0.01)
    p.add_argument("--device", default="cpu", help="cpu or cuda")
    p.add_argument("--mode", default="forward", choices=["forward", "gradient"])
    p.add_argument("--dtype", default="float32", choices=["float32", "float64"], help="DebyeCalculator only")
    p.add_argument("--pair-sum", default="grid", choices=["grid", "direct"], help="DebyeCalculator only")
    p.add_argument("--spread-step-factor", type=float, help="DebyeCalculator only: grid spacing 2^floor(log2(f/q_max))")
    p.add_argument("--bins-resolution", type=float, default=1.0, help="AES-Debye only")
    args = p.parse_args()

    record = {
        "key": f"{args.suite}/{args.program}/{args.variant}/{args.structure}", "suite": args.suite, "program": args.program,
        "variant": args.variant, "structure": args.structure, "atoms": args.atoms, "status": "ok", "samples": [],
        "threads": args.threads, "device": args.device, "mode": args.mode,
    }
    code = 0
    try:
        elements, coords = read_xyz(args.xyz)
        q = args.q_step*np.arange(1, args.n_q + 1)
        setup = setup_debyecalculator if args.program == "debyecalculator" else setup_aesdebye
        run, info = setup(args, elements, coords, q)
        record.update(info)

        t0 = time.perf_counter()
        I = run()
        first = time.perf_counter() - t0
        if args.reps == 0:
            record["samples"] = [first]
            record["single_run"] = True
        elif first > args.budget:
            record["samples"] = [first]
            record["status"] = "over budget"
            code = 3
        else:
            for _ in range(args.reps):
                t0 = time.perf_counter()
                I = run()
                record["samples"].append(time.perf_counter() - t0)
        if args.iq and I is not None:
            os.makedirs(os.path.dirname(args.iq), exist_ok=True)
            np.savetxt(args.iq, np.column_stack([q, np.asarray(I, dtype=np.float64)]), fmt="%.17g")
    except Exception as e:
        traceback.print_exc()
        record["status"] = f"failed: {type(e).__name__}: {e}"
        code = 1

    with open(args.results, "a") as f:
        f.write(json.dumps(record) + "\n")
    print(json.dumps(record))
    sys.exit(code)


if __name__ == "__main__":
    main()
