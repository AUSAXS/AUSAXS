---
title: GPU acceleration - AUSAXS
description: Install the optional AUSAXS GPU backend with ausaxs gpu install, either directly or from a manually downloaded release, and choose which GPU to use.
---

# GPU acceleration

AUSAXS runs on the CPU out of the box. To use a GPU, install the GPU backend, a separate download that AUSAXS loads when it is present. Prebuilt backends are available for Linux (x86-64), Windows (x86-64) and macOS (Apple silicon).

> [!NOTE]
> If you installed AUSAXS with `pip`, the `ausaxs` command may not be on your path. In that case, replace `ausaxs` with `python -m pyausaxs` in every command on this page, e.g. `python -m pyausaxs gpu install`.

## Installing the backend
Run:
```bash
ausaxs gpu install
```
This downloads the latest backend for your platform from the [AUSAXS-GPU releases](https://github.com/AUSAXS/AUSAXS-GPU/releases), lists the GPUs it finds, and asks you which one to use. Once a GPU is chosen, every AUSAXS tool uses it by default. Run the command again to choose a different GPU.

To skip the question, name the GPU with `--device`, either by its identifier as listed, or just by its API, e.g. `--device cuda` for the first CUDA GPU.

## Installing from a downloaded file
If the machine has no internet access, or you want a specific version, download the backend yourself from the [AUSAXS-GPU releases](https://github.com/AUSAXS/AUSAXS-GPU/releases):

| Platform               | File                         |
|------------------------|------------------------------|
| Linux (x86-64)         | `ausaxs-gpu-linux-x64.zip`   |
| Windows (x86-64)       | `ausaxs-gpu-windows-x64.zip` |
| macOS (Apple silicon)  | `ausaxs-gpu-macos-arm64.zip` |

Then install it from the file:
```bash
ausaxs gpu install --from <path to zip file>
```
`--from` also accepts a URL, and replaces any backend already installed. Add `--offline` to make sure no network requests are made.

## Managing the backend

| Command               | Description |
|-----------------------|-------------|
|`ausaxs gpu status`    | Show the installed backend, the chosen GPU, and whether it is used by default. |
|`ausaxs gpu disable`   | Stop using the GPU by default, without removing the backend. |
|`ausaxs gpu enable`    | Use the GPU by default again. |
|`ausaxs gpu remove`    | Remove the backend. |

For a single run, the `--gpu` and `--no-gpu` flags of `ausaxs fit`, `ausaxs em` and `ausaxs rigidbody` override the default.

If the backend does not find your GPU, see the [troubleshooting](troubleshooting.md) page.

## Next steps
You are ready to calculate your first profile. The [tutorials](tutorials.md) walk through complete analyses of real data, and the [SAXS](saxs.md) page explains what goes into the calculation.
