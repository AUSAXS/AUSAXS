---
title: GPU acceleration - AUSAXS
description: Install the optional AUSAXS GPU driver with ausaxs gpu install, either directly or from a manually downloaded release, and choose which GPU to use.
---

# GPU acceleration

AUSAXS runs on the CPU by default to keep the package size small. To use a GPU, you need to download and install a separate ~300MB AUSAXS GPU driver, which links AUSAXS to your graphics card. The driver supports NVIDIA and AMD GPUs on Linux (x86-64) and Windows (x86-64).

??? question "Is your platform or GPU not supported?"
    macOS, ARM systems and Intel GPUs are not supported yet, but AUSAXS still runs on the CPU on these systems. Adding support is currently a low priority. If you would like GPU acceleration on one of them, please [open an issue](https://github.com/AUSAXS/AUSAXS-GPU/issues/new) with your operating system, CPU architecture and GPU model, or 👍 an existing one. Requests from actual users are what moves a platform up the list!

## Installing the driver
Run:
```bash
ausaxs gpu install
```
This downloads the latest driver for your platform from the [AUSAXS-GPU releases](https://github.com/AUSAXS/AUSAXS-GPU/releases), lists the GPUs it finds, and asks you which one to use. Once a GPU is chosen, every AUSAXS tool uses it by default. Run the command again to choose a different GPU (without redownloading).

To skip the question, name the GPU with `--device`, either by its identifier as listed, or just by its API, e.g. `--device cuda` for the first CUDA GPU.

## Installing from a downloaded file
If the machine has no internet access, or you want a specific version, download the driver yourself from the [AUSAXS-GPU releases](https://github.com/AUSAXS/AUSAXS-GPU/releases):

| Platform               | File                         |
|------------------------|------------------------------|
| Linux (x86-64)         | `ausaxs-gpu-linux-x64.zip`   |
| Windows (x86-64)       | `ausaxs-gpu-windows-x64.zip` |
| macOS (Apple silicon)  | Still a work in progress, sorry! |

Then install it from the file:
```bash
ausaxs gpu install --from <path to zip file>
```
`--from` also accepts a URL, and replaces any driver already installed. Add `--offline` to make sure no network requests are made.

## Managing the driver

| Command               | Description |
|-----------------------|-------------|
|`ausaxs gpu status`    | Show the installed driver, the chosen GPU, and whether it is used by default. |
|`ausaxs gpu disable`   | Stop using the GPU by default, without removing the driver. |
|`ausaxs gpu enable`    | Use the GPU by default again. |
|`ausaxs gpu remove`    | Remove the driver. |

For a single run, the `--gpu` and `--no-gpu` flags of `ausaxs fit`, `ausaxs em` and `ausaxs rigidbody` override the default.

If the driver does not find your GPU, see the [troubleshooting](troubleshooting.md) page.

## Next steps
You are ready to calculate your first profile. The [tutorials](tutorials.md) walk through complete analyses of real data, and the [SAXS](saxs.md) page explains what goes into the calculation.
