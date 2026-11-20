---
title: Troubleshooting - AUSAXS
description: "Fix common AUSAXS problems on Windows, macOS and Linux: unknown residues, the q unit of data files, the GPU driver and installation errors."
---

# Troubleshooting

Please refer to the following if you encounter any issues.

> [!TIP]
> Can't find your issue? Please create a [discussion thread](https://github.com/AUSAXS/AUSAXS/discussions) with more information so we can help get it solved.

## Installing

??? question "Terminal says `ausaxs` is not found?"

    The folder pip installed `ausaxs` into is not on your `PATH`. The [installation](installation.md) page shows how to add it. Until then, `python -m pyausaxs` works in place of `ausaxs` in every command, e.g. `python -m pyausaxs --version`.

??? question "Windows: why do I get SmartScreen warnings?"

    The Microsoft SmartScreen filter may warn you about potential security risks when running AUSAXS for the first time. You should be able to ignore this warning by clicking "More info" and then "Run anyway". Removing this notice from our end requires a yearly fee to sign our program, which we are currently unable to support.

    **Warning!** You should only ignore SmartScreen warnings for software you trust! AUSAXS is completely open source, so if you have any security concerns, you are free to check the code.

??? question "macOS: why does the program crash immediately upon opening?"

    The binaries are built on a recent macOS version, and may not run on older ones.

    In case they do not work on your system, you can try to compile the software yourself as described [here](compiling.md#macos).

??? question "Linux: why can't I run the downloaded executable?"

    Usually there are no issues on this platform. Perhaps you simply lack the permissions to execute the program? Try marking it as executable: `chmod +x <program>`

## Reading structures and data

??? question "Unknown atom or unknown residue?"

    AUSAXS adds the hydrogens missing from a structure file itself, and needs the definition of each residue to know how many to add. Definitions are taken from the PDB's [Chemical Component Dictionary](https://www.wwpdb.org/data/ccd), and downloaded the first time a residue is met. An atom or residue it cannot find stops the run with an error naming it. Depending on the cause:

    - **Force-field names.** Files from simulations often use force-field names for protonation states, such as `HID`, `HSE` or `CYX`, and `SOL` for water. Add `hydrogens --force-field-names` to read them by their force-field meaning.
    - **No internet access.** With `--offline`, or on a machine without internet access, residues not met before cannot be downloaded. Run once with internet access, or leave out `--offline`.
    - **All hydrogens present.** If the file already contains every hydrogen, add `hydrogens --keep`. The residue definitions are then not needed.
    - **Truly unknown.** As a last resort, `--allow-unknown-atoms` and `--allow-unknown-residues` continue anyway. No hydrogens are added to these atoms, so their scattering is underestimated.

??? question "Is my data read in the right q unit?"

    $q$ may be given in either Å⁻¹ or nm⁻¹. If the header of the data file states its unit, as `[A]` or `[nm]`, that unit is used. Otherwise it is guessed from the values: data reaching above $q = 1$ is assumed to be in nm⁻¹, data starting below $q = 10^{-3}$ in Å⁻¹, and anything in between is assumed to be in Å⁻¹ with a warning. A message is printed whenever the unit is guessed or the data is converted from nm⁻¹.

    If the guess is wrong, give the unit explicitly with `data --unit A` or `data --unit nm`. The data in `ausaxs.fit` is always in Å⁻¹.

## GPU

If AUSAXS cannot use the GPU, it says why and continues on the CPU. The results are the same either way.

??? question "No usable GPU backend, or no usable device?"

    The driver package was installed, but no GPU it can use was found. If the machine has a supported GPU, it is usually missing the graphics driver of its vendor:

    - **NVIDIA:** install the proprietary NVIDIA driver. If `nvidia-smi` lists your GPU, the driver is working.
    - **AMD on Linux:** the GPU is reached through `/dev/kfd` and `/dev/dri`. Your user must be allowed to access them, which usually means being a member of the `render` and `video` groups.
    - **AMD on Windows:** install the AMD Software (Adrenalin) driver.

    Intel GPUs and Apple silicon are not supported yet; see [GPU acceleration](gpu.md).

??? question "The chosen GPU is not present?"

    The GPU chosen with `ausaxs gpu install` could not be found, so the GPU is not used. This happens if the GPU or its driver was changed, or if the same installation is shared between machines with different GPUs. The message lists the GPUs that are present. Run `ausaxs gpu install` again to choose one of them; nothing is downloaded again.

??? question "`ausaxs gpu status` says the installation is incomplete?"

    The installation was interrupted before the driver could be tested. Run `ausaxs gpu install` again.

??? question "The GPU backend failed while trying to ...?"

    The GPU stopped working in the middle of a calculation, so it is not used for the rest of the run, and the calculation continues on the CPU. If this happens again, please [open an issue](https://github.com/AUSAXS/AUSAXS-GPU/issues/new) with the full message, your operating system and your GPU.

??? question "The GPU backend could not be loaded, and has been removed again?"

    The downloaded driver could not be loaded on this machine, so it was removed again and nothing was changed. Please [open an issue](https://github.com/AUSAXS/AUSAXS-GPU/issues/new) with the full message, your operating system and your GPU.

