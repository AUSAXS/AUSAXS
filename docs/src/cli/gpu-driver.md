---
title: Managing the GPU driver - AUSAXS
description: Check, enable, disable or remove the AUSAXS GPU driver with the ausaxs gpu commands, and override the GPU setting for a single run.
---

# Managing the GPU driver

This page lists the commands for managing an installed AUSAXS GPU driver. To install it in the first place, see [GPU acceleration](../getting-started/gpu.md).

| Command               | Description |
|-----------------------|-------------|
|`ausaxs gpu status`    | Show the installed driver, the chosen GPU, and whether it is used by default. |
|`ausaxs gpu disable`   | Stop using the GPU by default, without removing the driver. |
|`ausaxs gpu enable`    | Use the GPU by default again. |
|`ausaxs gpu remove`    | Remove the driver. |

For a single run, the `--gpu` and `--no-gpu` flags of `ausaxs fit`, `ausaxs em` and `ausaxs rigidbody` override the default.
