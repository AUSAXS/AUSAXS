# Managing the GPU driver

This page lists the commands for managing an installed AUSAXS GPU driver. To install it in the first place, see [GPU acceleration](https://ausaxs.github.io/getting-started/gpu/).

| Command               | Description |
|-----------------------|-------------|
|`ausaxs gpu status`    | Show the installed driver, the chosen GPU, and whether it is used by default. |
|`ausaxs gpu disable`   | Stop using the GPU by default, without removing the driver. |
|`ausaxs gpu enable`    | Use the GPU by default again. |
|`ausaxs gpu remove`    | Remove the driver. |

For a single run, the `--gpu` and `--no-gpu` flags override the default.
