// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/api_helper.h>

/**
 * The main "ausaxs <tool> [options]" entrypoint. argv[0] is the program name and argv[1] the tool.
 * Returns the tool's exit code. Safe to call repeatedly in the same process.
 */
extern "C" API int cli_ausaxs(int argc, char const *argv[]);

/**
 * List an additional tool in the help text of cli_ausaxs.
 * Meant for wrappers that implement extra tools on top of the backend, and so intercept them before calling cli_ausaxs.
 * Adding a tool with an already listed name replaces its description.
 */
extern "C" API void cli_ausaxs_add_tool(char const* name, char const* description);
