// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/api_helper.h>

/**
 * The main "ausaxs <tool> [options]" entrypoint. argv[0] is the program name and argv[1] the tool.
 */
extern "C" API int cli_ausaxs(int argc, char const *argv[]);

/**
 * List an additional tool in the help text of cli_ausaxs.
 */
extern "C" API void cli_ausaxs_add_tool(char const* name, char const* description);
