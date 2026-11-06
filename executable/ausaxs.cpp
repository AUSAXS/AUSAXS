// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/cli/cli_ausaxs.h>

#include <iostream>

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    return cli_ausaxs(argc, argv);
}
