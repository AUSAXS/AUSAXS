// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/cli/cli_saxs_fitter.h>
#include <api/cli/cli_em_fitter.h>
#include <api/cli/cli_rigidbody.h>
#include <utility/Console.h>

#include <unordered_map>

using namespace ausaxs;

namespace {
    enum class Tool {Fit, EM, Rigidbody};
    std::unordered_map<std::string, Tool> tool_map {
        {"fit", Tool::Fit},
        {"saxs_fitter", Tool::Fit},
        {"em", Tool::EM},
        {"em_fitter", Tool::EM},
        {"rigidbody", Tool::Rigidbody},
        {"rigidbody_optimizer", Tool::Rigidbody}
    };
}

int main(int argc, char const *argv[]) {
    auto print_help = [&] () {
        console::print_text(
            "\nUsage: ausaxs <tool> [options]"
            "\n\nAvailable tools:"
            "\n  fit        - Fit SAXS data to a structure"
            "\n  em         - Fit EM map to SAXS data"
            "\n  rigidbody  - Rigid-body optimization"
            "\n\nFor tool-specific help:"
            "\n  ausaxs <tool> --help"            
        );
    };

    if (argc < 2) {
        print_help();
        return 1;
    }

    std::string tool = argv[1];
    if (!tool_map.contains(tool)) {
        console::print_warning("Unknown tool: " + tool);
        print_help();
        return 1;
    }
    switch (tool_map[tool]) {
        case Tool::Fit:
            return cli_saxs_fitter(argc-1, argv+1);
        case Tool::EM:
            return cli_em_fitter(argc-1, argv+1);
        case Tool::Rigidbody:
            return cli_rigidbody(argc-1, argv+1);
    }
}