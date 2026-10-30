// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/SaveElement.h>

#include <fitter/FitReporter.h>
#include <fitter/FitResult.h>
#include <io/detail/trajectory/XYZWriter.h>
#include <plots/PlotDataset.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/controller/IController.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <settings/GeneralSettings.h>

#include <cassert>
#include <unordered_map>
#include <utility>

using namespace ausaxs::rigidbody::sequencer;

namespace {
    // one '%' counter per path template, so "every_%.pdb" and "best_%.pdb" are numbered independently, while two
    // elements saving to the same template continue one series instead of overwriting each other's files
    std::unordered_map<std::string, int> counters;
    std::unordered_map<std::string, ausaxs::io::detail::xyz::XYZWriter> writers;

    bool is_supported_format(const std::string& ext) {
        return ext == ".pdb" || ext == ".fit" || ext == ".txt" || ext == ".png" || ext == ".xyz";
    }
}

SaveElement::SaveElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, io::File path) 
    : owner(owner), path(std::move(path))
{}

SaveElement::~SaveElement() {
    reset_statics();
}

void SaveElement::reset_statics() {
    counters.clear();
    writers.clear();
}

void SaveElement::run() {
    // find and replace '%' with the current counter value
    auto insert_counter = [this] (io::File file) {
        if (auto pos = file.stem().find('%'); pos != std::string::npos) {
            file.stem().replace(pos, 1, std::to_string(++counters[path.path()]));
        }
        return file;
    };

    // after a rejected step the fitter still holds the rejected candidate, so it is brought up to the molecule first.
    auto fit_current_state = [this] () {
        auto* controller = owner->_get_rigidbody()->controller.get();
        controller->update_fitter();
        return controller->get_fitter()->unconstrained_fit();
    };

    // PDB
    if (const auto& ext = path.extension(); ext == ".pdb") {
        owner->_get_molecule()->save(insert_counter(path));
    }

    // FIT
    else if (ext == ".fit") {
        // the same file the fitting CLIs write: all named columns, residuals included
        auto result = fit_current_state();
        result->curves.save(
            insert_counter(path),
            "chi2=" + std::to_string(result->fval/result->dof) + " dof=" + std::to_string(result->dof)
        );
    }

    // TXT - the fit report
    else if (ext == ".txt") {
        auto result = fit_current_state();
        fitter::FitReporter::save(result.get(), insert_counter(path));
    }

    // PNG
    else if (ext == ".png") {
        auto result = fit_current_state();
        plots::PlotDataset plot;
        plot.plot_residuals(
            result->curves.select_columns({0, 1, 2, 3}),
            plots::PlotOptions(
                style::draw::errors, {
                    {"color", style::color::black}, {"logx", true}, {"logy", true}, {"xlabel", "q [$\\AA$]"}, {"ylabel", "$I(q)$"}, {"zorder", -1},
                    {"legend", "chi2=" + std::to_string(result->fval/result->dof) + ", dof=" + std::to_string(result->dof)},
                    {"title", "Iteration " + std::to_string(ausaxs::rigidbody::sequencer::LoopElement::_get_current_iteration())}
                }
            )
        );
        plot.save(insert_counter(path));
    }

    // XYZ
    else if (ext == ".xyz") {
        auto p = path.path(); 
        if (!writers.contains(p)) {
            writers.emplace(p, path);
        }
        writers.at(p).write_frame(owner->_get_molecule());
    } else {
        assert(false && "SaveElement::run: unreachable, the format is checked when parsed.");
    }
}

std::vector<std::string> SaveElement::_valid_arguments() {
    return {};
}

InlineSignature SaveElement::_valid_inline_arguments() {
    return {.names = {"path"}, .min = 1, .max = 1};
}

// save [path] - resolved relative to the output folder
std::unique_ptr<GenericElement> SaveElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) { // NOLINT
    io::File path = settings::general::output + args.inlined[0];
    if (!is_supported_format(path.extension())) { // checked here, since at run time the whole optimization has already been done
        throw except::parse_error("save", args.inlined, "Unknown file format \"" + path.extension() + "\". Valid formats are .pdb, .fit, .txt, .png and .xyz.");
    }
    return std::make_unique<SaveElement>(owner, std::move(path));
}