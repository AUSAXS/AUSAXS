// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/SaveElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/Rigidbody.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <dataset/SimpleDataset.h>
#include <plots/PlotDataset.h>
#include <settings/GeneralSettings.h>
#include <io/detail/trajectory/XYZWriter.h>

#include <unordered_map>

using namespace ausaxs::rigidbody::sequencer;

namespace {
    int pdb_counter = 0;
    int fit_counter = 0;
    int png_counter = 0;
    std::unordered_map<std::string, ausaxs::io::detail::xyz::XYZWriter> writers;
}

SaveElement::SaveElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, const io::File& path) : LoopElementCallback(owner), path(path) {}
SaveElement::~SaveElement() {
    reset_statics();
}

void SaveElement::reset_statics() {
    pdb_counter = 0;
    fit_counter = 0;
    png_counter = 0;
    writers.clear();
}

void SaveElement::run() {
    // find and replace '%' with the current counter value
    auto insert_counter = [] (io::File file, int& count) {
        if (auto pos = file.stem().find('%'); pos != std::string::npos) {
            file.stem().replace(pos, 1, std::to_string(++count));
        }
        return file;
    };

    // PDB
    if (const auto& ext = path.extension(); ext == ".pdb") {
        owner->_get_molecule()->save(insert_counter(path, pdb_counter));
    }

    // FIT
    else if (ext == ".fit") {
        auto result = owner->_get_rigidbody()->controller->get_fitter()->fit();
        result->curves.select_columns({0, 1, 2, 3}).save(
            insert_counter(path, ++fit_counter),
            "chi2=" + std::to_string(result->fval) + ", dof=" + std::to_string(result->dof)
        );
    }

    // PNG
    else if (ext == ".png") {
        auto result = owner->_get_rigidbody()->controller->get_fitter()->fit();
        plots::PlotDataset plot;
        plot.plot_residuals(
            result->curves.select_columns({0, 1, 2, 3}),
            plots::PlotOptions(
                style::draw::errors, {
                    {"color", style::color::black}, {"logx", true}, {"logy", true}, {"xlabel", "q [$\\AA$]"}, {"ylabel", "$I(q)$"}, {"zorder", -1},
                    {"legend", "chi2=" + std::to_string(result->fval/result->dof) + ", dof=" + std::to_string(result->dof)},
                    {"title", "Iteration " + std::to_string(owner->_get_current_iteration())}
                }
            )
        );
        plot.save(insert_counter(path, png_counter));
    }

    // XYZ
    else if (ext == ".xyz") {
        auto p = path.path(); 
        if (!writers.contains(p)) {
            writers.emplace(p, path);
        }
        writers.at(p).write_frame(owner->_get_molecule());
    } else {
        throw std::runtime_error("SaveElement::run: Unknown file format: \"" + ext + "\"");
    }
}

std::vector<std::string> SaveElement::_valid_arguments() {
    return {};
}

InlineSignature SaveElement::_valid_inline_arguments() {
    return {.names = {"path"}, .min = 1, .max = 1};
}

// save [path] - resolved relative to the output folder
std::unique_ptr<GenericElement> SaveElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    return std::make_unique<SaveElement>(owner, settings::general::output + args.inlined[0]);
}