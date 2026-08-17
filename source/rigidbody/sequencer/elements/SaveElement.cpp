// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/SaveElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/elements/setup/SetupElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/BodyNameRegistry.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/continuation/ContinuationState.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/Rigidbody.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <dataset/SimpleDataset.h>
#include <plots/PlotDataset.h>
#include <settings/GeneralSettings.h>
#include <io/detail/trajectory/XYZWriter.h>
#include <utility/StringUtils.h>

#include <unordered_map>

using namespace ausaxs::rigidbody::sequencer;

namespace {
    // Output state that outlives a single sequencer run. It has to: a '%' counter numbers files across an entire
    // script, and a trajectory writer must keep one stream open across every frame it is handed. Both were function-
    // local statics, which is fine for a one-shot CLI invocation but wrong for a long-lived host (the GUI) that runs
    // several refinements in one process: counters kept climbing and a trajectory writer kept pointing at the previous
    // run's file. They are gathered here so SaveElement::_reset_output_state() can clear them between runs.
    struct OutputState {
        int pdb_counter = 0;
        int fit_counter = 0;
        int png_counter = 0;
        std::unordered_map<std::string, ausaxs::io::detail::xyz::XYZWriter> trajectory_writers;
    };

    OutputState& output_state() {
        static OutputState state;
        return state;
    }
}

SaveElement::SaveElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, const io::File& path) : LoopElementCallback(owner), path(path) {}
SaveElement::~SaveElement() = default;

void SaveElement::_reset_output_state() {
    output_state() = OutputState{};
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
        owner->_get_molecule()->save(insert_counter(path, output_state().pdb_counter));
    }

    // CONTINUE - the full refinement state, so a follow-up run can resume from exactly this pose
    else if (utility::to_lowercase(ext) == continuation::continuation_extension) {
        auto names = owner->_get_sequencer()->setup()._body_name_registry().base_body_names();
        continuation::write_continuation_state(path, *owner->_get_rigidbody(), names);
    }

    // FIT
    else if (ext == ".fit") {
        int& counter = output_state().fit_counter;
        auto result = owner->_get_rigidbody()->controller->get_fitter()->fit();
        result->curves.select_columns({0, 1, 2, 3}).save(
            insert_counter(path, ++counter),
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
        plot.save(insert_counter(path, output_state().png_counter));
    }

    // XYZ
    else if (ext == ".xyz") {
        auto& writers = output_state().trajectory_writers;
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

std::unique_ptr<GenericElement> SaveElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("save", "Unexpected named argument.");}
    if (args.inlined.size() != 1) {throw except::parse_error("save", "Expected only a single inline argument.");}
    return std::make_unique<SaveElement>(owner, settings::general::output + args.inlined[0]);
}