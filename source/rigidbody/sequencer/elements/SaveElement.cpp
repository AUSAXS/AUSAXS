// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/SaveElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/Rigidbody.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <dataset/SimpleDataset.h>
#include <plots/PlotDataset.h>
#include <settings/GeneralSettings.h>
#include <io/detail/trajectory/XYZWriter.h>

#include <unordered_map>

using namespace ausaxs::rigidbody::sequencer;

SaveElement::SaveElement(observer_ptr<rigidbody::sequencer::LoopElement> owner, const io::File& path) : LoopElementCallback(owner), path(path) {}
SaveElement::~SaveElement() = default;

void SaveElement::run() {
    // find and replace '%' with the current counter value
    auto insert_counter = [] (io::File file, int count) {
        if (auto pos = file.stem().find('%'); pos != std::string::npos) {
            file.stem().replace(pos, 1, std::to_string(count));
        }
        return file;
    };

    // PDB
    if (const auto& ext = path.extension(); ext == ".pdb") {
        static int counter = 0;
        owner->_get_molecule()->save(insert_counter(path, ++counter));
    }

    // FIT
    else if (ext == ".fit") {
        static int counter = 0;
        auto result = owner->_get_rigidbody()->controller->get_fitter()->fit();
        result->curves.select_columns({0, 1, 2, 3}).save(
            insert_counter(path, ++counter),
            "chi2=" + std::to_string(result->fval) + ", dof=" + std::to_string(result->dof)
        );
    }

    // XYZ
    else if (ext == ".xyz") {
        static std::unordered_map<std::string, io::detail::xyz::XYZWriter> writers;
        auto p = path.path(); 
        if (!writers.contains(p)) {
            writers.emplace(p, path);
        }
        writers.at(p).write_frame(owner->_get_molecule());
    } else {
        throw std::runtime_error("SaveElement::run: Unknown file format: \"" + ext + "\"");
    }
}
