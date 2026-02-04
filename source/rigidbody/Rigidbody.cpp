// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/parameters/OptimizableSymmetryStorage.h>
#include <rigidbody/controller/ControllerFactory.h>
#include <rigidbody/detail/Conformation.h>
#include <fitter/SmartFitter.h>
#include <fitter/FitResult.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <grid/Grid.h>
#include <utility/Logging.h>

using namespace ausaxs::rigidbody;

Rigidbody::~Rigidbody() = default;

Rigidbody::Rigidbody(data::Molecule&& _molecule) : molecule(std::move(_molecule)) {
    // Convert all symmetry storages to OptimizableSymmetryStorage for parameter optimization
    for (int i = 0; i < static_cast<int>(molecule.size_body()); ++i) {
        auto& body = molecule.get_body(i);
        if (auto obj = dynamic_cast<symmetry::OptimizableSymmetryStorage*>(body.symmetry().get_obj()); !obj) {
            body.symmetry().set_obj(std::make_unique<symmetry::OptimizableSymmetryStorage>(std::move(*body.symmetry().get_obj())));
        }
    }

    molecule.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramManagerMT);
    controller = factory::create_controller(this);
    body_selector = factory::create_selection_strategy(this);
    transformer = factory::create_transform_strategy(this);
    parameter_generator = factory::create_parameter_strategy(
        this, 
        settings::rigidbody::iterations,
        5,
        std::numbers::pi/3
    );
    constraints = std::make_unique<constraints::ConstraintManager>(this);
    conformation = std::make_unique<rigidbody::detail::Conformation>(this);
}

Rigidbody::Rigidbody(Rigidbody&& other) = default;
Rigidbody& Rigidbody::operator=(Rigidbody&& other) = default;

void Rigidbody::refresh_grid() {
    auto grid = molecule.get_grid();
    std::pair<Vector3<double>, Vector3<double>> bounds;
    bounds.first = Vector3<double>{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    bounds.second = Vector3<double>{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
    
    for (const auto& body : molecule.get_bodies()) {
        auto [min, max] = grid::Grid::bounding_box(body.get_atoms());
        auto w = body.get_waters();
        if (w.has_value()) {
            auto [min1, max1] = grid::Grid::bounding_box(w.value().get());
            for (int i = 0; i < 3; ++i) {
                min[i] = std::min(min[i], min1[i]);
                max[i] = std::max(max[i], max1[i]);
            }
        }
        
        // Account for symmetry bodies by computing their expected bounds
        for (std::size_t j = 0; j < body.size_symmetry(); ++j) {
            auto& sym = body.symmetry().get(j);
            auto cm = body.get_cm(false);
            
            for (int rep = 1; rep <= sym.repetitions; ++rep) {
                auto transform = sym.template get_transform<double>(cm, rep);
                
                // Transform the 8 corners of the bounding box to get symmetry-transformed bounds
                auto [body_min, body_max] = grid::Grid::bounding_box(body.get_atoms());
                std::vector<Vector3<double>> corners = {
                    {body_min.x(), body_min.y(), body_min.z()},
                    {body_max.x(), body_min.y(), body_min.z()},
                    {body_min.x(), body_max.y(), body_min.z()},
                    {body_max.x(), body_max.y(), body_min.z()},
                    {body_min.x(), body_min.y(), body_max.z()},
                    {body_max.x(), body_min.y(), body_max.z()},
                    {body_min.x(), body_max.y(), body_max.z()},
                    {body_max.x(), body_max.y(), body_max.z()}
                };
                
                for (const auto& corner : corners) {
                    auto transformed = transform(corner);
                    min[0] = std::min(min[0], transformed.x());
                    min[1] = std::min(min[1], transformed.y());
                    min[2] = std::min(min[2], transformed.z());
                    max[0] = std::max(max[0], transformed.x());
                    max[1] = std::max(max[1], transformed.y());
                    max[2] = std::max(max[2], transformed.z());
                }
            }
        }
        
        for (int i = 0; i < 3; ++i) {
            bounds.first[i] = std::min(bounds.first[i], min[i]);
            bounds.second[i] = std::max(bounds.second[i], max[i]);
        }
    }

    auto grid_bounds = grid->get_axes();
    if (grid_bounds.x.min < bounds.first.x() && grid_bounds.x.max > bounds.second.x() &&
        grid_bounds.y.min < bounds.first.y() && grid_bounds.y.max > bounds.second.y() &&
        grid_bounds.z.min < bounds.first.z() && grid_bounds.z.max > bounds.second.z()) {
        return;
    }

    logging::log("Rigidbody::refresh_grid: Refreshing grid to fit current conformation.");
    molecule.clear_grid();
    molecule.create_grid();
}