#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
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
    conformation = std::make_unique<rigidbody::detail::Conformation>(molecule.get_bodies());
}

Rigidbody::Rigidbody(Rigidbody&& other) = default;
Rigidbody& Rigidbody::operator=(Rigidbody&& other) = default;

void Rigidbody::refresh_grid() {
    auto grid = molecule.get_grid();
    std::pair<Vector3<double>, Vector3<double>> bounds;
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