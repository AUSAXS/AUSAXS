#pragma once

#include <grid/Grid.h>
#include <grid/detail/GridSurfaceDetection.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/GridSettings.h>

using namespace ausaxs;

class GridDebug : public grid::Grid {
    public: 
        using Grid::Grid;

		double get_atomic_radius(form_factor::form_factor_t) const override {return ra;}
		double get_hydration_radius() const override {return rh;}
        void set_atomic_radius(double ra) {this->ra = ra;}
        void set_hydration_radius(double rh) {this->rh = rh;}
		grid::detail::GridExcludedVolume generate_excluded_volume(bool determine_surface) override;

        static void generate_debug_grid(data::Molecule& protein) {
            settings::grid::min_bins = 20;
            auto grid = std::make_unique<GridDebug>(protein.get_bodies());
            grid->set_atomic_radius(0);
            protein.set_grid(std::move(grid));
        }

        std::vector<Vector3<double>> exv = {
            {0, 0, 0}, 
            { 1, 1, 1}, { 1, 1, -1}, { 1, -1, 1}, { 1, -1, -1}, 
            {-1, 1, 1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}
        };
    
    private:
        double ra = 0, rh = 0;
};

inline grid::detail::GridExcludedVolume GridDebug::generate_excluded_volume(bool determine_surface) {
    grid::detail::GridExcludedVolume vol;
    vol.interior = {exv[0]};
    vol.surface = {exv[1], exv[2], exv[3], exv[4], exv[5], exv[6], exv[7], exv[8]};

    if (!determine_surface) {
        vol.interior.insert(vol.interior.end(), vol.surface.begin(), vol.surface.end());
        vol.surface.clear();
    }

    return vol;
}