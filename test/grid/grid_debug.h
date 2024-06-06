#pragma once

#include <grid/Grid.h>
#include <grid/detail/GridSurfaceDetection.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <settings/GridSettings.h>

class GridDebug : public grid::Grid {
    public: 
        using Grid::Grid;

		double get_atomic_radius(constants::atom_t) const override {return ra;}
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
    
    private:
        double ra = 0, rh = 0;
};

inline grid::detail::GridExcludedVolume GridDebug::generate_excluded_volume(bool) {
    grid::detail::GridExcludedVolume vol;
    vol.interior = {{0, 0, 0}};
    vol.surface = {
        { 1, 1, 1}, { 1, 1, -1}, { 1, -1, 1}, { 1, -1, -1}, 
        {-1, 1, 1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}
    };
    return vol;
}