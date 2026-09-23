#pragma once

#include <data/Body.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
#include <grid/detail/GridSurfaceDetection.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <settings/GridSettings.h>

#include <cmath>

using namespace ausaxs;

class GridDebug : public grid::Grid {
    public: 
        using Grid::Grid;

		double get_atomic_radius(form_factor::form_factor_t /*atom*/) const override {return ra;}
		double get_hydration_radius() const override {return rh;}
        void set_atomic_radius(double ra) {this->ra = ra;}
        void set_hydration_radius(double rh) {this->rh = rh;}
		grid::exv::GridExcludedVolume generate_excluded_volume() override;

        static void generate_debug_grid(data::Molecule& protein) {
            settings::grid::min_bins = 20;
            auto grid = std::make_unique<GridDebug>(protein.get_bodies());
            grid->set_atomic_radius(0);
            protein.set_grid(std::move(grid));
        }

        inline static std::vector<Vector3<double>> exv = {
            {0, 0, 0}, 
            { 1, 1, 1}, { 1, 1, -1}, { 1, -1, 1}, { 1, -1, -1}, 
            {-1, 1, 1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}
        };

        /**
         * @brief Build an excluded volume from the given debug points, including the lattice sites the transform needs.
         *        The debug points all lie on a unit lattice within [-1, 1], so they are shifted by one to make the sites non-negative.
         */
        static grid::exv::GridExcludedVolume make_exv(std::vector<Vector3<double>> interior, std::vector<Vector3<double>> surface) {
            auto to_sites = [] (const std::vector<Vector3<double>>& points) {
                std::vector<Vector3<int>> sites;
                sites.reserve(points.size());
                for (const auto& p : points) {
                    sites.emplace_back(static_cast<int>(std::lround(p.x()))+1, static_cast<int>(std::lround(p.y()))+1, static_cast<int>(std::lround(p.z()))+1);
                }
                return sites;
            };

            grid::exv::GridExcludedVolume vol;
            vol.interior_sites = to_sites(interior);
            vol.surface_sites = to_sites(surface);
            vol.interior = std::move(interior);
            vol.surface = std::move(surface);
            vol.spacing = 1;
            return vol;
        }

    private:
        double ra = 0, rh = 0;
};

/**
 * @brief Debug version of the HistogramManagerMTFFGrid class, which uses a predictable excluded volume.
 */
 template<bool vbw>
 class DebugHistogramManagerMTFFGrid : public hist::HistogramManagerMTFFGrid<vbw> {
    public:
        using hist::HistogramManagerMTFFGrid<vbw>::HistogramManagerMTFFGrid;

        grid::exv::GridExcludedVolume get_exv() const override {
            return GridDebug::make_exv(GridDebug::exv, {});
        }
};

/**
 * @brief Debug version of the HistogramManagerMTFFGridScalableExv class, which uses a predictable excluded volume.
 */
template<bool vbw>
class DebugHistogramManagerMTFFGridScalableExv : public hist::HistogramManagerMTFFGridScalableExv<vbw> {
    public:
        using hist::HistogramManagerMTFFGridScalableExv<vbw>::HistogramManagerMTFFGridScalableExv;

        grid::exv::GridExcludedVolume get_exv() const override {
            return GridDebug::make_exv(GridDebug::exv, {});
        }
};

/**
 * @brief Debug version of the HistogramManagerMTFFGridSurface class, which uses a predictable excluded volume.
 */
template<bool vbw>
class DebugHistogramManagerMTFFGridSurface : public hist::HistogramManagerMTFFGridSurface<vbw> {
    public:
        using hist::HistogramManagerMTFFGridSurface<vbw>::HistogramManagerMTFFGridSurface;

        grid::exv::GridExcludedVolume get_exv() const override {
            return GridDebug::make_exv({GridDebug::exv[0]}, std::vector<Vector3<double>>(GridDebug::exv.begin()+1, GridDebug::exv.end()));
        }
};

inline grid::exv::GridExcludedVolume GridDebug::generate_excluded_volume() {
    auto res = Grid::generate_excluded_volume();

    if (!res.has_surface()) {return make_exv(exv, {});}
    return make_exv({exv[0]}, std::vector<Vector3<double>>(exv.begin()+1, exv.end()));
}