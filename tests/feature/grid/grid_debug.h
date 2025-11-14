#pragma once

#include <grid/Grid.h>
#include <grid/detail/GridSurfaceDetection.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <settings/GridSettings.h>

using namespace ausaxs;

class GridDebug : public grid::Grid {
    public: 
        using Grid::Grid;

		double get_atomic_radius(form_factor::form_factor_t) const override {return ra;}
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
            return {
                {
                    GridDebug::exv[0], 
                    GridDebug::exv[1], GridDebug::exv[2], GridDebug::exv[3], GridDebug::exv[4], 
                    GridDebug::exv[5], GridDebug::exv[6], GridDebug::exv[7], GridDebug::exv[8]
                }, 
                {}
            };
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
            return {
                {
                    GridDebug::exv[0], 
                    GridDebug::exv[1], GridDebug::exv[2], GridDebug::exv[3], GridDebug::exv[4], 
                    GridDebug::exv[5], GridDebug::exv[6], GridDebug::exv[7], GridDebug::exv[8]
                }, 
                {}
            };
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
            return {
                {GridDebug::exv[0]}, 
                {
                    GridDebug::exv[1], GridDebug::exv[2], GridDebug::exv[3], GridDebug::exv[4], 
                    GridDebug::exv[5], GridDebug::exv[6], GridDebug::exv[7], GridDebug::exv[8]
                }
            };
        }
};

inline grid::exv::GridExcludedVolume GridDebug::generate_excluded_volume() {
    auto res = Grid::generate_excluded_volume();

    grid::exv::GridExcludedVolume vol;
    vol.interior = {exv[0]};
    vol.surface = {exv[1], exv[2], exv[3], exv[4], exv[5], exv[6], exv[7], exv[8]};

    if (!res.has_surface()) {
        vol.interior.insert(vol.interior.end(), vol.surface.begin(), vol.surface.end());
        vol.surface.clear();
    }

    return vol;
}