#include <grid/exv/RawGridExv.h>

using namespace ausaxs::grid::exv;

GridExcludedVolume RawGridExv::create(observer_ptr<grid::Grid> grid) {
    std::vector<Vector3<double>> atoms;
    atoms.reserve(grid->get_volume());

    static detail::State exv_state = 
        detail::State::VOLUME   | 
        detail::State::A_CENTER | 
        detail::State::A_AREA   | 
        detail::State::W_CENTER | 
        detail::State::W_AREA
    ;

    int stride = std::max(1., std::round(settings::grid::exv::width/settings::grid::cell_width));
    int buffer = std::max(1., std::round(std::max(settings::grid::min_exv_radius, 2.)/settings::grid::cell_width));
    const auto& axes = grid->get_axes();

    auto[imin, imax] = grid->bounding_box_index();
    for (int i = std::max(imin.x()-buffer, 0); i < std::min<int>(imax.x()+buffer+1, axes.x.bins); i+=stride) {
        for (int j = std::max(imin.y()-buffer, 0); j < std::min<int>(imax.y()+buffer+1, axes.y.bins); j+=stride) {
            for (int k = std::max(imin.z()-buffer, 0); k < std::min<int>(imax.z()+buffer+1, axes.z.bins); k+=stride) {
                if (grid->index(i, j, k) & exv_state) {
                    atoms.emplace_back(grid->to_xyz(i, j, k));
                }
            }
        }
    }

    std::cout << "RawGridExv::create: added " << atoms.size() << " / " << grid->get_volume() << " atoms to the excluded volume." << std::endl;

    return GridExcludedVolume{std::move(atoms), {}};
}