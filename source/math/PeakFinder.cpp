/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <math/PeakFinder.h>
#include <algorithm>
#include <stdexcept>

using namespace ausaxs;

/*
    Alternative idea:
        Connect all local maxima with a line. This is essentially what we're doing now anyway, except we're doing it in a more complicated roundabout way.
*/

namespace peak_finder {struct Limit {double min = 0, max = 0;};}
std::vector<unsigned int> math::find_minima(const std::vector<double>& x, const std::vector<double>& y, unsigned int min_spacing, double min_prominence) {
	if (x.size() != y.size()) {throw std::invalid_argument("math::find_minima: x and y must have the same size.");}
    if (x.size() < 3) {return {};}
	unsigned int size = x.size();

    /**
     * @brief Get the linear equation of the line between the bounds of a local minimum
     * 
     * @return std::pair<double, double> The first element is the slope, the second is the y-intercept
     */
    auto bounds_eq = [&] (const peak_finder::Limit& bound) {
        double x1 = x[bound.min];
        double x2 = x[bound.max];
        double y1 = bound.min == 0      && y[bound.min] < y[bound.max] ? y[bound.max] : y[bound.min]; // account for the fact that we might be at the edge of the array
        double y2 = bound.max == size-1 && y[bound.max] < y[bound.min] ? y[bound.min] : y[bound.max]; // account for the fact that we might be at the edge of the array
        double a = (y2-y1)/(x2-x1);
        double b = y1 - a*x1;
        return std::make_pair(a, b);
    };

    /**
     * @brief Calculates the prominence of a local minimum
     * 
     * @param bound The bounds of the local minimum
     * @param xmin The x value of the local minimum
     * @param ymin The y value of the local minimum
     */
    auto calc_prominence = [&] (const peak_finder::Limit& bound, double xmin, double ymin) {
        // we want to intrapolate a line between the bounds and calculate the difference between it and the minimum
        auto [a, b] = bounds_eq(bound);
        double prominence = (a*xmin + b) - ymin;
        return prominence;
    };

    /**
     * @brief Relax the bounds of a local minimum to make its slope smaller.
     * 
     * @param bound The bounds of the local minimum.
     * @param index The index of the local minimum.
     */
    auto relax_bound = [&] (peak_finder::Limit& bound, unsigned int index) {
        if (bound.max - bound.min < 5) {return;}
        if (bound.min == 0 || bound.max == size-1) {return;}

        // determine which direction we have to move
        double slope = std::abs((y[bound.max] - y[bound.min]) / (x[bound.max] - x[bound.min]));
        double slope_left = std::abs((y[bound.max] - y[bound.min+1]) / (x[bound.max] - x[bound.min+1]));
        double slope_right = std::abs((y[bound.max-1] - y[bound.min]) / (x[bound.max-1] - x[bound.min]));

        // std::cout << "MOVING BOUNDS" << std::endl;
        // std::cout << "\tSlope is " << slope << std::endl;
        // std::cout << "\tSlope left is " << slope_left << std::endl;
        // std::cout << "\tSlope right is " << slope_right << std::endl;

        if (slope_left < slope_right) {
            // move the left bound
            bound.min++;
            do {
                // std::cout << "\tMoving left bound" << std::endl;
                // std::cout << "\t\tSlope is " << slope << std::endl;
                // std::cout << "\t\tSlope left is " << slope_left << std::endl;
                slope = slope_left*1.05;
                bound.min++;
                if (bound.min == index) {
                    bound.min--;
                    break;
                }
                slope_left = std::abs((y[bound.max] - y[bound.min]) / (x[bound.max] - x[bound.min]));
            }
            while (slope_left < slope);
            bound.min--;
        } else {
            // move the right bound
            bound.max--;
            do {
                // std::cout << "\tMoving right bound" << std::endl;
                // std::cout << "\t\tSlope is " << slope << std::endl;
                // std::cout << "\t\tSlope right is " << slope_right << std::endl;
                slope = slope_right*1.05;
                bound.max--;
                if (bound.max == index) {
                    bound.max++;
                    break;
                }
                slope_right = std::abs((y[bound.max] - y[bound.min]) / (x[bound.max] - x[bound.min]));

            }
            while (slope_right < slope);
            bound.max++;
        }
    };

    //######################################################//
    //###                FIND ALL MINIMA                 ###//
    //######################################################//
    std::vector<unsigned int> local_minima;
    {
        // special cases: first and last point
        if (y[0] < y[1]) {local_minima.push_back(0);}
        for (unsigned int i = 1; i < size-1; i++) {
            if (y[i] < y[i-1] && y[i] < y[i+1]) {
                local_minima.push_back(i);
            }
        }
        if (y[size-1] < y[size-2]) {local_minima.push_back(size-1);}
        if (local_minima.empty()) {return {};}
    }

    //######################################################//
    //###                ESTIMATE BOUNDS                 ###//
    //######################################################//
    std::vector<peak_finder::Limit> local_minima_bounds;
    local_minima_bounds.reserve(local_minima.size());
    {
        for (int index : local_minima) {
            peak_finder::Limit bounds;

            // first go left
            {
                int decreasing_left = 0;              // how many points in a row are lower than the maximal found point
                int left = std::max<int>(index-1, 0); // the index of the point we are currently looking at
                double max_diff = 0;                  // the highest point we have found so far
                while (decreasing_left < 3 && left >= 0) {
                    double diff = y[left] - y[index];
                    if (diff < max_diff*math::detail::min_slope) { // we want at least a min_slope increase for every point
                        decreasing_left++;
                    } else {
                        decreasing_left = 0;
                    }
                    max_diff = std::max(max_diff, diff);
                    left--;
                }
                bounds.min = std::max<int>(left+decreasing_left+1, 0);
            }

            // then go right
            {
                int decreasing_right = 0;              	    // how many points in a row are lower than the maximal found point
                int right = std::min<int>(index+1, size-1);	// the index of the point we are currently looking at
                double max_diff = 0;                        // the highest point we have found so far
                while (decreasing_right < 3 && right < int(size)) {
                    double diff = y[right] - y[index];
                    if (diff < max_diff*math::detail::min_slope) { // we want at least a 0% increase for every point
                        decreasing_right++;
                    } else {
                        decreasing_right = 0;
                    }
                    max_diff = std::max(max_diff, diff);
                    right++;
                }
                bounds.max = std::min<int>(right-decreasing_right-1, size-1);
            }

            // check for intersections
            {
                // structured binding not allowed due to Clang bug
                auto tmp = bounds_eq(bounds);
                auto& a = tmp.first;
                auto& b = tmp.second;
                auto fx = [&] (double x) {return a*x + b;};

                // check for intersection from the left
                for (int i = index-1; bounds.min < i; --i) {
                    if (fx(x[i]) < y[i]) {
                        bounds.min = i;
                        break;
                    }
                }

                // check for intersection from the right
                for (int i = index+1; i < bounds.max; ++i) {
                    if (fx(x[i]) < y[i]) {
                        bounds.max = i;
                        break;
                    }
                }
            }

            local_minima_bounds.push_back(std::move(bounds));
			#if DEBUG_PLOT
                for (unsigned int i = 0; i < local_minima.size(); i++) {
                    Limit& bound = local_minima_bounds[i];
                    unsigned int index = local_minima[i];

                    // bounds
    				SimpleDataset dummy1({x[bounds.min], x[bounds.max]}, {y[bounds.min], y[bounds.max]});
                    plot.plot(dummy1, plots::PlotOptions({{"color", style::color::green}, {"lw", 0.5}, {"zorder", 0}}));

                    // prominence
                    SimpleDataset dummy2({x[index], x[index]}, {y[index], y[index] + calc_prominence(bound, x[index], y[index])});
                    plot.plot(dummy2, plots::PlotOptions({{"color", style::color::green}, {"lw", 0.5}, {"zorder", 0}}));
                }
			#endif
        }
    }

    //######################################################//
    //###            MERGE OVERLAPPING BOUNDS            ###//
    //######################################################//
    if (local_minima_bounds.size() > 1) {
        std::vector<peak_finder::Limit> merged_bounds;
        std::vector<unsigned int> merged_minima;
        merged_bounds.reserve(local_minima.size());
        for (unsigned int i = 0; i < local_minima_bounds.size(); i++) {
            peak_finder::Limit bounds = local_minima_bounds[i];

            // go through all bounds and merge with the next one if they overlap
            unsigned int merge_count = 0;
            for (; i+merge_count+1 < local_minima_bounds.size(); ++merge_count) {
                if (bounds.max <= local_minima_bounds[i+merge_count+1].min) {break;}
                bounds.max = local_minima_bounds[i+merge_count+1].max;
            }

            // no merging done, just add the bounds & continue
            if (merge_count == 0) {
                merged_bounds.push_back(std::move(bounds));
                merged_minima.push_back(local_minima[i]);
                continue;
            }

            // go through all merged minima and find the one with the lowest y value
            unsigned int merged_minima_index = local_minima[i];
            for (unsigned int j = 1; j <= merge_count; j++) {
                if (y[local_minima[i+j]] < y[merged_minima_index]) {
                    merged_minima_index = local_minima[i+j];
                }
            }

            i += merge_count;
            merged_bounds.push_back(std::move(bounds));
            merged_minima.push_back(merged_minima_index);
        }

        if (!merged_bounds.empty()) {
            local_minima = std::move(merged_minima);
            local_minima_bounds = std::move(merged_bounds);
        }
    }

    //######################################################//
    //###                  RELAX BOUNDS                  ###//
    //######################################################//
    auto original_bounds = local_minima_bounds;
    for (unsigned int i = 0; i < local_minima.size(); ++i) {
        relax_bound(local_minima_bounds[i], local_minima[i]);
    }

    //######################################################//
    //###        FILTER PROMINENCE & MIN_SPACING         ###//
    //######################################################//
    if (0 < min_prominence) {
        // calculate all prominences
        std::vector<double> prominences(local_minima.size());
        for (unsigned int i = 0; i < local_minima.size(); i++) {
            prominences[i] = calc_prominence(local_minima_bounds[i], x[local_minima[i]], y[local_minima[i]]);
        }

        // update minimum prominence
        min_prominence *= *std::max_element(prominences.begin(), prominences.end());

        // filter out all local minima with a prominence smaller than the minimum prominence
        std::vector<unsigned int> filtered_minima;
        std::vector<peak_finder::Limit> filtered_bounds;
        for (unsigned int i = 0; i < local_minima.size(); ++i) {
            // if the prominence is smaller than the minimum prominence, we remove it
            if (prominences[i] < min_prominence) {

                // check if we can merge the bounds with a neighbouring minima
                // first check previous bound
                if (i != 0 && local_minima_bounds[i-1].max == local_minima_bounds[i].min) {
                    peak_finder::Limit new_bounds;
                    // restore original bounds and merge
                    new_bounds.min = original_bounds[i-1].min;
                    new_bounds.max = original_bounds[i].max;
                    unsigned int new_minima = y[local_minima[i-1]] < y[local_minima[i]] ? local_minima[i-1] : local_minima[i];
                    relax_bound(new_bounds, new_minima);

                    // recalculate prominence
                    double new_prominence = calc_prominence(new_bounds, x[new_minima], y[new_minima]);

                    // if the old prominence was better, discard the merge
                    if (new_prominence <= prominences[i-1]) {
                        continue;
                    }

                    // otherwise update the bounds and prominence
                    local_minima[i-1] = new_minima;
                    local_minima_bounds[i-1] = new_bounds;

                    // if the previous bound was not kept
                    if (prominences[i-1] < min_prominence) {
                        // check if the new prominence is large enough
                        if (min_prominence < new_prominence) {
                            filtered_minima.push_back(local_minima[i-1]);
                            filtered_bounds.push_back(local_minima_bounds[i-1]);
                        }
                    }

                    // else we just modify what we already have
                    else {
                        filtered_bounds.back() = local_minima_bounds[i-1];
                        filtered_minima.back() = local_minima[i-1];
                    }
                }

                // check next bound
                if (i != size-1 && local_minima_bounds[i].max == local_minima_bounds[i+1].min) {
                    peak_finder::Limit new_bounds;
                    // restore original bounds and merge
                    new_bounds.min = original_bounds[i].min;
                    new_bounds.max = original_bounds[i+1].max;
                    unsigned int new_minima = y[local_minima[i+1]] < y[local_minima[i]] ? local_minima[i+1] : local_minima[i];
                    relax_bound(new_bounds, new_minima);

                    // recalculate prominence
                    double new_prominence = calc_prominence(new_bounds, x[new_minima], y[new_minima]);
                    
                    // if the old prominence was better, discard the merge
                    if (new_prominence <= prominences[i+1]) {
                        continue;
                    }

                    // otherwise update the bounds and prominence
                    local_minima[i+1] = new_minima;
                    local_minima_bounds[i+1] = new_bounds;
                    prominences[i+1] = new_prominence;
                }
            } 
            
            // else prominence is larger than minimum prominence, so we keep it
            else {
                filtered_minima.push_back(local_minima[i]);
                filtered_bounds.push_back(local_minima_bounds[i]);
            }
        }

        local_minima = std::move(filtered_minima);
        local_minima_bounds = std::move(filtered_bounds);
    }

    // now we want to filter out the ones that are too close to each other
    if (0 != min_spacing) {
        std::vector<unsigned int> filtered_minima = {local_minima.front()};
        for (int i : local_minima) {
            if (min_spacing <= i - filtered_minima.back()) {
                filtered_minima.push_back(i);
            } else {
                if (y[i] < y[filtered_minima.back()]) {
                    filtered_minima.back() = i;
                }
            }
        }
        local_minima = std::move(filtered_minima);
    }

    //? TODO: perform a Gaussian fit of each minima to check quality?
    return local_minima;
}