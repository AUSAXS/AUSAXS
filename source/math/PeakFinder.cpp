#include <math/PeakFinder.h>
#include <plots/PlotDataset.h>
#include <dataset/SimpleDataset.h>
#include <utility/Exceptions.h>

#include <algorithm>

#define DEBUG_PLOT true
#define DEBUG_OUTPUT true

std::vector<unsigned int> math::find_minima(const std::vector<double>& x, std::vector<double> y, unsigned int min_spacing, double min_prominence) {
	if (x.size() != y.size()) {throw except::invalid_argument("math::find_minima: x and y must have the same size.");}
    if (x.size() < 3) {return {};}
	unsigned int size = x.size();

	#ifdef DEBUG_PLOT
		plots::PlotDataset plot;
		{
			SimpleDataset data(x, y);
			plot.plot(data);
		}
	#endif

    // first we want to find all local minima
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

	#ifdef DEBUG_OUTPUT
		std::cout << "Found " << local_minima.size() << " local minima." << std::endl;
		for (int i : local_minima) {
			std::cout << "\t(" << x[i] << ", " << y[i] << ")" << std::endl;
		}
	#endif

    // now we want to filter out the ones that are too close to each other
    if (min_spacing > 0) {
        std::vector<unsigned int> filtered_minima;
        for (int i : local_minima) {
            if (filtered_minima.empty()) {
                filtered_minima.push_back(i);
            } else {
                if (min_spacing <= i - filtered_minima.back()) {
                    filtered_minima.push_back(i);
                }
            }
        }
        local_minima = std::move(filtered_minima);
    }

	#ifdef DEBUG_OUTPUT
		std::cout << "After filtering spacing " << local_minima.size() << " local minima." << std::endl;
		for (int i : local_minima) {
			std::cout << "\t(" << x[i] << ", " << y[i] << ")" << std::endl;
		}
	#endif

    // estimate the bounds of each minima
    std::vector<Limit> local_minima_bounds;
    local_minima_bounds.reserve(local_minima.size());
    {
        for (int i : local_minima) {
            Limit bounds;

            // first go left
            {
                int decreasing_left = 0;            // how many points in a row are lower than the maximal found point
                int left = std::max<int>(i-1, 0);   // the index of the point we are currently looking at
                double max_diff = 0;                // the highest point we have found so far
                // std::cout << "checking minima (" << x(i) << ", " << y(i) << ")" << std::endl;
                while (decreasing_left < 3 && left >= 0) {
                    // std::cout << "\tleft: (" << x(left) << ", " << y(left) << ")" << std::endl;
                    double diff = y[left] - y[i];
                    if (diff < max_diff*1.05) {     // we want at least a 5% increase for every point
                        decreasing_left++;
                        // std::cout << "\t\tnot accepted. decreasing left:" << decreasing_left << std::endl;
                    } else {
                        decreasing_left = 0;
                        // std::cout << "\t\taccepted." << std::endl;
                    }
                    max_diff = std::max(max_diff, diff);
                    left--;
                }
                bounds.min = std::max<int>(left+decreasing_left+1, 0);
            }

            // then go right
            {
                int decreasing_right = 0;              	// how many points in a row are lower than the maximal found point
                int right = std::min<int>(i+1, size-1);	// the index of the point we are currently looking at
                double max_diff = 0;                    // the highest point we have found so far
                while (decreasing_right < 3 && right < size) {
                    // std::cout << "\tright: (" << x(right) << ", " << y(right) << ")" << std::endl;
                    double diff = y[right] - y[i];
                    if (diff < max_diff*1.05) {             // we want at least a 5% increase for every point
                        // std::cout << "\t\tnot accepted. decreasing right:" << decreasing_right << std::endl;
                        decreasing_right++;
                    } else {
                        decreasing_right = 0;
                        // std::cout << "\t\taccepted." << std::endl;
                    }
                    max_diff = std::max(max_diff, diff);
                    right++;
                }
                bounds.max = std::min<int>(right-decreasing_right-1, size-1);
            }

			// #ifdef DEBUG_PLOT
			// 	SimpleDataset dummy({x[bounds.min], x[bounds.max]}, {y[bounds.min], y[bounds.max]});
			// 	dummy.add_plot_options({{"color", style::color::green}, {"ls", style::line::dashed}});
			// 	plot.plot(dummy);
			// #endif
			#ifdef DEBUG_OUTPUT
				std::cout << "\tBounds: (" << x[bounds.min] << ", " << y[bounds.min] << ") -> (" << x[bounds.max] << ", " << y[bounds.max] << ")" << std::endl;
			#endif
            if (bounds.min == 0 && y[0] < y[bounds.max]) {y[0] = y[bounds.max];}
            else if (bounds.max == size-1 && y[size-1] < y[bounds.min]) {y[size-1] = y[bounds.min];}
            local_minima_bounds.push_back(std::move(bounds));
        }

        // merge overlapping bounds
        if (local_minima_bounds.size() > 1) {
            std::vector<Limit> merged_bounds;
            std::vector<unsigned int> merged_minima;
            merged_bounds.reserve(local_minima.size());
            for (unsigned int i = 0; i < local_minima_bounds.size()-1; i++) {
                Limit bounds = local_minima_bounds[i];

                // go through all bounds and merge with the next one if they overlap
                unsigned int merge_count = 1;
                for (; i+merge_count < local_minima_bounds.size(); ++merge_count) {
                    if (bounds.max <= local_minima_bounds[i+merge_count].min) {break;}
					#ifdef DEBUG_OUTPUT
            	        std::cout << "Merging bounds (" << bounds.min << ", " << bounds.max << ") and (" << local_minima_bounds[i+merge_count].min << ", " << local_minima_bounds[i+merge_count].max << "). ";
					#endif

                    bounds.max = local_minima_bounds[i+merge_count].max;

					#ifdef DEBUG_OUTPUT
            	        std::cout << "New bounds: (" << bounds.min << ", " << bounds.max << ")" << std::endl;
					#endif
                }
				merge_count--;
                if (merge_count == 0) {continue;}

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

		#ifdef DEBUG_OUTPUT
			std::cout << "After merging bounds " << local_minima.size() << " local minima." << std::endl;
			for (auto bound : local_minima_bounds) {
				std::cout << "\tBounds: (" << x[bound.min] << ", " << y[bound.min] << ") -> (" << x[bound.max] << ", " << y[bound.max] << ")" << std::endl;
			}
		#endif

        // we have to move one end of the bounds to make them more horizontal
        auto relax_bound = [&] (Limit& bound, unsigned int index) {
            if (bound.max - bound.min < 5) {return;}
            if (bound.min == 0 || bound.max == size-1) {return;}

            // determine which direction we have to move
            double slope = std::abs((y[bound.max] - y[bound.min]) / (x[bound.max] - x[bound.min]));
            double slope_left = std::abs((y[bound.max] - y[bound.min+1]) / (x[bound.max] - x[bound.min+1]));
            double slope_right = std::abs((y[bound.max-1] - y[bound.min]) / (x[bound.max-1] - x[bound.min]));

            std::cout << "MOVING BOUNDS" << std::endl;
            std::cout << "\tSlope is " << slope << std::endl;
            std::cout << "\tSlope left is " << slope_left << std::endl;
            std::cout << "\tSlope right is " << slope_right << std::endl;

            if (slope_left < slope_right) {
                // move the left bound
                bound.min++;
                do {
                    std::cout << "\tMoving left bound" << std::endl;
                    std::cout << "\t\tSlope is " << slope << std::endl;
                    std::cout << "\t\tSlope left is " << slope_left << std::endl;
                    slope = slope_left*1.05;
                    bound.min++;
                    if (bound.min == index) {
                        bound.min--;
                        break;
                    }
                    slope_left = std::abs((y[bound.max] - y[bound.min]) / (x[bound.max] - x[bound.min]));
                }
                while (slope_left < slope);
            } else {
                // move the right bound
                bound.max--;
                do {
                    std::cout << "\tMoving right bound" << std::endl;
                    std::cout << "\t\tSlope is " << slope << std::endl;
                    std::cout << "\t\tSlope right is " << slope_right << std::endl;
                    slope = slope_right*1.05;
                    bound.max--;
                    if (bound.max == index) {
                        bound.max++;
                        break;
                    }
                    slope_right = std::abs((y[bound.max] - y[bound.min]) / (x[bound.max] - x[bound.min]));

                }
                while (slope_right < slope);
            }
        };

        auto original_bounds = local_minima_bounds;
        for (unsigned int i = 0; i < local_minima.size(); ++i) {
            relax_bound(local_minima_bounds[i], local_minima[i]);
        }

        // now we want to filter based on the prominence. 
        auto calc_prominence = [&] (unsigned int index) {
            // we want to intrapolate a line between the bounds and calculate the difference between it and the minimum
            double x1 = x[local_minima_bounds[index].min];
            double x2 = x[local_minima_bounds[index].max];
            double y1 = y[local_minima_bounds[index].min];
            double y2 = y[local_minima_bounds[index].max];
            double a = (y2-y1)/(x2-x1);
            double b = y1 - a*x1;
            double prominence = std::abs(y[local_minima[index]] - (a*x[local_minima[index]] + b));
            return prominence;
        };

        if (min_prominence > 0) {
            // calculate all prominences
            std::vector<double> prominences(local_minima.size());
            for (unsigned int i = 0; i < local_minima.size(); i++) {
                prominences[i] = calc_prominence(i);
                std::cout << "prominence " << i << " is " << prominences[i] << std::endl;
            }

            // update minimum prominence
            min_prominence *= *std::max_element(prominences.begin(), prominences.end());
            std::cout << "Minimum prominence is " << min_prominence << std::endl;

            // filter out all local minima with a prominence smaller than the minimum prominence
            std::vector<unsigned int> filtered_minima;
            std::vector<Limit> filtered_bounds;
            for (unsigned int i = 0; i < local_minima.size(); i++) {
                // if the prominence is smaller than the minimum prominence, we remove it
                if (prominences[i] < min_prominence) {
                    // check if we can merge the bounds with a neighbouring minima

                    // check previous bound
                    if (i != 0 && local_minima_bounds[i-1].max == local_minima_bounds[i].min) {
                        std::cout << "Merging bounds " << i-1 << " and " << i << std::endl;
                        // restore original bounds and merge
                        local_minima_bounds[i-1].min = original_bounds[i-1].min;
                        local_minima_bounds[i-1].max = original_bounds[i].max;
                        relax_bound(local_minima_bounds[i-1], local_minima[i-1]);
                        filtered_bounds.back() = local_minima_bounds[i-1];
                        filtered_minima.back() = std::min(local_minima[i-1], local_minima[i]);

                        // recalculate prominence
                        double prev_prominence = prominences[i-1];
                        prominences[i-1] = calc_prominence(i-1);
                        std::cout << "\tNew prominence is " << prominences[i-1] << std::endl;

                        // if the prominence is now large enough, we keep it
                        if (prev_prominence < min_prominence && min_prominence < prominences[i-1]) {
                            filtered_minima.push_back(local_minima[i-1]);
                            filtered_bounds.push_back(local_minima_bounds[i-1]);
                        }
                    }

                    // check next bound
                    if (i != size-1 && local_minima_bounds[i].max == local_minima_bounds[i+1].min) {
                        std::cout << "Merging bounds " << i+1 << " and " << i << std::endl;
                        local_minima_bounds[i+1].min = original_bounds[i].min;
                        local_minima_bounds[i+1].max = original_bounds[i+1].max;
                        relax_bound(local_minima_bounds[i+1], local_minima[i+1]);
                        local_minima[i+1] = std::min(local_minima[i], local_minima[i+1]);

                        // recalculate prominence
                        prominences[i+1] = calc_prominence(i+1);
                        std::cout << "\tNew prominence is " << prominences[i+1] << std::endl;
                    }
                } else {
                    filtered_minima.push_back(local_minima[i]);
                    filtered_bounds.push_back(local_minima_bounds[i]);
                }
            }

            local_minima = std::move(filtered_minima);
            local_minima_bounds = std::move(filtered_bounds);
        }

		#ifdef DEBUG_OUTPUT
			std::cout << "After filtering prominence " << local_minima.size() << " local minima." << std::endl;
			for (int i : local_minima) {
				std::cout << "\t(" << x[i] << ", " << y[i] << ")" << std::endl;
			}
		#endif
		#ifdef DEBUG_PLOT
            for (unsigned int i = 0; i < local_minima.size(); i++) {
                Limit& bound = local_minima_bounds[i];

                // bounds
                SimpleDataset dummy1({x[bound.min], x[bound.max]}, {y[bound.min], y[bound.max]});
                dummy1.add_plot_options({{"color", style::color::red}, {"ls", style::line::dashed}});
                plot.plot(dummy1);

                // prominence
                SimpleDataset dummy2({x[local_minima[i]], x[local_minima[i]]}, {y[local_minima[i]], y[local_minima[i]] + calc_prominence(i)});
                dummy2.add_plot_options({{"color", style::color::red}, {"ls", style::line::dashed}});
                plot.plot(dummy2);
            }
	        plot.save("temp/temp/local_minima.png");
		#endif

        //? TODO: perform a Gaussian fit of each minima to check quality?
        return local_minima;
    }
}