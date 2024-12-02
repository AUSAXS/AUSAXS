/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/ImageStack.h>
#include <settings/All.h>
#include <plots/All.h>
#include <mini/All.h>
#include <fitter/LinearFitter.h>
#include <fitter/SmartFitter.h>
#include <mini/detail/Parameter.h>
#include <em/detail/ExtendedLandscape.h>
#include <em/manager/ProteinManager.h>
#include <data/Molecule.h>
#include <utility/Console.h>
#include <utility/Limit.h>
#include <utility/Utility.h>
#include <constants/Constants.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <settings/EMSettings.h>
#include <settings/HistogramSettings.h>
#include <hydrate/generation/RadialHydration.h>
#include <math/Vector3.h>

#include <fstream>
#include <cassert>

using namespace ausaxs;
using namespace ausaxs::em;
using namespace ausaxs::fitter;

ImageStack::ImageStack(const io::ExistingFile& file) : ImageStackBase(file) {}

ImageStack::ImageStack(const std::vector<Image>& images) : ImageStackBase(images) {}

ImageStack::~ImageStack() = default;

double ImageStack::get_mass(double cutoff) const {
    auto p = get_protein_manager()->get_protein(cutoff);
    p->clear_grid();
    return p->get_excluded_volume_mass()/1e3;
}

// std::unique_ptr<EMFitResult> ImageStack::fit(std::unique_ptr<hist::ICompositeDistanceHistogram> h) {
//     Limit lim = {from_level(settings::em::alpha_levels.min), from_level(settings::em::alpha_levels.max)};
//     mini::Parameter param("cutoff", lim.center(), lim);
//     return fit(std::move(h), param);
// }

// std::unique_ptr<EMFitResult> ImageStack::fit(std::unique_ptr<hist::ICompositeDistanceHistogram> h, mini::Parameter& param) {
//     if (!param.has_bounds()) {return fit(std::move(h));} // ensure parameter bounds are present

//     SmartFitter fitter(get_data(), std::move(h));
//     std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(std::move(h), limit) : std::make_shared<LinearFitter>(std::move(h), limit);
//     return fit_helper(std::move(fitter), param);
// }

std::unique_ptr<EMFitResult> ImageStack::fit(const io::ExistingFile& file) {
    Limit lim = {from_level(settings::em::alpha_levels.min), from_level(settings::em::alpha_levels.max)};
    mini::Parameter param("cutoff", lim.center(), lim);
    return fit(file, param);
}

std::unique_ptr<EMFitResult> ImageStack::fit(const io::ExistingFile& file, mini::Parameter& param) {
    if (!param.has_bounds()) {return fit(file);} // ensure parameter bounds are present
    std::unique_ptr<SmartFitter> fitter = std::make_unique<SmartFitter>(file);
    return fit_helper(std::move(fitter), param);
}

std::shared_ptr<FitResult> last_fit; //? not the prettiest option, but it works for now
std::function<double(std::vector<double>)> ImageStack::prepare_function(std::shared_ptr<SmartFitter> _fitter) {
    // convert the calculated intensities to absolute scale
    // utility::print_warning("Warning in ImageStack::prepare_function: Not using absolute scale.");
    // auto protein = phm->get_protein(1);
    // double c = settings::em::concentration;                                // concentration
    // double m = protein->get_absolute_mass()*constants::unit::mg;          // mass
    // double DrhoV2 = std::pow(protein->get_relative_charge(), 2);          // charge
    // double re2 = pow(constants::radius::electron*constants::unit::cm, 2); // squared scattering length
    // double I0 = DrhoV2*re2*c/m;
    // fitter.normalize_intensity(I0);

    // stored vars for optimization
    static double last_c;
    static unsigned int counter;
    last_c = 5;
    counter = 0;

    // fitter is captured by value to guarantee its lifetime will be the same as the lambda
    // 'this' is ok since prepare_function is private and thus only used within the class itself
    hydrate::RadialHydration::set_noise_generator([] () {return Vector3<double>{0, 0, 0};}); // ensure hydration shell is deterministic
    return [this, fitter = std::move(_fitter)] (const std::vector<double>& params) -> double {
        auto p = get_protein_manager()->get_protein(params[0]);
        if (settings::em::hydrate) {
            p->clear_grid();                // clear grid from previous iteration
            p->generate_new_hydration();    // generate a new hydration layer

            // pointer cast is ok since the type should always be HydrationFitter when hydration is enabled
            fitter->set_guess({mini::Parameter{constants::fit::to_string(constants::fit::Parameters::SCALING_WATER), last_c, {0, 200}}});
            fitter->set_algorithm(mini::algorithm::SCAN);
            fitter->set_model(p->get_histogram());

            auto mass = p->get_excluded_volume_mass()/1e3;                                                // mass in kDa
            last_fit = fitter->fit();                                                                     // do the fit
            water_factors.push_back(last_fit->get_parameter(constants::fit::Parameters::SCALING_WATER));  // record c value
            last_c = last_fit->get_parameter(constants::fit::Parameters::SCALING_WATER).value;            // update c for next iteration
            evals.push_back(detail::ExtendedLandscape(params[0], mass, p->get_volume_grid(), std::move(last_fit->evaluated_points)));  // record evaluated points
        } else {
            p->clear_grid();                                    // clear grid from previous iteration
            auto mass = p->get_excluded_volume_mass()/1e3;      // mass in kDa
            fitter->set_model(p->get_histogram());
            last_fit = fitter->fit();
            evals.push_back(detail::ExtendedLandscape(params[0], mass, p->get_volume_grid(), std::move(last_fit->evaluated_points)));  // record evaluated points
        }

        double val = last_fit->fval;
        progress.notify(counter++);
        if (settings::fit::verbose) {
            std::cout << "Step " << utility::print_element(counter, 4) << ": Evaluated cutoff value " << utility::print_element(params[0], 8) << " with chi2 " << utility::print_element(val, 8) << std::flush << "\r";
        }
        return val;
    }; 
}

std::unique_ptr<EMFitResult> ImageStack::fit_helper(std::shared_ptr<SmartFitter> fitter, mini::Parameter& param) {
    //##########################################################//
    //###                       SETUP                        ###//
    //##########################################################//
    if (settings::em::plot_landscapes && settings::em::hydrate) {
        fitter->set_algorithm(mini::algorithm::SCAN);
    }

    update_charge_levels(*param.bounds);
    set_minimum_bounds(param.bounds->min);
    auto func = prepare_function(fitter);
    mini::Landscape evals; // since we'll be using multiple minimizers, we'll need to store the evaluated points manually
    unsigned int dof = fitter->dof()-1; // minus one because we're also fitting the cutoff

    if (settings::general::verbose) {
        std::cout << "The mass range [" << std::left << std::setw(8) << get_mass(param.bounds->min) 
                                << ", " << std::left << std::setw(8) << get_mass(param.bounds->max) << "] kDa will be scanned." << std::endl;
    }

    EMFitResult::EMFitInfo plots;

    //##########################################################//
    //###                DETERMINE LANDSCAPE                 ###//
    //##########################################################//
    mini::LimitedScan minimizer(func, param, settings::fit::max_iterations);
    minimizer.set_limit(5, true);
    SimpleDataset chi2_data;
    {
        auto l = minimizer.landscape(settings::fit::max_iterations);
        evals.append(l);
        chi2_data = l.as_dataset();
    }

    chi2_data.sort_x();
    auto min_abs = chi2_data.find_minimum();
    std::cout << "Minimum at " << min_abs.x << " with chi2 " << min_abs.y << std::endl;

    //##########################################################//
    //### CHECK LANDSCAPE IS OK FOR AVERAGING & INTERPLATION ###//
    //##########################################################//
    chi2_data.limit_y(0, min_abs.y*5);  // focus on the area near the absolute minimum
    if (chi2_data.size_rows() < 10) {       // if we have too few points after imposing the limit, we must sample some more
        Limit bounds;                       // first we determine the bounds of the area we want to sample
        if (chi2_data.size_rows() < 3) {    // if we only have one or two points, sample the area between the neighbouring points
            double s = (param.bounds->max - param.bounds->min)/settings::fit::max_iterations;
            bounds = {min_abs.x - s, min_abs.x + s};
        }
        else { // otherwise just use the new bounds of the limited landscape
            bounds = chi2_data.span_x();
        }

        // prepare a new minimizer with the new bounds
        console::print_warning("Function is varying strongly. Sampling more points around the minimum.");
        mini::LimitedScan mini2(func, mini::Parameter("cutoff", bounds), settings::fit::max_iterations/4);
        {
            auto l = mini2.landscape(settings::fit::max_iterations/2);
            evals.append(l);
            chi2_data = l.as_dataset();
        }
        chi2_data.sort_x();
        min_abs = chi2_data.find_minimum();
        std::cout << "New minimum at " << min_abs.x << " with chi2 " << min_abs.y << std::endl;
        chi2_data.limit_y(0, min_abs.y*5);

        if (chi2_data.size_rows() < 10) {
            throw except::unexpected("ImageStack::fit: Could not sample enough points around the minimum. Function varies too much.");
        }
    }

    //##########################################################//
    //###         AVERAGE & INTERPLATE MORE POINTS           ###//
    //##########################################################//
    Dataset data_avg_int; // cutoff, chi2, mass
    {
        auto ra = chi2_data.rolling_average(7).interpolate(5); // impose a moving average filter
        data_avg_int = Dataset(ra.size(), 3);
        data_avg_int.set_col_names({"cutoff", "chi2", "mass"});
        data_avg_int.col("cutoff") = ra.x();
        data_avg_int.col("chi2") = ra.y();
    }

    double spacing = data_avg_int.x(1)-data_avg_int.x(0); 
    auto minima = data_avg_int.find_minima(static_cast<int>(0.1*data_avg_int.size()), 0.1); // find all minima. they should be fairly spaced out (10% seems reasonable?)
    {   // find the absolute minimum in the smoothed landscape
        auto tmp = data_avg_int.find_minimum(1);
        if (tmp[1] < min_abs.y) {
            std::cout << "Warning: Absolute minimum in the smoothed landscape is different from the original minimum. Using the smoothed minimum." << std::endl;
            std::cout << "\tOriginal: " << min_abs.x << " " << min_abs.y << std::endl;
            std::cout << "\tSmoothed: " << tmp[0] << " " << tmp[1] << std::endl;
            min_abs = {tmp[0], tmp[1], tmp[2]};
        }
    }

    // prepare the mass axis
    if (settings::em::mass_axis) {
        Dataset mass_data(this->evals.size(), 2);
        for (unsigned int i = 0; i < this->evals.size(); ++i) {
            mass_data[i] = {this->evals[i].cutoff, this->evals[i].mass};
        }
        mass_data.sort_x();
        data_avg_int.col("mass") = mass_data.interpolate(data_avg_int.x()).y();
    }

    { // remove minima that are too far away from the absolute minimum
        std::vector<unsigned int> to_keep;
        for (auto m : minima) {
            if (data_avg_int.y(m) < min_abs.y*2) {to_keep.push_back(m);}
        }
        minima = std::move(to_keep);
    }
    
    // update our parameter since we interpolated more points
    param.guess = min_abs.x;
    param.bounds = Limit(min_abs.x-3*spacing, min_abs.x+3*spacing); // uncertainty is 3*spacing between points

    // save .pdb structures of the other minima
    if (settings::em::save_pdb && 1 < minima.size()) {
        unsigned int enumerate = 0;
        std::string info;
        for (auto m : minima) {
            if (data_avg_int.x(m) == min_abs.x) {continue;}
            auto temp_protein = get_protein_manager()->get_protein(data_avg_int.x(m));
            if (settings::em::hydrate) {
                temp_protein->clear_grid();
                temp_protein->generate_new_hydration();
            }
            temp_protein->save(settings::general::output + "models/model_" + std::to_string(++enumerate) + ".pdb");
            info += "Model " + std::to_string(enumerate) + ": (σ, χ²) = " + std::to_string(to_level(data_avg_int.x(m))) + " " + std::to_string(data_avg_int.y(m)) + "\n";
            if (settings::em::mass_axis) {
                info += "  Estimated mass = " + std::to_string(data_avg_int.col("mass")[m]) + " kDa\n";
            }
        }
        std::ofstream out(settings::general::output + "models/info.txt");
        out << info;
    }

    if (settings::general::generate_plots) {
        { // make a nice plot of the landscape within some range of the minimum; this is often nicer to look at than the full landscape due to the reduced y-range
            // plot the minimum in blue
            SimpleDataset p_min, chi2_copy = chi2_data, avg_copy = data_avg_int;
            for (auto m : minima) {
                // if the minimum is too close to the absolute minimum & the absolute minimum is lower, plot the absolute minimum instead
                if (std::abs(data_avg_int.x(m) - min_abs.x) < from_level(0.5) && min_abs.y < data_avg_int.y(m)) {
                    p_min.push_back(to_level(min_abs.x), min_abs.y/dof);
                    continue;
                }
                p_min.push_back(to_level(data_avg_int.x(m)), data_avg_int.y(m)/dof);
            }

            // convert cutoff to std levels & normalize chi2
            for (unsigned int i = 0; i < chi2_copy.size(); ++i) {
                chi2_copy.x(i) = to_level(chi2_copy.x(i));
                chi2_copy.y(i) /= dof;
            }
            for (unsigned int i = 0; i < avg_copy.size(); ++i) {
                avg_copy.x(i) = to_level(avg_copy.x(i));
                avg_copy.y(i) /= dof;
            }
            plots.chi2_limited = chi2_copy;

            // prepare rest of the plot
            plots::PlotDataset plot(avg_copy, plots::PlotOptions(style::draw::line, {{"color", style::color::red}, {"xlabel", "cutoff [$\\sigma$]"}, {"ylabel", "$\\chi_r^2$"}}));
            plot.plot(chi2_copy, plots::PlotOptions(style::draw::points, {}));
            plot.plot(p_min, plots::PlotOptions(style::draw::points, {{"color", style::color::blue}, {"s", 12}}));
            plot.save(settings::general::output + "chi2_evaluated_points_limited." + settings::plots::format);

            if (settings::em::mass_axis) {
                // create chi2 / mass dataset
                SimpleDataset mass_avg_copy(data_avg_int.col("mass"), data_avg_int.col("chi2")/dof);

                SimpleDataset mass_p_min;
                for (auto m : minima) {
                    // if the minimum is too close to the absolute minimum & the absolute minimum is lower, plot the absolute minimum instead
                    if (std::abs(data_avg_int.x(m) - min_abs.x) < from_level(0.5) && min_abs.y < data_avg_int.y(m)) {
                        mass_p_min.push_back(data_avg_int.interpolate_x(min_abs.x, 2), min_abs.y/dof);
                        continue;
                    }
                    mass_p_min.push_back(data_avg_int.col("mass")[m], data_avg_int.col("chi2")[m]/dof);
                }
                plots.mass_limited = mass_avg_copy;

                // make the plot
                plots::PlotDataset plot_mass(mass_avg_copy, plots::PlotOptions(style::draw::line, {{"color", style::color::red}, {"xlabel", "mass [kDa]"}, {"ylabel", "$\\chi_r^2$"}}));
                plot_mass.plot(SimpleDataset(data_avg_int.interpolate(chi2_data.x()).col(2), chi2_copy.y()), plots::PlotOptions(style::draw::points, {}));
                plot_mass.plot(mass_p_min, plots::PlotOptions(style::draw::points, {{"color", style::color::blue}, {"s", 12}}));
                plot_mass.save(settings::general::output + "chi2_evaluated_points_limited_mass." + settings::plots::format);
            }

            if (settings::em::hydrate && settings::general::supplementary_plots) {
                plots.water_factors = get_fitted_water_factors_dataset();
                plots::PlotDataset::quick_plot(
                    plots.water_factors,
                    plots::PlotOptions(style::draw::points, {{"xlabel", "Iteration"}, {"ylabel", "Scaling factor"}}),
                    settings::general::output + "water_factors." + settings::plots::format
                );
            }
        }

        { // plot all evaluated points
            { // chi2 landscape
                auto l = evals.as_dataset();
                l.sort_x();
                for (unsigned int i = 0; i < l.size(); ++i) {
                    l.x(i) = to_level(l.x(i));
                    l.y(i) /= dof;
                }
                plots.chi2_full = l;

                plots::PlotDataset::quick_plot(
                    l, 
                    plots::PlotOptions(style::draw::points, {{"xlabel", "cutoff [$\\sigma$]"}, {"ylabel", "$\\chi_r^2$"}}), 
                    settings::general::output + "chi2_evaluated_points_full." + settings::plots::format
                );
            }

            // volume as a function of cutoff
            if (settings::general::supplementary_plots) {
                SimpleDataset volume_data(this->evals.size()); 
                for (unsigned int i = 0; i < this->evals.size(); ++i) {
                    volume_data.x(i) = to_level(this->evals[i].cutoff);
                    volume_data.y(i) = this->evals[i].mass;
                }
                volume_data.sort_x();
                plots.volume = volume_data;

                plots::PlotDataset::quick_plot(
                    volume_data, 
                    plots::PlotOptions(style::draw::points, {{"xlabel", "cutoff [$\\sigma$]"}, {"ylabel", "Volume [Å³]"}, {"title", "Volume as a function of cutoff"}}), 
                    settings::general::output + "volume." + settings::plots::format
                );
            }

            // plot with mass axis
            // if (settings::em::mass_axis) {
            //     plots::PlotDataset::quick_plot(
            //         mass_data,
            //         plots::PlotOptions(style::draw::points, {{"xlabel", "mass [kDa]"}, {"ylabel", "$\\chi_r^2$"}}),
            //         settings::general::output + "chi2_evaluated_points_full_mass." + settings::plots::format
            //     );
            // }
        }
    }

    //##########################################################//
    //###             EXPLORE AREA AROUND MINIMUM            ###//
    //##########################################################//
    // if hydration is enabled, the chi2 will oscillate heavily around the minimum
    // we therefore want to sample the area near the minimum to get an average
    mini::Result res;
    update_charge_levels({from_level(to_level(min_abs.x)-0.5), from_level(to_level(min_abs.x)+0.5)});
    if (settings::em::hydrate) {
        // reset evaluated points
        this->evals.clear();

        // sample the area around the minimum
        mini::MinimumExplorer explorer(func, param, settings::fit::max_iterations);
        res = explorer.minimize();

        // check if we found a better absolute minima
        auto explored_points = explorer.get_evaluated_points().as_dataset();
        if (auto new_min = explored_points.find_minimum(); new_min.y < min_abs.y) {
            min_abs = new_min;
        }

        // plot evaluated points near the minimum
        if (settings::general::generate_plots) {
            explored_points.y() = explored_points.y()/dof;

            // calculate the mean & standard deviation of the sampled points
            double mu = explored_points.mean();
            double sigma = explored_points.std();

            // plot the starting point in blue
            SimpleDataset p_start;
            p_start.push_back(min_abs.x, min_abs.y/dof);
            plots.chi2_minimum = explored_points;

            // do the actual plotting
            plots::PlotDataset(explored_points, plots::PlotOptions(style::draw::points, {{"xlabel", "cutoff [$\\sigma$]"}, {"ylabel", "$\\chi_r^2$"}}))
                .hline(mu, plots::PlotOptions(style::draw::line, {{"color", style::color::red}}))
                .hline(mu+sigma, plots::PlotOptions(style::draw::line, {{"color", style::color::red}, {"linestyle", "--"}}))
                .hline(mu-sigma, plots::PlotOptions(style::draw::line, {{"color", style::color::red}, {"linestyle", "--"}}))
                .plot(p_start, plots::PlotOptions(style::draw::points, {{"color", style::color::blue}, {"s", 9}}))
            .save(settings::general::output + "chi2_near_minimum." + settings::plots::format);

            // make mass version
            if (settings::em::mass_axis) {
                Dataset mass_cutoff(0, 2);
                // this->evals also records the masses, so the last explored_points.size() entries are the ones we want
                for (int i = static_cast<int>(this->evals.size() - explored_points.size()); i < static_cast<int>(this->evals.size()); ++i) {
                    mass_cutoff.push_back({this->evals[i].cutoff, this->evals[i].mass});
                }
                mass_cutoff.sort_x();

                // create chi2 / mass dataset
                explored_points.x() = mass_cutoff.y();

                // interpolate start point
                p_start.x(0) = mass_cutoff.interpolate_x(p_start.x(0), 1);

                plots.mass_minimum = explored_points;

                // make the plot
                plots::PlotDataset(explored_points, plots::PlotOptions(style::draw::points, {{"xlabel", "mass [kDa]"}, {"ylabel", "$\\chi_r^2$"}}))
                    .hline(mu, plots::PlotOptions(style::draw::line, {{"color", style::color::red}}))
                    .hline(mu+sigma, plots::PlotOptions(style::draw::line, {{"color", style::color::red}, {"linestyle", "--"}}))
                    .hline(mu-sigma, plots::PlotOptions(style::draw::line, {{"color", style::color::red}, {"linestyle", "--"}}))
                    .plot(p_start, plots::PlotOptions(style::draw::points, {{"color", style::color::blue}, {"s", 9}}))
                .save(settings::general::output + "chi2_near_minimum_mass." + settings::plots::format);
            }
        }
    } 
    
    // otherwise do a quick fit to ensure we're at the very bottom of the valley
    else {
        mini::Golden golden(func, param);
        res = golden.minimize();
        if (res.fval < min_abs.y) {
            min_abs = golden.get_evaluated_points().as_dataset().find_minimum();
        }
    }

    // Make 3D landscape plot
    if (settings::em::plot_landscapes && settings::em::hydrate) {
        mini::Landscape l;
        l.evals.reserve(1000);
        for (int i = 0; i < static_cast<int>(this->evals.size()); ++i) {
            for (int j = 0; j < static_cast<int>(this->evals[i].strip.evals.size()); ++j) {
                double x = this->evals[i].cutoff;
                double y = this->evals[i].strip.evals[j].vals.front();
                double z = this->evals[i].strip.evals[j].fval;
                l.evals.push_back(mini::Evaluation({x, y}, z));
            }
        }
        plots::PlotLandscape::quick_plot(l, plots::PlotOptions({{"xlabel", "cutoff"}, {"ylabel", "c"}, {"zlabel", "$\\chi^2$"}}), settings::general::output + "chi2_data." + settings::plots::format);
    }

    // update the fitter with the optimal cutoff, such that the returned fit is actually the best one
    double fval = func({min_abs.x});
    assert(std::abs(fval - min_abs.y) < 1e-6 && "ImageStack::fit: The minimum found by the minimizer does not match the minimum found in the dataset.");

    std::unique_ptr<fitter::EMFitResult> emfit = std::make_unique<EMFitResult>(res, fval, dof+3); // +3 because they'll be subtracted again by the add_fit call
    {
        auto data = fitter->get_data();
        emfit->set_data_curves(
            data.x(), 
            data.y(), 
            data.yerr(), 
            fitter->get_model_curve({last_fit->get_parameter(constants::fit::Parameters::SCALING_WATER)}), 
            fitter->get_residuals({last_fit->get_parameter(constants::fit::Parameters::SCALING_WATER)})
        );
    }
    emfit->add_fit(last_fit.get(), true);
    emfit->fevals = evals.evals.size();
    emfit->em_info = std::move(plots);
    emfit->evaluated_points = std::move(evals);
    emfit->level = to_level(min_abs.x);
    if (settings::em::mass_axis) {emfit->mass = data_avg_int.interpolate_x(min_abs.x, 2);}
    if (settings::em::save_pdb) {
        auto temp_protein = get_protein_manager()->get_protein(min_abs.x);
        if (settings::em::hydrate) {
            temp_protein->clear_grid();
            temp_protein->generate_new_hydration();
        }
        temp_protein->save(settings::general::output + "model.pdb");}
    return emfit;
}

// mini::Landscape ImageStack::cutoff_scan(const Axis& points, const io::ExistingFile& file) {
//     std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(file) : std::make_shared<LinearFitter>(file);
//     return cutoff_scan_helper(points, std::move(fitter));
// }

// mini::Landscape ImageStack::cutoff_scan(unsigned int points, const io::ExistingFile& file) {
//     Axis axis(from_level(settings::em::alpha_levels.min), from_level(settings::em::alpha_levels.max), points);
//     return cutoff_scan(axis, file);
// }

// mini::Landscape ImageStack::cutoff_scan(const Axis& points, std::unique_ptr<hist::ICompositeDistanceHistogram> h) {
//     auto limit = Limit(settings::axes::qmin, settings::axes::qmax);
//     std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(std::move(h), limit) : std::make_shared<LinearFitter>(std::move(h), limit);
//     return cutoff_scan_helper(points, std::move(fitter));
// }

// mini::Landscape ImageStack::cutoff_scan(unsigned int points, std::unique_ptr<hist::ICompositeDistanceHistogram> h) {
//     Axis axis(from_level(settings::em::alpha_levels.min), from_level(settings::em::alpha_levels.max), points);
//     return cutoff_scan(axis, std::move(h));
// }

// std::pair<EMFitResult, mini::Landscape> ImageStack::cutoff_scan_fit(unsigned int points, std::unique_ptr<hist::ICompositeDistanceHistogram> h) {
//     Axis axis(from_level(settings::em::alpha_levels.min), from_level(settings::em::alpha_levels.max), points);
//     return cutoff_scan_fit(axis, std::move(h));
// }

// std::pair<EMFitResult, mini::Landscape> ImageStack::cutoff_scan_fit(const Axis& points, const io::ExistingFile& file) {
//     std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(file) : std::make_shared<LinearFitter>(file);    
//     return cutoff_scan_fit_helper(points, std::move(fitter));
// }

// std::pair<EMFitResult, mini::Landscape> ImageStack::cutoff_scan_fit(unsigned int points, const io::ExistingFile& file) {
//     Axis axis(from_level(settings::em::alpha_levels.min), from_level(settings::em::alpha_levels.max), points);
//     std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(file) : std::make_shared<LinearFitter>(file);    
//     return cutoff_scan_fit_helper(axis, std::move(fitter));
// }

// mini::Landscape ImageStack::cutoff_scan_helper(const Axis& points, std::shared_ptr<LinearFitter> fitter) {
//     update_charge_levels(points.limits());
//     set_minimum_bounds(points.min);
//     auto func = prepare_function(std::move(fitter));

//     mini::Golden minimizer(std::move(func), mini::Parameter{"cutoff", points.limits()});
//     return minimizer.landscape(points.bins);
// }

// std::pair<EMFitResult, mini::Landscape> ImageStack::cutoff_scan_fit(const Axis& points, std::unique_ptr<hist::ICompositeDistanceHistogram> h) {
//     auto limit = Limit(settings::axes::qmin, settings::axes::qmax);
//     std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(std::move(h), limit) : std::make_shared<LinearFitter>(std::move(h), limit);
//     return cutoff_scan_fit_helper(points, std::move(fitter));
// }

// std::pair<EMFitResult, mini::Landscape> ImageStack::cutoff_scan_fit_helper(const Axis& points, std::shared_ptr<LinearFitter> fitter) {
//     update_charge_levels(points.limits());
//     set_minimum_bounds(points.min);
//     auto func = prepare_function(fitter);

//     // cutoff scan
//     mini::Golden minimizer(std::move(func), mini::Parameter{"cutoff", points.limits()});
//     mini::Landscape landscape = minimizer.landscape(points.bins);

//     // fit
//     double l = from_level(1);
//     Limit limit(settings::em::alpha_levels.min*l, settings::em::alpha_levels.max*l);
//     minimizer.clear_parameters();
//     minimizer.add_parameter({"cutoff", limit.center(), limit});
//     auto res = minimizer.minimize();

//     EMFitResult emfit(fitter.get(), res, res.fval);
//     emfit.evaluated_points = minimizer.get_evaluated_points();

//     return {emfit, landscape};
// }

const std::vector<mini::FittedParameter>& ImageStack::get_fitted_water_factors() const {
    return water_factors;
}

SimpleDataset ImageStack::get_fitted_water_factors_dataset() const {
    std::vector<double> x(water_factors.size()), y(water_factors.size());
    for (unsigned int i = 0; i < water_factors.size(); i++) {
        x[i] = i;
        y[i] = water_factors[i].value;
    }
    return SimpleDataset(x, y, "Iteration", "Scaling factor");
}

void ImageStack::update_charge_levels(const Limit& limit) const noexcept {
    std::vector<double> levels;
    for (unsigned int i = 0; i < settings::em::charge_levels; i++) {
        levels.push_back(limit.min + i*limit.span()/settings::em::charge_levels);
    }
    get_protein_manager()->set_charge_levels(levels);
}