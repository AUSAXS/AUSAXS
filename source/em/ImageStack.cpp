#include <em/ImageStack.h>
#include <settings/EMSettings.h>
#include <settings/PlotSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/ProteinSettings.h>
#include <settings/FitSettings.h>
#include <utility/Console.h>
#include <plots/all.h>
#include <fitter/LinearFitter.h>
#include <fitter/HydrationFitter.h>
#include <mini/all.h>
#include <em/detail/ExtendedLandscape.h>
#include <em/manager/ProteinManager.h>
#include <data/Protein.h>

using namespace em;
using namespace fitter;

std::shared_ptr<EMFit> ImageStack::fit(const hist::ScatteringHistogram& h) {
    Limit lim = {from_level(settings::em::alpha_levels.min), from_level(settings::em::alpha_levels.max)};
    mini::Parameter param("cutoff", lim.center(), lim);
    return fit(h, param);
}

std::shared_ptr<EMFit> ImageStack::fit(const hist::ScatteringHistogram& h, mini::Parameter param) {
    if (!param.has_bounds()) {return fit(h);} // ensure parameter bounds are present

    std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(h, get_limits()) : std::make_shared<LinearFitter>(h, get_limits());
    return fit_helper(fitter, param);
}

std::shared_ptr<EMFit> ImageStack::fit(const io::ExistingFile& file) {
    Limit lim = {from_level(settings::em::alpha_levels.min), from_level(settings::em::alpha_levels.max)};
    mini::Parameter param("cutoff", lim.center(), lim);
    return fit(file, param);
}

std::shared_ptr<EMFit> ImageStack::fit(const io::ExistingFile& file, mini::Parameter param) {
    if (!param.has_bounds()) {return fit(file);} // ensure parameter bounds are present
    std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(file) : std::make_shared<LinearFitter>(file);
    return fit_helper(fitter, param);
}

std::shared_ptr<EMFit> ImageStack::fit_helper(std::shared_ptr<LinearFitter> fitter, mini::Parameter param) {
    update_charge_levels(*param.bounds);
    set_minimum_bounds(param.bounds->min);
    auto f = prepare_function(fitter);
    mini::Landscape evals; // since we'll be using multiple minimizers, we'll need to store the evaluated points manually

    //##########################################################//
    //###                DETERMINE LANDSCAPE                 ###//
    //##########################################################//
    mini::LimitedScan minimizer(f, param, settings::fit::max_iterations);
    minimizer.set_limit(fitter->dof()*20);
    SimpleDataset chi2_landscape;
    {
        auto l = minimizer.landscape(settings::fit::max_iterations);
        evals.append(l);
        chi2_landscape = l.as_dataset();
    }

    chi2_landscape.sort_x();
    auto min = chi2_landscape.find_minimum();

    //##########################################################//
    //### CHECK LANDSCAPE IS OK FOR AVERAGING & INTERPLATION ###//
    //##########################################################//
    chi2_landscape.limit_y(0, min.y*5);  // focus on the area near the minimum
    if (chi2_landscape.size() < 10) {    // if we have too few points after imposing the limit, we must sample some more
        Limit bounds;                    // first we determine the bounds of the area we want to sample
        if (chi2_landscape.size() < 3) { // if we only have one or two points, sample the area between the neighbouring points
            double s = (param.bounds->max - param.bounds->min)/settings::fit::max_iterations;
            bounds = {min.x - s, min.x + s};
        }
        else { // otherwise just use the new bounds of the limited landscape
            bounds = chi2_landscape.span_x();
        }

        // prepare a new minimizer with the new bounds
        console::print_warning("Function is varying strongly. Sampling more points around the minimum.");
        mini::LimitedScan mini2(f, mini::Parameter("cutoff", bounds), settings::fit::max_iterations/4);
        {
            auto l = mini2.landscape(settings::fit::max_iterations/2);
            evals.append(l);
            chi2_landscape = l.as_dataset();
        }
        chi2_landscape.sort_x();
        min = chi2_landscape.find_minimum();
        chi2_landscape.limit_y(0, min.y*5);

        if (chi2_landscape.size() < 10) {
            throw except::unexpected("ImageStack::fit: Could not sample enough points around the minimum. Function varies too much.");
        }
    }

    //##########################################################//
    //###         AVERAGE & INTERPLATE MORE POINTS           ###//
    //##########################################################//
    SimpleDataset avg = chi2_landscape.rolling_average(7);  // impose a moving average filter 
    avg.interpolate(5);                                     // interpolate more points

    min = avg.find_minimum();
    double spacing = avg.x(1)-avg.x(0); 
    param.guess = min.x;
    param.bounds = Limit(min.x-3*spacing, min.x+3*spacing); // uncertainty is 3*spacing between points

    // Plot all evaluated points
    if (settings::general::supplementary_plots) {
        { // Plot evaluated points around the minimum
            // plot the minimum in blue
            SimpleDataset p_min, chi2_copy = chi2_landscape;
            p_min.push_back(min.x, min.y);
            p_min.add_plot_options(style::draw::points, {{"color", style::color::blue}, {"s", 9}});

            // prepare rest of the plot
            avg.add_plot_options(style::draw::line, {{"color", style::color::red}, {"xlabel", "cutoff"}, {"ylabel", "$\\chi^2$"}});
            plots::PlotDataset plot(avg);
            chi2_copy.add_plot_options(style::draw::points);
            plot.plot(chi2_copy);
            plot.plot(p_min);
            plot.save(settings::general::output + "chi2_evaluated_points_limited." + settings::plots::format);
        }

        { // Plot all evaluated points
            auto l = evals.as_dataset();
            l.sort_x();
            l.add_plot_options(style::draw::points, {{"xlabel", "cutoff"}, {"ylabel", "$\\chi^2$"}});
            plots::PlotDataset::quick_plot(l, settings::general::output + "chi2_evaluated_points_full." + settings::plots::format);
        }
    }

    //##########################################################//
    //###             EXPLORE AREA AROUND MINIMUM            ###//
    //##########################################################//
    // if hydration is enabled, the chi2 will oscillate heavily around the minimum
    // we therefore want to sample the area near the minimum to get an average
    mini::Result res;
    if (settings::em::hydrate) {
        // sample the area around the minimum
        mini::MinimumExplorer explorer(f, param, settings::fit::max_iterations);
        res = explorer.minimize();

        SimpleDataset area;
        {
            auto l = explorer.landscape();
            evals.append(l);
            area = l.as_dataset();
        }

        // Plot evaluated points near the minimum
        if (settings::general::supplementary_plots) {
            // calculate the mean & standard deviation of the sampled points
            double mu = area.mean();
            double sigma = area.std();

            // plot horizontal lines at the mean and mean +/- sigma
            auto xspan = area.span_x();
            SimpleDataset l({xspan.min, xspan.max}, {mu, mu});
            SimpleDataset lp({xspan.min, xspan.max}, {mu+sigma, mu+sigma});
            SimpleDataset lm({xspan.min, xspan.max}, {mu-sigma, mu-sigma});
            l.add_plot_options(style::draw::line, {{"color", style::color::red}});
            lp.add_plot_options(style::draw::line, {{"color", style::color::red}, {"linestyle", "--"}});
            lm.add_plot_options(style::draw::line, {{"color", style::color::red}, {"linestyle", "--"}});

            // plot the starting point in blue
            SimpleDataset p_start;
            p_start.push_back(min.x, min.y);
            p_start.add_plot_options(style::draw::points, {{"color", style::color::blue}, {"s", 9}});

            // do the actual plotting
            area.add_plot_options(style::draw::points, {{"xlabel", "cutoff"}, {"ylabel", "$\\chi^2$"}});
            plots::PlotDataset plot(area);
            plot.plot(l);
            plot.plot(lm);
            plot.plot(lp);
            plot.plot(p_start);
            plot.save(settings::general::output + "chi2_near_minimum." + settings::plots::format);
        }
    } 
    
    // otherwise do a quick fit to ensure we're at the very bottom of the valley
    else {
        mini::Golden golden(f, param);
        res = golden.minimize();
        evals.append(golden.get_evaluated_points());
    }

    // Make 3D landscape plot
    if (settings::general::supplementary_plots && settings::em::hydrate) {
        mini::Landscape l;
        l.evals.reserve(1000);
        for (unsigned int i = 0; i < this->evals.size(); i++) {
            for (unsigned int j = 0; j < this->evals[i].strip.evals.size(); j++) {
                double x = this->evals[i].cutoff;
                double y = this->evals[i].strip.evals[j].vals.front();
                double z = this->evals[i].strip.evals[j].fval;
                l.evals.push_back(mini::Evaluation({x, y}, z));
            }
        }
        l.add_plot_options({{"xlabel", "cutoff"}, {"ylabel", "c"}, {"zlabel", "$\\chi^2$"}});
        plots::PlotLandscape::quick_plot(l, settings::general::output + "chi2_landscape." + settings::plots::format);
    }

    // update the fitter with the optimal cutoff, such that the returned fit is actually the best one
    min = evals.as_dataset().find_minimum();
    f({min.x});

    std::shared_ptr<fitter::EMFit> emfit = std::make_shared<EMFit>(*fitter, res, res.fval);
    emfit->evaluated_points = evals;
    emfit->fevals = evals.evals.size();
    emfit->level = to_level(min.x);
    if (settings::em::save_pdb) {get_protein_manager()->get_protein()->save(settings::general::output + "model.pdb");}
    return emfit;
}

std::function<double(std::vector<double>)> ImageStack::prepare_function(std::shared_ptr<LinearFitter> fitter) {
    // convert the calculated intensities to absolute scale
    // utility::print_warning("Warning in ImageStack::prepare_function: Not using absolute scale.");
    // auto protein = phm->get_protein(1);
    // double c = settings::em::concentration;                                // concentration
    // double m = protein->get_absolute_mass()*constants::unit::mg;          // mass
    // double DrhoV2 = std::pow(protein->get_relative_charge(), 2);          // charge
    // double re2 = pow(constants::radius::electron*constants::unit::cm, 2); // squared scattering length
    // double I0 = DrhoV2*re2*c/m;
    // fitter.normalize_intensity(I0);

    // fit function
    settings::protein::center = false;   // do not center the protein - this may cause issues
    if (settings::em::plot_landscapes && settings::em::hydrate) {
        std::static_pointer_cast<HydrationFitter>(fitter)->set_algorithm(mini::type::SCAN); 
    }
    
     // fitter is captured by value to guarantee its lifetime will be the same as the lambda
     // 'this' is ok since prepare_function is private and thus only used within the class itself
    std::function<double(std::vector<double>)> chi2 = [this, fitter] (std::vector<double> params) {
        static unsigned int counter = 0;
        static double last_c = 5;
        auto p = get_protein_manager()->get_protein(params[0]);

        std::shared_ptr<Fit> fit;
        if (settings::em::hydrate) {
            p->clear_grid();                // clear grid from previous iteration
            p->generate_new_hydration();    // generate a new hydration layer

            // pointer cast is ok since the type should always be HydrationFitter when hydration is enabled
            std::static_pointer_cast<HydrationFitter>(fitter)->set_guess(mini::Parameter{"c", last_c, {0, 200}}); 
            fitter->set_scattering_hist(p->get_histogram());

            fit = fitter->fit();                                                            // do the fit
            water_factors.push_back(fit->get_parameter("c"));                               // record c value
            last_c = fit->get_parameter("c").value;                                         // update c for next iteration
            evals.push_back(detail::ExtendedLandscape(params[0], fit->evaluated_points));   // record evaluated points
        } else {
            fitter->set_scattering_hist(p->get_histogram());
            fit = fitter->fit();
        }

        double val = fit->fval;
        if (settings::fit::verbose) {
            std::cout << "Step " << utility::print_element(counter++, 5) << ": Evaluated cutoff value " << utility::print_element(params[0], 12) << " with chi2 " << utility::print_element(val, 12) << std::flush << "\r";
        }
        return val;
    }; 
    return chi2;
}

mini::Landscape ImageStack::cutoff_scan(const Axis& points, const io::ExistingFile& file) {
    std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(file) : std::make_shared<LinearFitter>(file);
    return cutoff_scan_helper(points, fitter);
}

mini::Landscape ImageStack::cutoff_scan(unsigned int points, const io::ExistingFile& file) {
    Axis axis(points, from_level(settings::em::alpha_levels.min), from_level(settings::em::alpha_levels.max));
    return cutoff_scan(axis, file);
}

mini::Landscape ImageStack::cutoff_scan(const Axis& points, const hist::ScatteringHistogram& h) {
    std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(h, get_limits()) : std::make_shared<LinearFitter>(h, get_limits());
    return cutoff_scan_helper(points, fitter);
}

mini::Landscape ImageStack::cutoff_scan(unsigned int points, const hist::ScatteringHistogram& h) {
    Axis axis(points, from_level(settings::em::alpha_levels.min), from_level(settings::em::alpha_levels.max));
    return cutoff_scan(axis, h);
}

std::pair<EMFit, mini::Landscape> ImageStack::cutoff_scan_fit(unsigned int points, const hist::ScatteringHistogram& h) {
    Axis axis(points, from_level(settings::em::alpha_levels.min), from_level(settings::em::alpha_levels.max));
    return cutoff_scan_fit(axis, h);
}

std::pair<EMFit, mini::Landscape> ImageStack::cutoff_scan_fit(const Axis& points, std::string file) {
    std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(file) : std::make_shared<LinearFitter>(file);    
    return cutoff_scan_fit_helper(points, fitter);
}

std::pair<EMFit, mini::Landscape> ImageStack::cutoff_scan_fit(unsigned int points, std::string file) {
    Axis axis(points, from_level(settings::em::alpha_levels.min), from_level(settings::em::alpha_levels.max));
    std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(file) : std::make_shared<LinearFitter>(file);    
    return cutoff_scan_fit_helper(axis, fitter);
}

mini::Landscape ImageStack::cutoff_scan_helper(const Axis& points, std::shared_ptr<LinearFitter> fitter) {
    update_charge_levels(points.limits());
    set_minimum_bounds(points.min);
    auto func = prepare_function(fitter);

    mini::Golden minimizer(func, mini::Parameter{"cutoff", points.limits()});
    return minimizer.landscape(points.bins);
}

std::pair<EMFit, mini::Landscape> ImageStack::cutoff_scan_fit(const Axis& points, const hist::ScatteringHistogram& h) {
    std::shared_ptr<LinearFitter> fitter = settings::em::hydrate ? std::make_shared<HydrationFitter>(h, get_limits()) : std::make_shared<LinearFitter>(h, get_limits());
    return cutoff_scan_fit_helper(points, fitter);
}

std::pair<EMFit, mini::Landscape> ImageStack::cutoff_scan_fit_helper(const Axis& points, std::shared_ptr<LinearFitter> fitter) {
    update_charge_levels(points.limits());
    set_minimum_bounds(points.min);
    auto func = prepare_function(fitter);

    // cutoff scan
    mini::Golden minimizer(func, mini::Parameter{"cutoff", points.limits()});
    mini::Landscape landscape = minimizer.landscape(points.bins);

    // fit
    double l = from_level(1);
    Limit limit(settings::em::alpha_levels.min*l, settings::em::alpha_levels.max*l);
    minimizer.clear_parameters();
    minimizer.add_parameter({"cutoff", limit.center(), limit});
    auto res = minimizer.minimize();

    EMFit emfit(*fitter, res, res.fval);
    emfit.evaluated_points = minimizer.get_evaluated_points();

    return {emfit, landscape};
}

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