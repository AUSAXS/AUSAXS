#include <em/ImageStack.h>
#include <em/PartialHistogramManager.h>
#include <fitter/IntensityFitter.h>
#include <plots/all.h>

using namespace em;

std::shared_ptr<EMFit> ImageStack::fit(const hist::ScatteringHistogram& h) {
    Limit lim = {from_level(setting::em::alpha_levels.min), from_level(setting::em::alpha_levels.max)};
    mini::Parameter param("cutoff", lim.center(), lim);
    return fit(h, param);
}

std::shared_ptr<EMFit> ImageStack::fit(const hist::ScatteringHistogram& h, mini::Parameter param) {
    if (!param.has_bounds()) {return fit(h);} // ensure parameter bounds are present

    std::shared_ptr<SimpleIntensityFitter> fitter = setting::em::hydrate ? std::make_shared<IntensityFitter>(h, get_limits()) : std::make_shared<SimpleIntensityFitter>(h, get_limits());
    return fit_helper(fitter, param);
}

std::shared_ptr<EMFit> ImageStack::fit(std::string file) {
    Limit lim = {from_level(setting::em::alpha_levels.min), from_level(setting::em::alpha_levels.max)};
    mini::Parameter param("cutoff", lim.center(), lim);
    return fit(file, param);
}

std::shared_ptr<EMFit> ImageStack::fit(std::string file, mini::Parameter param) {
    if (!param.has_bounds()) {return fit(file);} // ensure parameter bounds are present
    std::shared_ptr<SimpleIntensityFitter> fitter = setting::em::hydrate ? std::make_shared<IntensityFitter>(file) : std::make_shared<SimpleIntensityFitter>(file);
    return fit_helper(fitter, param);
}

std::shared_ptr<EMFit> ImageStack::fit_helper(std::shared_ptr<SimpleIntensityFitter> fitter, mini::Parameter param) {
    update_charge_levels(*param.bounds);
    set_minimum_bounds(param.bounds->min);
    auto f = prepare_function(fitter);
    mini::Landscape evals; // since we'll be using multiple minimizers, we'll need to store the evaluated points manually

    //##########################################################//
    //###                DETERMINE LANDSCAPE                 ###//
    //##########################################################//
    mini::LimitedScan minimizer(f, param, setting::em::evals);
    minimizer.set_limit(fitter->dof()*20);
    SimpleDataset d;
    {
        auto l = minimizer.landscape(setting::em::evals);
        evals.append(l);
        d = l.as_dataset();
    }

    d.sort_x();
    auto min = d.find_minimum();

    //##########################################################//
    //### CHECK LANDSCAPE IS OK FOR AVERAGING & INTERPLATION ###//
    //##########################################################//
    d.limit_y(0, min.y*5); // focus on the area near the minimum
    if (d.size() < 10) {    // if we have too few points after imposing the limit, we must sample some more
        Limit bounds;       // first we determine the bounds of the area we want to sample
        if (d.size() < 3) { // if we only have one or two points, sample the area between the neighbouring points
            double s = (param.bounds->max - param.bounds->min)/setting::em::evals;
            bounds = {min.x - s, min.x + s};
        }
        else { // otherwise just use the new bounds of the limited landscape
            bounds = d.span_x();
        }

        // prepare a new minimizer with the new bounds
        utility::print_warning("Function is varying strongly. Sampling more points around the minimum.");
        mini::LimitedScan mini2(f, mini::Parameter("cutoff", bounds), setting::em::evals/4);
        {
            auto l = mini2.landscape(setting::em::evals/2);
            evals.append(l);
            d = l.as_dataset();
        }
        d.sort_x();
        min = d.find_minimum();
        d.limit_y(0, min.y*5);

        if (d.size() < 10) {
            throw except::unexpected("ImageStack::fit: Could not sample enough points around the minimum. Function varies too much.");
        }
    }

    //##########################################################//
    //###         AVERAGE & INTERPLATE MORE POINTS           ###//
    //##########################################################//
    SimpleDataset avg = d.rolling_average(7); // impose a moving average filter 
    avg.interpolate(5);                       // interpolate more points

    min = avg.find_minimum();
    double spacing = avg.x(1)-avg.x(0); 
    param.guess = min.x;
    param.bounds = Limit(min.x-3*spacing, min.x+3*spacing); // uncertainty is 3*spacing between points

    // optional plot
    if (setting::plot::em::additional_plots) {
        // plot the minimum in blue
        SimpleDataset p_start;
        p_start.push_back(min.x, min.y);
        p_start.add_plot_options(style::draw::points, {{"color", style::color::blue}, {"s", 9}});

        avg.add_plot_options(style::draw::line, {{"color", style::color::red}, {"xlabel", "cutoff"}, {"ylabel", "$\\chi^2$"}});
        plots::PlotDataset plot(avg);
        d.add_plot_options(style::draw::points);
        plot.plot(d);
        plot.plot(p_start);
        plot.save(setting::plot::path + "chi2_evaluated_points." + setting::plot::format);
    }

    //##########################################################//
    //###             EXPLORE AREA AROUND MINIMUM            ###//
    //##########################################################//
    // if hydration is enabled, the chi2 will oscillate heavily around the minimum
    // we therefore want to sample the area near the minimum to get an average
    mini::Result res;
    if (setting::em::hydrate) {
        // sample the area around the minimum
        mini::MinimumExplorer explorer(f, param, setting::em::evals);
        res = explorer.minimize();

        SimpleDataset area;
        {
            auto l = explorer.landscape();
            evals.append(l);
            area = l.as_dataset();
        }
        if (setting::plot::em::additional_plots) {
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
            plot.save(setting::plot::path + "chi2_near_minimum." + setting::plot::format);
        }
    } 
    
    // otherwise do a quick fit to ensure we're at the very bottom of the valley
    else {
        mini::Golden golden(f, param);
        res = golden.minimize();
        evals.append(golden.get_evaluated_points());
    }

    // make landscape plot
    if (setting::plot::em::additional_plots && setting::em::hydrate) {
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
        plots::PlotLandscape::quick_plot(l, setting::plot::path + "chi2_landscape." + setting::plot::format);
    }

    // update the fitter with the optimal cutoff, such that the returned fit is actually the best one
    min = evals.as_dataset().find_minimum();
    f({min.x});

    std::shared_ptr<EMFit> emfit = std::make_shared<EMFit>(*fitter, res, min.y);
    emfit->evaluated_points = evals;
    emfit->fevals = evals.evals.size();
    emfit->level = to_level(min.x);
    if (setting::em::save_pdb) {get_histogram_manager()->get_protein()->save(setting::plot::path + "model.pdb");}
    return emfit;
}

std::function<double(std::vector<double>)> ImageStack::prepare_function(std::shared_ptr<SimpleIntensityFitter> fitter) {
    // convert the calculated intensities to absolute scale
    // utility::print_warning("Warning in ImageStack::prepare_function: Not using absolute scale.");
    // auto protein = phm->get_protein(1);
    // double c = setting::em::concentration;                                // concentration
    // double m = protein->get_absolute_mass()*constants::unit::mg;          // mass
    // double DrhoV2 = std::pow(protein->get_relative_charge(), 2);          // charge
    // double re2 = pow(constants::radius::electron*constants::unit::cm, 2); // squared scattering length
    // double I0 = DrhoV2*re2*c/m;
    // fitter.normalize_intensity(I0);

    // fit function
    static unsigned int counter;        // counter for the number of fevals. must be static to avoid going out of scope
    counter = 0;                        // reset counter to 0 every time prepare_function is called
    setting::protein::center = false;   // do not center the protein - this may cause issues
    setting::grid::percent_water = 0.05;
    if (setting::plot::em::landscape && setting::em::hydrate) {
        std::static_pointer_cast<IntensityFitter>(fitter)->set_algorithm(mini::type::SCAN); 
    }
    
     // fitter is captured by value to guarantee its lifetime will be the same as the lambda
     // 'this' is ok since prepare_function is private and thus only used within the class itself
    std::function<double(std::vector<double>)> chi2 = [this, fitter] (std::vector<double> params) {
        double last_c = 5; // arbitrary start value
        auto p = get_histogram_manager()->get_protein(params[0]);

        std::shared_ptr<Fit> fit;
        if (setting::em::hydrate) {
            p->clear_grid();                // clear grid from previous iteration
            p->generate_new_hydration();    // generate a new hydration layer

            // pointer cast is ok since the type should always be IntensityFitter when hydration is enabled
            std::static_pointer_cast<IntensityFitter>(fitter)->set_guess(mini::Parameter{"c", last_c, {0, 200}}); 
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
        if (setting::fit::verbose) {
            std::cout << "Step " << utility::print_element(counter++, 5) << ": Evaluated cutoff value " << utility::print_element(params[0], 12) << " with chi2 " << utility::print_element(val, 12) << std::flush << "\r";
        }
        return val;
    }; 
    return chi2;
}

mini::Landscape ImageStack::cutoff_scan(const Axis& points, std::string file) {
    std::shared_ptr<SimpleIntensityFitter> fitter = setting::em::hydrate ? std::make_shared<IntensityFitter>(file) : std::make_shared<SimpleIntensityFitter>(file);
    return cutoff_scan_helper(points, fitter);
}

mini::Landscape ImageStack::cutoff_scan(unsigned int points, std::string file) {
    Axis axis(points, from_level(setting::em::alpha_levels.min), from_level(setting::em::alpha_levels.max));
    return cutoff_scan(axis, file);
}

mini::Landscape ImageStack::cutoff_scan(const Axis& points, const hist::ScatteringHistogram& h) {
    std::shared_ptr<SimpleIntensityFitter> fitter = setting::em::hydrate ? std::make_shared<IntensityFitter>(h, get_limits()) : std::make_shared<SimpleIntensityFitter>(h, get_limits());
    return cutoff_scan_helper(points, fitter);
}

mini::Landscape ImageStack::cutoff_scan(unsigned int points, const hist::ScatteringHistogram& h) {
    Axis axis(points, from_level(setting::em::alpha_levels.min), from_level(setting::em::alpha_levels.max));
    return cutoff_scan(axis, h);
}

std::pair<EMFit, mini::Landscape> ImageStack::cutoff_scan_fit(unsigned int points, const hist::ScatteringHistogram& h) {
    Axis axis(points, from_level(setting::em::alpha_levels.min), from_level(setting::em::alpha_levels.max));
    return cutoff_scan_fit(axis, h);
}

std::pair<EMFit, mini::Landscape> ImageStack::cutoff_scan_fit(const Axis& points, std::string file) {
    std::shared_ptr<SimpleIntensityFitter> fitter = setting::em::hydrate ? std::make_shared<IntensityFitter>(file) : std::make_shared<SimpleIntensityFitter>(file);    
    return cutoff_scan_fit_helper(points, fitter);
}

std::pair<EMFit, mini::Landscape> ImageStack::cutoff_scan_fit(unsigned int points, std::string file) {
    Axis axis(points, from_level(setting::em::alpha_levels.min), from_level(setting::em::alpha_levels.max));
    std::shared_ptr<SimpleIntensityFitter> fitter = setting::em::hydrate ? std::make_shared<IntensityFitter>(file) : std::make_shared<SimpleIntensityFitter>(file);    
    return cutoff_scan_fit_helper(axis, fitter);
}

mini::Landscape ImageStack::cutoff_scan_helper(const Axis& points, std::shared_ptr<SimpleIntensityFitter> fitter) {
    update_charge_levels(points.limits());
    set_minimum_bounds(points.min);
    auto func = prepare_function(fitter);

    mini::Golden minimizer(func, mini::Parameter{"cutoff", points.limits()});
    return minimizer.landscape(points.bins);
}

std::pair<EMFit, mini::Landscape> ImageStack::cutoff_scan_fit(const Axis& points, const hist::ScatteringHistogram& h) {
    std::shared_ptr<SimpleIntensityFitter> fitter = setting::em::hydrate ? std::make_shared<IntensityFitter>(h, get_limits()) : std::make_shared<SimpleIntensityFitter>(h, get_limits());
    return cutoff_scan_fit_helper(points, fitter);
}

std::pair<EMFit, mini::Landscape> ImageStack::cutoff_scan_fit_helper(const Axis& points, std::shared_ptr<SimpleIntensityFitter> fitter) {
    update_charge_levels(points.limits());
    set_minimum_bounds(points.min);
    auto func = prepare_function(fitter);

    // cutoff scan
    mini::Golden minimizer(func, mini::Parameter{"cutoff", points.limits()});
    mini::Landscape landscape = minimizer.landscape(points.bins);

    // fit
    double l = from_level(1);
    Limit limit(setting::em::alpha_levels.min*l, setting::em::alpha_levels.max*l);
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

void ImageStack::update_charge_levels(Limit limit) const noexcept {
    std::vector<double> levels;
    for (unsigned int i = 0; i < setting::em::charge_levels; i++) {
        levels.push_back(limit.min + i*limit.span()/setting::em::charge_levels);
    }
    get_histogram_manager()->set_charge_levels(levels);
}