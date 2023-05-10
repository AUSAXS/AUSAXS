#include <mini/MinimumExplorer.h>
#include <mini/detail/Parameter.h>
#include <mini/detail/FittedParameter.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

using namespace mini;

MinimumExplorer::MinimumExplorer(double(&func)(std::vector<double>), unsigned int evals) : Minimizer(func) {
    set_max_evals(evals);
}

MinimumExplorer::MinimumExplorer(std::function<double(std::vector<double>)> func, unsigned int evals) : Minimizer(func) {
    set_max_evals(evals);
}

MinimumExplorer::MinimumExplorer(double(&func)(std::vector<double>), const Parameter& param, unsigned int evals) : Minimizer(func) {
    set_max_evals(evals);
    add_parameter(param);
}

MinimumExplorer::MinimumExplorer(std::function<double(std::vector<double>)> func, const Parameter& param, unsigned int evals) : Minimizer(func) {
    set_max_evals(evals);
    add_parameter(param);
}

mini::Landscape MinimumExplorer::landscape(unsigned int evals) {
    if (parameters.empty()) {throw except::bad_order("MinimumExplorer::landscape: No parameters were supplied.");}
    if (!evaluations.evals.empty()) {return evaluations;} // if the minimizer has already been called, we can just reuse its result

    const Parameter& param = parameters[0];
    double xmin = *param.guess;
    double xmid = xmin;
    double x = xmin;
    double fmin = function({xmin});

    //###############################################################//
    //###        DETERMINE SPACING BETWEEN EVALUATIONS            ###//
    //###############################################################//
    // we want to find the smallest spacing that still changes the function value
    double spacing = 1e-5;
    unsigned int vchanges = 0;  // we want at least 2 changes for a decent estimate
    unsigned int counter = 0;   // counter to prevent infinite loop
    double factor = 2;          // step scaling factor
    record_evaluations(false);  // disable recording while determining spacing

    // reset the spacing
    auto reset_spacing = [&] () {
        spacing = 1e-5;
        if (param.has_bounds()) {
            spacing = std::min(spacing, (param.bounds->max - param.bounds->min)/evals);
        }
    };

    // update the spacing
    double fprev = fmin;
    auto update_spacing = [&] (double x) {
        double f = function({x});

        // check if the function value changed
        if (1e-6 < std::abs(f - fprev)) {
            vchanges++;
            if (2 < vchanges) {return true;} // stop after fval changed twice
            fprev = f;
            factor = 1.3;                    // scale slower since we're close to the final step size
        } else {
            spacing *= factor;
        }
        counter++;
        return false;
    };

    reset_spacing();
    while (counter < 20) { 
        x += spacing;                   // move to the right
        if (update_spacing(x)) {break;} // check if we have enough changes
        x = xmid;                       // move back to the middle
    }

    // if counter is 20 we couldn't determine the spacing by going to the right, so we try again but now going left
    if (counter == 20) {
        counter = 0;
        x = xmid;
        reset_spacing();
        while (counter < 20) {
            x -= spacing;                   // move to the left
            if (update_spacing(x)) {break;} // check if we have enough changes
            x = xmid;                       // move back to the middle
        }
    }

    // if counter is still 20 going left didn't work either. the function is probably ill-defined
    if (counter == 20) {
        throw except::bad_order("MinimumExplorer::landscape: Could not determine spacing for landscape.");
    }

    spacing /= 2; // step size is twice the distance between fval changes after ending the earlier loop


    //###############################################################//
    //###                 ESTIMATE MINIMUM VALUE                  ###//
    //###############################################################//
    record_evaluations(true); // start recording again

    // go three steps to the left
    x = xmid;               // go back to the middle
    fprev = fmin;           // keep track of last value
    counter = 0;            // reset counter
    unsigned int iter = 0;  // keep track of how many iterations we've done
    double start_space = spacing;
    for (int i = 0; i < 4; i++) {
        x -= spacing;
        double f = function({x});

        // check if the function value actually changed
        if (std::abs(fprev - f) < 1e-6) {
            // if not, refine the spacing and try again
            iter++;
            x += spacing;
            spacing *= 1.3;
            i--;
            evaluations.evals.pop_back();

            // if we've tried too many times going to the left doesn't work
            if (iter == 10) {
                spacing = start_space;
                counter = 4; // mark left side as being finished
                break;
            }
        }

        // check if this is a new minimum
        if (f < fmin) {
            fmin = f;
            xmin = x;
            continue;
        } 

        // check if this value is higher than the previous one
        if (fprev < f) {
            fprev = f;
            counter++;
        }
    }
    // if counter == 4 the function is monotonically increasing to the left, and we shouldn't explore it further
    bool left = !(counter == 4);

    // go three steps to the right
    x = xmid;       // go back to the middle
    fprev = fmin;   // reset fprev
    counter = 0;    // reset counter
    iter = 0;       // reset iter
    start_space = spacing;
    for (int i = 0; i < 4; i++) {
        x += spacing;
        double f = function({x});

        // check if the function value actually changed
        if (std::abs(fprev - f) < 1e-6) {
            // if not, refine the spacing and try again
            iter++;
            x -= spacing;
            spacing *= 1.3;
            i--;
            evaluations.evals.pop_back();

            // if we've tried too many times going to the right doesn't work
            if (iter == 10) {
                spacing = start_space;
                counter = 4; // mark right side as being finished
                break;
            }
        }

        // check if this is a new minimum
        if (f < fmin) {
            fmin = f;
            xmin = x;
            continue;
        }

        // check if this value is higher than the previous one
        if (fprev < f) {
            fprev = f;
            counter++;
        }
    }
    // if counter == 4 the function is monotonically increasing to the right, and we shouldn't explore it further
    bool right = !(counter == 4);

    // we now change tactics: instead of requiring 3 monotonic increases in fval before stopping, we now just want it to be higher than the mean four times in a row
    auto points = get_evaluated_points().as_dataset();
    double mu = points.mean();

    if (right) {
        // now go the remaining steps to the right, terminating if four consecutive evals are all above the mean
        counter = 0;
        x = xmid + 4*spacing;   // start four steps to the right of the middle
        unsigned int above = 0; // number of consecutive points higher than the mean
        while (above < 4 && counter++ < (evals-9)/2) {
            x += spacing;
            double f = function({x});
            if (f < fmin) {
                fmin = f;
                xmin = x;
            }

            if (mu < f) {
                above++;
            } else {
                above = 0;
            }
        }
    }

    if (left) {
        // repeat for left-steps
        counter = 0;
        unsigned int above = 0;
        x = xmid - 4*spacing;   // start four steps to the left of the middle
        while (above < 4 && counter++ < (evals-9)/2) {
            x -= spacing;
            double f = function({x});
            if (f < fmin) {
                fmin = f;
                xmin = x;
            }

            if (mu < f) {
                above++;
            } else {
                above = 0;
            }
        }
    }

    return get_evaluated_points();
}

Result MinimumExplorer::minimize_override() {
    auto l = landscape(max_evals).as_dataset();
    auto min = l.find_minimum();
    FittedParameter p(parameters[0], min.x, l.span_x() - min.x);
    return Result(p, l.mean(), fevals);
}

void MinimumExplorer::add_parameter(const Parameter& param) {
    if (!param.has_guess()) {throw except::invalid_operation("MinimumExplorer::add_parameter: Guess value must be supplied.");}
    if (!parameters.empty()) {throw except::invalid_operation("MinimumExplorer::add_parameter: This minimizer only supports 1D problems.");}
    parameters.push_back(param);
}