#include <utility/Dataset.h>
#include <utility/Settings.h>

#include <fstream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

void Dataset::rebin() noexcept {
    Dataset rebinned; // rebinned dataset

    for (unsigned int i = 0; i < size(); i++) {
        // determine how many data points to fold into one
        unsigned int fold;
        if (0.1 < data[i].x) {fold = 8;}
        else if (0.06 < data[i].x) {fold = 4;}
        else if (0.03 < data[i].x) {fold = 2;}
        else {fold = 1;}

        std::cout << "now folding " << i << " to " << i + fold << std::endl;

        // loop over each data point to be folded
        double siginv = 0, sumw = 0, qsum = 0;
        unsigned int ss = 0;
        for (; ss < fold; ss++) {
            std::cout << "checkpoint1" << std::endl;
            if (i == size()) {break;}
            std::cout << "checkpoint1" << std::endl;
            siginv += (std::pow(data[i].yerr, -2));
            std::cout << "checkpoint1" << std::endl;
            sumw += data[i].y/(std::pow(data[i].yerr, 2));
            std::cout << "checkpoint1" << std::endl;
            qsum += data[i].x;
            std::cout << "checkpoint1" << std::endl;
            i++;
        }

        // average their values into a single new one
        double q = qsum/ss;
        double I = sumw/siginv;
        double Ierr = std::pow(siginv, -0.5);
        push_back(ErrorPoint2D(q, I, 0, Ierr));
    }
    rebinned.save("temp/dataset/test.dat");
    *this = rebinned;
}

void SAXSDataset::load(const std::string file) {
    // check if file was succesfully opened
    std::ifstream input(file);
    if (!input.is_open()) {throw std::ios_base::failure("Error in IntensityFitter::read: Could not open file \"" + file + "\"");}

    std::string line; // placeholder for the current line
    while(getline(input, line)) {
        if (line[0] == ' ') {line = line.substr(1);} // fix leading space
        std::vector<std::string> tokens;
        boost::split(tokens, line, boost::is_any_of(" ,\t")); // spaces, commas, and tabs can all be used as separators (but not a mix of them)

        // determine if we are in some sort of header
        if (tokens.size() < 2 || tokens.size() > 4) {continue;} // too many separators
        bool skip = false;
        for (unsigned int i = 0; i < tokens.size(); i++) { // check if they are numbers
            if (!tokens[i].empty() && tokens[i].find_first_not_of("0123456789-.Ee") != std::string::npos) {skip = true;}
        }
        if (skip) {continue;}

        // now we are most likely beyond any headers
        double _q, _I, _sigma;
        _q = std::stod(tokens[0]); // we know for sure that the strings are convertible to numbers (boost check)
        _I = std::stod(tokens[1]);
        if (_q > 10) {continue;} // probably not a q-value if it's larger than 10

        // check user-defined limits
        if (_q < setting::fit::q_low) {continue;}
        if (_q > setting::fit::q_high) {continue;}

        // add the values to our vectors
        ErrorPoint2D point(_q, _I, 0, 0);

        // if x | y | yerr
        if (tokens.size() == 3) {
            point.yerr = std::stod(tokens[2]);
        }

        data.push_back(point);
    }
    input.close();
}