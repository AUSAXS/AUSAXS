#include <crystal/miller/AllMillers.h>
#include <crystal/Fval.h>

std::vector<crystal::Miller> crystal::AllMillers::generate() const {
    std::vector<Miller> millers;
    for (int h = 0; h <= this->h; h++) {
        for (int k = -this->k; k <= this->k; k++) {
            for (int l = -this->l; l <= this->l; l++) {
                Miller m(h, k, l);
                double q = crystal::Fval::Q(m).norm();
                if (q < settings::crystal::max_q) {
                    millers.push_back(Miller(h, k, l));
                }
            }
        }
    }
    return millers;
}
