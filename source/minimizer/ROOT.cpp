#include <minimizer/ROOT.h>

mini::ROOT::ROOT(std::string package, std::string algorithm) {}

mini::ROOT::ROOT(std::string package, std::string algorithm, double(&f)(double*), Parameter param = Parameter()) {}

mini::Result mini::ROOT::minimize() const {}

void mini::ROOT::add_parameter(const Parameter& param) {}

Dataset mini::ROOT::landscape(unsigned int evals = 100) const {}