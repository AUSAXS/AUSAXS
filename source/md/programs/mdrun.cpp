/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/programs/mdrun.h>
#include <md/programs/mdrun/MDExecutor.h>
#include <md/utility/Exceptions.h>

#include <stdexcept>

using namespace ausaxs;
using namespace ausaxs::md;

mdrun::mdrun() {
    cmd.append("mdrun -v -cpi");
}

mdrun::mdrun(const TPRFile& tpr) : mdrun() {
    input(tpr);
}

mdrun& mdrun::input(const TPRFile& tpr) {
    this->tpr = tpr;
    options.push_back(std::make_shared<shell::Argument>("-s", tpr.absolute_path()));
    return *this;
}

mdrun& mdrun::output(const io::Folder& folder, const std::string& prefix) {
    this->folder = folder;
    options.push_back(std::make_shared<shell::Argument>("-deffnm", folder.absolute_path() + prefix));
    return *this;
}

mdrun& mdrun::jobname(const std::string& name) {
    this->name = name;
    return *this;
}

mdrun& mdrun::plumed(const io::File& plumed_input) {
    options.push_back(std::make_shared<shell::Argument>("-plumed", plumed_input.absolute_path()));
    return *this;
}

mdrun& mdrun::multidir(const std::vector<io::Folder>& dirs) {
    if (dirs.empty()) {throw std::invalid_argument("mdrun::multidir: no directories provided");}

    // Derive shared base folder (parent of the first replica directory).
    std::string p = dirs[0].absolute_path();
    auto slash = p.rfind('/');
    multidir_base = io::Folder((slash == std::string::npos) ? p : p.substr(0, slash));
    this->folder  = multidir_base;

    // Build space-separated list of replica directory paths for -multidir.
    std::string multi_value;
    for (const auto& d : dirs) {
        if (!multi_value.empty()) {multi_value += " ";}
        multi_value += d.absolute_path();
    }
    options.push_back(std::make_shared<shell::Option>("-multidir", multi_value));
    options.push_back(std::make_shared<shell::Argument>("-deffnm", "prod"));
    return *this;
}

std::unique_ptr<Executor<MDRunResult>> mdrun::run(std::unique_ptr<executor::type> executor) {
    return executor->md_runner(folder, command().get());
}

std::unique_ptr<Executor<MultiMDRunResult>> mdrun::run_multi(std::unique_ptr<executor::type> executor) {
    if (multidir_base.empty()) {throw std::logic_error("mdrun::run_multi: call multidir() before run_multi()");}
    return executor->multi_md_runner(multidir_base, command().get());
}

void mdrun::validate() const {
    // For multi-replica runs multidir_base is set instead of tpr.
    if (!multidir_base.empty()) {return;}
    if (!tpr.exists()) {
        throw except::missing_option("mdrun: No tpr file specified");
    }
}
