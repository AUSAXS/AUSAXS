/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/programs/mdrun.h>
#include <md/programs/mdrun/MDExecutor.h>
#include <md/utility/Exceptions.h>

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

std::unique_ptr<Executor<MDRunResult>> mdrun::run(std::unique_ptr<executor::type> executor) {
    return executor->md_runner(folder, command().get());
}

void mdrun::validate() const {
    if (!tpr.exists()) {
        throw except::missing_option("mdrun: No tpr file specified");
    }
}