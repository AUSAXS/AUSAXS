/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/programs/mdrun.h>
#include <md/programs/mdrun/Execution.h>
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
    options.push_back(std::make_shared<shell::Argument>("-s", tpr));
    return *this;
}

mdrun& mdrun::output(const io::Folder& folder, const std::string& prefix) {
    this->folder = folder;
    options.push_back(std::make_shared<shell::Argument>("-deffnm", folder.str() + prefix));
    return *this;
}

mdrun& mdrun::output(const io::Folder& folder) {
    this->folder = folder;
    options.push_back(std::make_shared<shell::Argument>("-multidir", folder));
    return *this;
}

mdrun& mdrun::jobname(const std::string& name) {
    this->name = name;
    return *this;
}

std::unique_ptr<shell::Jobscript<MDRunResult>> mdrun::run(RunLocation where, std::string jobscript) {
    switch (where) {
        case RunLocation::local: {
            return std::make_unique<LocalExecution<MDRunResult>>([*this]() {auto tmp = *this; return tmp.execute();}, folder);
        }
        case RunLocation::lusi: {
            cmd.append("-nt 12 -nice 19 -pin on -pinstride 1 -pinoffset 0 -gpu_id 0");
            return std::make_unique<LocalExecution<MDRunResult>>([*this](){auto tmp = *this; return tmp.execute();}, folder);
        }
        case RunLocation::smaug: {
            cmd.append("-ntmpi 1 -nt $cpupergpu -cpi -stepout 5000 -maxh 48 >& md.lis");
            return std::make_unique<SmaugExecution<MDRunResult>>(command().get(), folder);
        }
        default: 
            throw except::invalid_argument("mdrun: Unknown location.");
    }
}

void mdrun::validate() const {
    if (!tpr.exists()) {
        throw except::missing_option("mdrun: No tpr file specified");
    }
}