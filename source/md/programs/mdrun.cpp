#include <md/programs/mdrun.h>
#include <md/programs/mdrun/Execution.h>
#include <md/utility/Exceptions.h>

using namespace gmx;

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

mdrun& mdrun::output(const Folder& folder, const std::string& prefix) {
    this->folder = folder;
    options.push_back(std::make_shared<shell::Argument>("-deffnm", folder + prefix));
    return *this;
}

mdrun& mdrun::output(const Folder& folder) {
    this->folder = folder;
    options.push_back(std::make_shared<shell::Argument>("-multidir", folder));
    return *this;
}

mdrun& mdrun::jobname(const std::string& name) {
    this->name = name;
    return *this;
}

std::unique_ptr<shell::Jobscript<MDRunResult>> mdrun::run(location where, std::string jobscript) {
    switch (where) {
        case location::local: {
            return std::make_unique<LocalExecution<MDRunResult>>([*this]() {auto tmp = *this; return tmp.execute();}, folder);
        }
        case location::lucy: {
            cmd.append("-nt 12 -nice 19 -pin on -pinstride 1 -pinoffset 0 -gpu_id 0");
            return std::make_unique<LocalExecution<MDRunResult>>([*this](){auto tmp = *this; return tmp.execute();}, folder);
        }
        case location::smaug: {
            cmd.append("-ntmpi 1 -nt $cpupergpu -cpi -stepout 5000");
            return std::make_unique<SmaugExecution<MDRunResult>>(tpr, folder, name, jobscript);
        }
        default: 
            throw except::invalid_argument("mdrun: Unknown location.");
    }
}

void mdrun::validate() const {
    if (tpr.empty()) {
        throw except::missing_option("mdrun: No tpr file specified");
    }
}