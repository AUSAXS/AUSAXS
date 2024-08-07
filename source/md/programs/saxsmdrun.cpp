/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/programs/saxsmdrun.h>
#include <md/programs/mdrun/Execution.h>
#include <md/utility/Exceptions.h>

using namespace ausaxs;
using namespace ausaxs::md;

saxsmdrun::saxsmdrun() {
    cmd.append("mdrun -v -cpi");
}

saxsmdrun::saxsmdrun(const TPRFile& mol, const TPRFile& buf) : saxsmdrun() {
    input(mol, buf);
}

saxsmdrun& saxsmdrun::input(const TPRFile& mol, const TPRFile& buf) {
    this->moltpr = mol;
    this->buftpr = buf;
    options.push_back(std::make_shared<shell::Argument>("-s", mol.absolute_path()));
    options.push_back(std::make_shared<shell::Argument>("-sw", buf.absolute_path()));
    return *this;
}

saxsmdrun& saxsmdrun::rerun(const XTCFile& mol, const XTCFile& buf) {
    options.push_back(std::make_shared<shell::Argument>("-rerun", mol.absolute_path()));
    options.push_back(std::make_shared<shell::Argument>("-fw", buf.absolute_path()));
    return *this;
}

saxsmdrun& saxsmdrun::output(const io::Folder& folder, const std::string&) {
    this->folder = folder;
    options.push_back(std::make_shared<shell::Argument>("-multidir", folder));
    options.push_back(std::make_shared<shell::Argument>("-ow", folder));
    return *this;
}

saxsmdrun& saxsmdrun::jobname(const std::string& name) {
    this->name = name;
    return *this;
}

saxsmdrun& saxsmdrun::env_var(const std::string& var, const std::string& value) {
    _export += "export " + var + "=" + value + "; ";
    return *this;
}

std::unique_ptr<shell::Jobscript<SAXSRunResult>> saxsmdrun::run(location where, std::string jobscript) {
    switch (where) {
        case location::local: {
            cmd.prepend(_export);
            return std::make_unique<LocalExecution<SAXSRunResult>>([*this](){auto tmp = *this; return tmp.execute();}, folder);
        }
        case location::lusi: {
            cmd.append("-nt 12 -nice 19 -pin on -pinstride 1 -pinoffset 0 -gpu_id 0");
            cmd.prepend(_export + "cd " + folder + ";");
            return std::make_unique<LocalExecution<SAXSRunResult>>([*this](){auto tmp = *this; return tmp.execute();}, folder);
        }
        case location::smaug: {
            std::string args = "";
            for (auto& option : options) {
                if (option->name == "-s") {
                    args.append("-tpr " + option->value + " ");
                } else if (option->name == "-rerun") {
                    args.append("-xtc " + option->value + " ");
                } else if (option->name == "-sw") {
                    args.append(option->get() + " ");
                } else if (option->name == "-fw") {
                    args.append(option->get() + " ");
                }
            }
            cmd.append("-ntmpi 1 -nt $cpupergpu -cpi -stepout 5000");
            return std::make_unique<SmaugExecution<SAXSRunResult>>(args, _export, folder, name, jobscript);
        }
        default: 
            throw except::invalid_argument("saxsmdrun: Unknown location.");
    }
}

void saxsmdrun::validate() const {
    if (!moltpr.exists() || !buftpr.exists()) {
        throw except::missing_option("saxsmdrun: No tpr file specified");
    }
}