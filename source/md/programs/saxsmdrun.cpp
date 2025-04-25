/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <md/programs/saxsmdrun.h>
#include <md/programs/mdrun/MDExecutor.h>
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

saxsmdrun& saxsmdrun::output(const io::Folder& folder, const std::string& prefix) {
    this->folder = folder;
    options.push_back(std::make_shared<shell::Argument>("-deffnm", folder.absolute_path() + prefix));
    return *this;
}

saxsmdrun& saxsmdrun::jobname(const std::string& name) {
    this->name = name;
    return *this;
}

saxsmdrun& saxsmdrun::env_var(const std::string& var, const std::string& value) {
    cmd.prepend("export " + var + "=" + value + "; ");
    return *this;
}

std::unique_ptr<Executor<SAXSRunResult>> saxsmdrun::run(std::unique_ptr<executor::type> executor) {
    return executor->saxs_runner(folder, command().get());
}

void saxsmdrun::validate() const {
    if (!moltpr.exists() || !buftpr.exists()) {
        throw except::missing_option("saxsmdrun: No tpr file specified");
    }
}