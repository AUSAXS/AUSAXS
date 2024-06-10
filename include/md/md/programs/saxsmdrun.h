#pragma once

#include <md/programs/gmx.h>
#include <md/programs/mdrun/MDRunResult.h>
#include <md/programs/mdrun/Execution.h>
#include <md/utility/files/all.h>

namespace gmx {
    class saxsmdrun : private gmx {
        public: 
            saxsmdrun();
            saxsmdrun(const TPRFile& moltpr, const TPRFile& buftpr);
            saxsmdrun& input(const TPRFile& moltpr, const TPRFile& buftpr);
            saxsmdrun& output(const Folder& folder, const std::string& prefix);
            saxsmdrun& rerun(const XTCFile& mol, const XTCFile& buf);
            saxsmdrun& env_var(const std::string& var, const std::string& value);
            std::unique_ptr<shell::Jobscript<SAXSRunResult>> run(location where, std::string jobscript = "");
            saxsmdrun& jobname(const std::string& name);

        private: 
            TPRFile moltpr, buftpr;
            std::string folder, _export;
            std::string name = "saxsrun";

            void validate() const override;
    };
}