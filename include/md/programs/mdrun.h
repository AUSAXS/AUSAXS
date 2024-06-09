#pragma once

#include <programs/gmx.h>
#include <programs/mdrun/MDRunResult.h>
#include <programs/mdrun/Execution.h>
#include <utility/files/all.h>

namespace gmx {
    class mdrun : private gmx {
        public: 
            mdrun();
            mdrun(const TPRFile& tpr);
            mdrun& input(const TPRFile& tpr);
            mdrun& output(const Folder& folder, const std::string& prefix);
            mdrun& output(const Folder& folder);
            mdrun& jobname(const std::string& name);
            std::unique_ptr<shell::Jobscript<MDRunResult>> run(location where, std::string jobscript = "");

        private: 
            TPRFile tpr;
            std::string folder;
            std::string name = "run";

            void validate() const override;
    };
}