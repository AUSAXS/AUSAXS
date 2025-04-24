#pragma once

#include <md/programs/gmx.h>
#include <md/programs/mdrun/MDRunResult.h>
#include <md/programs/mdrun/MDExecutor.h>
#include <md/utility/files/all.h>

namespace ausaxs::md {
    class mdrun : private gmx {
        public: 
            mdrun();
            mdrun(const TPRFile& tpr);
            mdrun& input(const TPRFile& tpr);
            mdrun& output(const io::Folder& folder, const std::string& prefix);
            mdrun& output(const io::Folder& folder);
            mdrun& jobname(const std::string& name);
            std::unique_ptr<Executor<MDRunResult>> run(std::unique_ptr<executor::type> executor);

        private: 
            TPRFile tpr;
            std::string folder;
            std::string name = "run";

            void validate() const override;
    };
}