#pragma once

#include <md/programs/gmx.h>
#include <md/programs/mdrun/MDRunResult.h>
#include <md/programs/mdrun/MDExecutor.h>
#include <md/utility/files/all.h>

#include <vector>

namespace ausaxs::md {
    class mdrun : private gmx {
        public: 
            mdrun();
            mdrun(const TPRFile& tpr);
            mdrun& input(const TPRFile& tpr);
            mdrun& output(const io::Folder& folder, const std::string& prefix);
            mdrun& jobname(const std::string& name);

            /**
             * @brief Add a PLUMED input file (passed as -plumed to mdrun).
             */
            mdrun& plumed(const io::File& plumed_input);

            /**
             * @brief Configure a multi-replica run (-multidir).
             *
             * Replaces output(): the command will use -multidir <dirs...> -deffnm prod.
             * The base folder (parent of the first replica directory) is stored
             * for use by run_multi().
             */
            mdrun& multidir(const std::vector<io::Folder>& dirs);

            std::unique_ptr<Executor<MDRunResult>> run(std::unique_ptr<executor::type> executor);

            /**
             * @brief Launch a multi-replica run and return an executor handle.
             *
             * Must be called after multidir(). The returned Executor's result()
             * yields a MultiMDRunResult with individual per-replica outputs.
             */
            std::unique_ptr<Executor<MultiMDRunResult>> run_multi(std::unique_ptr<executor::type> executor);

        private: 
            TPRFile tpr;
            io::Folder folder;
            io::Folder multidir_base;
            std::string name = "run";

            void validate() const override;
    };
}
