#pragma once

#include <md/shell/Jobscript.h>
#include <md/utility/files/all.h>
#include <md/programs/mdrun/MDRunResult.h>

#include <functional>

namespace ausaxs::md {
    enum class RunLocation {
        local,
        lusi,
        smaug
    };

    /**
     * @brief A Jobscript that does not execute anything. 
     * 
     * This can be used when the execution of an mdrun must be skipped, and only existing results are needed.
     */
    template<typename T>
    class NoExecution : public shell::Jobscript<T> {
        public: 
            NoExecution(const std::string& folder) : folder(folder) {}
            virtual ~NoExecution() = default;
            void create() override {}
            void submit() override {}
            void wait() override {}
            T result() override {return T(folder);}
        
        private:
            std::string folder;
    };

    /**
     * @brief A Jobscript that executes the mdrun locally.
     */
    template<typename T>
    class LocalExecution : public shell::Jobscript<T> {
        public: 
            LocalExecution(std::function<std::string()> func, const std::string& folder) : func(std::move(func)), folder(folder) {}
            ~LocalExecution() override = default;
            void create() override {}

            void submit() override {
                if (submitted) {return;}
                func();
                submitted = true;
            }

            void wait() override {
                if (!submitted) {submit();}
            }

            T result() override {
                if (!submitted) {submit();}
                return T(folder);
            }

        private:
            bool submitted = false;
            std::function<std::string()> func;
            io::Folder folder; 
    };

    /**
     * @brief A Jobscript that executes the mdrun on the Smaug cluster.
     */
    template<typename T>
    class SmaugExecution : public shell::Slurmscript<T> {
        public: 
            ~SmaugExecution() override = default;

            SmaugExecution(std::string_view run_cmd, const io::Folder& folder) : folder(folder) {
                this->filename = folder + "/job.sh";
                std::string out = 
                    "#!/bin/bash\n"
                    "#SBATCH -o " + folder.absolute_path() + "/job.out\n"
                    "#SBATCH -e " + folder.absolute_path() + "/job.err\n"
                    "#SBATCH --time 24:00:00\n"
                    "#SBATCH --nodes 1\n"
                    "#SBATCH --gres gpu:1\n"
                    "#SBATCH --ntasks-per-node 12\n"
                    "#SBATCH --exclude fang[51,53]\n\n"
                    "cpu_count=$(grep -c ^processor /proc/cpuinfo)\n"
                    "gpu_count=$(nvidia-smi -L | wc -l)\n"
                    "threads=$[cpu_count/gpu_count]\n\n"
                    "module use /data/shared/spack/0.21.1+240303/modules\n"
                    "module add gromacs-swaxs\n\n"
                    "cd " + folder.absolute_path() + "\n"
                    "" + std::string(run_cmd) + " -ntmpi 1 -ntomp $threads -maxh 48 -stepout 5000 >& md.lis\n";
                io::File(folder + "job.sh").create(out);
            }

            T result() override {return T(folder);}

        private: 
            std::string folder;
    };
}
