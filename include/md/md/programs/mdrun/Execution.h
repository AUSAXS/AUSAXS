#pragma once

#include <md/shell/Jobscript.h>
#include <md/utility/files/all.h>
#include <md/programs/mdrun/MDRunResult.h>
#include <md/gmx/Settings.h>

#include <functional>

namespace md {
    enum class location {
        local,
        lucy,
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
            LocalExecution(std::function<std::string()> func, const std::string& folder) : func(func), folder(folder) {}
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
            std::string folder; 
    };

    /**
     * @brief A Jobscript that executes the mdrun on the Smaug cluster.
     */
    template<typename T>
    class SmaugExecution : public shell::Slurmscript<T> {
        public: 
            ~SmaugExecution() override = default;

            SmaugExecution(const TPRFile& tpr, const std::string& folder, const std::string& name, const SHFile& jobscript) : folder(folder) {
                this->filename = tpr.parent_path() + "/job.sh";
                this->cmd.set("cd " + tpr.parent_path() + "; " + jobscript.path + " -j " + name + " -version release-2021.swaxs -gpu-gener any+ampere -tpr " + tpr.absolute());
            }

            SmaugExecution(const std::string& args, const std::string& _export, const std::string& folder, const std::string& name, const SHFile& jobscript) : folder(folder) {
                this->filename = folder + "/job.sh";
                this->cmd.set(_export + "cd " + folder + "; " + jobscript.path + " -j " + name + " -version release-2021.swaxs -gpu-gener any+ampere " + args);
            }

            T result() override {return T(folder);}

        private: 
            std::string folder;
    };
}
