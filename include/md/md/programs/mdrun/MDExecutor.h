#pragma once

#include <md/utility/files/SHFile.h>
#include <md/programs/mdrun/MDRunResult.h>
#include <md/shell/Command.h>
#include <md/utility/Exceptions.h>

#include <fstream>
#include <chrono>
#include <ctime>
#include <thread>

namespace ausaxs::md {
    template<typename T>
    concept MDResult = std::constructible_from<T, const io::Folder&>;

    struct ExecutorBase {
        virtual ~ExecutorBase() = default;
        virtual void submit() = 0;
        virtual void wait() = 0;
    };

    template<MDResult T>
    struct Executor : public ExecutorBase {
        Executor(const io::Folder& output) : output(output) {}
        virtual ~Executor() = default;

        T result() {
            wait();
            return T(output);
        }

        io::Folder output;
    };

    template<MDResult T>
    struct NoExecutor : public Executor<T> {
        NoExecutor(const io::Folder& output) : Executor<T>(output) {}
        void submit() override {}
        void wait() override {}
    };

    template<MDResult T>
    struct LocalExecutor : public Executor<T> {
        LocalExecutor(const io::Folder& output, std::string_view run_cmd) : Executor<T>(output), cmd(run_cmd){
            cmd.prepend("cd " + output.absolute_path() + "; ");
        }

        void submit() override {
            auto result = cmd.execute();
        }

        void wait() override {
            if (!submitted) {submit();}
        }

        bool submitted = false;
        shell::Command cmd;
    };

    template<MDResult T>
    struct TemplateExecutor : public Executor<T>, SHFile {
        TemplateExecutor(const io::Folder& folder) : Executor<T>(folder) {}
        TemplateExecutor(const io::Folder& folder, const SHFile& template_file, std::string_view run_cmd) : Executor<T>(folder) {
            if (!template_file.exists()) {
                throw except::io_error("TemplateRunner: The template file does not exist.");
            }

            SHFile::operator=(folder + "/job.sh");
            std::string out;
            {
                std::ifstream in(template_file.absolute_path());
                std::string line;
                bool found = false;
                while (std::getline(in, line) && !found) {
                    if (line.starts_with("$(mdrun_cmd)")) {
                        line.replace(line.find("$(mdrun_cmd)"), 12, "cd " + folder.absolute_path() + "; " + std::string(run_cmd));
                        found = true;
                    }
                    out += line + "\n";
                }
            }
            SHFile::create(out);
        }

        void submit() override {
            if (!SHFile::exists()) {throw except::io_error("TemplateRunner::submit: The script does not exist.");}
            shell::Command("sh " + this->absolute_path()).execute();
        }

        void wait() override {
            if (!submitted) {submit();}
        }

        bool submitted = false;
    };

    template<MDResult T>
    struct SlurmRunner : public TemplateExecutor<T> {
        SlurmRunner(const io::Folder& folder, const SHFile& template_file, std::string_view run_cmd, std::string_view jobname) : TemplateExecutor<T>(folder) {
            if (!template_file.exists()) {
                throw except::io_error("TemplateRunner: The template file does not exist.");
            }

            SHFile::operator=(folder + "/job.sh");
            std::string out;
            {
                std::ifstream in(template_file.absolute_path());
                std::string line;
                bool found = false;
                while (std::getline(in, line) && !found) {
                    if (line.starts_with("$(mdrun_cmd)")) {
                        line.replace(line.find("$(mdrun_cmd)"), 12, "cd " + folder.absolute_path() + "; " + std::string(run_cmd));
                        found = true;
                    } else if (line.find("$(jobname)") != std::string::npos) {
                        line.replace(line.find("$(jobname)"), 10, jobname);
                    }
                    out += line + "\n";
                }
            }
            SHFile::create(out);
        }

        void submit() override {
            if (!SHFile::exists()) {throw except::io_error("TemplateRunner::submit: The script does not exist.");}
            auto result = shell::Command("sbatch " + this->absolute_path()).execute();
            if (result.exit_code != 0) {
                throw std::runtime_error("TemplateRunner::submit: Error submitting job: \"" + result.out + "\"");
            }

            for (auto& c : result.out) {
                if (isdigit(c)) {
                    jobid += c;
                }
            }
            std::cout << "Submitted job " << jobid << std::endl;
        }

        void wait() override {
            int sleep_time = 5;
            if (jobid.empty()) {submit();}
            while (true) {
                auto result = shell::Command("squeue").execute();
                if (result.out.find(jobid) == std::string::npos) {
                    break;
                }

                char timestamp[50]{ 0 };
                std::time_t time = std::time(nullptr);
                std::strftime(timestamp, 30, "[%H:%M:%S]", std::localtime(&time));

                // print current hour, minute and second
                std::cout << "\r" << timestamp << ": Waiting for job " << jobid << " to finish..." << std::flush;
                std::this_thread::sleep_for(std::chrono::seconds(sleep_time));
                sleep_time = std::min(sleep_time * 2, 120);
            }
        }

        std::string jobid;
    };

    namespace executor {
        enum class id {
            local,
            none,
            templated,
            slurm
        };

        struct type {
            virtual std::unique_ptr<Executor<MDRunResult>> md_runner(const io::Folder& output, std::string_view run_cmd) = 0;
            virtual std::unique_ptr<Executor<SAXSRunResult>> saxs_runner(const io::Folder& output, std::string_view run_cmd) = 0;
        };

        struct local : type {
            static std::unique_ptr<type> construct() {return std::make_unique<local>();}
            std::unique_ptr<Executor<MDRunResult>> md_runner(const io::Folder& output, std::string_view run_cmd) override {
                return std::make_unique<LocalExecutor<MDRunResult>>(output, run_cmd);
            }
            std::unique_ptr<Executor<SAXSRunResult>> saxs_runner(const io::Folder& output, std::string_view run_cmd) override {
                return std::make_unique<LocalExecutor<SAXSRunResult>>(output, run_cmd);
            }
        };

        struct none : type {
            static std::unique_ptr<type> construct() {return std::make_unique<none>();}
            std::unique_ptr<Executor<MDRunResult>> md_runner(const io::Folder& output, std::string_view) override {
                return std::make_unique<NoExecutor<MDRunResult>>(output);
            }
            std::unique_ptr<Executor<SAXSRunResult>> saxs_runner(const io::Folder& output, std::string_view) override {
                return std::make_unique<NoExecutor<SAXSRunResult>>(output);
            }
        };

        struct templated : type {
            templated(const SHFile& template_file) : template_file(template_file) {}
            static std::unique_ptr<type> construct(const SHFile& template_file) {return std::make_unique<templated>(template_file);}
            std::unique_ptr<Executor<MDRunResult>> md_runner(const io::Folder& output, std::string_view run_cmd) override {
                return std::make_unique<TemplateExecutor<MDRunResult>>(output, template_file, run_cmd);
            }
            std::unique_ptr<Executor<SAXSRunResult>> saxs_runner(const io::Folder& output, std::string_view run_cmd) override {
                return std::make_unique<TemplateExecutor<SAXSRunResult>>(output, template_file, run_cmd);
            }
            SHFile template_file;
        };

        struct slurm : type {
            slurm(const SHFile& template_file, std::string_view jobname) : template_file(template_file), jobname(jobname) {}
            static std::unique_ptr<type> construct(const SHFile& template_file, std::string_view jobname) {return std::make_unique<slurm>(template_file, jobname);}
            std::unique_ptr<Executor<MDRunResult>> md_runner(const io::Folder& output, std::string_view run_cmd) override {
                return std::make_unique<SlurmRunner<MDRunResult>>(output, template_file, run_cmd, jobname);
            }
            std::unique_ptr<Executor<SAXSRunResult>> saxs_runner(const io::Folder& output, std::string_view run_cmd) override {
                return std::make_unique<SlurmRunner<SAXSRunResult>>(output, template_file, run_cmd, jobname);
            }
            SHFile template_file;
            std::string jobname;
        };
    }
}
