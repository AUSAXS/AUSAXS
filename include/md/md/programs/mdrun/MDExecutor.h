#pragma once

#include <md/utility/files/SHFile.h>
#include <md/programs/mdrun/MDRunResult.h>
#include <md/shell/Command.h>
#include <md/utility/Exceptions.h>
#include <settings/MDSettings.h>

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
            cmd.prepend(gmx::env_string());
        }

        void submit() override {
            auto result = cmd.execute();
            if (result.exit_code != 0) {
                throw std::runtime_error("LocalExecutor: command failed (exit " + std::to_string(result.exit_code) + "): " + cmd.get());
            }
            submitted = true;
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
                    if (line.find("$(mdrun_cmd)") != std::string::npos) {
                        out += gmx::env_string() + "\n"; // add env
                        line.replace(line.find("$(mdrun_cmd)"), 12, "cd " + folder.absolute_path() + "; " + std::string(run_cmd));
                        found = true;
                    } else if (line.find("$(output)") != std::string::npos) {
                        line.replace(line.find("$(output)"), 9, folder.absolute_path());
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

    /**
     * @brief Like TemplateExecutor but also substitutes $(n_replicas) in the
     *        script, enabling HPC templates to use it in SBATCH headers
     *        (e.g. #SBATCH --ntasks=$(n_replicas)).
     */
    template<MDResult T>
    struct MultiTemplateExecutor : public TemplateExecutor<T> {
        MultiTemplateExecutor(const io::Folder& folder, const SHFile& template_file, std::string_view run_cmd, int n_replicas) : TemplateExecutor<T>(folder) {
            if (!template_file.exists()) {
                throw except::io_error("MultiTemplateExecutor: The template file does not exist.");
            }

            SHFile::operator=(folder + "/job.sh");
            std::string out;
            {
                std::ifstream in(template_file.absolute_path());
                std::string line;
                bool found = false;
                while (std::getline(in, line) && !found) {
                    if (line.find("$(mdrun_cmd)") != std::string::npos) {
                        out += gmx::env_string() + "\n"; // add env
                        line.replace(line.find("$(mdrun_cmd)"), 12, "cd " + folder.absolute_path() + "; " + std::string(run_cmd));
                        found = true;
                    } else if (line.find("$(output)") != std::string::npos) {
                        line.replace(line.find("$(output)"), 9, folder.absolute_path());
                    } else if (line.find("$(n_replicas)") != std::string::npos) {
                        line.replace(line.find("$(n_replicas)"), 13, std::to_string(n_replicas));
                    }
                    out += line + "\n";
                }
            }
            SHFile::create(out);
        }
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
                    if (line.find("$(mdrun_cmd)") != std::string::npos) {
                        out += gmx::env_string() + "\n"; // add env
                        line.replace(line.find("$(mdrun_cmd)"), 12, "cd " + folder.absolute_path() + "; " + std::string(run_cmd));
                        found = true;
                    } else if (line.find("$(jobname)") != std::string::npos) {
                        line.replace(line.find("$(jobname)"), 10, jobname);
                    } else if (line.find("$(output)") != std::string::npos) {
                        line.replace(line.find("$(output)"), 9, folder.absolute_path());
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

    /**
     * @brief Like SlurmRunner but also substitutes $(n_replicas) everywhere in
     *        the template script (e.g. for #SBATCH --ntasks=$(n_replicas)).
     */
    template<MDResult T>
    struct MultiSlurmRunner : public TemplateExecutor<T> {
        MultiSlurmRunner(const io::Folder& folder, const SHFile& template_file,
                         std::string_view run_cmd, int n_replicas, std::string_view jobname)
            : TemplateExecutor<T>(folder)
        {
            if (!template_file.exists()) {
                throw except::io_error("MultiSlurmRunner: The template file does not exist.");
            }

            SHFile::operator=(folder + "/job.sh");
            std::string out;
            {
                std::ifstream in(template_file.absolute_path());
                std::string line;
                bool found = false;
                while (std::getline(in, line) && !found) {
                    if (line.find("$(mdrun_cmd)") != std::string::npos) {
                        out += gmx::env_string() + "\n"; // add env
                        line.replace(line.find("$(mdrun_cmd)"), 12, "cd " + folder.absolute_path() + "; " + std::string(run_cmd));
                        found = true;
                    } else if (line.find("$(jobname)") != std::string::npos) {
                        line.replace(line.find("$(jobname)"), 10, jobname);
                    } else if (line.find("$(output)") != std::string::npos) {
                        line.replace(line.find("$(output)"), 9, folder.absolute_path()); 
                    } else if (line.find("$(n_replicas)") != std::string::npos) {
                        line.replace(line.find("$(n_replicas)"), 13, std::to_string(n_replicas));
                    }
                    out += line + "\n";
                }
            }
            SHFile::create(out);
        }

        void submit() override {
            if (!SHFile::exists()) {throw except::io_error("MultiSlurmRunner::submit: The script does not exist.");}
            auto result = shell::Command("sbatch " + this->absolute_path()).execute();
            if (result.exit_code != 0) {
                throw std::runtime_error("MultiSlurmRunner::submit: Error submitting job: \"" + result.out + "\"");
            }
            for (auto& c : result.out) {
                if (isdigit(c)) {jobid += c;}
            }
            std::cout << "Submitted job " << jobid << std::endl;
        }

        void wait() override {
            int sleep_time = 5;
            if (jobid.empty()) {submit();}
            while (true) {
                auto result = shell::Command("squeue").execute();
                if (result.out.find(jobid) == std::string::npos) {break;}
                char timestamp[50]{0};
                std::time_t time = std::time(nullptr);
                std::strftime(timestamp, 30, "[%H:%M:%S]", std::localtime(&time));
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
            virtual std::unique_ptr<Executor<MultiMDRunResult>> multi_md_runner(const io::Folder& output, std::string_view run_cmd) = 0;
        };

        struct local : type {
            static std::unique_ptr<type> construct() {return std::make_unique<local>();}
            std::unique_ptr<Executor<MDRunResult>> md_runner(const io::Folder& output, std::string_view run_cmd) override {
                return std::make_unique<LocalExecutor<MDRunResult>>(output, run_cmd);
            }
            std::unique_ptr<Executor<SAXSRunResult>> saxs_runner(const io::Folder& output, std::string_view run_cmd) override {
                return std::make_unique<LocalExecutor<SAXSRunResult>>(output, run_cmd);
            }
            std::unique_ptr<Executor<MultiMDRunResult>> multi_md_runner(const io::Folder& output, std::string_view run_cmd) override {
                return std::make_unique<LocalExecutor<MultiMDRunResult>>(output, run_cmd);
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
            std::unique_ptr<Executor<MultiMDRunResult>> multi_md_runner(const io::Folder& output, std::string_view) override {
                return std::make_unique<NoExecutor<MultiMDRunResult>>(output);
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
            std::unique_ptr<Executor<MultiMDRunResult>> multi_md_runner(const io::Folder& output, std::string_view run_cmd) override {
                return std::make_unique<TemplateExecutor<MultiMDRunResult>>(output, template_file, run_cmd);
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
            std::unique_ptr<Executor<MultiMDRunResult>> multi_md_runner(const io::Folder& output, std::string_view run_cmd) override {
                return std::make_unique<SlurmRunner<MultiMDRunResult>>(output, template_file, run_cmd, jobname);
            }
            SHFile template_file;
            std::string jobname;
        };

        /**
         * @brief Executor for multi-replica (PLUMED -multidir) runs.
         *
         * Wraps local, template, or Slurm execution and automatically:
         *  - Prepends "mpirun -np N" to the mdrun command in multi_md_runner.
         *  - Substitutes $(n_replicas) throughout HPC template scripts so that
         *    headers like "#SBATCH --ntasks=$(n_replicas)" work as expected.
         *
         * Usage:
         *   executor::multi::local(N)
         *   executor::multi::templated(N, template_file)
         *   executor::multi::slurm(N, template_file, jobname)
         */
        struct multi : type {
            enum class backend { local_be, templated_be, slurm_be };

            multi(int n_replicas, backend be, SHFile template_file = {}, std::string jobname = {})
                : n_replicas(n_replicas), be(be), template_file(std::move(template_file)), jobname(std::move(jobname)) 
            {}

            static std::unique_ptr<type> local(int n_replicas) {
                return std::make_unique<multi>(n_replicas, backend::local_be);
            }
            static std::unique_ptr<type> templated(int n_replicas, const SHFile& template_file) {
                return std::make_unique<multi>(n_replicas, backend::templated_be, template_file);
            }
            static std::unique_ptr<type> slurm(int n_replicas, const SHFile& template_file, std::string_view jobname) {
                return std::make_unique<multi>(n_replicas, backend::slurm_be, template_file, std::string(jobname));
            }

            // Single-replica paths pass straight through to local/template/slurm without mpirun.
            std::unique_ptr<Executor<MDRunResult>> md_runner(const io::Folder& output, std::string_view run_cmd) override {
                switch (be) {
                    case backend::local_be:     return std::make_unique<LocalExecutor<MDRunResult>>(output, run_cmd);
                    case backend::templated_be: return std::make_unique<TemplateExecutor<MDRunResult>>(output, template_file, run_cmd);
                    case backend::slurm_be:     return std::make_unique<SlurmRunner<MDRunResult>>(output, template_file, run_cmd, jobname);
                }
                throw std::logic_error("executor::multi::md_runner: unknown backend");
            }
            std::unique_ptr<Executor<SAXSRunResult>> saxs_runner(const io::Folder& output, std::string_view run_cmd) override {
                switch (be) {
                    case backend::local_be:     return std::make_unique<LocalExecutor<SAXSRunResult>>(output, run_cmd);
                    case backend::templated_be: return std::make_unique<TemplateExecutor<SAXSRunResult>>(output, template_file, run_cmd);
                    case backend::slurm_be:     return std::make_unique<SlurmRunner<SAXSRunResult>>(output, template_file, run_cmd, jobname);
                }
                throw std::logic_error("executor::multi::saxs_runner: unknown backend");
            }

            // Multi-replica path: prefix PLUMED_KERNEL env var (if configured) and mpirun -np N.
            std::unique_ptr<Executor<MultiMDRunResult>> multi_md_runner(const io::Folder& output, std::string_view run_cmd) override {
                std::string mpi_cmd;
                if (!settings::md::plumed_kernel.empty()) {
                    mpi_cmd = "PLUMED_KERNEL=" + settings::md::plumed_kernel + " ";
                }
                mpi_cmd += "mpirun -np " + std::to_string(n_replicas) + " " + std::string(run_cmd);
                switch (be) {
                    case backend::local_be:     return std::make_unique<LocalExecutor<MultiMDRunResult>>(output, mpi_cmd);
                    case backend::templated_be: return std::make_unique<MultiTemplateExecutor<MultiMDRunResult>>(output, template_file, mpi_cmd, n_replicas);
                    case backend::slurm_be:     return std::make_unique<MultiSlurmRunner<MultiMDRunResult>>(output, template_file, mpi_cmd, n_replicas, jobname);
                }
                throw std::logic_error("executor::multi::multi_md_runner: unknown backend");
            }

            int n_replicas;
            backend be;
            SHFile template_file;
            std::string jobname;
        };
    }
}
