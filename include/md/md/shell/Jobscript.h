#pragma once

#include <md/shell/Command.h>
#include <md/utility/Exceptions.h>

#include <time.h>
#include <thread>
#include <chrono>

namespace shell {
    template<typename T>
    class Jobscript {
        public: 
            /**
             * @brief Create a new jobscript without submitting it.
             */
            virtual void create() {
                cmd.execute();
                created = true;
            }

            /**
             * @brief Submit the jobscript to the scheduler.
             */
            virtual void submit() = 0;

            /**
             * @brief Wait for the jobscript to finish. 
             */
            virtual void wait() = 0;

            virtual T result() = 0;

        protected:
            Command cmd;
            std::string filename;
            bool created = false;
    };

    template<typename T>
    class Slurmscript : public Jobscript<T> {
        public: 
            void submit() override {
                if (!this->created) {this->create();}
                if (this->filename.empty()) {throw except::unexpected("Jobscript: no file to submit!");}
                auto result = shell::Command("sbatch " + this->filename).execute();
                if (result.exit_code != 0) {
                    throw std::runtime_error("Error submitting job: " + result.out);
                }

                for (auto& c : result.out) {
                    if (isdigit(c)) {
                        jobid += c;
                    }
                }

                std::cout << "Submitted job " << jobid << std::endl;
            }

            /**
             * @brief Blocking wait for the job to finish. 
             * 
             * This is done by checking the output of `squeue` every 2 minutes. 
             */
            void wait() override {
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
                    std::this_thread::sleep_for(std::chrono::seconds(120));
                }
            }
        
        private: 
            std::string jobid;
    };
}