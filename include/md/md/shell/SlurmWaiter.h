#pragma once

#include <md/shell/Command.h>

#include <mutex>
#include <condition_variable>
#include <queue>
#include <cstdlib>
#include <chrono>

namespace shell {
    class SlurmWait {
        public:
            SlurmWait(const std::string& job_id) :
                job_id(job_id), done(false), interrupted(false) {}

            void wait() {
                std::unique_lock<std::mutex> lock(mtx);
                while (!done && !interrupted) {
                    // Check the status of the Slurm job
                    Command cmd("squeue -h -j " + job_id + " -t PD,R,CD,S,ST -o %T");
                    if (cmd.execute().out.find("COMPLETED") != std::string::npos) {
                        done = true;
                    }

                    // If the job is not done and not interrupted, wait for 2 minutes
                    if (!done && !interrupted) {
                        cv.wait_for(lock, std::chrono::minutes(2));
                    }
                }
            }

            void interrupt() {
                std::unique_lock<std::mutex> lock(mtx);
                interrupted = true;
                cv.notify_all();
            }

            std::queue<double> get_results() {
                std::unique_lock<std::mutex> lock(mtx);
                return results;
            }

        private:
            std::string job_id;
            bool done;
            bool interrupted;
            std::mutex mtx;
            std::condition_variable cv;
            std::queue<double> results;
    };
}