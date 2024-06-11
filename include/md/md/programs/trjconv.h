#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>

namespace md {
    class trjconv : private gmx {
        public: 
            trjconv() {
                cmd.append("trjconv");
            }

            trjconv(const XTCFile& xtc) : trjconv() {
                input(xtc);
            }

            trjconv& input(const XTCFile& xtc) {
                options.push_back(std::make_shared<shell::Argument>("-f", xtc));
                return *this;
            }

            trjconv& output(const XTCFile& xtc) {
                this->xtc = xtc;
                options.push_back(std::make_shared<shell::Argument>("-o", xtc));
                return *this;
            }

            trjconv& startframe(int frame) {
                options.push_back(std::make_shared<shell::Argument>("-b", frame));
                return *this;
            }

            trjconv& endframe(int frame) {
                options.push_back(std::make_shared<shell::Argument>("-e", frame));
                return *this;
            }

            trjconv& skip_every_n_frame(unsigned int n) {
                options.push_back(std::make_shared<shell::Argument>("-skip", n));
                return *this;
            }

            trjconv& center() {
                options.push_back(std::make_shared<shell::Flag>("-center"));
                return *this;
            }

            trjconv& boxcenter(const std::string& boxcenter) {
                options.push_back(std::make_shared<shell::Argument>("-boxcenter", boxcenter));
                return *this;
            }

            trjconv& group(const std::string& group) {
                cmd.prepend("echo " + group + " 0 |");
                return *this;
            }

            trjconv& pbc(const std::string& pbc) {
                options.push_back(std::make_shared<shell::Argument>("-pbc", pbc));
                return *this;
            }

            trjconv& ur(const std::string& ur) {
                options.push_back(std::make_shared<shell::Argument>("-ur", ur));
                return *this;
            }

            trjconv& index(const NDXFile& ndx) {
                options.push_back(std::make_shared<shell::Argument>("-n", ndx));
                return *this;
            }

            trjconv& runfile(const TPRFile& tpr) {
                options.push_back(std::make_shared<shell::Argument>("-s", tpr));
                return *this;
            }

            std::tuple<XTCFile> run() {
                gmx::execute();
                return std::make_tuple(xtc);
            }

        private: 
            XTCFile xtc;
            void validate() const override {}
    };
}