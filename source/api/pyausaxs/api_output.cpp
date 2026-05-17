// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_output.h>

#include <memory>
#include <iostream>
#include <streambuf>
#include <string>

namespace {
    /**
     * A std::streambuf that forwards output to a C callback, flushing on each
     * newline and on sync() (called by std::endl and std::flush).
     */
    class CCallbackStreamBuf : public std::streambuf {
        public:
            explicit CCallbackStreamBuf(ausaxs_output_cb cb, std::streambuf* original)
                : cb_(cb), original_(original) {}

            std::streambuf* original() const { return original_; }

        protected:
            int overflow(int c) override {
                if (c != EOF) {
                    char ch = static_cast<char>(c);
                    xsputn(&ch, 1);
                }
                return c;
            }

            std::streamsize xsputn(const char* s, std::streamsize n) override {
                for (std::streamsize i = 0; i < n; ++i) {
                    buf_ += s[i];
                    if (s[i] == '\n') {
                        cb_(buf_.c_str(), static_cast<int>(buf_.size()));
                        buf_.clear();
                    }
                }
                return n;
            }

            int sync() override {
                if (!buf_.empty()) {
                    cb_(buf_.c_str(), static_cast<int>(buf_.size()));
                    buf_.clear();
                }
                return 0;
            }

        private:
            ausaxs_output_cb cb_;
            std::string buf_;
            std::streambuf* original_;
    };

    std::unique_ptr<CCallbackStreamBuf> g_current = nullptr;
}

void set_output_callback(ausaxs_output_cb cb) {
    // Always restore any previously installed callback first.
    if (g_current) {
        std::cout.rdbuf(g_current->original());
        g_current = nullptr;
    }
    if (cb) {
        g_current = std::make_unique<CCallbackStreamBuf>(cb, std::cout.rdbuf());
        std::cout.rdbuf(g_current.get());
    }
}

void reset_output_callback() {
    set_output_callback(nullptr);
}