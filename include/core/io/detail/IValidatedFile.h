#pragma once

#include <io/File.h>
#include <utility/observer_ptr.h>

namespace io::detail {
    namespace {
        template<typename T> concept file_validation_t = requires(observer_ptr<File> file) {
            {T::validate(file)};
        };
    }

    /**
     * @brief Virtual interface for a file that must be validated upon construction.
     */
    template<file_validation_t F>
    class IValidatedFile : public File {
        public:
            IValidatedFile() = default;
            IValidatedFile(const IValidatedFile&) = default;
            IValidatedFile(IValidatedFile&&) noexcept = default;
            IValidatedFile& operator=(const IValidatedFile&) = default;
            IValidatedFile& operator=(IValidatedFile&&) noexcept = default;
            virtual ~IValidatedFile() = default;

            IValidatedFile(const File& file) : File(file) {F::validate(this);}
            IValidatedFile(File&& file) : File(std::move(file)) {F::validate(this);}

            template<::detail::string_like T>
            IValidatedFile(const T& path) : IValidatedFile(std::string_view(path)) {}
            IValidatedFile(std::string_view path) : File(path) {F::validate(this);}

            template<::detail::string_like T>
            IValidatedFile& operator=(const T& path) {return *this = std::string_view(path);}
            IValidatedFile& operator=(std::string_view path) {File::operator=(path); F::validate(this); return *this;}
    };
}