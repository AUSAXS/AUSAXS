#include <dataset/DatasetFactory.h>

#include <dataset/detail/DatasetConstructor.h>
#include <dataset/detail/DATConstructor.h>
#include <dataset/detail/XVGConstructor.h>
#include <utility/Utility.h>
#include <utility/Constants.h>

std::shared_ptr<Dataset> factory::DatasetFactory::construct(std::string filename) {
    std::unique_ptr<detail::DatasetConstructor> constructor;
    auto ext = utility::extension(filename);
    if (ext == ".dat" || ext == ".txt" || ext == ".rsr") {
        constructor = std::make_unique<detail::DATConstructor>();
    } else if (ext == ".xvg") {
        constructor = std::make_unique<detail::XVGConstructor>();
    } else {
        throw except::invalid_operation("factory::create: Unknown file extension \"" + ext + "\".");
    }
    return constructor->construct(filename);
}