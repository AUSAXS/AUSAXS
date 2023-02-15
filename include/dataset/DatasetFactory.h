#pragma once

#include <utility/Dataset.h>
#include <utility/SimpleDataset.h>
#include <utility/Dataset2D.h>

namespace factory {
    static std::shared_ptr<SimpleDataset> create(std::string filename) {
        std::shared_ptr<Dataset> data = std::make_shared<Dataset>(filename);
        if (data->M == 3) {
            return std::static_pointer_cast<SimpleDataset>(data);
        } else if (data->M == 4) {
            return std::static_pointer_cast<Dataset2D>(std::move(data));
        } else {
            throw except::invalid_operation("factory::create: Dataset has wrong number of columns.");
        }
    }
}