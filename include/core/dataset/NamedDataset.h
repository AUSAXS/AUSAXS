#include <dataset/NamedWrapper.h>
#include <dataset/Dataset.h>
#include <dataset/Dataset2D.h>
#include <dataset/SimpleDataset.h>

namespace ausaxs {
    using NamedDataset = NamedWrapper<Dataset>;
    using NamedSimpleDataset = NamedWrapper<Dataset>;
    using NamedDataset2D = NamedWrapper<Dataset2D>;
}