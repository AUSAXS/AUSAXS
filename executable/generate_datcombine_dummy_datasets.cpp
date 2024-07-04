#include <dataset/SimpleDataset.h>

int main(int argc, char const *argv[]) {
    SimpleDataset data("data/SASDJG5/SASDJG5.dat");
    data.save("output/datcombine/1.dat");

    for (auto& y : data.y()) {
        y = y + 1e-9;
    }
    data.save("output/datcombine/2.dat");
}