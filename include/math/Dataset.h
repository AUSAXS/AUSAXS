#include <math/Matrix.h>
#include <plots/PlotOptions.h>

class Dataset : Matrix<double> {
    public: 
        Dataset() noexcept : Matrix() {}

};

class SAXSDataset : Matrix<double> {
    public: 
        /**
         * @brief Construct a new empty dataset.
         */
        SAXSDataset() noexcept : Matrix() {
            initialize();
        }

        /**
         * @brief Construct a new dataset with x, y, and yerr values. 
         */
        SAXSDataset(std::vector<double> x, std::vector<double> y, std::vector<double> yerr) : Matrix(x.size(), 3) {
            for (unsigned int i = 0; i < x.size(); i++) {
                row(i) = {x[i], y[i], yerr[i]};
            }
            initialize();
        }

        /**
         * @brief Destructor.
         */
        ~SAXSDataset() override = default;


        const ConstColumn<double> x() const {return ConstColumn<double>(data, N, M, 0);}
        Column<double> x() {return Column<double>(data, N, M, 0);}

        const ConstColumn<double> y() const {return ConstColumn<double>(data, N, M, 1);}
        Column<double> y() {return Column<double>(data, N, M, 1);}

        const ConstColumn<double> yerr() const {return ConstColumn<double>(data, N, M, 2);}
        Column<double> yerr() {return Column<double>(data, N, M, 2);}

        /**
         * @brief Add a new point at the end of the dataset.
         */
        void push_back(double x, double y, double yerr) {
            extend(1);
            row(N) = {x, y, yerr};
        }

    private: 
        plots::PlotOptions options;

        void initialize() {
            options.xlabel = "q";
            options.ylabel = "I";
        }
};

void test() {
    SAXSDataset data; 
    auto y = data.y();
}