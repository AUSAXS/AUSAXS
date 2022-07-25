#pragma once

#include <vector>

template<typename T>
class View {
    static_assert(std::is_base_of<detail::IPoint, T>::value, "T must be a subclass of IPoint");

    public: 
        View(std::vector<T>& data) : data(data) {}
        virtual double& operator[](unsigned int i) = 0;
        virtual const double& operator[](unsigned int i) const = 0;
        unsigned int size() const noexcept {return data.size();}

    protected: 
        std::vector<T>& data;
};

template<typename T>
class ErrorView : public View<T> {
    using View<T>::View;
    static_assert(std::is_base_of<detail::IPointError, T>::value, "T must be a subclass of IPointError");
};

template<typename T>
class XView : public View<T> {
    static_assert(T::dim() >= 1, "Dataset1D can only be used with 1D data or higher.");

    public: 
        XView(std::vector<T>& data) : View<T>(data) {}

        /**
         * @brief Get a value from this view. 
         *        Complexity is O(1).
         */
        double& operator[](unsigned int i) override {return this->data[i].x;}

        /**
         * @brief Get a value from this view. 
         *        Complexity is O(1).
         */
        const double& operator[](unsigned int i) const override {return this->data[i].x;}

        /**
         * @brief Assign a vector to this view. It will be truncated to the size of the vector.
         *        Complexity is O(n).
         */
        XView& operator=(std::initializer_list<double> list) {
            this->data.resize(list.size());
            std::transform(this->begin(), this->end(), list.begin(), this->begin(), [](T& p, double x) {p.x = x;});
            return *this;
        }

        /**
         * @brief Get a vector representation of this view.
         *        Complexity is O(n).
         */
        std::vector<double> as_vector() const {
            std::vector<double> v(this->size());
            for (unsigned int i = 0; i < this->size(); i++) {
                v[i] = this->data[i].x;
            }
            return v;
        }

        /**
         * @brief Get a vector representation of this view.
         *        Complexity is O(n).
         */
        operator std::vector<double>() const {
            return as_vector();
        }
};

template<typename T>
class XErrorView : public ErrorView<T> {
    static_assert(T::dim() >= 1, "Dataset1D can only be used with 1D data or higher.");

    public: 
        XErrorView(std::vector<T>& data) : ErrorView<T>(data) {}

        /**
         * @brief Get a value from this view. 
         *        Complexity is O(1).
         */
        double& operator[](unsigned int i) override {return this->data[i].xerr;}

        /**
         * @brief Get a value from this view. 
         *        Complexity is O(1).
         */
        const double& operator[](unsigned int i) const override {return this->data[i].xerr;}

        /**
         * @brief Assign a vector to this view. It will be truncated to the size of the vector.
         *        Complexity is O(n).
         */
        XErrorView& operator=(std::initializer_list<double> list) {
            this->data.resize(list.size());
            std::transform(this->begin(), this->end(), list.begin(), this->begin(), [](T& p, double xerr) {p.xerr = xerr;});
            return *this;
        }

        /**
         * @brief Get a vector representation of this view.
         *        Complexity is O(n).
         */
        std::vector<double> as_vector() const {
            std::vector<double> v(this->size());
            for (unsigned int i = 0; i < this->size(); i++) {
                v[i] = this->data[i].xerr;
            }
            return v;
        }

        /**
         * @brief Get a vector representation of this view.
         *        Complexity is O(n).
         */
        operator std::vector<double>() const {
            return as_vector();
        }
};

template<typename T>
class YView : public View<T> {
    static_assert(T::dim() >= 2, "Dataset1D can only be used with 2D data or higher.");

    public: 
        YView(std::vector<T>& data) : View<T>(data) {}

        /**
         * @brief Get a value from this view. 
         *        Complexity is O(1).
         */
        double& operator[](unsigned int i) override {return this->data[i].y;}

        /**
         * @brief Get a value from this view. 
         *        Complexity is O(1).
         */
        const double& operator[](unsigned int i) const override {return this->data[i].y;}

        /**
         * @brief Assign a vector to this view. It will be truncated to the size of the vector.
         *        Complexity is O(n).
         */
        YView& operator=(std::initializer_list<double> list) {
            this->data.resize(list.size());
            std::transform(this->begin(), this->end(), list.begin(), this->begin(), [](T& p, double y) {p.y = y;});
            return *this;
        }

        /**
         * @brief Get a vector representation of this view.
         *        Complexity is O(n).
         */
        std::vector<double> as_vector() const {
            std::vector<double> v(this->size());
            for (unsigned int i = 0; i < this->size(); i++) {
                v[i] = this->data[i].y;
            }
            return v;
        }

        /**
         * @brief Get a vector representation of this view.
         *        Complexity is O(n).
         */
        operator std::vector<double>() const {
            return as_vector();
        }
};

template<typename T>
class YErrorView : public ErrorView<T> {
    static_assert(T::dim() >= 2, "Dataset1D can only be used with 1D data or higher.");

    public: 
        YErrorView(std::vector<T>& data) : ErrorView<T>(data) {}

        /**
         * @brief Get a value from this view. 
         *        Complexity is O(1).
         */
        double& operator[](unsigned int i) override {return this->data[i].yerr;}

        /**
         * @brief Get a value from this view. 
         *        Complexity is O(1).
         */
        const double& operator[](unsigned int i) const override {return this->data[i].yerr;}

        /**
         * @brief Assign a vector to this view. It will be truncated to the size of the vector.
         *        Complexity is O(n).
         */
        YErrorView& operator=(std::initializer_list<double> list) {
            this->data.resize(list.size());
            std::transform(this->begin(), this->end(), list.begin(), this->begin(), [](T& p, double yerr) {p.yerr = yerr;});
            return *this;
        }

        /**
         * @brief Get a vector representation of this view.
         */
        std::vector<double> as_vector() const {
            std::vector<double> v(this->size());
            for (unsigned int i = 0; i < this->size(); i++) {
                v[i] = this->data[i].yerr;
            }
            return v;
        }

        /**
         * @brief Get a vector representation of this view.
         */
        operator std::vector<double>() const {
            return as_vector();
        }
};

template<typename T>
class ConstView {
    static_assert(std::is_base_of<detail::IPoint, T>::value, "T must be a subclass of IPoint");

    public: 
        ConstView(const std::vector<T>& data) : data(data) {}
        virtual const double& operator[](unsigned int i) const = 0;
        unsigned int size() const noexcept {return data.size();}

    protected: 
        const std::vector<T>& data;
};

template<typename T>
class ConstErrorView : public ConstView<T> {
    using ConstView<T>::ConstView;
    static_assert(std::is_base_of<detail::IPointError, T>::value, "T must be a subclass of IPointError");
};

template<typename T>
class XConstView : public ConstView<T> {
    static_assert(T::dim() >= 1, "Dataset1D can only be used with 1D data or higher.");

    public: 
        XConstView(const std::vector<T>& data) : View<T>(data) {}

        /**
         * @brief Get a value from this view. 
         *        Complexity is O(1).
         */
        const double& operator[](unsigned int i) const override {return this->data[i].x;}

        /**
         * @brief Get a vector representation of this view.
         *        Complexity is O(n).
         */
        std::vector<double> as_vector() const {
            std::vector<double> v(this->size());
            for (unsigned int i = 0; i < this->size(); i++) {
                v[i] = this->data[i].x;
            }
            return v;
        }

        /**
         * @brief Get a vector representation of this view.
         *        Complexity is O(n).
         */
        operator std::vector<double>() const {
            return as_vector();
        }
};

template<typename T>
class XConstErrorView : public ConstErrorView<T> {
    static_assert(T::dim() >= 1, "Dataset1D can only be used with 1D data or higher.");

    public: 
        XConstErrorView(const std::vector<T>& data) : ErrorView<T>(data) {}

        /**
         * @brief Get a value from this view. 
         *        Complexity is O(1).
         */
        const double& operator[](unsigned int i) const override {return this->data[i].xerr;}

        /**
         * @brief Get a vector representation of this view.
         *        Complexity is O(n).
         */
        std::vector<double> as_vector() const {
            std::vector<double> v(this->size());
            for (unsigned int i = 0; i < this->size(); i++) {
                v[i] = this->data[i].xerr;
            }
            return v;
        }

        /**
         * @brief Get a vector representation of this view.
         *        Complexity is O(n).
         */
        operator std::vector<double>() const {
            return as_vector();
        }
};

template<typename T>
class YConstView : public ConstView<T> {
    static_assert(T::dim() >= 2, "Dataset1D can only be used with 2D data or higher.");

    public: 
        YConstView(const std::vector<T>& data) : View<T>(data) {}

        /**
         * @brief Get a value from this view. 
         *        Complexity is O(1).
         */
        const double& operator[](unsigned int i) const override {return this->data[i].y;}

        /**
         * @brief Get a vector representation of this view.
         *        Complexity is O(n).
         */
        std::vector<double> as_vector() const {
            std::vector<double> v(this->size());
            for (unsigned int i = 0; i < this->size(); i++) {
                v[i] = this->data[i].y;
            }
            return v;
        }

        /**
         * @brief Get a vector representation of this view.
         *        Complexity is O(n).
         */
        operator std::vector<double>() const {
            return as_vector();
        }
};

template<typename T>
class YConstErrorView : public ConstErrorView<T> {
    static_assert(T::dim() >= 2, "Dataset1D can only be used with 1D data or higher.");

    public: 
        YConstErrorView(const std::vector<T>& data) : ErrorView<T>(data) {}

        /**
         * @brief Get a value from this view. 
         *        Complexity is O(1).
         */
        const double& operator[](unsigned int i) const override {return this->data[i].yerr;}

        /**
         * @brief Get a vector representation of this view.
         */
        std::vector<double> as_vector() const {
            std::vector<double> v(this->size());
            for (unsigned int i = 0; i < this->size(); i++) {
                v[i] = this->data[i].yerr;
            }
            return v;
        }

        /**
         * @brief Get a vector representation of this view.
         */
        operator std::vector<double>() const {
            return as_vector();
        }
};