#pragma once

namespace hist::detail {
    /**
     * @brief A small struct which keeps track of which distances has been added to it. 
     *        This is useful for later using weighted bins.
     */
    template<typename T>
    class WeightedEntry {
        public:
            WeightedEntry() = default;
            WeightedEntry(T value) : distance(0), value(value) {}

            void add(float distance, T value) {
                this->distance += distance;
                this->value += value;
                count++;
            }

            T distance;
            T value;
            T count;
            // union {
            //     struct {
            //         T distance;
            //         T value;
            //     };
            //     T data[2];
            // }
    };
}