#include <utility/Observable.h>

using namespace utility;

template<typename T>
void Observable<T>::attach_observer(observer_ptr<Observer<T>> o) {
    observers.push_back(o);
}

template<typename T>
void Observable<T>::detach_observer(observer_ptr<Observer<T>> o) {
    observers.erase(std::remove(observers.begin(), observers.end(), &o), observers.end());
}

template<typename T>
void Observable<T>::notify(const T& value) {
    for (auto& o : observers) {
        o->notify(value);
    }
}