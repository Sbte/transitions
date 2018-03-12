#ifndef TIMESTEPPER_HPP
#define TIMESTEPPER_HPP

#include <vector>
#include <functional>
#include <random>
#include <algorithm>

template<class T>
class TimeStepper
{
    std::function<T(T const &)> F_;
    std::function<T(T const &)> G_;

    std::function<double ()> generator_;

    double rho_;
    double cdist_;

public:
    TimeStepper(std::function<T(T const &)> F);
    TimeStepper(std::function<T(T const &)> F,
                std::function<T(T const &)> G);

    T transient(T x, double dt, double tmax);
    T stoch_time_step(T const &x, double dt);
    T stoch_transient(T x, double dt, double tmax);
};

template<class T>
TimeStepper<T>::TimeStepper(std::function<T(T const &)> F)
    :
    F_(F)
{}

template<class T>
TimeStepper<T>::TimeStepper(std::function<T(T const &)> F,
                            std::function<T(T const &)> G)
    :
    F_(F),
    G_(G)
{
    std::random_device rd;
    std::default_random_engine engine(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);
    generator_ = std::bind(distribution, std::ref(engine));
}

template<class T>
T TimeStepper<T>::stoch_transient(T x, double dt, double tmax)
{
    int end = tmax / dt;
    for (int i = 0; i < end; i++)
        x = stoch_time_step(x, dt);
    return x;
}

#endif
