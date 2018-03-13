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
    std::function<T(T const &)> dist_fun_;

    std::function<double ()> generator_;

    double rho_;

public:
    TimeStepper(std::function<T(T const &)> F);
    TimeStepper(std::function<T(T const &)> F,
                std::function<T(T const &)> G,
                std::function<T(T const &)> dist_fun,
                double rho);

    T transient(T x, double dt, double tmax);
    T stoch_time_step(T const &x, double dt);
    T stoch_transient(T x, double dt, double tmax);
    double stoch_transient_max_distance(
        T x, double dt, double tmax, double max_distance);
};

template<class T>
TimeStepper<T>::TimeStepper(std::function<T(T const &)> F)
    :
    F_(F)
{}

template<class T>
TimeStepper<T>::TimeStepper(std::function<T(T const &)> F,
                            std::function<T(T const &)> G,
                            std::function<T(T const &)> dist_fun,
                            double rho)
    :
    F_(F),
    G_(G),
    dist_fun_(dist_fun),
    rho_(rho)
{
    std::random_device rd;
    std::default_random_engine engine(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);
    generator_ = std::bind(distribution, std::ref(engine));
}

template<class T>
T TimeStepper<T>::stoch_transient(T x, double dt, double tmax)
{
    for (double t = 0; t < tmax; t += dt)
        x = stoch_time_step(x, dt);
    return x;
}

template<class T>
double TimeStepper<T>::stoch_transient_max_distance(
    T x, double dt, double tmax, double max_distance)
{
    const double lim = max_distance - rho_;
    for (double t = 0; t < tmax; t += dt)
    {
        x = stoch_time_step(x, dt);
        if (dist_fun_(x) > lim)
            return t;
    }
    return -1;
}

#endif
