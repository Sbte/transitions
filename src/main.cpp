#include <vector>
#include <iostream>
#include <cmath>

#include "timestepper.hpp"

int main()
{
    auto F = [](double const &x) {
        return x - x * x * x;
    };
    auto G = [](double const &x) {
        return sqrt(0.1);
    };

    TimeStepper<double> timestepper(F, G);

    double tmax = 2.0;
    double dt = 0.01;

    double x = timestepper.transient(-1, dt, tmax);
    std::cout << x << std::endl;
    x = timestepper.stoch_transient(-1, dt, tmax);
    std::cout << x << std::endl;
}
