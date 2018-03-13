#include <vector>
#include <iostream>
#include <cmath>

#include "timestepper.hpp"

int main()
{
    constexpr auto F = [](double const &x) {
        return x - x * x * x;
    };
    constexpr auto G = [](double const &x) {
        constexpr double s = sqrt(0.1);
        return s;
    };
    auto dist_fun = [](double const &x) {
        const double y = (x - 1) / 2;
        const double z = (x + 1) / 2;
        constexpr double f1 = 0.25;
        constexpr double f2 = 0.75;
        return f1 - f1 * exp(-8.0 * z * z) + f2 * exp(-8.0 * y * y);
    };

    TimeStepper<double> timestepper(F, G, dist_fun, 0.05);

    constexpr double tmax = 10.0;
    constexpr double dt = 0.001;

    double x = timestepper.transient(-1, dt, tmax);
    std::cout << x << std::endl;
    x = timestepper.stoch_transient(-1, dt, tmax);
    std::cout << x << std::endl;

    int ntrans = 0;
    constexpr int numexp = 100000;
    for (int i = 0; i < numexp; i++)
    {
        double t = timestepper.stoch_transient_max_distance(-1, dt, tmax, 1);
        if (t > 0)
            ntrans++;
    }

    std::cout << "Transition probability direct: " << (double)ntrans / (double)numexp << std::endl;
}
