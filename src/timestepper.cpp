#include "timestepper.hpp"

#include <iostream>

extern "C"
{
    void daxpy_(int *n, double *a, double const *x, int *incx,
                double *y, int *incy);
    void dgemv_(char *transa, int *m, int *n,
                double *alpha, double const *a, int *lda,
                double const *x, int *incx, double *beta,
                double *y, int *incy);
    double dnrm2_(int *n, double const *x, int *incx);
}

namespace BlasWrapper
{

void DAXPY(int n, double a, double const *x, double *y)
{
    int inc = 1;
    daxpy_(&n, &a, x, &inc, y, &inc);
}

void DAXPY(double a, std::vector<double> const &x, std::vector<double> &y)
{
    int n = x.size();
    DAXPY(n, a, &x[0], &y[0]);
}

double DNRM2(int n, double const *x, int incx)
{
    return dnrm2_(&n, x, &incx);
}

double DNRM2(std::vector<double> const &x)
{
    return DNRM2(x.size(), &x[0], 1);
}

void DGEMV(char trans, int m, int n,
           double alpha, double const *a, int lda,
           double const *x, int incx, double beta,
           double *y, int incy)
{
    dgemv_(&trans, &m, &n, &alpha,
           a, &lda, x, &incx, &beta, y, &incy);
}

void DGEMV(std::vector<double> const &A, std::vector<double> const &x, std::vector<double> &y)
{
    int n = x.size();
    int m = A.size() / n;
    DGEMV('N', m, n, 1.0, &A[0], m, &x[0], 1, 0.0, &y[0], 1);
}

void DGEMV(std::vector<double> const &A, std::vector<double> const &x,
           double beta, std::vector<double> &y)
{
    int n = x.size();
    int m = A.size() / n;
    DGEMV('N', m, n, 1.0, &A[0], m, &x[0], 1, beta, &y[0], 1);
}

}

void dump(std::vector<double> const &x)
{
    for (double i: x)
        std::cout << i << ", ";
    std::cout << std::endl;
}

template<>
std::vector<double> TimeStepper<std::vector<double> >::transient(
    std::vector<double> x,
    double dt, double tmax)
{
    int end = tmax / dt;
    for (int i = 0; i < end; i++)
    {
        std::vector<double> y = F_(x);
        BlasWrapper::DAXPY(dt, y, x);
    }
    return x;
}

template<>
std::vector<double> TimeStepper<std::vector<double> >::stoch_time_step(
    std::vector<double> const &x,
    double dt)
{
    std::vector<double> Gx = G_(x);
    std::vector<double> Fx = F_(x);

    int n = x.size();
    int m = Gx.size() / n;

    std::vector<double> pert(m);
    std::generate(pert.begin(), pert.end(), generator_);

    std::vector<double> y(x);
    BlasWrapper::DAXPY(dt, Fx, y);
    BlasWrapper::DGEMV(Gx, pert, 1.0, y);

    return y;
}

template<>
std::vector<double> TimeStepper<std::vector<double> >::stoch_transient(
    std::vector<double> x,
    double dt, double tmax)
{
    int end = tmax / dt;
    for (int i = 0; i < end; i++)
    {
        x = stoch_time_step(x, dt);
    }
    return x;
}

template<>
double TimeStepper<double>::transient(
    double x,
    double dt, double tmax)
{
    int end = tmax / dt;
    for (int i = 0; i < end; i++)
        x = F_(x) * dt + x;
    return x;
}

template<>
double TimeStepper<double>::stoch_time_step(
    double const &x,
    double dt)
{
    double Gx = G_(x);
    double Fx = F_(x);
    double pert = sqrt(dt) * generator_();

    return  Fx * dt + Gx * pert;
}

template<>
double TimeStepper<double>::stoch_transient(
    double x,
    double dt, double tmax)
{
    int end = tmax / dt;
    for (int i = 0; i < end; i++)
        x = stoch_time_step(x, dt);
    return x;
}
