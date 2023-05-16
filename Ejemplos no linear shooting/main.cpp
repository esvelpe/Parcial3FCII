#include <iostream>
#include <iomanip>
#include <math.h>
// #include "LinearShooting.cpp"
#include "ShootingClass.h"
using namespace std;

double p(double x)
{
    if (x != 0)
    {
        return -2 / x;
    }
    else
    {
        exit(1);
    }
}

double q(double x)
{
    if (x != 0)
    {
        return 2 / (x * x);
    }
    else
    {
        exit(1);
    }
}

double r(double x)
{
    if (x != 0)
    {
        return sin(log(x)) / (x * x);
    }
    else
    {
        exit(1);
    }
}

double fun(double x, double y, double yp)
{
    return -yp*yp - y + log(x);
}

double fy(double x, double y, double yp)
{
    return -1;
}

double fyp(double x, double y, double yp)
{
    return -2*yp;
}

int main()
{
    // double (*f[3])(double) = {p, q, r};
    //  double w[11][2];
    //  LinearShooting(1.0, 2.0, 1.0, 2.0, 10, f, w);
    double (*f_no[3])(double, double, double) = {fun, fy, fyp};

    // Shooting EDOb(1.0, 2.0, 1.0, 2.0, 10, f);
    // double **w = EDOb.linear_solutions();

    Shooting EDO_No_Lin(1, 2, 0.0, log(2), 20, f_no, 1e-5, 10);
    double **w_no_lin = EDO_No_Lin.no_linear_solutions();

    cout << "# x" << setw(12) << "y" << endl;
    for (int i = 0; i < 21; i++)
    {
        cout << w_no_lin[0][i] << setw(12) << w_no_lin[1][i] << endl;
    }
    return 0;
}
