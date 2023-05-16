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
    return (32.0 + 2.0 * pow(x, 3) - y * yp) / 8.0;
}

double fy(double x, double y, double yp)
{
    return -yp / 8.0;
}

double fyp(double x, double y, double yp)
{
    return -y / 8.0;
}

double real_lineal(double x)
{
    double c2 = (8.0 - 12.0 * sin(log(2)) - 4.0 * cos(log(2))) / 70.0;
    double c1 = 11.0 / 10.0 - c2;
    return c1 * x + c2 * pow(x, -2) - 3 * sin(log(x)) / 10.0 - cos(log(x)) / 10.0;
}

double real_no_lineal(double x)
{
    return pow(x, 2) + 16.0 / x;
}

int main()
{
    double (*f[3])(double) = {p, q, r};
    //  double w[11][2];
    //  LinearShooting(1.0, 2.0, 1.0, 2.0, 10, f, w);
    double (*f_no[3])(double, double, double) = {fun, fy, fyp};

    Shooting EDOb(1.0, 2.0, 1.0, 2.0, 10, f);
    double **w = EDOb.linear_solutions();
    vector<double> real_lin;
    for (int i = 0; i < EDOb.getN(); i++)
    {
        real_lin.push_back(real_lineal(w[0][i]));
    }
    EDOb.printTable(real_lin);

    Shooting EDO_No_Lin(1.0, 3.0, 17.0, 43.0 / 3.0, 20, f_no, 1e-5, 10);
    double **w_no_lin = EDO_No_Lin.no_linear_solutions();
    vector<double> real_no_lin;
    for (int i = 0; i < EDO_No_Lin.getN(); i++)
    {
        real_no_lin.push_back(real_no_lineal(w_no_lin[0][i]));
    }
    EDO_No_Lin.printTable(real_no_lin);

    // cout << "# x" << setw(12) << "y" << endl;
    // for (int i = 0; i < 11; i++)
    // {
    //     cout << w[0][i] << setw(12) << w[1][i] << endl;
    // }
    return 0;
}