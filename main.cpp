#include <iostream>
#include <iomanip>
#include <math.h>
#include "LinearShooting.cpp"
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

int main()
{
    double (*f[3])(double) = {p, q, r};
    double w[11][2];

    LinearShooting(1.0, 2.0, 1.0, 2.0, 10, f, w);

    cout << "# x" << setw(8) << "y" << endl;
    for (int i = 0; i < 11; i++)
    {
        cout << w[i][0] << setw(8) << w[i][1] << endl;
    }
    return 0;
}