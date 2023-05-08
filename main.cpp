#include <iostream>
#include <iomanip>
#include <math.h>
//#include "LinearShooting.cpp"
#include "LinearShootingClass.h"
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
    //double w[11][2];
    //LinearShooting(1.0, 2.0, 1.0, 2.0, 10, f, w);
    
    LinearShooting EDOb(1.0, 2.0, 1.0, 2.0, 10, f);
    double **w = EDOb.solutions();
    
    cout << "# x" << setw(8) << "y" << endl;
    for (int i = 0; i < 11; i++)
    {
        cout << w[0][i] << setw(8) << w[1][i] << endl;
    }
    return 0;
}