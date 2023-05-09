#include <iostream>
#include <iomanip>
#include <math.h>
//#include "LinearShooting.cpp"
#include "LinearShootingClass.h"
using namespace std;

double p(double x)
{
    if (x >= 0)
    {
        return 1;
    }
    else
    {
        exit(1);
    }
}

double q(double x)
{
    if (x >= 0)
    {
        return 2;
    }
    else
    {
        exit(1);
    }
}

double r(double x)
{
    if (x >= 0)
    {
        return cos(x);
    }
    else
    {
        exit(1);
    }
}

int main()
{
    double (*f[3])(double) = {p, q, r};
    
    LinearShooting EDOb(0.0, 3.14159265/2, -0.3, -0.1, 20, f); //double a, double b, double alpha, double beta, unsigned int N, double (*g[3])(double)
    double **w = EDOb.solutions();
    
    cout << "# x" << " " << "y" << endl;
    for (int i = 0; i < 21; i++)
    {
        cout << w[0][i] << " " << w[1][i] << endl;
    }
    return 0;
}