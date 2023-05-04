#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

double p(double x, double (*f)(double x))
{
    return ((f)(x));
}

void LinearShooting(double a, double b, double alpha, double beta, int N, double (*f[3])(double), double w[][2])
{
    double h = (b - a) / (double)N;
    double u1[N] = {alpha};
    double u2[N] = {0.0};
    double v1[N] = {0.0};
    double v2[N] = {1.0};
    double x, w1, w2;
    // Definicion de valores de RK para y
    double k1_1, k2_1, k3_1, k4_1;
    double k1_2, k2_2, k3_2, k4_2;
    // Definicion de valores de RK para y'
    double kp1_1, kp2_1, kp3_1, kp4_1;
    double kp1_2, kp2_2, kp3_2, kp4_2;

    w[0][0] = a;
    w[0][1] = alpha;

    for (int i = 1; i < N + 1; i++)
    {
        x = a + i * h;
        w[i][0] = x;

        // RK1 para K
        k1_1 = h * u2[i - 1];
        k1_2 = h * (u2[i - 1] * (*f[0])(x) + u1[i - 1] * (*f[1])(x) + (*f[2])(x));
        // RK2
        k2_1 = h * (u2[i - 1] + 0.5 * k1_2);
        k2_2 = h * ((*f[0])(x + 0.5 * h) * (u2[i - 1] + 0.5 * k1_2) + (*f[1])(x + 0.5 * h) * (u1[i - 1] + 0.5 * k1_1) + (*f[2])(x + 0.5 * h));
        // RK3
        k3_1 = h * (u2[i - 1] + 0.5 * k2_2);
        k3_2 = h * ((*f[0])(x + 0.5 * h) * (u2[i - 1] + 0.5 * k2_2) + (*f[1])(x + 0.5 * h) * (u1[i - 1] + 0.5 * k2_1) + (*f[2])(x + 0.5 * h));
        // RK4
        k4_1 = h * (u2[i - 1] + k3_2);
        k4_2 = h * ((*f[0])(x + h) * (u2[i - 1] + k3_2) + (*f[1])(x + h) * (u1[i - 1] + k3_1) + (*f[2])(x + h));

        // Redefinimos u
        u1[i] = u1[i - 1] + (k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1) / 6.0;
        u2[i] = u2[i - 1] + (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2) / 6.0;

        // RK1 para kp
        kp1_1 = h * v2[i - 1];
        kp1_2 = h * (v2[i - 1] * (*f[0])(x) + v1[i - 1] * (*f[1])(x));
        // RK2
        kp2_1 = h * (v2[i - 1] + 0.5 * kp1_2);
        kp2_2 = h * ((*f[0])(x + 0.5 * h) * (v2[i - 1] + 0.5 * kp1_2) + (*f[1])(x + 0.5 * h) * (v1[i - 1] + 0.5 * kp1_1));
        // RK3
        kp3_1 = h * (v2[i - 1] + 0.5 * kp2_2);
        kp3_2 = h * ((*f[0])(x + 0.5 * h) * (v2[i - 1] + 0.5 * kp2_2) + (*f[1])(x + 0.5 * h) * (v1[i - 1] + 0.5 * kp2_1));
        // RK4
        kp4_1 = h * (v2[i - 1] + kp3_2);
        kp4_2 = h * ((*f[0])(x + h) * (v2[i - 1] + kp3_2) + (*f[1])(x + h) * (v1[i - 1] + kp3_1));

        // Redefinimos v
        v1[i] = v1[i - 1] + (kp1_1 + 2 * kp2_1 + 2 * kp3_1 + kp4_1) / 6.0;
        v2[i] = v2[i - 1] + (kp1_2 + 2 * kp2_2 + 2 * kp3_2 + kp4_2) / 6.0;
    }

    
    w1 = alpha;
    w2 = (beta - u1[N]) / v1[N];

    for (int i = 1; i < N+1; i++)
    {
        w[i][1] = u1[i] + w2 * v1[i];
    }
}