#include "ShootingClass.h"
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

Shooting::Shooting(double a, double b, double alpha, double beta, unsigned int N, double (*g[3])(double))
{
    // Inicializa las varibales
    this->a = a;
    this->b = b;
    this->alpha = alpha;
    this->beta = beta;
    this->N = N;
    this->h = (b - a) / (double)N;
    *f_lineal = *g;

    w = new double *[2];
    w[0] = new double[N + 1];
    w[1] = new double[N + 1]; // Inicializa cada una de las columnas con un array de N+1 entradas
    // u1[N] = {alpha};
    // u2[N] = {0.0};
    // v1[N] = {0.0};
    // v2[N] = {1.0};
    u1 = new vector<double>(N);
    u1->push_back(alpha);
    u2 = new vector<double>(N);
    u2->push_back(0.0);
    v1 = new vector<double>(N);
    v1->push_back(0.0);
    v2 = new vector<double>(N);
    v2->push_back(1.0);
}

Shooting::Shooting(double a, double b, double alpha, double beta, unsigned int N, double (*g[3])(double, double, double), double TOL, int N_max)
{
    this->a = a;
    this->b = b;
    this->alpha = alpha;
    this->beta = beta;
    this->N = N;
    this->h = (b - a) / (double)N;
    this->TOL = TOL;
    this->N_max = N_max;
    this->TK = (beta - alpha) / (b - a);
    //*f_no_lineal = *g;
    f_no_lineal[0] = g[0];
    f_no_lineal[1] = g[1];
    f_no_lineal[2] = g[2];

    w = new double *[2];
    w[0] = new double[N + 1];
    w[1] = new double[N + 1]; // Inicializa cada una de las columnas con un array de N+1 entradas
    // u1[N] = {alpha};
    // u2[N] = {TK};
    // v1[N] = {0.0};
    // v2[N] = {1.0};
    u1 = new vector<double>(N+1);
    u1->at(0) = alpha;
    u2 = new vector<double>(N+1);
    u2->at(0) = 0.0;
    w1 = 0.0;
    w2 = 1.0;
    // v1 = new vector<double>(N);
    // v1->push_back(0.0);
    // v2 = new vector<double>(N);
    // v2->push_back(1.0);
}

double **Shooting::linear_solutions()
{

    for (int i = 1; i < N + 1; i++)
    {
        x = a + i * h;
        w[0][i] = x;

        // RK1 para K
        k1_1 = h * u2->at(i - 1);
        k1_2 = h * (u2->at(i - 1) * (*f_lineal[0])(x) + u1->at(i - 1) * (*f_lineal[1])(x) + (*f_lineal[2])(x));
        // RK2
        k2_1 = h * (u2->at(i - 1) + 0.5 * k1_2);
        k2_2 = h * ((*f_lineal[0])(x + 0.5 * h) * (u2->at(i - 1) + 0.5 * k1_2) + (*f_lineal[1])(x + 0.5 * h) * (u1->at(i - 1) + 0.5 * k1_1) + (*f_lineal[2])(x + 0.5 * h));
        // RK3
        k3_1 = h * (u2->at(i - 1) + 0.5 * k2_2);
        k3_2 = h * ((*f_lineal[0])(x + 0.5 * h) * (u2->at(i - 1) + 0.5 * k2_2) + (*f_lineal[1])(x + 0.5 * h) * (u1->at(i - 1) + 0.5 * k2_1) + (*f_lineal[2])(x + 0.5 * h));
        // RK4
        k4_1 = h * (u2->at(i - 1) + k3_2);
        k4_2 = h * ((*f_lineal[0])(x + h) * (u2->at(i - 1) + k3_2) + (*f_lineal[1])(x + h) * (u1->at(i - 1) + k3_1) + (*f_lineal[2])(x + h));

        // Redefinimos u
        u1->at(i) = u1->at(i - 1) + (k1_1 + 2 * k2_1 + 2 * k3_1 + k4_1) / 6.0;
        u2->at(i) = u2->at(i - 1) + (k1_2 + 2 * k2_2 + 2 * k3_2 + k4_2) / 6.0;

        // RK1 para kp
        kp1_1 = h * v2->at(i - 1);
        kp1_2 = h * (v2->at(i - 1) * (*f_lineal[0])(x) + v1->at(i - 1) * (*f_lineal[1])(x));
        // RK2
        kp2_1 = h * (v2->at(i - 1) + 0.5 * kp1_2);
        kp2_2 = h * ((*f_lineal[0])(x + 0.5 * h) * (v2->at(i - 1) + 0.5 * kp1_2) + (*f_lineal[1])(x + 0.5 * h) * (v1->at(i - 1) + 0.5 * kp1_1));
        // RK3
        kp3_1 = h * (v2->at(i - 1) + 0.5 * kp2_2);
        kp3_2 = h * ((*f_lineal[0])(x + 0.5 * h) * (v2->at(i - 1) + 0.5 * kp2_2) + (*f_lineal[1])(x + 0.5 * h) * (v1->at(i - 1) + 0.5 * kp2_1));
        // RK4
        kp4_1 = h * (v2->at(i - 1) + kp3_2);
        kp4_2 = h * ((*f_lineal[0])(x + h) * (v2->at(i - 1) + kp3_2) + (*f_lineal[1])(x + h) * (v1->at(i - 1) + kp3_1));

        // Redefinimos v
        v1->at(i) = v1->at(i - 1) + (kp1_1 + 2 * kp2_1 + 2 * kp3_1 + kp4_1) / 6.0;
        v2->at(i) = v2->at(i - 1) + (kp1_2 + 2 * kp2_2 + 2 * kp3_2 + kp4_2) / 6.0;
    }

    w1 = alpha;
    w2 = (beta - u1->at(N)) / v1->at(N);

    for (int i = 1; i < N + 1; i++)
    {
        w[i][1] = u1->at(i) + w2 * v1->at(i);
    }

    return w;
}

double **Shooting::no_linear_solutions()
{
    int k = 1;
    while (k <= N_max)
    {
        u1->at(0) = alpha;
        u2->at(0) = TK;
        w1 = 0.0;
        w2 = 1.0;

        w[0][0] = a;
        for (int i = 1; i < N+1; i++)
        {
            x = a + (i)*h;
            w[0][i] = x;

            k1_1 = h * u2->at(i - 1);
            k1_2 = h * (*f_no_lineal[0])(x, u1->at(i - 1), u2->at(i - 1));
            k2_1 = h * (u2->at(i - 1) + 0.5 * k1_2);
            k2_2 = h * (*f_no_lineal[0])(x + 0.5 * h, u1->at(i - 1) + 0.5 * k1_1, u2->at(i - 1) + 0.5 * k1_2);
            k3_1 = h * (u2->at(i - 1) + 0.5 * k2_2);
            k3_2 = h * (*f_no_lineal[0])(x + 0.5 * h, u1->at(i - 1) + 0.5 * k2_1, u2->at(i - 1) + 0.5 * k2_2);
            k4_1 = h * (u2->at(i - 1) + k3_2);
            k4_2 = h * (*f_no_lineal[0])(x + 0.5 * h, u1->at(i - 1) + k3_1, u2->at(i - 1) + k3_2);

            u1->at(i) = u1->at(i - 1) + (k1_1 + 2.0 * k2_1 + 2.0 * k3_1 + k4_1) / 6.0;
            u2->at(i) = u2->at(i - 1) + (k1_2 + 2.0 * k2_2 + 2.0 * k3_2 + k4_2) / 6.0;

            kp1_1 = h * w2;
            kp1_2 = h * ((*f_no_lineal[1])(x, u1->at(i - 1), u2->at(i - 1)) * w1 + (*f_no_lineal[2])(x, u1->at(i - 1), u2->at(i - 1)) * w2);
            kp2_1 = h * (w2 + 0.5 * kp1_2);
            kp2_2 = h * ((*f_no_lineal[1])(x + 0.5 * h, u1->at(i - 1), u2->at(i - 1)) * (w1 + 0.5 * kp1_1) + (*f_no_lineal[2])(x + 0.5 * h, u1->at(i - 1), u2->at(i - 1)) * (w2 + 0.5 * kp1_2));
            kp3_1 = h * (w2 + 0.5 * kp2_2);
            kp3_2 = h * ((*f_no_lineal[1])(x + 0.5 * h, u1->at(i - 1), u2->at(i - 1)) * (w1 + 0.5 * kp2_1) + (*f_no_lineal[2])(x + 0.5 * h, u1->at(i - 1), u2->at(i - 1)) * (w2 + 0.5 * kp2_2));
            kp4_1 = h * (w2 + kp3_2);
            kp4_2 = h * ((*f_no_lineal[1])(x + 0.5 * h, u1->at(i - 1), u2->at(i - 1)) * (w1 + kp3_1) + (*f_no_lineal[2])(x + 0.5 * h, u1->at(i - 1), u2->at(i - 1)) * (w2 + kp3_2));
            w1 += (kp1_1 + 2.0 * kp2_1 + 2.0 * kp3_1 + kp4_1) / 6.0;
            w2 += (kp1_2 + 2.0 * kp2_2 + 2.0 * kp3_2 + kp4_2) / 6.0;
        }

        if (abs(u1->at(N) - beta) <= TOL)
        {
            for (int i = 0; i < N+1; i++)
            {
                w[1][i] = u1->at(i);
            }
            return w;
            break;
        }

        TK = TK - (u1->at(N) - beta) / w1;
        k++;
    }
    cout << "Se excedió el número máximo de iteraciones" << endl;
    return w;
}
