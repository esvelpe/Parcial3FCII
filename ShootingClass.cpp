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
    this->f_lineal[0] = g[0];
    this->f_lineal[1] = g[1];
    this->f_lineal[2] = g[2];
    this->w1 = 0.0;
    this->w2 = 0.0;

    // Inicializa la matriz de la solución
    w = new double *[2];
    // Inicializa cada una de las columnas con un array de N+1 entradas
    w[0] = new double[N + 1];
    w[1] = new double[N + 1];

    u1 = new vector<double>(N + 1);
    u1->at(0) = alpha;
    u2 = new vector<double>(N + 1);
    u2->at(0) = 0.0;
    v1 = new vector<double>(N + 1);
    v1->at(0) = 0.0;
    v2 = new vector<double>(N + 1);
    v2->at(0) = 1.0;
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
    this->f_no_lineal[0] = g[0];
    this->f_no_lineal[1] = g[1];
    this->f_no_lineal[2] = g[2];

    w = new double *[2];
    // Inicializa cada una de las columnas con un array de N+1 entradas
    w[0] = new double[N + 1];
    w[1] = new double[N + 1];

    u1 = new vector<double>(N + 1);
    u1->at(0) = alpha;
    u2 = new vector<double>(N + 1);
    u2->at(0) = 0.0;
    this->w1 = 0.0;
    this->w2 = 1.0;
}

Shooting::~Shooting(){
	delete [] w;
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
    w[0][0] = a;
    for (int i = 0; i < N + 1; i++)
    {
        w[1][i] = u1->at(i) + w2 * v1->at(i);
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
        for (int i = 1; i < N + 1; i++)
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
            for (int i = 0; i < N + 1; i++)
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

void Shooting::setX(double x)
{
    this->x = x;
};
double Shooting::getX() const
{
    return this->x;
};
void Shooting::setA(double a)
{
    this->a = a;
};
double Shooting::getA() const
{
    return this->a;
};
void Shooting::setB(double b)
{
    this->b = b;
};
double Shooting::getB() const
{
    return this->b;
};
void Shooting::setAlpha(double alpha)
{
    this->alpha = alpha;
};
double Shooting::getAlpha() const
{
    return this->alpha;
};
void Shooting::setBeta(double beta)
{
    this->beta = beta;
};
double Shooting::getBeta() const
{
    return this->beta;
};
void Shooting::setN(unsigned int N)
{
    this->N = N;
};
unsigned int Shooting::getN() const
{
    return this->N;
};
void Shooting::setH()
{
    this->h = (getB() - getA()) / (double)getN();
};
double Shooting::getH() const
{
    return this->h;
};
void Shooting::setTOL(double TOL)
{
    this->TOL = TOL;
};
double Shooting::getTOL() const
{
    return this->TOL;
};
void Shooting::setN_max(unsigned int N_max)
{
    this->N_max = N_max;
};
unsigned int Shooting::getN_max() const
{
    return this->N_max;
};
void Shooting::setTK()
{
    this->TK = (getBeta() - getAlpha()) / (getB() - getA());
};
double Shooting::getTK() const
{
    return this->TK;
};
void Shooting::printTable(const vector<double> &vector) const
{
    cout << "x_i" << setw(15) << "w(x_i)" << setw(15) << "y(x_i)" << setw(25) << "|w(x_i)-y(x_i)|" << endl;
    for (int i = 0; i < getN(); i++)
    {
        cout << this->w[0][i] << setw(15) << this->w[1][i] << setw(15) << vector.at(i) << setw(25) << abs(this->w[1][i] - vector.at(i)) << endl;
    }
};
