#include "LinearShootingClass.h"

LinearShooting::LinearShooting(double a, double b, double alpha, double beta, unsigned int N, double (*g[3])(double)){
    //Inicializa las varibales
    this->a = a; this->b = b;
    this->alpha = alpha; this->beta = beta;
    this->N = N; this->h = (b - a) / (double)N;

    //Inicializa las funciones P, Q y R
    f[0] = g[0];
    f[1] = g[1];
    f[2] = g[2];

    //Inicializa cada una de las columnas con un array de N+1 entradas
    w = new double*[2];
    w[0] = new double[N+1]; w[1] = new double[N+1]; 
    
    //Inicializa los vectores u1, u2, v1, v2
    u1 = new double[N]; u1[0] = alpha; 
    u2 = new double[N]; u2[0] = 0.0;

    v1 = new double[N]; v1[0] = 0; 
    v2 = new double[N]; v2[0] = 1.0;

}

LinearShooting::~LinearShooting(){
    delete[] w[0]; delete[] w[1]; delete[] w;
    delete[] u1; delete[] u2;
    delete[] v1; delete[] v2;
}

//Doble puntero porque retorna un arreglo de arreglos (Matriz)
double **LinearShooting::solutions()
{

    w[0][0] = a;
    w[1][0] = alpha;

        for (int i = 1; i < N + 1; i++)
        {

            x = a + i * h;
            w[0][i] = x;


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
            w[1][i] = u1[i] + w2 * v1[i];
        }

        return w;
}