#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

class LinearShooting{
    public:
    
    LinearShooting(double , double , double , double , unsigned int , double (**)(double));
    ~LinearShooting();
    double ** solutions();

    private:

    double x{0};
    double a{0}, b{0}, alpha{0}, beta{0};
    unsigned int N{0}; double h{0};
    
    double w1{0}, w2{0};
    double ** w; 

    double *u1 = NULL; double *u2 = NULL;
    double *v1 = NULL; double *v2 = NULL;

    double (*f[3])(double);

    // Definicion de valores de RK para y
    double k1_1{0.0}, k2_1{0.0}, k3_1{0.0}, k4_1{0.0};
    double k1_2{0.0}, k2_2{0.0}, k3_2{0.0}, k4_2{0.0};
    // Definicion de valores de RK para y'
    double kp1_1{0.0}, kp2_1{0.0}, kp3_1{0.0}, kp4_1{0.0};
    double kp1_2{0.0}, kp2_2{0.0}, kp3_2{0.0}, kp4_2{0.0};
};