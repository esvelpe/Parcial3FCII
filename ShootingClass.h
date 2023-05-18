#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

class Shooting
{
public:
    Shooting(double, double, double, double, unsigned int, double (**)(double));
    Shooting(double, double, double, double, unsigned int, double (**)(double, double, double), double, int);
    ~Shooting();
    double **linear_solutions();
    double **no_linear_solutions();
    void setX(double);
    double getX() const;
    void setA(double);
    double getA() const;
    void setB(double);
    double getB() const;
    void setAlpha(double);
    double getAlpha() const;
    void setBeta(double);
    double getBeta() const;
    void setN(unsigned int);
    unsigned int getN() const;
    void setH();
    double getH() const;
    void setTOL(double);
    double getTOL() const;
    void setN_max(unsigned int);
    unsigned int getN_max() const;
    void setTK();
    double getTK() const;
    void printTable(const vector<double> &) const;

private:
    double x{0};
    double a{0}, b{0}, alpha{0}, beta{0};
    unsigned int N{0};
    double h{0};
    double TOL{0};
    double N_max{0};
    double TK{0};

    double w1{0}, w2{0};
    double **w;

    vector<double> *u1 = NULL;
    vector<double> *u2 = NULL;
    vector<double> *v1 = NULL;
    vector<double> *v2 = NULL;

    double (*f_lineal[3])(double);
    double (*f_no_lineal[3])(double, double, double);

    // Definicion de valores de RK para y
    double k1_1{0.0}, k2_1{0.0}, k3_1{0.0}, k4_1{0.0};
    double k1_2{0.0}, k2_2{0.0}, k3_2{0.0}, k4_2{0.0};
    // Definicion de valores de RK para y'
    double kp1_1{0.0}, kp2_1{0.0}, kp3_1{0.0}, kp4_1{0.0};
    double kp1_2{0.0}, kp2_2{0.0}, kp3_2{0.0}, kp4_2{0.0};
};
