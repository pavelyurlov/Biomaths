#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <vector>

#include <eigen3/Eigen/Dense>

const double PI = 3.141592653589793238462643383279502884197169;

typedef Eigen::ArrayXd Xarr;
typedef Eigen::ArrayXXd XXarr;


class Answer {
 public:
    int N1;
    int N2;
    Eigen::ArrayXXd D11;
    Eigen::ArrayXXd D12;
    Eigen::ArrayXXd D22;
    int iter;
};

Answer solver(
        double N1, double N2,
        Eigen::ArrayXXd D11, Eigen::ArrayXXd D12, Eigen::ArrayXXd D22,
        Eigen::ArrayXd w11, Eigen::ArrayXd w12, Eigen::ArrayXd w21, Eigen::ArrayXd w22,
        double d11, double d12, double d21, double d22,
        Eigen::ArrayXd m1, Eigen::ArrayXd m2,
        double b1, double b2,
        double d1, double d2,
        double h, double A, double al, int N, int dim) {
    double eps = 1e-8, eps2 = 1e-4;
    int max_iter = 500;
    double y11_old = 1, y12_old = 1, y21_old = 1, y22_old = 1, mistake = 1;
    double n1_old = 100000, n2_old = 100000, mistake2 = 1000;

    Eigen::ArrayXd r = Eigen::ArrayXd::LinSpaced(A / (N - 1), 0, A);
    Eigen::ArrayXd k = PI / A * r;

    /*
    TODO:
    ht - Hankel transform of order 0
    iht - inverse Hankel transform of order 0
    [~, I] = ht(r, r, k);
    [~, J] = iht(r, k, r);
    */

    int iter = 0;
    while (mistake > eps && iter < max_iter && mistake2 > eps2) {
        D12 = d12_calc(D12, D11, D22, m1, m2, w11, w12, w21,
                b1, b2, d1, d2, d11, d12, d21, d22,
                N1, N2, al, A, N, dim, I, J);
        Xarr temp = n_calc(D11, D12, D22, w11, w12, w21, w22,
                d11, d12, d21, d22, b1, b2, d1, d2, h, A, N, dim);
        N1 = temp(0); N2 = temp(1);
        D11 = d11_calc(D11, D12, m1, w11, w12, b1, d1,
                d11, d12, N1, N2, al, A, N, dim, I, J);
        D22 = d22_calc(D22, D12, m2, w22, w21, b2, d2,
                d22, d21, N2, N1, al, A, N, dim, I, J);
        double y11 = y_calc(w11, D11, d11, h, N, A, dim);
        double y12 = y_calc(w12, D12, d12, h, N, A, dim);
        double y21 = y_calc(w21, D12, d21, h, N, A, dim);
        double y22 = y_calc(w22, D22, d22, h, N, A, dim);
        mistake = abs(y11 - y11_old) / y11 +
                abs(y12 - y12_old) / y12 +
                abs(y21 - y21_old) / y21 +
                abs(y22 - y22_old) / y22;
        y11_old = y11;
        y12_old = y12;
        y21_old = y21;
        y22_old = y22;
        ++iter;
        temp = n_calc(D11, D12, D22, w11, w12, w21, w22,
                d11, d12, d21, d22, b1, b2, d1, d2, h, A, N, dim);
        N1 = temp(0); N2 = temp(1);
        mistake2 = abs(N1 - n1_old) + abs(N2 - n2_old);
        n1_old = N1;
        n2_old = N2;
    }

    Answer res;
    res.N1 = N1;
    res.N2 = N2;
    res.D11 = D11;
    res.D12 = D12;
    res.D22 = D22;
    res.iter = iter;

    return res;
}
