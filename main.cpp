#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <vector>


#include <eigen3/Eigen/Dense>


// using namespace std;
// using namespace Eigen;


const double PI = 3.141592653589793238462643383279502884197169;
const double ONE_OVER_SQRT_2PI = 0.39894228040143267793994605993438;

inline double square(double x) {
    return x * x;
}

Eigen::ArrayXd
    my_vec_normpdf(const Eigen::ArrayXd& vec, double mu, double sigma) {
        int n = vec.size();
        Eigen::ArrayXd res(n);
        for (int i = 0; i < n; ++i) {
            res(i) = (ONE_OVER_SQRT_2PI/sigma)*
                exp(-0.5 * square((vec(i) - mu) / sigma));
        }
        return res;
}

/* TODO: answer structure */


int main() {
    int N = 512;
    int A = 2;
    double h = (0. + A) / (N - 1);
    Eigen::ArrayXd r(N);
    r = Eigen::ArrayXd::LinSpaced(h, 0, A);
    double sm1 = 0.04, sm2 = 0.06;
    double b1 = 0.4, b2 = 0.4;
    double d1 = 0.2, d2 = 0.2;
    double d11 = 0.001, d12 = 0.001, d21 = 0.001, d22 = 0.001;
    double sw11 = 0.04, sw22 = 0.04;
    double sw12 = 0.04, sw21 = 0.04;

    Eigen::ArrayXd m1 = b1 * my_vec_normpdf(r, 0, sm1);
    Eigen::ArrayXd m2 = b2 * my_vec_normpdf(r, 0, sm2);

    Eigen::ArrayXd w11 = d11 * my_vec_normpdf(r, 0, sw11);
    Eigen::ArrayXd w12 = d12 * my_vec_normpdf(r, 0, sw12);
    Eigen::ArrayXd w21 = d21 * my_vec_normpdf(r, 0, sw21);
    Eigen::ArrayXd w22 = d22 * my_vec_normpdf(r, 0, sw22);

    double al = 0.8;
    int N1 = 0;
    int N2 = 0;

    double high = 0.15, low = 0.001;
    double step = (high - low) / 99;
    Eigen::ArrayXd sw1 = Eigen::ArrayXd::LinSpaced(step, low, high);
    Eigen::ArrayXd sw2 = Eigen::ArrayXd::LinSpaced(step, low, high);

    Eigen::ArrayXXd Z = Eigen::ArrayXXd::Zero(sw1.size(), r.size());

    for (int i = 0; i < sw1.size(); ++i) {
        for (int j = 0; j < sw2.size(); ++j) {
            if (i + j == 100 - 2) {
                w11 = d11 * my_vec_normpdf(r, 0, sw1(i));
                w12 = d12 * my_vec_normpdf(r, 0, sw2(j));
                w21 = d21 * my_vec_normpdf(r, 0, sw2(j));
                w22 = d22 * my_vec_normpdf(r, 0, sw1(i));

                Eigen::ArrayXXd D11 = Eigen::ArrayXXd::Zero(1, N);
                Eigen::ArrayXXd D12 = Eigen::ArrayXXd::Zero(1, N);
                Eigen::ArrayXXd D22 = Eigen::ArrayXXd::Zero(1, N);
                N1 = 0;
                N2 = 0;

                answer = solver(N1, N2, D11, D12, D22,
                        w11, w12, w21, w22,
                        d11, d12, d21, d22,
                        m1, m2, b1, b2, d1, d2,
                        h, A, al, N, 3);
            }
        }
    }

    /* write to file */

    return 0;
}
