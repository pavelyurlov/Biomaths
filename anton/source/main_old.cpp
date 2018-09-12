#include <iomanip>
#include <chrono>
#include <fstream>

#include <initializer_list>
#include <cstdlib>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <vector>


#define EIGEN_FFTW_DEFAULT

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Householder>

#include "convolution.h"


const double PI = 3.141592653589793238462643383279502884197169;
const double ONE_OVER_SQRT_2PI = 0.39894228040143267793994605993438;



typedef Eigen::ArrayXd arr;
typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;



inline double square(double x) {
    // std::cout << __PRETTY_FUNCTION__ << std::endl;

    return x * x;
}

double normpdf(double x, double mu, double sigma) {
    // std::cout << __PRETTY_FUNCTION__ << std::endl;

    return (ONE_OVER_SQRT_2PI / sigma) *
        std::exp(-0.5 * square((x - mu) / sigma));
}

vec normpdf(const vec& x, double mu, double sigma) {
    // std::cout << __PRETTY_FUNCTION__ << std::endl;

    size_t n = x.size();
    vec res(n);
    for (size_t i = 0; i < n; ++i) {
        res(i) = normpdf(x(i), mu, sigma);
    }
    return res;
}

/*
 * TODO:
 * normpdf for 2 and 3 dimensions
 */


class InputSet {
 public:
    double N1, N2;
    vec D11, D12, D22;
    vec w11, w12, w21, w22;
    double d11, d12, d21, d22;
    vec m1, m2;
    double b1, b2;
    double d1, d2;
    double h, A, al;
    size_t N;
    size_t dim;
};

class OutputSet {
 public:
    double N1, N2;
    vec D11, D12, D22;
    int iter;
};


double y_calc(const vec&, const vec&, double, double, size_t, double, size_t);
vec n_calc(InputSet&);
vec d11_calc(InputSet&);
vec d12_calc(InputSet&);
vec d22_calc(InputSet&);


OutputSet solver(
        InputSet& s) {
    /*
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    */

    double eps = 1e-8, eps2 = 1e-4;
    int max_iter = 500;
    double y11_old = 1, y12_old = 1, y21_old = 1, y22_old = 1, mistake = 1;
    double n1_old = 100000, n2_old = 100000, mistake2 = 1000;

    Eigen::ArrayXd r = Eigen::ArrayXd::LinSpaced(s.N, 0, s.A);
    Eigen::ArrayXd k = PI / s.A * r;

    /*
    TODO:
    ht - Hankel transform of order 0
    iht - inverse Hankel transform of order 0
    [~, I] = ht(r, r, k);
    [~, J] = iht(r, k, r);
    */

    int iter = 0;
    while (mistake > eps && iter < max_iter && mistake2 > eps2) {
        s.D12 = d12_calc(s);
        vec temp = n_calc(s);
        s.N1 = temp(0); s.N2 = temp(1);
        s.D11 = d11_calc(s);
        s.D22 = d22_calc(s);
        double y11 = y_calc(s.w11, s.D11, s.d11, s.h, s.N, s.A, s.dim);
        double y12 = y_calc(s.w12, s.D12, s.d12, s.h, s.N, s.A, s.dim);
        double y21 = y_calc(s.w21, s.D12, s.d21, s.h, s.N, s.A, s.dim);
        double y22 = y_calc(s.w22, s.D22, s.d22, s.h, s.N, s.A, s.dim);
        mistake = std::abs(y11 - y11_old) / y11 +
                std::abs(y12 - y12_old) / y12 +
                std::abs(y21 - y21_old) / y21 +
                std::abs(y22 - y22_old) / y22;
        y11_old = y11;
        y12_old = y12;
        y21_old = y21;
        y22_old = y22;
        ++iter;
        temp = n_calc(s);
        s.N1 = temp(0); s.N2 = temp(1);
        mistake2 = std::abs(s.N1 - n1_old) + std::abs(s.N2 - n2_old);
        n1_old = s.N1;
        n2_old = s.N2;
    }

    OutputSet res;
    res.N1 = s.N1;
    res.N2 = s.N2;
    res.D11 = s.D11;
    res.D12 = s.D12;
    res.D22 = s.D22;
    res.iter = iter;

    return res;
}

vec n_calc(InputSet& s) {


    vec left(2);
    left << s.b1 - s.d1, s.b2 - s.d2;

    double y11 = y_calc(s.w11, s.D11, s.d11, s.h, s.N, s.A, s.dim);
    double y12 = y_calc(s.w12, s.D12, s.d12, s.h, s.N, s.A, s.dim);
    double y21 = y_calc(s.w21, s.D12, s.d21, s.h, s.N, s.A, s.dim);
    double y22 = y_calc(s.w22, s.D22, s.d22, s.h, s.N, s.A, s.dim);


    mat right(2, 2);
    right << y11, y12, y21, y22;

    mat right_inv(2, 2);
    right_inv = right.inverse();

    return right_inv * left;
}


double y_calc(
        const vec& w,
        const vec& C,
        double d, double h, size_t N, double A, size_t dim) {


    vec r;
    switch(dim) {
        case 1:
            return h * w.cwiseProduct(C).sum() + d;

        case 2:
            r.setLinSpaced(N, 0, A);
            return 2 * PI * h * w.cwiseProduct(C.cwiseProduct(r)).sum() + d;

        case 3:
            r.setLinSpaced(N, 0, A);
            return 4 * PI * h * w.cwiseProduct(C.cwiseProduct(r.cwiseProduct(r))).sum() + d;

        default:
            throw std::runtime_error("Unsupported dimension\n");
    }
}




vec d11_calc(InputSet& s) {


    vec ones = vec::Ones(s.N);
    vec left = ((1. - s.al / 2) * s.b1 + s.al / 2 *
            (s.d1 + s.N1 * s.d11 + s.N2 * s.d12)) * ones + s.w11;
    vec right = conv_dim(s.m1, s.D11, s.N, s.A, s.dim) -
        s.al / 2 * s.N1 * ((s.D11 + 2 * ones).cwiseProduct(
                conv_dim(s.w11, s.D11, s.N, s.A, s.dim)) +
                conv_dim((s.w11).cwiseProduct(s.D11), s.D11, s.N, s.A, s.dim)) -
        s.al / 2 * s.N2 * ((s.D11 + 2 * ones).cwiseProduct(
                conv_dim(s.w12, s.D12, s.N, s.A, s.dim)) +
                conv_dim((s.w12).cwiseProduct(s.D12), s.D12, s.N, s.A, s.dim)) +
        (1. / s.N1) * s.m1 - s.w11;
    // right = right + (1. / s.N1) * s.m1 - s.w11;
    return right.cwiseQuotient(left);
}



vec d12_calc(InputSet& s) {


    const vec ones = vec::Ones(s.N);
    vec left = ((1. - s.al / 2.) * (s.b1 + s.b2) + s.al / 2. *
        (s.d1 + s.d2 + s.N1 * (s.d11 + s.d21) + s.N2 * (s.d12 + s.d22))) * ones +
        s.w12 + s.w21;
    vec right = conv_dim(s.m1 + s.m2, s.D12, s.N, s.A, s.dim) -
        s.al / 2. * s.N1 * (
            (s.D12 + 2. * ones).cwiseProduct(
                conv_dim(s.w11, s.D12, s.N, s.A, s.dim) +
                conv_dim(s.w21, s.D11, s.N, s.A, s.dim)
            ) +
                conv_dim((s.w21).cwiseProduct(s.D12), s.D11, s.N, s.A, s.dim) +
                conv_dim((s.w11).cwiseProduct(s.D11), s.D12, s.N, s.A, s.dim)
        ) -
    // right -=
        s.al / 2. * s.N2 * (
            (s.D12 + 2. * ones).cwiseProduct(
                conv_dim(s.w12, s.D22, s.N, s.A, s.dim) +
                conv_dim(s.w22, s.D12, s.N, s.A, s.dim)
            ) +
                conv_dim((s.w22).cwiseProduct(s.D22), s.D12, s.N, s.A, s.dim) +
                conv_dim((s.w12).cwiseProduct(s.D12), s.D22, s.N, s.A, s.dim)
        ) - s.w12 - s.w21;
    return right.cwiseQuotient(left);
}


/* D12 == D21 */
vec d22_calc(InputSet& s) {

    vec ones = vec::Ones(s.N);
    vec left = ((1 - s.al / 2) * s.b2 + s.al / 2 *
            (s.d2 + s.N1 * s.d21 + s.N2 * s.d22)) * ones + s.w22;
    vec right = conv_dim(s.m2, s.D22, s.N, s.A, s.dim) -
        s.al / 2 * s.N2 * ((s.D22 + 2 * ones).cwiseProduct(
                conv_dim(s.w22, s.D22, s.N, s.A, s.dim)) +
                conv_dim((s.w22).cwiseProduct(s.D22), s.D22, s.N, s.A, s.dim)) -
        s.al / 2 * s.N1 * ((s.D22 + 2 * ones).cwiseProduct(
                conv_dim(s.w21, s.D12, s.N, s.A, s.dim)) +
                conv_dim((s.w21).cwiseProduct(s.D12), s.D12, s.N, s.A, s.dim)) +
        (1. / s.N2) * s.m2 - s.w22;
    // right = right + (1. / s.N2) * s.m2 - s.w22;
    return right.cwiseQuotient(left);
}




void plot(const vec &x, std::initializer_list<vec> y) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    FILE *out = fopen("data", "w");
    for (ssize_t i = 0; i < x.size(); ++i) {
        fprintf(out, "%lf ", x(i));
        for (const auto &v : y) {
            fprintf(out, "%lf ", v(i));
        }
        fprintf(out, "\n");
    }
    fclose(out);
}



int main() {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    InputSet s;
    s.dim = 1;
    s.N = 512;
    s.A = 2.;
    vec r;
    r.setLinSpaced(s.N, 0, s.A);
    s.h = r(1) - r(0);
    double sm1 = 0.04, sm2 = 0.1;
    s.b1 = 0.4; s.b2 = 0.4;
    s.d1 = 0.2; s.d2 = 0.2;
    s.d11 = s.d21 = s.d22 = 0.001; s.d12 = 0.0005;
    double sw11 = 0.04, sw12 = 0.04, sw21 = 0.04, sw22 = 0.04;

    s.m1 = s.b1 * normpdf(r, 0, sm1);
    s.m2 = s.b2 * normpdf(r, 0, sm2);

    s.w11 = s.d11 * normpdf(r, 0, sw11);

    s.w21 = s.d21 * normpdf(r, 0, sw21);
    s.w22 = s.d22 * normpdf(r, 0, sw22);

    s.al = 0.4;
    s.N1 = 0.;
    s.N2 = 0.;

    vec d12;
    d12.setLinSpaced(1000, 0, 0.0015);

    vec N1_ans = vec::Zero(d12.size());
    vec N2_ans = vec::Zero(d12.size());

    auto start = std::chrono::steady_clock::now();
    double total_time = 0.;

    for (ssize_t i = 0; i < d12.size(); ++i) {
        s.w12 = d12(i) * normpdf(r, 0, sw12);
        s.D11 = s.D12 = s.D22 = vec::Zero(s.N);
        s.d12 = d12(i);
        s.N1 = s.N2 = 0.;

        OutputSet res = solver(s);

        N1_ans(i) = res.N1;
        N2_ans(i) = res.N2;

        auto step = std::chrono::steady_clock::now();
        double time = std::chrono::duration_cast<std::chrono::milliseconds>(step - start).count() / 1000.;
        total_time += time;
        start = step;
        std::cout << i << " " << time << " " << total_time << std::endl;
    }

    plot(d12, {N1_ans, N2_ans});

    return 0;
}


