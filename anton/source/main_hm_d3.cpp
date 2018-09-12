/* Part of this is forked from https://github.com/GadS06/NumericalMethod
 */


#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <iostream>


/* TODO:
 * download gsl (GNU Scientific Library)
 * and find out what specialfunctions.h is

#include <specialfunctions.h>
*/
#include <gsl/gsl_dht.h>


#include <fftw3.h>



#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Householder>
#include <eigen3/unsupported/Eigen/FFT>





typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;

typedef double num;
typedef unsigned int uint;


vec hankel_2D_gsl(vec f, double A)
{

	static size_t the_size = 0; // /2
	static gsl_dht *the_dht;
	static double *the_in;
	static double *the_out;

	vec res(f.size());

   size_t fsize = (size_t) f.size();
	if (the_size != 0 && the_size != fsize)
	{
		gsl_dht_free(the_dht);
		the_size = 0;
		delete the_in;
		delete the_out;
	}

	if (the_size == 0)
	{
		the_size = f.size();
		the_dht = gsl_dht_new(the_size / 2 + the_size % 2, 0, A);
		the_in = new double[the_size / 2 + the_size % 2];
		the_out = new double[the_size / 2 + the_size % 2];
	}

	for (size_t i = 0; i < the_size / 2 + the_size % 2; i++)
	{
		the_in[the_size / 2 + the_size % 2 - i - 1] = f(i);
	}

	gsl_dht_apply(the_dht, the_in, the_out);

	for (size_t i = 0; i < the_size / 2 + the_size % 2; i++)
	{
		res(the_size / 2 + the_size % 2 - i - 1) = res(the_size / 2 + i) = the_out[i];
	}
	return res;
}

// двумерное преобразование Фурье радиально-симметричных функций
// середина входного вектора считается точкой ноль
vec conv_radial_2D(const vec& a, const vec& b, double A)
{


	vec res;
	// res = \frac{ 1 }{2\pi} H[(2\pi) ^ 2 H[f] \cdot H[g]]
	// so, we need H
	res = hankel_2D_gsl((hankel_2D_gsl(a, A)).cwiseProduct(hankel_2D_gsl(b, A)) * 4 * M_PI * M_PI, A) * (1 / (2 * M_PI));
	return res;
}


vec conv_fourier_1D(const vec &f, const vec &g) {

    /* Victor's version: */

    assert(f.size() == g.size());
    ssize_t n = f.size() * 2 - 1;
    vec res(n);
    for (ssize_t k = 0; k < n; ++k) {
        for (
                ssize_t j = std::max(ssize_t(0), k - f.size() + 1);
                j <= std::min(k, f.size() - 1);
                ++j) {
            res[k] = f[j] * g[k - j];
        }
    }
    res = res.segment((n - g.size()) / 2, f.size());
    return res;



    // std::cout << __PRETTY_FUNCTION__ << std::endl;

    /*
    int const nf = f.size();
    int const ng = g.size();
    int const n  = nf + ng - 1;
    vec out(n);
    for(auto i(0); i < n; ++i) {
        int const jmn = (i >= ng - 1)? i - (ng - 1) : 0;
        int const jmx = (i <  nf - 1)? i            : nf - 1;
        for(auto j(jmn); j <= jmx; ++j) {
            out[i] = (f[j] * g[i - j]);
        }
    }
    return out.segment((n - nf) / 2, nf);
    */


/*
    Eigen::VectorXcd ca, cb;
    vec res;
    static Eigen::FFT<double> fft;
    fft.fwd(ca, a);
    fft.fwd(cb, b);
    ca = ca.cwiseProduct(cb);
    fft.inv(res, ca);
    return res;
    */
}



vec conv_dim(const vec& a, const vec& b, size_t N, double A, size_t dim) {
    // std::cout << __PRETTY_FUNCTION__ << std::endl;


    vec r, ar;
    switch (dim) {
        case 1:
            return conv_fourier_1D(a, b);

        case 2:
		      return conv_radial_2D(a, b, A);

        case 3:
            r = vec::LinSpaced(N, 0, A);
            ar = a.cwiseProduct(r);
            return 4 * M_PI * conv_fourier_1D(ar, b);

        default:
            throw std::runtime_error("Unsupported dimension\n");
	}
}





/*
 *
 *
 *
 *
 *
*/

/*
The following code was previously in main.cpp
*/

/*
 *
 *
 *
 *
 *
*/





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




void plot_dim2(const vec &x, std::initializer_list<vec> y) {
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


void surf_plot(const vec &x, const vec &y, const mat &z, const std::string &fname) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    FILE *out = fopen(fname.data(), "w");

    fprintf(out, "%zd", x.size());
    for (ssize_t i = 0; i < x.size(); ++i) {
        fprintf(out, " %lf", x(i));
    }
    fprintf(out, "\n");

    for (ssize_t i = 0; i < y.size(); ++i) {
        fprintf(out, "%lf", y(i));
        for (ssize_t j = 0; j < x.size(); ++j) {
            fprintf(out, " %lf", z(j, i));
        }
        fprintf(out, "\n");
    }
    fclose(out);
}




void ccto_d1() {
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

    plot_dim2(d12, {N1_ans, N2_ans});
}



void ccto_d2() {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    InputSet s;
    s.dim = 2;
    s.N = 512;
    s.A = 2.;
    vec r;
    r.setLinSpaced(s.N, 0, s.A);
    s.h = r(1) - r(0);
    double sm1 = 0.04;
    s.b1 = 0.4; s.b2 = 0.4;
    s.d1 = 0.2; s.d2 = 0.2;
    s.d11 = s.d21 = s.d22 = 0.001; s.d12 = 0.0005;
    double sw11 = 0.04, sw12 = 0.04, sw21 = 0.04, sw22 = 0.04;

    s.m1 = s.b1 * normpdf(r, 0, sm1);

    s.w11 = s.d11 * normpdf(r, 0, sw11);

    s.w21 = s.d21 * normpdf(r, 0, sw21);
    s.w22 = s.d22 * normpdf(r, 0, sw22);

    s.al = 0.4;
    s.N1 = 0.;
    s.N2 = 0.;

    vec d12;
    d12.setLinSpaced(100, 0, 0.0015);
    vec sm2;
    sm2.setLinSpaced(100, 0.0001, 0.2);

    mat N1_ans = mat::Zero(d12.size(), sm2.size());
    mat N2_ans = mat::Zero(d12.size(), sm2.size());

    auto start = std::chrono::steady_clock::now();
    double total_time = 0.;

    for (ssize_t i = 0; i < d12.size(); ++i) {
        for (ssize_t j = 0; j < sm2.size(); ++j) {
            s.m2 = s.b2 * normpdf(r, 0, sm2(j));
            s.w12 = d12(i) * normpdf(r, 0, sw12);

            s.D11 = s.D12 = s.D22 = vec::Zero(s.N);
            s.d12 = d12(i);
            s.N1 = s.N2 = 0.;

            OutputSet res = solver(s);

            N1_ans(i, j) = res.N1;
            N2_ans(i, j) = res.N2;

            auto step = std::chrono::steady_clock::now();
            double time = std::chrono::duration_cast<std::chrono::milliseconds>(step - start).count() / 1000.;
            total_time += time;
            start = step;
            std::cout << i << " " << j << " " << time << " " << total_time << std::endl;
        }
    }

    surf_plot(d12, sm2, N1_ans, "N1.data");
    surf_plot(d12, sm2, N2_ans, "N2.data");
}

void ccto_d3() {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    InputSet s;
    s.dim = 3;
    s.N = 512;
    s.A = 2.;
    vec r;
    r.setLinSpaced(s.N, 0, s.A);
    s.h = r(1) - r(0);
    double sm1 = 0.04;
    s.b1 = 0.4; s.b2 = 0.4;
    s.d1 = 0.2; s.d2 = 0.2;
    s.d11 = s.d21 = s.d22 = 0.001; s.d12 = 0.0005;
    double sw11 = 0.04, sw12 = 0.04, sw21 = 0.04, sw22 = 0.04;

    s.m1 = s.b1 * normpdf(r, 0, sm1);

    s.w11 = s.d11 * normpdf(r, 0, sw11);

    s.w21 = s.d21 * normpdf(r, 0, sw21);
    s.w22 = s.d22 * normpdf(r, 0, sw22);

    s.al = 0.4;
    s.N1 = 0.;
    s.N2 = 0.;

    vec d12;
    d12.setLinSpaced(100, 0, 0.0015);
    vec sm2;
    sm2.setLinSpaced(100, 0.001, 0.2);

    mat N1_ans = mat::Zero(d12.size(), sm2.size());
    mat N2_ans = mat::Zero(d12.size(), sm2.size());

    auto start = std::chrono::steady_clock::now();
    double total_time = 0.;

    for (ssize_t i = 0; i < d12.size(); ++i) {
        for (ssize_t j = 0; j < sm2.size(); ++j) {
            s.m2 = s.b2 * normpdf(r, 0, sm2(j));
            s.w12 = d12(i) * normpdf(r, 0, sw12);

            s.D11 = s.D12 = s.D22 = vec::Zero(s.N);
            s.d12 = d12(i);
            s.N1 = s.N2 = 0.;

            OutputSet res = solver(s);

            N1_ans(i, j) = res.N1;
            N2_ans(i, j) = res.N2;

            auto step = std::chrono::steady_clock::now();
            double time = std::chrono::duration_cast<std::chrono::milliseconds>(step - start).count() / 1000.;
            total_time += time;
            start = step;
            std::cout << i << " " << j << " " << time << " " << total_time << std::endl;
        }
    }

    surf_plot(d12, sm2, N1_ans, "N1.data");
    surf_plot(d12, sm2, N2_ans, "N2.data");
}


void hm_d1() {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    InputSet s;
    s.N = 512;
    s.dim = 1;
    s.A = 2.;
    vec r;
    r.setLinSpaced(s.N, 0, s.A);
    s.h = r(1) - r(0);
    double sm1 = 0.04, sm2 = 0.06;
    s.b1 = 0.4; s.b2 = 0.4;
    s.d1 = 0.2; s.d2 = 0.2;
    s.d11 = s.d12 = s.d21 = s.d22 = 0.001;
    double sw11 = 0.04, sw12 = 0.04, sw21 = 0.04, sw22 = 0.04;

    s.m1 = s.b1 * normpdf(r, 0, sm1);
    s.m2 = s.b2 * normpdf(r, 0, sm2);

    s.w11 = s.d11 * normpdf(r, 0, sw11);
    s.w12 = s.d12 * normpdf(r, 0, sw12);
    s.w21 = s.d21 * normpdf(r, 0, sw21);
    s.w22 = s.d22 * normpdf(r, 0, sw22);

    s.al = 0.4;
    s.N1 = 0.;
    s.N2 = 0.;

    vec sw1;
    sw1.setLinSpaced(100, 0.001, 0.15);
    vec sw2;
    sw2.setLinSpaced(100, 0.001, 0.15);

    mat N1_ans = mat::Zero(sw1.size(), sw2.size());
    mat N2_ans = mat::Zero(sw1.size(), sw2.size());

    auto start = std::chrono::steady_clock::now();
    double total_time = 0.;

    for (ssize_t i = 0; i < sw1.size(); ++i) {
        for (ssize_t j = 0; j < sw2.size(); ++j) {
            s.w11 = s.d11 * normpdf(r, 0, sw1(i));
            s.w12 = s.d12 * normpdf(r, 0, sw2(j));
            s.w21 = s.d21 * normpdf(r, 0, sw2(j));
            s.w22 = s.d22 * normpdf(r, 0, sw1(i));

            s.D11 = s.D12 = s.D22 = vec::Zero(s.N);
            s.N1 = s.N2 = 0.;

            OutputSet res = solver(s);

            N1_ans(i, j) = res.N1;
            N2_ans(i, j) = res.N2;

            auto step = std::chrono::steady_clock::now();
            double time = std::chrono::duration_cast<std::chrono::milliseconds>(step - start).count() / 1000.;
            total_time += time;
            start = step;
            std::cout << i << " " << j << " " << time << " " << total_time << std::endl;
        }
    }

    surf_plot(sw1, sw2, N1_ans, "N1.data");
    surf_plot(sw1, sw2, N2_ans, "N2.data");
}



void hm_d2() {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    InputSet s;
    s.N = 512;
    s.dim = 2;
    s.A = 2.;
    vec r;
    r.setLinSpaced(s.N, 0, s.A);
    s.h = r(1) - r(0);
    double sm1 = 0.04, sm2 = 0.06;
    s.b1 = 0.4; s.b2 = 0.4;
    s.d1 = 0.2; s.d2 = 0.2;
    s.d11 = s.d12 = s.d21 = s.d22 = 0.001;
    double sw11 = 0.04, sw12 = 0.04, sw21 = 0.04, sw22 = 0.04;

    s.m1 = s.b1 * normpdf(r, 0, sm1);
    s.m2 = s.b2 * normpdf(r, 0, sm2);

    s.w11 = s.d11 * normpdf(r, 0, sw11);
    s.w12 = s.d12 * normpdf(r, 0, sw12);
    s.w21 = s.d21 * normpdf(r, 0, sw21);
    s.w22 = s.d22 * normpdf(r, 0, sw22);

    s.al = 0.4;
    s.N1 = 0.;
    s.N2 = 0.;

    vec sw1;
    sw1.setLinSpaced(100, 0.001, 0.15);
    vec sw2;
    sw2.setLinSpaced(100, 0.001, 0.15);

    mat N1_ans = mat::Zero(sw1.size(), sw2.size());
    mat N2_ans = mat::Zero(sw1.size(), sw2.size());

    auto start = std::chrono::steady_clock::now();
    double total_time = 0.;

    for (ssize_t i = 0; i < sw1.size(); ++i) {
        for (ssize_t j = 0; j < sw2.size(); ++j) {
            s.w11 = s.d11 * normpdf(r, 0, sw1(i));
            s.w12 = s.d12 * normpdf(r, 0, sw2(j));
            s.w21 = s.d21 * normpdf(r, 0, sw2(j));
            s.w22 = s.d22 * normpdf(r, 0, sw1(i));

            s.D11 = s.D12 = s.D22 = vec::Zero(s.N);
            s.N1 = s.N2 = 0.;

            OutputSet res = solver(s);

            N1_ans(i, j) = res.N1;
            N2_ans(i, j) = res.N2;

            auto step = std::chrono::steady_clock::now();
            double time = std::chrono::duration_cast<std::chrono::milliseconds>(step - start).count() / 1000.;
            total_time += time;
            start = step;
            std::cout << i << " " << j << " " << time << " " << total_time << std::endl;
        }
    }

    surf_plot(sw1, sw2, N1_ans, "N1.data");
    surf_plot(sw1, sw2, N2_ans, "N2.data");
}



void hm_d3() {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    InputSet s;
    s.N = 512;
    s.dim = 2;
    s.A = 2.;
    vec r;
    r.setLinSpaced(s.N, 0, s.A);
    s.h = r(1) - r(0);
    double sm1 = 0.04, sm2 = 0.06;
    s.b1 = 0.4; s.b2 = 0.4;
    s.d1 = 0.2; s.d2 = 0.2;
    s.d11 = s.d12 = s.d21 = s.d22 = 0.001;
    double sw11 = 0.04, sw12 = 0.04, sw21 = 0.04, sw22 = 0.04;

    s.m1 = s.b1 * normpdf(r, 0, sm1);
    s.m2 = s.b2 * normpdf(r, 0, sm2);

    s.w11 = s.d11 * normpdf(r, 0, sw11);
    s.w12 = s.d12 * normpdf(r, 0, sw12);
    s.w21 = s.d21 * normpdf(r, 0, sw21);
    s.w22 = s.d22 * normpdf(r, 0, sw22);

    s.al = 0.4;
    s.N1 = 0.;
    s.N2 = 0.;

    vec sw1;
    sw1.setLinSpaced(100, 0.001, 0.15);
    vec sw2;
    sw2.setLinSpaced(100, 0.001, 0.15);

    mat N1_ans = mat::Zero(sw1.size(), sw2.size());
    mat N2_ans = mat::Zero(sw1.size(), sw2.size());

    auto start = std::chrono::steady_clock::now();
    double total_time = 0.;

    for (ssize_t i = 0; i < sw1.size(); ++i) {
        for (ssize_t j = 0; j < sw2.size(); ++j) {
            s.w11 = s.d11 * normpdf(r, 0, sw1(i));
            s.w12 = s.d12 * normpdf(r, 0, sw2(j));
            s.w21 = s.d21 * normpdf(r, 0, sw2(j));
            s.w22 = s.d22 * normpdf(r, 0, sw1(i));

            s.D11 = s.D12 = s.D22 = vec::Zero(s.N);
            s.N1 = s.N2 = 0.;

            OutputSet res = solver(s);

            N1_ans(i, j) = res.N1;
            N2_ans(i, j) = res.N2;

            auto step = std::chrono::steady_clock::now();
            double time = std::chrono::duration_cast<std::chrono::milliseconds>(step - start).count() / 1000.;
            total_time += time;
            start = step;
            std::cout << i << " " << j << " " << time << " " << total_time << std::endl;
        }
    }

    surf_plot(sw1, sw2, N1_ans, "N1.data");
    surf_plot(sw1, sw2, N2_ans, "N2.data");
}

int main() {
    hm_d3();
    return 0;
}

