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
