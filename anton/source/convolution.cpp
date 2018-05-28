#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>


#include <specialfunctions.h>
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

	if (the_size != 0 && the_size != f.size())
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
		res[the_size / 2 + the_size % 2 - i - 1] = res[the_size / 2 + i] = the_out[i];
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


vec conv_fourier_1D(const vec &a, const vec &b) {
    Eigen::VectorXcd ca, cb;
    vec res;
    static Eigen::FFT<double> fft;
    fft.fwd(ca, a);
    fft.fwd(cb, b);
    ca = ca.cwiseProduct(cb);
    fft.inv(res, ca);
    return res;
}



vec conv_dim(const vec& a, const vec& b, size_t N, double A, size_t dim) {
    switch (dim) {
        case 1:
            return conv_fourier_1D(a, b);

        case 2:
		      return conv_radial_2D(a, b, A);

        case 3:
            vec r = vec::LinSpaced(N, 0, A);
            vec ar = a.cwiseProduct(r);
            return 4 * M_PI * conv_fourier_1D(ar, b);

        default:
            throw std::runtime_error("Unsupported dimension\n");
	}
}
