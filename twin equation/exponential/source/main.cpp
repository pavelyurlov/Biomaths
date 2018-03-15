#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <vector>


#include <eigen3/Eigen/Dense>


using namespace std;
using namespace Eigen;



// const double PI = 3.141592653589793238462643383279502884197169;


double m(double sigma, double x) {
    return exp(-sigma * abs(x));
}


double w(double A, double rho, double x) {
    return A * exp(-rho * abs(x));
}

double c_s(double A, double rho, double sigma, size_t s) {
    double product = pow(A, s + 1.);
    if (s % 2 == 0) {
        product = -product;
    }
    for (size_t j = 0; j <= s; ++j) {
        product *= (1. - sigma * sigma / (rho * rho * (j + 1.) * (j + 1.)));
    }
    return product;
}

double exact(double A, double rho, double sigma, size_t limit, double x) {
    double result = 1.;
    for (size_t s = 0; s < limit; ++s) {
        double c_tmp = c_s(A, rho, sigma, s);
        double power = -(s + 1.) * rho * abs(x);
        double e_tmp = exp(power);
        result += c_tmp * e_tmp;
    }
    return result;
}

double C(double A, double rho, double sigma, size_t limit) {
    double result = 0.;
    for (size_t s = 0; s < limit; ++s) {
        result -= c_s(A, rho, sigma, s) * (s + 1.) * rho * sigma /
            (rho * rho * (s + 1.) * (s + 1.) - sigma * sigma);
    }
    cout << "C(A, rho, sigma, limit) = " << result << endl;
    return result;
}


double xi(double R, size_t N, size_t i) {
    return -R + R / N + 2. * R * i / N;
}


bool compare_functions(size_t N, const VectorXd& f, const VectorXd& g, double epsilon) {
    for (size_t i = 0; i < N; ++i) {
        if (abs(f[i] - g[i]) >= epsilon) {
            return false;
        }
    }
    return true;
}

VectorXd neumann(const MatrixXd& K, const VectorXd& f, size_t N, double epsilon, size_t *iters) {
    VectorXd P_prev(N), P_cur(N);
    P_prev = f;
    size_t i = 0;
    while (1) {
        P_cur = K * P_prev + f;
        ++i;
        if (i % 10 == 0) {
            if (compare_functions(N, P_cur, P_prev, epsilon)) {
                break;
            }
        }
        P_prev = P_cur;
    }
    *iters = i;
    return P_cur;
}




int main() {
    double A, rho, sigma;
    size_t limit = 1000;
    cout << "Enter A, rho, such that |A| < 1, 2 > rho > 0," <<
        " (2 - rho) and 2 / rho are not natural numbers\n";
    cin >> A >> rho;
    sigma = 2.;

    cout << "N points in [-R, R]\n";
    cout << "Enter N, R, epsilon\n";
    size_t N;
    double R, epsilon;
    cin >> N >> R >> epsilon;


    auto start = std::chrono::steady_clock::now();

    VectorXd f(N);
    double C_0 = C(A, rho, sigma, limit);
    for (size_t i = 0; i < N; ++i) {
        double x = xi(R, N, i);
        f[i] = (C_0 * m(sigma, x) - w(A, rho, x)) / (1 + w(A, rho, x));
        // f[i] = (C_0 * m(sigma, x) + 2. / sigma) / (1. + w(A, rho, x)) - 1.;
        // f[i] = (C_0 * m(sigma, x)) / (1. + w(A, rho, x));
    }


    MatrixXd K(N, N);
    double delta = (double) R / (double) N;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            double xi_i = xi(R, N, i);
            double xi_j = xi(R, N, j);
            /*
            if (i == j) {
                K(i, j) = delta * m(sigma, xi_i - xi_j) / (1. + w(A, rho, xi_i));
            } else {
                K(i, j) = 0;
            }
            */

            /*
            trapezoidal rule
            */
            K(i, j) = delta * (m(sigma, xi_i + xi_j - delta) +
                    m(sigma, xi_i + xi_j + delta)) / (1. + w(A, rho, xi_i));
        }
    }

    size_t iters = 0;
    VectorXd my_solution = neumann(K, f, N, epsilon, &iters);
    cout << "Number of iterations: " << iters << endl;

    for (size_t i = 0; i < N; ++i) {
        my_solution[i] += 1.;
    }


    auto finish = std::chrono::steady_clock::now();

    double execution_time = chrono::duration_cast<chrono::milliseconds>(finish - start).count() / 1000.0;
    cout << "Execution time: " << execution_time << " seconds." << endl;

    cout << "my[N / 2] =  " << my_solution[N / 2] << endl;
    cout << "my[0] = " << my_solution[0] << endl;

    ofstream out("data");
    out << "#x\t\texact solution\tmy solution\terror in %" << endl;
    for (size_t i = 0; i < N; ++i) {
        double exact_solution = exact(A, rho, sigma, limit, xi(R, N, i));
        if (i == 0) {
            cout << "exact[0] = " << exact_solution << endl;
        }
        if (i == N / 2) {
            cout << "exact[N / 2] =  " << exact_solution << endl;
        }
        double error = 100. * abs(my_solution[i] / exact_solution - 1.);
        out << xi(R, N, i) << setprecision(10) << "\t"
            << exact_solution << setprecision(10) << "\t"
            << my_solution[i] << setprecision(10) << "\t"
            << error << setprecision(10) << endl;
    }
    out.close();


    out.open("graph");
    out << "set term png size 1920, 1080" << endl;
    out << "set output \"graph.png\"" << endl;
    out << "set title \"2 graphs. Exponential kernel. " << N << " points from " <<
        -R << " to " << R << ". " << iters << " iterations. A = " <<
        A << ", rho = " << rho << " sigma = " <<sigma << ". Execution time: " <<
        execution_time << " seconds.\"\n";
    out << "plot \"data\" using 1:3 w l title \"My\", \
        \"data\" using 1:2 w l title \"Exact\"" << endl;
    out.close();

    out.open("error");
    out << "set term png size 1920, 1080" << endl;
    out << "set output \"error.png\"" << endl;
    out << "set title \"Error. Exponential kernel. " << N << " points from " <<
        -R << " to " << R << ". " << iters << " iterations. A = " << A <<
        ", rho = " << rho << " sigma = " << sigma << ". Execution time: " <<
        execution_time << " seconds.\"\n";
    out << "plot \"data\" u 1:4 w l title \"abs(my / exact - 1) * 100%\" axes x1y2 lt rgb \"red\"" << endl;
    out.close();
}



