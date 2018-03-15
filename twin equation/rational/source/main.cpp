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



#define PI 3.141592653589793238462643383279502884197169



double m(double x, double p) {
    return p / (PI  * (x * x + p * p));
}


double w(int n, double A, double x, double p) {
    return A / (x * x + (n + 1) * (n + 1) * p * p);
}


double S_k(int k, int n, double A, double p) {
    double result = 1.;
    for (int s = 1; s <= k; ++s) {
        double tmp = (double) ((n + 1) * (n + 1) - s * s);
        result *= tmp / (tmp + A / (p * p));
    }
    return result;
}


double C(int n, double A, double p) {
    double result = 0.;
    for (int k = 1; k <= n; ++k) {
        result += k * S_k(k, n, A, p) / ((n + 1) * (n + 1) - k * k);
    }
    result /= PI;
    result /= p;

    result += (n + 1.) * p * S_k(n, n, A, p) / (PI * A);
    result =  1. / result;
    return result;
}

double r_k(int k, int n, double A, double C_0, double p) {
    return p * C_0 * k * S_k(k, n, A, p);
}


double exact(int n, double A, double C_0, double x, double p) {
    double result = 0.;
    for (int k = 1; k <= n; ++k) {
        result += r_k(k, n, A, C_0, p) / (x * x + k * k * p * p);
    }
    result /= PI;
    result += 1.;
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
    double A, p;
    int n;
    p = 1;
    cout << "Enter A, n, such that A > 0, n \\in N\n";
    cin >> A >> n;

    cout << "N points in [-R, R]\n";
    cout << "Enter N, R, epsilon\n";
    size_t N;
    double R, epsilon;
    cin >> N >> R >> epsilon;


    auto start = std::chrono::steady_clock::now();

    VectorXd f(N);
    double C_0 = C(n, A, p);
    cout << "C = " << C_0 << endl;

    for (size_t i = 0; i < N; ++i) {
        double x = xi(R, N, i);
        f[i] = (C_0 * m(x, p) - w(n, A, x, p)) / (1 + w(n, A, x, p));
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
            K(i, j) = delta * (m(xi_i + xi_j - delta, p) +
                    m(xi_i + xi_j + delta, p)) / (1. + w(n, A, xi_i, p));
        }
    }

    size_t iters;
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
        double exact_solution = exact(n, A, C_0, xi(R, N, i), p);
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
    out << "set title \"2 graphs. Rational kernel. " << N << " points from " <<
        -R << " to " << R << ". " << iters << " iterations. A = " <<
        A << ", n = " << n << ". Execution time: " <<
        execution_time << " seconds.\"\n";
    out << "plot \"data\" using 1:3 w l title \"My\", \
        \"data\" using 1:2 w l title \"Exact\"" << endl;
    out.close();

    out.open("error");
    out << "set term png size 1920, 1080" << endl;
    out << "set output \"error.png\"" << endl;
    out << "set title \"Error. Rational kernel. " << N << " points from " <<
        -R << " to " << R << ". " << iters << " iterations. A = " << A <<
        ", n = " << n << ". Execution time: " <<
        execution_time << " seconds.\"\n";
    out << "plot \"data\" u 1:4 w l title \"|my / exact - 1| * 100%\" axes x1y2 lt rgb \"red\"" << endl;
    out.close();
}



