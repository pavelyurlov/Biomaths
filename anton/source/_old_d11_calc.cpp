#include <cmath>

#include <eigen3/Eigen/Dense>

typedef Eigen::ArrayXd Xarr;
typedef Eigen::ArrayXXd XXarr;

Eigen::ArrayXd d11_calc(
        XXarr D11, XXarr D12,
        Xarr m1,
        Xarr w11, Xarr w12,
        double b1, double d1, double d11, double d12,
        int N1, int N2, double al, double A, int N, int dim,
        ? I, ? J) {

    XXarr ones = Eigen::ArrayXXd::Constant(1, N, 1);
    Xarr left = ((1 - al / 2) * b1 + al / 2 *
            (d1 + N1 * d11 + N2 * d12)) * ones + w11;
    XXarr right = conv_dim(m1, D11, N, A, dim, I, J) -
        al / 2 * N1 * ((D11 + 2 * ones) *
                conv_dim(w11, D11, N, A, dim, I, J) +
                conv_dim(w11 * D11, D11, N, A, dim, I, J)) -
        al / 2 * N2 * ((D11 + 2 * ones) *
                conv_dim(w12, D12, N, A, dim, I, J) +
                conv_dim(w12 * D12, D12, N, A, dim, I, J));
    right = right + (1. / N1) * m1 - w11;
    return right / left;
}

