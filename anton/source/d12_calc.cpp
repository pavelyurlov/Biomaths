#include <cmath>

#include <eigen3/Eigen/Dense>

typedef Eigen::ArrayXd Xarr;
typedef Eigen::ArrayXXd XXarr;


Eigen::ArrayXd d12_calc(
        XXarr D12, XXarr D11, XXarr D22,
        Xarr m1, Xarr m2,
        Xarr w11, Xarr w12, Xarr w21, Xarr w22,
        double b1, double b2, double d1, double d2,
        double d11, double d12, double d21, double d22,
        int N1, int N2, double al, double A, int N, int dim,
        ? I, ? J) {

    XXarr ones = Eigen::ArrayXXd::Constant(1, N, 1);
    Xarr left = ((1. - al / 2) * (b1 + b2) + al / 2 *
            (d1 + d2 + N1 * (d11 + d21) + N2 * (d12 + d22))) * ones +
        w12 + w21;
    XXarr right = conv_dim(m1 + m2, D12, N, A, dim, I, J) -
        al / 2 * N1 * ((D12 + 2 * ones) *
                (conv_dim(w11, D12, N, A, dim, I, J) +
                conv_dim(w21, D11, N, A, dim, I, J)) +
                conv_dim(w21 * D12, D11, N, A, dim, I, J) +
                conv_dim(w11 * D11, D12, N, A, dim, I, J)) -
        al / 2 * N2 * ((D12 + 2 * ones) *
                (conv_dim(w12, D22, N, A, dim, I, J) +
                conv_dim(w22, D12, N, A, dim, I, J)) +
                conv_dim(w22 * D22, D12, N, A, dim, I, J) +
                conv_dim(w12 * D12, D22, N, A, dim, I, J));
    right = right - w12 - w21;
    return right / left;
}

