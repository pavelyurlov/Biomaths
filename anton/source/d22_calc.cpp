#include <cmath>

#include <eigen3/Eigen/Dense>

typedef Eigen::ArrayXd Xarr;
typedef Eigen::ArrayXXd XXarr;


Eigen::ArrayXd d22_calc(
        XXarr D22, XXarr D21,
        Xarr m2,
        Xarr w22, Xarr w21,
        double b2, double d2, double d22, double d21,
        int N2, int N1, double al, double A, int N, int dim,
        ? I, ? J) {

    XXarr ones = Eigen::ArrayXXd::Constant(1, N, 1);
    Xarr left = ((1 - al / 2) * b2 + al / 2 *
            (d2 + N1 * d21 + N2 * d22)) * ones + w22;
    XXarr right = conv_dim(m2, D22, N, A, dim, I, J) -
        al / 2 * N1 * ((D22 + 2 * ones) *
                conv_dim(w22, D22, N, A, dim, I, J) +
                conv_dim(w22 * D22, D22, N, A, dim, I, J)) -
        al / 2 * N2 * ((D22 + 2 * ones) *
                conv_dim(w21, D21, N, A, dim, I, J) +
                conv_dim(w21 * D21, D21, N, A, dim, I, J));
    right = right + (1. / N1) * m2 - w22;
    return right / left;
}

