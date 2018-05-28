#include <cmath>

#include <eigen3/Eigen/Dense>

typedef Eigen::ArrayXd Xarr;
typedef Eigen::ArrayXXd XXarr;

const double PI = 3.141592653589793238462643383279502884197169;


TODO:   convolution (conv()),
        types (???, ?),
        Hankel transform (ht, iht)


??? conv_dim(
        ?arr u, ?arr v, int N, double A, int dim, ? I, ? J) {
    if (dim == 1) {
        return conv(u, v, 'same');
    }
    if (dim == 2) {
        Xarr r = Xarr::LinSpaced(A / (N - 1), 0, A);
        Xarr k = PI / A * r;
        return iht(ht(u, r, k, I) * ht(v, r, k, I), k, r, J);
    }
    // if (dim == 3)
    Xarr r = Xarr::LinSpaced(A / (N - 1), 0, A);
    return 4 * PI * conv(u * r, v, 'same');
}
