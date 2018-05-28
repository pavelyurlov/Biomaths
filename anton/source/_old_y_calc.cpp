#include <cmath>

#include <eigen3/Eigen/Dense>

typedef Eigen::ArrayXd Xarr;
typedef Eigen::ArrayXXd XXarr;

const double PI = 3.141592653589793238462643383279502884197169;

double y_calc(
        Xarr w, Xarr C, double d, double h, int N, double A, int dim) {
    if (dim == 1) {
        return h * (w * C).sum() + d;
    }
    if (dim == 2) {
        Xarr r = Xarr::LinSpaced(A / (N - 1), 0, A);
        return 2 * PI * h * (w * C * r).sum() + d;
    }
    // if (dim == 3)
    Xarr r = Xarr::LinSpaced(A / (N - 1), 0, A);
    return 4 * PI * (h * w * C * r * r).sum() + d;
}
