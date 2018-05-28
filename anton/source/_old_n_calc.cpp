#include <cmath>

#include <eigen3/Eigen/Dense>

typedef Eigen::ArrayXd Xarr;
typedef Eigen::ArrayXXd XXarr;

const double PI = 3.141592653589793238462643383279502884197169;

Xarr n_calc(
        XXarr D11, XXarr D12, XXarr D22,
        Xarr w11, Xarr w12, Xarr w21, Xarr w22,
        double d11, double d12, double d21, double d22,
        double b1, double b2, double d1, double d2,
        double h, int N, double A, int dim) {
    Xarr left(2);
    left << b1 - d1, b2 - d2;
    double y11 = y_calc(w11, D11, d11, h, N, A, dim);
    double y12 = y_calc(w12, D12, d12, h, N, A, dim);
    double y21 = y_calc(w21, D12, d21, h, N, A, dim);
    double y22 = y_calc(w22, D22, d22, h, N, A, dim);

    XXarr right(2, 2);
    right << y11, y12, y21, y22;
    return right.matrix().Eigen::PartialPivLU().solve(left.matrix());
}
