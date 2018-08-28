#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>

int N = 100;
using namespace std;
typedef vector<double> vec;
typedef vector<vector<double> > mat;


double my_abs(double a) {
    return a;
}


void surf_plot(const vec &x, const vec &y, const mat &z, const std::string &fname) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    FILE *out = fopen(fname.data(), "w");

    fprintf(out, "%zd", x.size());
    for (ssize_t i = 0; i < N; ++i) {
        fprintf(out, " %lf", x[i]);
    }
    fprintf(out, "\n");

    for (ssize_t i = 0; i < N; ++i) {
        fprintf(out, "%lf", y[i]);
        for (ssize_t j = 0; j < N; ++j) {
            fprintf(out, " %lf", z[j][i]);
        }
        fprintf(out, "\n");
    }
    fclose(out);
}


    /*
    FILE *data_my_1 = fopen(argv[1], "r");
    FILE *data_my_2 = fopen(argv[2], "r");
    FILE *data_anton_1 = fopen(argv[3], "r");
    FILE *data_anton_2 = fopen(argv[4], "r");
    */

int main(int argc, char** argv) {
    if (argc != 11) {
        cout << "need exactly 10 arguments\n";
        return 1;
    }

    /*
    vec d12_vec;
    d12_vec.setLinSpaced(100, 0, 0.0015);
    vec sm2_vec;
    sm2_vec.setLinSpaced(100, 0.0001, 0.2);
    */

    double low1, low2, high1, high2, step1, step2;
    sscanf(argv[5], "%lf", &low1);
    sscanf(argv[6], "%lf", &high1);
    step1 = (high1 - low1) / (N - 1);
    sscanf(argv[7], "%lf", &low2);
    sscanf(argv[8], "%lf", &high2);
    step2 = (high2 - low2) / (N - 1);
    vec first(N);
    vec second(N);
    for (int i = 0; i < N; ++i) {
        first[i] = low1 + i * step1;
        second[i] = low2 + i * step2;
    }
    mat my_1(N, vector<double>(N));
    mat my_2(N, vector<double>(N));
    mat anton_1(N, vector<double>(N));
    mat anton_2(N, vector<double>(N));

    fstream fs_my_1;
    fstream fs_my_2;
    fstream fs_anton_1;
    fstream fs_anton_2;

    fs_my_1.open(argv[1], fstream::in);
    fs_my_2.open(argv[2], fstream::in);
    fs_anton_1.open(argv[3], fstream::in);
    fs_anton_2.open(argv[4], fstream::in);

    // reading data
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fs_my_1 >> my_1[i][j];
            fs_my_2 >> my_2[i][j];
            fs_anton_1 >> anton_1[i][j];
            fs_anton_2 >> anton_2[i][j];
        }
    }

    // compare
    mat error_1(N, vector<double>(N));
    mat error_2(N, vector<double>(N));
    int count1 = 0, count2 = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double e;
            e = error_1[i][j] = my_abs((my_1[i][j] - anton_1[i][j]) / anton_1[i][j]) * 100;
            if (e > 1 || e < -1) {
                cout << i << ' ' << j << " : " << e << endl;
                ++count1;
            }
            e = error_2[i][j] = my_abs((my_2[i][j] - anton_2[i][j]) / anton_2[i][j]) * 100;
            if (e > 1 || e < -1) {
                cout << i << ' ' << j << " : " << e << endl;
                ++count2;
            }
        }
    }
    cout << "count1 = " << count1 << endl;
    cout << "count2 = " << count2 << endl;
    cout << "total = " << count1 + count2 << endl;

    int dim;
    sscanf(argv[10], "%d", &dim);
    cout << "dim = " << dim << endl;

    char str1[50], str2[50];
    string ccto, hm;
    ccto = "ccto";
    hm = "hm";
    if (argv[9] == ccto) {
        cout << ccto << endl;
        sprintf(str1, "cctoD%d_E1.data", dim);
        sprintf(str2, "cctoD%d_E2.data", dim);
    } else if (argv[9] == hm) {
        cout << hm << endl;
        sprintf(str1, "hmD%d_E1.data", dim);
        sprintf(str2, "hmD%d_E2.data", dim);
    } else {
        cout << "neither ccto, nor hm\n";
        return 1;
    }

    // plot error
    surf_plot(first, second, error_1, str1);
    surf_plot(first, second, error_2, str2);

    return 0;
}
