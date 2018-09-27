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
    return a < 0 ? -a : a;
}

void surf_plot(const vec &x, const vec &y, const mat &z, const string &fname) {
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

void write_plot_file(char *type, const int dim,
        char *filename1, char *filename2) {
    char plotfilename[50];
    sprintf(plotfilename, "plot_%s_d%d", type, dim);

    FILE *f = fopen(plotfilename, "w");

    fprintf(f, "set term png size 1920, 1080\n");
    fprintf(f, "set output \"abs_c++_vs_matlab_%s_d%d.png\"\n", type, dim);
    fprintf(f, "set title \"Absolute error C++ vs Matlab (C++ - Matlab): %s D%d\"\n", type, dim);
    fprintf(f, "set view 30, 45, 1, 1\n");
    fprintf(f, "splot \'%s\' nonuniform matrix w l", filename1);
    fprintf(f, ", \'%s\' nonuniform matrix w l\n", filename2);

    fclose(f);
}






int main(int argc, char** argv) {
    if (argc != 5) {
        cerr << "need exactly 4 arguments: 'hm' or 'ccto', dim (1, 2 or 3) and 2 folders: C++ and Matlab\n";
        return 1;
    }

    int dim;
    sscanf(argv[2], "%d", &dim);
    cout << "dim = " << dim << endl;

    char outname1[50], outname2[50];
    char incpp1[50], incpp2[50], inmatlab1[50], inmatlab2[50];

    string ccto = "ccto", hm = "hm";
    if (argv[1] == ccto) {
        cout << ccto << endl;
    } else if (argv[1] == hm) {
        cout << hm << endl;
    } else {
        cerr << "neither ccto, nor hm\n";
        return 1;
    }

    sprintf(outname1, "abs_c++_vs_matlab_%sD%d_N1.data", argv[1], dim);
    sprintf(outname2, "abs_c++_vs_matlab_%sD%d_N2.data", argv[1], dim);
    sprintf(incpp1, "%s/c++_%sD%d_N1.data", argv[3], argv[1], dim);
    sprintf(incpp2, "%s/c++_%sD%d_N2.data", argv[3], argv[1], dim);
    sprintf(inmatlab1, "%s/matlab_%sD%d_N1.txt", argv[4], argv[1], dim);
    sprintf(inmatlab2, "%s/matlab_%sD%d_N2.txt", argv[4], argv[1], dim);

    cerr << "Files:" << endl << " " << incpp1 << endl << " " << incpp2 << endl
        << " " << inmatlab1 << endl << " " << inmatlab2 << endl;


    vec first(N);
    vec second(N);

    mat my_1(N, vector<double>(N));
    mat my_2(N, vector<double>(N));
    mat anton_1(N, vector<double>(N));
    mat anton_2(N, vector<double>(N));

    fstream fs_my_1;
    fstream fs_my_2;
    fstream fs_anton_1;
    fstream fs_anton_2;

    fs_my_1.open(incpp1, fstream::in);
    fs_my_2.open(incpp2, fstream::in);
    fs_anton_1.open(inmatlab1, fstream::in);
    fs_anton_2.open(inmatlab2, fstream::in);

    // reading c++ data
    double num;
    fs_my_1 >> num;
    fs_my_2 >> num;
    for (int i = 0; i < N; ++i) {
        fs_my_1 >> first[i];
        fs_my_2 >> first[i];
    }
    for (int i = 0; i < N; ++i) {
        fs_my_1 >> second[i];
        fs_my_2 >> second[i];

        for (int j = 0; j < N; ++j) {
            fs_my_1 >> my_1[j][i];
            fs_my_2 >> my_2[j][i];
        }
    }


    // reading matlab data
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            char num[100];
            fs_anton_1 >> anton_1[i][j];
            sprintf(num, "%lf", anton_1[i][j]);
            sscanf(num, "%lf", &anton_1[i][j]);

            fs_anton_2 >> anton_2[i][j];
            sprintf(num, "%lf", anton_2[i][j]);
            sscanf(num, "%lf", &anton_2[i][j]);
        }
    }

    // compare
    mat abs_error_1(N, vector<double>(N));
    mat abs_error_2(N, vector<double>(N));

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            abs_error_1[i][j] = my_1[i][j] - anton_1[i][j];
            abs_error_2[i][j] = my_2[i][j] - anton_2[i][j];
        }
    }


    // plot absolute error
    surf_plot(first, second, abs_error_1, outname1);
    surf_plot(first, second, abs_error_2, outname2);

    write_plot_file(argv[1], dim, outname1, outname2);

    return 0;
}
