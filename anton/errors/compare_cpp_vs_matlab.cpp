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
    fprintf(f, "set output \"c++_vs_matlab_%s_d%d.png\"\n", type, dim);
    fprintf(f, "set title \"Error C++ vs Matlab: %s D%d, in %%\"\n", type, dim);
    fprintf(f, "set view 30, 45, 1, 1\n");
    fprintf(f, "splot \'%s\' nonuniform matrix w l", filename1);
    fprintf(f, ", \'%s\' nonuniform matrix w l\n", filename2);

    fclose(f);
}






int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "need exactly 2 arguments: 'hm' or 'ccto', and dim (1, 2 or 3)\n";
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

    sprintf(outname1, "c++_vs_matlab_%sD%d_E1.data", argv[1], dim);
    sprintf(outname2, "c++_vs_matlab_%sD%d_E2.data", argv[1], dim);
    sprintf(incpp1, "c++_%sD%d_N1.data", argv[1], dim);
    sprintf(incpp2, "c++_%sD%d_N2.data", argv[1], dim);
    sprintf(inmatlab1, "matlab_%sD%d_N1.txt", argv[1], dim);
    sprintf(inmatlab2, "matlab_%sD%d_N2.txt", argv[1], dim);

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
    mat error_1(N, vector<double>(N));
    mat error_2(N, vector<double>(N));
    int count1_percent = 0, count2_percent = 0;
    int count1_abs = 0, count2_abs = 0;

    double delta;
    cerr << "Enter delta (for %% error)\n";
    cin >> delta;
    cout << "delta = " << delta << endl;

    double difference;
    cerr << "Enter difference (for absolute error)\n";
    cin >> difference;
    cout << "difference = " << difference << endl;


    vector<pair<int, int> > mistakes1;
    mistakes1.reserve(N * N);
    vector<pair<int, int> > mistakes2;
    mistakes2.reserve(N * N);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double e;
            e = error_1[i][j] = (my_1[i][j] - anton_1[i][j]) / anton_1[i][j] * 100;
            if (my_abs(e) > delta) {
                cout << "N1:   " << i << ' ' << j << " : " << e << " : " << my_1[i][j] << " vs " << anton_1[i][j] << endl;
                ++count1_percent;
                if (my_abs(my_1[i][j] - anton_1[i][j]) > difference) {
                    ++count1_abs;
                    mistakes1.push_back(make_pair(i, j));
                }
            }
            e = error_2[i][j] = (my_2[i][j] - anton_2[i][j]) / anton_2[i][j] * 100;
            if (my_abs(e) > delta) {
                cout << "N2:   " << i << ' ' << j << " : " << e << " : " << my_2[i][j] << " vs " << anton_2[i][j] << endl;
                ++count2_percent;
                if (my_abs(my_2[i][j] - anton_2[i][j]) > difference) {
                    ++count2_abs;
                    mistakes2.push_back(make_pair(i, j));
                }
            }
        }
    }
    cout << "count1 in % = " << count1_percent << endl;
    cout << "count2 in % = " << count2_percent << endl;
    cout << "total in % = " << count1_percent + count2_percent << endl;

    cerr << "count1 in % = " << count1_percent << endl;
    cerr << "count2 in % = " << count2_percent << endl;
    cerr << "total in % = " << count1_percent + count2_percent << endl;



    cout << "\n\n\n\n\n";
    for (int i = 0; i < N; ++i) {
        cout << "*";
    }
    cout << "\n\n\n\n\n";

    cout << "\n\n\n\n\n";
    for (int i = 0; i < N; ++i) {
        cout << "*";
    }
    cout << "\n\n\n\n\n";

    cout << "\n\n\n\n\n";
    for (int i = 0; i < N; ++i) {
        cout << "*";
    }
    cout << "\n\n\n\n\n";

    cout << "\n\n\n\n\n";
    for (int i = 0; i < N; ++i) {
        cout << "*";
    }
    cout << "\n\n\n\n\n";

    cout << "\n\n\n\n\n";
    for (int i = 0; i < N; ++i) {
        cout << "*";
    }
    cout << "\n\n\n\n\n";



    for (auto p : mistakes1) {
        int i = p.first, j = p.second;
        double diff = my_1[i][j] - anton_1[i][j];
        cout << "N1:   " << i << ' ' << j << " : " << diff << " : " << my_1[i][j] << " vs " << anton_1[i][j] << endl;
    }

    for (auto p : mistakes2) {
        int i = p.first, j = p.second;
        double diff = my_2[i][j] - anton_2[i][j];
        cout << "N2:   " << i << ' ' << j << " : " << diff << " : " << my_2[i][j] << " vs " << anton_2[i][j] << endl;
    }

    cout << "Among them:" << endl;
    cout << "count1 in abs = " << count1_abs << endl;
    cout << "count2 in abs = " << count2_abs << endl;
    cout << "total in abs = " << count1_abs  + count2_abs << endl;

    cerr << "Among them:" << endl;
    cerr << "count1 in abs = " << count1_abs << endl;
    cerr << "count2 in abs = " << count2_abs << endl;
    cerr << "total in abs = " << count1_abs  + count2_abs << endl;




    // plot error
    surf_plot(first, second, error_1, outname1);
    surf_plot(first, second, error_2, outname2);

    write_plot_file(argv[1], dim, outname1, outname2);

    return 0;
}
