#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <cmath>

int N = 100;
using namespace std;


void myplot(const vector<int> &values, const vector<int> &counts1,
        const vector<int> &counts2, char *fname) {
    cerr << __PRETTY_FUNCTION__ << std::endl;

    for (int i = 0; i < 5; ++i) {
        fprintf(stdout, "%d %d %d\n", values[i], counts1[i], counts2[i]);
    }

    FILE *out = fopen(fname, "w");

    int size = values.size();

    for (ssize_t i = 0; i < size; ++i) {
        fprintf(out, "%d %d %d\n", values[i], counts1[i], counts2[i]);
    }

    fclose(out);
}


void write_plot_file(char *type, const int dim,
        char *filename) {
    char plotfilename[50];
    sprintf(plotfilename, "plot_stat_%s_d%d", type, dim);

    FILE *f = fopen(plotfilename, "w");

    fprintf(f, "set term png size 1920, 1080\n");
    fprintf(f, "set output \"stat_c++_vs_matlab_%s_d%d.png\"\n", type, dim);
    fprintf(f, "set title \"Statistics C++ vs Matlab (C++ - Matlab): %s D%d\"\n", type, dim);
    fprintf(f, "set view 30, 45, 1, 1\n");
    fprintf(f, "plot \'%s\' using 1:2 w l title \"N1\"", filename);
    fprintf(f, ", \'%s\' using 1:3 w l title \"N2\"\n", filename);

    fclose(f);
}





int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "need exactly 2 arguments: 'hm' or 'ccto' and dim (1, 2 or 3)\n";
        return 1;
    }

    int dim;
    sscanf(argv[2], "%d", &dim);
    cerr << "dim = " << dim << endl;

    char in1[50], in2[50];

    string ccto = "ccto", hm = "hm";
    if (argv[1] == ccto) {
        cerr << ccto << endl;
    } else if (argv[1] == hm) {
        cerr << hm << endl;
    } else {
        cerr << "neither ccto, nor hm\n";
        return 1;
    }

    sprintf(in1, "abs_c++_vs_matlab_%sD%d_N1.data", argv[1], dim);
    sprintf(in2, "abs_c++_vs_matlab_%sD%d_N2.data", argv[1], dim);

    cerr << "Files:" << endl << " " << in1 << endl << " " << in2 << endl;


    vector<double> first(N);
    vector<double> second(N);

    vector<vector<double> > N1(N, vector<double>(N));
    vector<vector<double> > N2(N, vector<double>(N));

    fstream fs1;
    fstream fs2;

    fs1.open(in1, fstream::in);
    fs2.open(in2, fstream::in);

    // read data
    double num;
    fs1 >> num;
    fs2 >> num;
    double max = -1000000000, min = 1000000000;
    for (int i = 0; i < N; ++i) {
        fs1 >> first[i];
        fs2 >> first[i];
    }
    for (int i = 0; i < N; ++i) {
        fs1 >> second[i];
        fs2 >> second[i];

        for (int j = 0; j < N; ++j) {
            fs1 >> N1[j][i];
            fs2 >> N2[j][i];

            if (max < N1[j][i]) max = N1[j][i];
            if (max < N2[j][i]) max = N2[j][i];

            if (min > N1[j][i]) min = N1[j][i];
            if (min > N2[j][i]) min = N2[j][i];
        }
    }
    cerr << "max = " << max << endl << "min = " << min << endl;


    // process

    int min_int = round(min), max_int = round(max);
    int steps = max_int - min_int + 1;

    vector<int> values;
    values.reserve(steps);
    for (int i = 0; i < steps; ++i) {
        values[i] = min_int + i;
    }

    vector<int> counts1, counts2;
    counts1.reserve(steps);
    counts2.reserve(steps);
    for (int i = 0; i < steps; ++i) {
        counts1[i] = counts2[i] = 0;
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // N1
            int pos1 = round(N1[i][j]) - min_int;
            counts1[pos1] += 1;

            // N2
            int pos2 = round(N2[i][j]) - min_int;
            counts2[pos2] += 1;
        }
    }


    // plot
    char outname[100];
    sprintf(outname, "stat_of_c++_min_matlab_%sD%d.txt", argv[1], dim);
    // myplot(values, counts1, counts2, outname);

    for (int i = 0; i < steps; ++i) {
        cout << values[i] << " " << counts1[i] << " " << counts2[i] << endl;
    }

    write_plot_file(argv[1], dim, outname);

    return 0;
}
