#ifndef SUPPORT_H
#define SUPPORT_H

#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>

using namespace std;

class AppContext {
public:
    int N_STEPS = 100;
    int N_PRINT = 1;
    bool FILE_PRINT = true;
    bool LOSE_HEAT = false;
    bool USE_LAPACK = true;

    int n;
    int N;

    int radius = -1;

    double D;
    double dx;
    double dt;
    double diff;

    double T1_heat = 0;
    double T2_heat = 300;

    int work_size;
    int work_solve_size;
    int *ipiv;
    double *work;
    double *work_solve;

    double **d_block;
    double **u_block;
    double **l_block;
    double **b_block;
    double **f_block;
    double **out_block;

    double **alfa_inv;
    double **gamma;

    long long tot_bytes;

    vector<pair<int, int>> indexes;
    int preset = 0;

    void setUpContext(int lapack_work_size) {
        // Method/function defined inside the class
        // allocare e calcolare diff, worksize, ecc
        //// calcolaer tot_bytes
        N = n * n;
        if (radius == -1) {
            radius = min(5, n / 2);
        }

        diff = D * dt / (dx * dx);
        ipiv = new int[n];
        work_size = std::max(lapack_work_size, n * n * 2);
        work_solve_size = n;
        work = new double[work_size];
        work_solve = new double[work_solve_size];

        d_block = new double *[n];
        u_block = new double *[n];
        l_block = new double *[n];
        b_block = new double *[n];
        f_block = new double *[n];
        out_block = new double *[n];

        alfa_inv = new double *[n];
        gamma = new double *[n];

        for (int i = 0; i < n; i++) {
            if (i == 0) {
                d_block[0] = new double[n * n * n]{};
                u_block[0] = new double[n * n]{};
                l_block[0] = new double[n * n]{};
                b_block[0] = new double[n * n]{};
                f_block[0] = new double[n * n]{};
                out_block[0] = new double[n * n]{};
                alfa_inv[0] = new double[n * n * (n + 1)]{};
                gamma[0] = alfa_inv[0] + n;
            } else {
                d_block[i] = d_block[0] + n * n * i;
                u_block[i] = u_block[0] + n * i;
                l_block[i] = l_block[0] + n * i;
                b_block[i] = b_block[0] + n * i;
                f_block[i] = f_block[0] + n * i;
                out_block[i] = out_block[0] + n * i;
                alfa_inv[i] = alfa_inv[0] + n * (n + 1) * i;
                gamma[i] = alfa_inv[i] + n;
            }
        }

        tot_bytes = (work_size + work_solve_size + n * n * 5 + n * n * n * 1 + n * n * (n + 1) + n * 8 + n) * 8;
        cout << "Per eseguire il programma servono circa " << tot_bytes / (1024.0 * 1024) << " MB" << endl;
    }

    void deallocate() {
        // dealloca solo lo zero
        for (int i = 0; i < 1; i++) {
            delete[] d_block[i];
            delete[] u_block[i];
            delete[] l_block[i];
            delete[] b_block[i];
            delete[] f_block[i];
            delete[] out_block[i];
            delete[] alfa_inv[i];
        }

        delete[] d_block;
        delete[] u_block;
        delete[] l_block;
        delete[] b_block;
        delete[] f_block;
        delete[] out_block;
        delete[] work;
        delete[] work_solve;

        delete[] gamma;
        delete[] alfa_inv;
        delete[] ipiv;
    }

    void reset() {
        for (int i = 0; i < n; i++) {
            // std::fill(d_block[i], d_block[i] + n * n, 0);
            // std::fill(u_block[i], u_block[i] + n, 0);
            // std::fill(l_block[i], l_block[i] + n, 0);
            std::fill(b_block[i], b_block[i] + n, T1_heat);
            std::fill(f_block[i], f_block[i] + n, 0);
            std::fill(out_block[i], out_block[i] + n, T1_heat);
            // std::fill(alfa_inv[i], alfa_inv[i] + n * (n + 1), 0);
        }

        std::fill(ipiv, ipiv + n, 0);
        std::fill(work, work + work_size, 0);
        std::fill(work_solve, work_solve + work_solve_size, 0);
    }
};

void print_matrix(string desc, int m, int n, double *a, int precision = 5) {
    cout << setprecision(precision) << fixed;
    // cout << endl;
    if (desc.length()) {
        cout << desc << endl;
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            cout << setw(precision + 5) << a[i * n + j];

        cout << endl;
    }
}

void save_snapshot(int index, int n, double **d_matrix, double **u_strip, double **l_strip,
                   double **b, double **f, double **out, double **gamma, double **alfa_inv) {
    ofstream fileg;
    int int_length = index == 0 ? 1 : int(log10(index)) + 1;
    string padding = std::string(3 - int_length, '0') + to_string(index);
    string nomeg = (string) "./matrices/matrix_" + padding + ".dat";
    fileg.open(nomeg, ios::out);
    fileg.precision(15);
    fileg << n << endl;

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n * n; k++) {
            fileg << d_matrix[i][k] << " ";
        }
    }
    fileg << endl;

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            fileg << u_strip[i][k] << " ";
        }
    }
    fileg << endl;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            fileg << l_strip[i][k] << " ";
        }
    }
    fileg << endl;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            fileg << b[i][k] << " ";
        }
    }
    fileg << endl;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            fileg << f[i][k] << " ";
        }
    }
    fileg << endl;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            fileg << out[i][k] << " ";
        }
    }
    fileg << endl;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n * n; k++) {
            fileg << gamma[i][k] << " ";
        }
    }
    fileg << endl;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n * n; k++) {
            fileg << alfa_inv[i][k] << " ";
        }
    }
}

chrono::steady_clock::time_point init_time;
void time_it(string str, bool start = false) {

    if (start) {
        init_time = chrono::steady_clock::now();
    } else {
        chrono::steady_clock::time_point now = chrono::steady_clock::now();
        cout << str << "  " << (chrono::duration_cast<chrono::microseconds>(now - init_time).count()) / 1000000.0 << endl;
        init_time = chrono::steady_clock::now();
    }
}

#endif