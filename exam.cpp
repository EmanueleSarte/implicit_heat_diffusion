#include "render.h"
#include "support.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;
extern "C" {
extern void dsytri_(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *info);
extern void dsytrf_(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
}

void print_matrix(string desc, int m, int n, double *a, int precision);
void save_snapshot(int index, int n, double **d_matrix, double **u_strip, double **l_strip,
                   double **b, double **f, double **out, double **gamma, double **alfa_inv);
void time_it(string str, bool start);

// ##########################

vector<pair<int, int>> get_circle_indexes(int n, int r, int temp) {
    int c = n / 2;
    int r2 = (r - 1) * (r - 1);

    vector<pair<int, int>> indexes;
    for (int k = c - r + 1; k < c + r; k++) {
        for (int j = c - r + 1; j < c + r; j++) {
            if ((k - c) * (k - c) + (j - c) * (j - c) <= r2) {
                indexes.push_back(make_pair(k * n + j, temp));
            }
        }
    }

    return indexes;
}

int get_lapack_work_size(int side) {
    int info;
    double workopt;
    int lwork = -1;
    dsytrf_((char *)"L", &side, NULL, &side, NULL, &workopt, &lwork, &info);
    int work_size = (int)workopt;
    if (info != 0) {
        cout << "errore get_lapack_work_size" << endl;
    }
    return work_size;
}

void get_inverse(int side, double *mat, double *work, int work_size, int *ipiv) {

    int info;
    dsytrf_((char *)"L", &side, mat, &side, ipiv, work, &work_size, &info);
    if (info != 0) {
        cout << "errore" << endl;
    }
    dsytri_((char *)"L", &side, mat, &side, ipiv, work, &info);
    if (info != 0) {
        cout << "errore" << endl;
    }
}

void get_inverse_manual_lite(int n, double *mat, double *work1, double *work2) {
    double *dec = work1;

    // print_matrix("mat:", n, n, mat);

    // mat è valorizzata soltanto nel triangolo sopra
    // bisogna tenerne conto qui
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            int cell = i * n + j;
            int cellT = j * n + i;
            for (int k = 0; k < j; k++) {
                sum += dec[i * n + k] * dec[j * n + k];
            }

            if (i == j) {
                dec[cell] = sqrt(mat[cell] - sum);
            } else {
                if (j < i) {
                    // Se sono nella parte non valorizzata, prendo il triangolo sopra
                    dec[cell] = (1.0 / dec[j * n + j] * (mat[cellT] - sum));

                } else {
                    dec[cell] = (1.0 / dec[j * n + j] * (mat[cell] - sum));
                }
            }
        }
    }

    // print_matrix("dec:", n, n, dec);

    double *temp = work2;
    fill(temp, temp + n * n, 0);
    for (int k = 0; k < n; k++) {
        for (int j = 0; j < n; j++) {
            int cell = k * n + j;
            for (int i = 0; i < k; i++) {
                temp[cell] += -1 * dec[k * n + i] * temp[i * n + j];
            }
            if (k == j) {
                temp[cell] += 1;
            }
            temp[cell] /= dec[k * n + k];
        }
    }

    // print_matrix("temp", n, n, temp);
    // Si può fare senza fill, ora probabilmente fa molte operazioni in più
    fill(mat, mat + n * n, 0);
    // print_matrix("mat", n, n, mat);
    for (int k = n - 1; k > -1; k--) {
        for (int j = 0; j < n; j++) {
            int cell = k * n + j;
            for (int i = k + 1; i < n; i++) {
                mat[cell] += -dec[i * n + k] * mat[i * n + j];
            }
            mat[cell] += temp[cell];
            mat[cell] = mat[cell] / dec[k * n + k];
        }
    }

    // print_matrix("out", n, n, mat);
}

void matrix_triup_mult_vector(int side, double *mat, double *b, double *out) {
    int n = side;
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += mat[j * n + i] * b[j];
        }
        for (int j = i; j < n; j++) {
            sum += mat[i * n + j] * b[j];
        }
        out[i] = sum;
    }
}

void matrix_triup_mult_vector_with_sum(int side, double *mat, double *b1, double *b2, double *out) {
    int n = side;
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += mat[j * n + i] * (b1[j] + b2[j]);
        }
        for (int j = i; j < n; j++) {
            sum += mat[i * n + j] * (b1[j] + b2[j]);
        }
        out[i] = sum;
    }
}

// u_block e l_block hanno n elementi ciascuno lungo n (d_block invece ha n elementi lunghi n*n)
void get_block_matrix_lite(int lil_side, double diff, double **d_block, double **u_block, double **l_block) {
    int n = lil_side;
    int N = n * n;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Crea le matrici diagonali
            int d_index = j * n + j;
            d_block[i][d_index] = 1 + 4 * diff;
            if (j != n - 1) {
                d_block[i][d_index + 1] = -diff;
            }
            if (j != 0) {
                d_block[i][d_index - 1] = -diff;
            }

            // Impostiamo le matrici laterali (u e b)
            u_block[i][j] = -diff;
            l_block[i][j] = -diff;
        }
    }
}

void set_internal_condition(int lil_side, double diff, double **d_block, double **u_block, double **l_block, double **b_block, double **f_block,
                            vector<pair<int, int>> conds) {

    int n = lil_side;
    for (pair<int, int> &cond : conds) {
        int pos = cond.first;
        int ind = pos % n;
        int big_ind = pos / n;
        int d_ind = ind * n + ind;

        d_block[big_ind][d_ind] = 1;
        l_block[big_ind][ind] = 0;
        u_block[big_ind][ind] = 0;

        if (ind != 0) {
            d_block[big_ind][d_ind - 1] = 0;
            d_block[big_ind][d_ind - n] = 0;
            // f_block[big_ind][ind - 1] += diff * cond.second;
        }
        if (ind != n - 1) {
            d_block[big_ind][d_ind + 1] = 0;
            d_block[big_ind][d_ind + n] = 0;
            // f_block[big_ind][ind + 1] += diff * cond.second;
        }
        if (big_ind != 0) {
            u_block[big_ind - 1][ind] = 0;
            // f_block[big_ind - 1][ind] += diff * cond.second;
        }
        if (big_ind != n - 1) {
            l_block[big_ind + 1][ind] = 0;
            // f_block[big_ind + 1][ind] += diff * cond.second;
        }

        b_block[big_ind][ind] = cond.second;
        // Togliamo eventuali valori di f nella posizione di questa condizione
        // perchè non vogliamo che evolva
        f_block[big_ind][ind] = 0;
    }

    for (pair<int, int> &cond : conds) {
        int pos = cond.first;
        int ind = pos % n;
        int big_ind = pos / n;
        int d_ind = ind * n + ind;
        if (ind != 0 && d_block[big_ind][d_ind - n - 1] != 1) {
            f_block[big_ind][ind - 1] += diff * cond.second;
        }
        if (ind != n - 1 && d_block[big_ind][d_ind + n + 1] != 1) {
            f_block[big_ind][ind + 1] += diff * cond.second;
        }
        if (big_ind != 0 && d_block[big_ind - 1][d_ind] != 1) {
            f_block[big_ind - 1][ind] += diff * cond.second;
        }
        if (big_ind != n - 1 && d_block[big_ind + 1][d_ind] != 1) {
            f_block[big_ind + 1][ind] += diff * cond.second;
        }
    }
}

void initialize_alfa_inverse(int lil_side, double **d_block, double **u_block,
                             double **l_block, int work_size, double *work, int *ipiv,
                             double **gamma, double **alfa_inv, bool use_lapack = true) {
    int n = lil_side;
    int N = n * n;

    // Metto in alfa_inv[0] la matrice alfa[0] che corrisponde a d_block[0]
    std::copy(d_block[0], d_block[0] + N, alfa_inv[0]);
    // print_matrix("alfa:", n, n, alfa_inv[0]);
    // L'inversa viene fatta soltanto nel triangolo sopra, in quanto la matrice di partenza è simmetrica
    if (use_lapack) {
        get_inverse(n, alfa_inv[0], work, work_size, ipiv);
    } else {
        get_inverse_manual_lite(n, alfa_inv[0], work, work + n * n);
    }
    // print_matrix("alfa:", n, n, alfa_inv[0]);

    for (int i = 1; i < n; i++) {
        // calcoliamo prima gamma[i-1] usando alfa_inv[i-1] e u_block[i-1] facendo attenzione che:
        // se c'è un elemento sulla diagonale che è 0, la matrice potrebbe non essere simmetrica
        // ma in quel caso tutti gli elementi sotto lo 0 (nella stessa colonna) sono zero
        // possiamo prendere quindi gli elementi a destra dello zero (sulla stessa riga) e metterli
        // al posto degli zero della colonna
        //
        // Dopo procediamo a calcolarci alfa[i] usando l'elemento di gamma[i-1] appena costruito
        // sempre stando attenti al fatto che gamma[i-1] è leggermente modificata
        for (int k = 0; k < n; k++) {
            // k è il numero della colonna in cui stiamo lavorando
            // mult è il numero per cui la colonna di alfa_inv[i-1] viene moltiplicata e messa in gamma[i-1]
            double mult = u_block[i - 1][k];
            if (mult != 0) {
                // Qua semplicemente calcoliamo gamma[i-1] tenendo conto di come sono salvate gamma e alfa_inv
                // Appena calcolata il valore di gamma[i-1] lo usiamo per calcolare alfa[i]
                // Stiamo ciclando sulle righe partendo dalla diagonale e andando verso l'ultima
                for (int j = k; j < n; j++) {
                    gamma[i - 1][j * n + k] = alfa_inv[i - 1][k * n + j] * mult;
                    alfa_inv[i][k * n + j] = d_block[i][k * n + j] - l_block[i][j] * gamma[i - 1][j * n + k];
                }
            } else {
                // se mult è zero, allora dobbiamo fare il trick descritto sopra
                // La diagonale di gamma[i-1] sarà 0, e quindi alfa_inv[i] sarà semplicemente d[i]
                // questo caso lo facciamo a parte rispetto al resto
                gamma[i - 1][k * n + k] = 0;
                alfa_inv[i][k * n + k] = d_block[i][k * n + k];

                // Qua facciamo come sopra soltanto che invece di moltiplicare il valore della cella per mult
                // lo moltiplichiamo per il moltiplicatore della colonna successiva in quanto stiamo prendendo
                // i valori a destra della diagonale che sarebbero moltiplicati per gli altri valori di u_block[i-1]
                // in questo caso alfa_inv è banale, ed il fatto che gamma non sia simmetrica non influisce sul calcolo di
                // alfa in quanto equivalente a: alfa[i] = d[i] - l[i] * alfa[i-1]^-1 * u[i-1] ed essendo l ed u
                // diagonali con un buco nella stessa riga, alfa ha un buco lunga quella riga e quella colonna
                for (int j = k + 1; j < n; j++) {
                    gamma[i - 1][j * n + k] = alfa_inv[i - 1][k * n + j] * u_block[i - 1][j];
                    alfa_inv[i][k * n + j] = d_block[i][k * n + j];
                }
            }

            // print_matrix("alfa_inv:", n, n, alfa_inv[i - 1]);
        }

        // print_matrix("gamma:", n, n, gamma[0]);
        // print_matrix("alfa:", n, n, alfa_inv[i]);
        // Abbiamo l'inversa di alfa nella tridiagonale inferiore
        // print_matrix("alfa_inv:", n, n, alfa_inv[i]);
        if (use_lapack) {
            get_inverse(n, alfa_inv[i], work, work_size, ipiv);
        } else {
            get_inverse_manual_lite(n, alfa_inv[i], work, work + n * n);
        }
        // print_matrix("alfa_inv:", n, n, alfa_inv[i]);
        // print_matrix("totale i-1 :", n + 1, n, alfa_inv[i - 1]);
    }
}

void solve_block_tridiag(int lil_side, double **u_block, double **l_block, double **b_block, double **f_block,
                         double **out_block, double *work, double **gamma, double **alfa_inv) {
    int n = lil_side;
    int N = n * n;

    // print_matrix("alfa", n, n, alfa_inv[0]);
    // Moltiplichiamo una matrice triangolare superiore con un vettore della forma b + f e mettiamo
    // il risultato in out, stiamo risolvendo il sistema avendo già le matrici pronte
    matrix_triup_mult_vector_with_sum(n, alfa_inv[0], b_block[0], f_block[0], out_block[0]);

    for (int i = 1; i < n; i++) {
        // Calcoliamo (b[i] + f[i]) - l[i] @ y[i - 1] sfruttando il fatto che l è diagonale
        for (int k = 0; k < n; k++) {
            work[k] = b_block[i][k] + f_block[i][k] - l_block[i][k] * out_block[i - 1][k];
        }
        // Ora moltiplichiamo alfa_inv[i] per il vettore calcolato prima
        // y[i] = alfa[i]^-1 @ (b[i] + f[i] - l[i] @ y[i - 1]) = alfa[i]^-1 @ work
        matrix_triup_mult_vector(n, alfa_inv[i], work, out_block[i]);
    }

    for (int i = n - 2; i > -1; i--) {
        // Calcoliamo out[i] = y[i] - gamma[i] @ out[i + 1] tenendo conto che gamma è leggermente modificata
        // per ogni colonna di gamma (l'idea è di partire dalla diagonale in modo da fare una moltiplicazione
        // colonna k per elemento k-esimo del vettore da moltiplicare che è uguale a x_k, dove x = x_1 + x_2 + ... + x_n)
        for (int k = 0; k < n; k++) {

            // Se la diagonale di gamma è zero, allora devo supporre che tutti gli elementi sotto valgono zero
            // e gli elementi a destra si trovano sotto la diagonale
            // Ma dato che stiamo facendo la moltiplicazione matrice per vettore muovendoci colonna per colonna
            // se la diagonale è zero, allora anche tutti quelli sotto lo sono, e quindi non devo ciclare su
            // quei elementi, perchè dovrebbero essere zero (ma non lo sono), invece quando devo utilizzare gli
            // elementi sulla destra della diagonale (che appunto sono salvati sotto) vado a prendere semplicemente
            // quelli sotto in maniera automatica (supponendo che la matrice sia simmetrica)
            if (gamma[i][k * n + k] != 0) {

                double mult = out_block[i + 1][k];
                for (int j = 0; j <= k; j++) {
                    out_block[i][j] -= gamma[i][k * n + j] * mult;
                }
                for (int j = k + 1; j < n; j++) {
                    out_block[i][j] -= gamma[i][j * n + k] * mult;
                }
            }
        }
    }
}

// ridurre d_block che è ancora n x n
// Siamo sicuri che sia meglio calcolarsi l'inversa completamente e non tenersi il risultato della risoluzione
// del sistema lineare? cioè se ho Ax = b, è vero che  b cambia sempre e quindi dovrei rifarlo, ma lapack lascia
// A fattorizzata, quindi si può sfruttare quello? Problema: Devo passare A completamente
int main(int argc, const char *argv[]) {

    // initScene(100);
    // renderScene(0, NULL);

    AppContext ctx;
    ctx.N_STEPS = 50;
    ctx.N_PRINT = 1;
    ctx.FILE_PRINT = true;
    ctx.LOSE_HEAT = false;
    ctx.USE_LAPACK = true;
    int n = 100;
    int N = n * n;

    ctx.n = n;
    ctx.T2_heat = 300;
    ctx.D = 3;
    ctx.dx = 0.6;
    ctx.dt = 1 / 25.0;

    for (int i = 0; i < argc; i++) {
        std::string s(argv[i]);
        if (s == "no-lapack") {
            ctx.USE_LAPACK = false;
            cout << "Non userò lapack" << endl;
        } else if (s == "lose_heat") {
            ctx.LOSE_HEAT = true;
            cout << "Le sorgenti all'interno della griglia non sono infinite" << endl;
        } else if (s == "-r" || s == "-radius") {
            ctx.radius = std::stoi(argv[i + 1], nullptr);
            i++;
        } else if (s == "-t" || s == "-temp" || s == "-temperature") {
            ctx.T2_heat = atof(argv[i + 1]);
            i++;
        } else if (s == "-D" || s == "-diff" || s == "-diffusivity") {
            ctx.D = atof(argv[i + 1]);
            i++;
        }
    }

    ctx.setUpContext(get_lapack_work_size(ctx.n));

    cout << "Diffusività (m^s/s): " << ctx.D << "; dx (m): " << ctx.dx << "; dt (s): " << ctx.dt << ";" << endl;
    cout << "Raggio di T2 (m): " << ctx.radius * ctx.dx << " (o " << ctx.radius << " blocchi); ";
    cout << "Temperatura di T2 (K): " << ctx.T2_heat << endl;

    // double diff = D * dt / (dx * dx); // Coefficiente che moltiplica le temp nella matrice
    // int work_size = get_lapack_work_size(n);
    // if (work_size < n * n * 2) {
    //     work_size = n * n * 2;
    // }
    // int work_solve_size = n;
    // int *ipiv = new int[n];
    // double *work = new double[work_size];
    // double *work_solve = new double[work_solve_size];

    // double **d_block = new double *[n];
    // double **u_block = new double *[n];
    // double **l_block = new double *[n];
    // double **b_block = new double *[n];
    // double **f_block = new double *[n];
    // double **out_block = new double *[n];
    // // Per risolvere a blocchi:
    // double **alfa_inv = new double *[n];
    // double **gamma = new double *[n];

    // for (int i = 0; i < n; i++) {
    //     d_block[i] = new double[n * n]{};
    //     u_block[i] = new double[n]{};
    //     l_block[i] = new double[n]{};
    //     b_block[i] = new double[n]{};
    //     f_block[i] = new double[n]{};
    //     out_block[i] = new double[n]{};

    //     alfa_inv[i] = new double[n * (n + 1)]{};
    //     gamma[i] = alfa_inv[i] + n;
    // }

    // // 5: u,l,b,f,out; 1: d; 1: alfa_inv; n * 8 puntatori a double; + ipiv (n)
    // long int tot_bytes = (work_size + work_solve_size + n * n * 5 + n * n * n * 1 + n * n * (n + 1) + n * 8 + n) * 8;
    // cout << "Per eseguire il programma servono circa " << tot_bytes / (1024.0 * 1024) << " MB" << endl;

    get_block_matrix_lite(n, ctx.diff, ctx.d_block, ctx.u_block, ctx.l_block);

    vector<pair<int, int>> indexes = get_circle_indexes(n, ctx.radius, ctx.T2_heat);
    // vector<pair<int, int>> indexes;
    indexes.push_back(make_pair(85, 500));

    // if (LOSE_HEAT) {
    //     for (auto &cond : indexes) {
    //         int pos = cond.first;
    //         int ind = pos % n;
    //         int big_ind = pos / n;
    //         b_block[big_ind][ind] = cond.second;
    //     }
    // } else {
    set_internal_condition(n, ctx.diff, ctx.d_block, ctx.u_block,
                           ctx.l_block, ctx.b_block, ctx.f_block, indexes);
    // }

    // double temp = 0;
    // for (int k = 0; k < n; k++) {
    //     for (int j = 0; j < n; j++) {
    //         int pos = k * n + j;
    //         if (pos < n || pos > n * (n - 1) || pos % n == 0 || pos % n == n - 1) {
    //             f_block[k][j] = temp * diff;
    //         }
    //         if (pos == 0 || pos == n - 1 || pos == n * n - n || pos == n * n - 1) {
    //             f_block[k][j] += temp * diff;
    //         }
    //     }
    // }

    // vector<pair<int, int>> indexes2;
    // indexes2.push_back(make_pair(11325, 1000));
    // set_internal_condition(n, diff, d_block, u_block, l_block, b_block, f_block, indexes2);

    // for (int i = 0; i < n; i++) {
    //     b_block[49][i] = 1000;
    //     b_block[50][i] = 1000;
    // }

    time_it("", true);
    initialize_alfa_inverse(n, ctx.d_block, ctx.u_block, ctx.l_block,
                            ctx.work_size, ctx.work, ctx.ipiv,
                            ctx.gamma, ctx.alfa_inv, ctx.USE_LAPACK);
    time_it("Tempo per creare le inverse (s): ");

    ///////////////////////////////////////////////////////RIMETTERLO?
    // for (int i = 0; i < n; i++) {
    //     delete[] d_block[i];
    // }
    // delete[] d_block;

    ofstream fileg;
    if (ctx.FILE_PRINT) {
        string nomeg = (string) "heat_2d.dat";
        fileg.open(nomeg, ios::out);
        fileg.precision(1);
        fileg.setf(ios::fixed);
        fileg << n << endl;
    }

    time_it("", true);
    int KSTPES = std::max(ctx.N_STEPS / 10, 1);
    for (int k = 0; k < ctx.N_STEPS; k++) {

        if (ctx.FILE_PRINT && (k % ctx.N_PRINT == 0)) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (k == 0) {
                        fileg << ctx.b_block[i][j] << " ";
                    } else {
                        fileg << ctx.out_block[i][j] << " ";
                    }
                }
            }
            // scrivere \n al posto di endl è decisamente più veloce
            fileg << '\n';
        }

        solve_block_tridiag(n, ctx.u_block, ctx.l_block, ctx.b_block,
                            ctx.f_block, ctx.out_block, ctx.work_solve,
                            ctx.gamma, ctx.alfa_inv);
        // save_snapshot(k, n, d_block, u_block, l_block, b_block, f_block, out_block, gamma, alfa_inv);

        if (k % KSTPES == 0) {
            cout << k << "/" << ctx.N_STEPS << endl;
        }

        for (int i = 0; i < n; i++) {
            double *tmp = ctx.out_block[i];
            ctx.out_block[i] = ctx.b_block[i];
            ctx.b_block[i] = tmp;
        }
    }

    fileg.close();

    // chrono::steady_clock::time_point end = chrono::steady_clock::now();
    time_it("Tempo per eseguire il programma (sec) =", false);
    // cout << "Tempo per eseguire il programma (sec) = " << (chrono::duration_cast<chrono::microseconds>(end - begin).count()) / 1000000.0 << endl;

    ctx.deallocate();

    return 0;
}