#include "render.h"
#include "support.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

// Defining lapack methods used
// DSYTRF - Computes the factorization of a real symmetric matrix A using the Bunch-Kaufman diagonal pivoting method
// https://www.netlib.org/lapack/explore-3.2-html/dsytrf.f.html or
// https://www.netlib.org/lapack/explore-html/d8/d0e/group__hetrf_ga431b081d6c9c48af82ec003a7d3070ff.html
//
// DSYTRI - Computes the actual inversion of a real symmetric indefinite using the factorization above
// https://www.netlib.org/lapack/explore-3.2-html/dsytri.f.html
// https://www.netlib.org/lapack/explore-html/da/dfa/group__hetri_gac3efdf326b7a18865d26c5801635beea.html

extern "C" {
extern void dsytri_(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *info);
extern void dsytrf_(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
}

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

// Return the maximum memory used by the routine dsytrf with a specific side
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

// Actually compute the inverse in place on 'mat'
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

// In case we don't want to use lapack there is a manual inversion algorithm
void get_inverse_manual_lite(int n, double *mat, double *work1, double *work2) {

    double *dec = work1;
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

// Calculate mat * b = out
// where 'mat' is an upper triangular matrix of side 'side', 'b' is a vector and 'out' the resulting vector
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

// Calculate mat * (b1+b2) = out
// where 'mat' is an upper triangular matrix with side 'side',
// 'b1' and 'b2' are two vectors and 'out' is the resulting vector 
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
    // int N = n * n;

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

// Set the boundary condition of the matrix
void set_boundary_condition(int lil_side, double diff, double **d_block, double **u_block, double **l_block, double **b_block, double **f_block,
                            double bound_temp) {
    int n = lil_side;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {

            if (i == 0 || k == 0 || i == n - 1 || k == n - 1) {
                f_block[i][k] += diff * bound_temp;
            }
            // Doppie condizioni sui quattro vertici
            if (i == 0 && (k == 0 || k == n - 1)) {
                f_block[i][k] += diff * bound_temp;
            } else if (i == n - 1 && (k == 0 || k == n - 1)) {
                f_block[i][k] += diff * bound_temp;
            }
        }
    }
}

void set_initial_configuration(int lil_side, int preset, double **b_block, double t2, double t1) {

    switch (preset) {
        case 4:
            for (int i = 0; i < lil_side; i++) {
                for (int k = 0; k < lil_side; k++) {
                    b_block[i][k] = (sin(k * 0.25) * 0.5 + 0.5) * t2;
                }
            }
            break;
        case 5: {
            double period = (1.0 / lil_side) * 2 * M_PI;
            for (int i = 0; i < lil_side; i++) {
                for (int k = 0; k < lil_side; k++) {
                    b_block[i][k] = (sin(k * period) * 0.5 + 0.5) * t2;
                }
            }
            break;
        }

        default:
            break;
    }
}

void initialize_alfa_inverse(int lil_side, double **d_block, double **u_block,
                             double **l_block, int work_size, double *work, int *ipiv,
                             double **gamma, double **alfa_inv, bool use_lapack) {
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

void add_preset(int preset, int n, double t2, vector<pair<int, int>> &indexes) {
    if (preset == 1) {
        for (int i = 0; i < 19; i++) {
            int step = n / 20;
            indexes.push_back(make_pair(3 * n + (i + 1) * step, t2 * 0.25));
        }
    } else if (preset == 2) {
        indexes.push_back(make_pair(5000 - 60 + n, 0));
        indexes.push_back(make_pair(5000 - 60 + n - 1, 0));
        indexes.push_back(make_pair(5000 - 60 + n - 2, 0));
        indexes.push_back(make_pair(5000 - 40 + n, 0));
        indexes.push_back(make_pair(5000 - 40 + n + 1, 0));
        indexes.push_back(make_pair(5000 - 40 + n + 2, 0));

        int hside = 10;
        int y1 = n * 0.35;
        int x1 = n * 0.5;
        int hole = 8;
        for (int i = 0; i < hside; i++) {
            int offset = y1 * n + x1 - hside;
            indexes.push_back(make_pair(offset - hole + i, 0));
            // cout << y1 * n - hside + i << "  ";
        }

        for (int i = 0; i < hside; i++) {
            int offset = y1 * n + x1;
            indexes.push_back(make_pair(offset + hole + i, 0));
            // cout << y1 * n + i << "  ";
        }
        cout << "\n";
    } else if (preset == 3) {
        indexes.push_back(make_pair(85, t2 * 2 / 3));
        indexes.push_back(make_pair(250, t2));
    }
}

void next_step(int frame, void *ctx_ptr, float *colors, double *delta_time, bool save_snap = false) {

    AppContext &ctx = *((AppContext *)ctx_ptr);
    if (frame != 0) {
        solve_block_tridiag(ctx.n, ctx.u_block, ctx.l_block, ctx.b_block,
                            ctx.f_block, ctx.out_block, ctx.work_solve,
                            ctx.gamma, ctx.alfa_inv);
    }

    if (colors != NULL) {
        double range = ctx.T2_heat - ctx.T1_heat + 0.01;
        for (int i = 0; i < ctx.n; i++) {
            for (int j = 0; j < ctx.n; j++) {
                colors[i * ctx.n + j] = (ctx.out_block[i][j] - ctx.T1_heat + 0.005) / range;
            }
        }
    }
    if (save_snap) {
        save_snapshot(frame, ctx.n, ctx.d_block, ctx.u_block, ctx.l_block, ctx.b_block,
                      ctx.f_block, ctx.out_block, ctx.gamma, ctx.alfa_inv);
    }

    if (frame != 0) {
        for (int i = 0; i < ctx.n; i++) {
            double *tmp = ctx.out_block[i];
            ctx.out_block[i] = ctx.b_block[i];
            ctx.b_block[i] = tmp;
        }
        *delta_time = ctx.dt;
    } else {
        *delta_time = 0;
    }
}

void reset_application(void *ctx_ptr) {
    AppContext &ctx = *((AppContext *)ctx_ptr);
    ctx.reset();
    int n = ctx.n;
    set_boundary_condition(n, ctx.diff, ctx.d_block, ctx.u_block,
                           ctx.l_block, ctx.b_block, ctx.f_block, ctx.T1_heat);
    if (ctx.LOSE_HEAT) {
        for (auto &cond : ctx.indexes) {
            int pos = cond.first;
            int ind = pos % n;
            int big_ind = pos / n;
            ctx.b_block[big_ind][ind] = cond.second;
        }
    } else {
        set_internal_condition(n, ctx.diff, ctx.d_block, ctx.u_block,
                               ctx.l_block, ctx.b_block, ctx.f_block, ctx.indexes);
    }

    set_initial_configuration(ctx.n, ctx.preset, ctx.b_block, ctx.T2_heat, ctx.T1_heat);
    for (int i = 0; i < n; i++) {
        copy(ctx.b_block[i], ctx.b_block[i] + n, ctx.out_block[i]);
    }
}

void on_canvas_click(void *ctx_ptr, int x, int y, int width, int height) {
    AppContext &ctx = *((AppContext *)ctx_ptr);
    int xi = x * 1.0 / width * ctx.n;
    int yi = y * 1.0 / height * ctx.n;
    int index = yi * ctx.n + xi;
    //usiamo il b block in quanto dopo next_step è tutto in bblock l'output
    double value = ctx.b_block[yi][xi]; // al contrario
    // double temp = value * (ctx.T2_heat - ctx.T1_heat) + ctx.T1_heat;
    double temp = value;

    std::cout << "Quadrato: " << xi << " " << yi << " (" << index << "),  Temp: " << temp << " K" << std::endl;
}

int main(int argc, char **argv) {

    AppContext ctx;
    ctx.N_STEPS = 10;
    ctx.N_PRINT = 1;
    ctx.FILE_PRINT = true;
    ctx.LOSE_HEAT = false;
    ctx.USE_LAPACK = true;
    int n = 100;
    int N = n * n;

    ctx.n = n;
    ctx.T1_heat = 0;
    ctx.T2_heat = 300;
    // Copper 111e-6, Oro: 127e-6, Acqua: 0.143e-6, Ferro: 23e-6
    ctx.D = 11e-6;
    ctx.dx = 0.1 / n; // 10 cm divisi in n pezzetti
    ctx.dt = 0.022;

    for (int i = 0; i < argc; i++) {
        std::string s(argv[i]);
        if (s == "no-lapack") {
            ctx.USE_LAPACK = false;
            cout << "Non userò lapack" << endl;
        } else if (s == "-lose_heat") {
            ctx.LOSE_HEAT = true;
            cout << "Le sorgenti all'interno della griglia non sono infinite" << endl;
        } else if (s == "-r" || s == "-radius") {
            ctx.radius = std::stoi(argv[i + 1], nullptr);
            i++;
        } else if (s == "-t2" || s == "-temp2" || s == "-temperature2") {
            ctx.T2_heat = atof(argv[i + 1]);
            i++;
        } else if (s == "-t1" || s == "-temp1" || s == "-temperature1") {
            ctx.T1_heat = atof(argv[i + 1]);
            i++;
        } else if (s == "-D" || s == "-diff" || s == "-diffusivity") {
            ctx.D = atof(argv[i + 1]);
            i++;
        } else if (s == "-time" || s == "-dt") {
            ctx.dt = atof(argv[i + 1]);
            i++;
        } else if (s == "-dx") {
            ctx.dx = atof(argv[i + 1]);
            i++;
        } else if (s == "-preset") {
            ctx.preset = atoi(argv[i + 1]);
            i++;
        }
    }

    ctx.setUpContext(get_lapack_work_size(ctx.n));

    cout << "Diffusività (m^2/s): " << ctx.D << "; dx (m): " << ctx.dx << "; dt (s): " << ctx.dt << ";" << endl;
    cout << "Raggio di T2 (m): " << ctx.radius * ctx.dx << " (o " << ctx.radius << " blocchi); ";
    cout << "Temperatura di T1 (K): " << ctx.T1_heat << "; ";
    cout << "Temperatura di T2 (K): " << ctx.T2_heat << endl;

    get_block_matrix_lite(n, ctx.diff, ctx.d_block, ctx.u_block, ctx.l_block);

    vector<pair<int, int>> indexes = get_circle_indexes(n, ctx.radius, ctx.T2_heat);
    add_preset(ctx.preset, n, ctx.T2_heat, indexes);

    ctx.indexes = indexes;

    reset_application((void *)&ctx);

    time_it("", true);
    initialize_alfa_inverse(n, ctx.d_block, ctx.u_block, ctx.l_block,
                            ctx.work_size, ctx.work, ctx.ipiv,
                            ctx.gamma, ctx.alfa_inv, ctx.USE_LAPACK);
    time_it("Tempo per creare le inverse (s): ");

    char str_info[128];
    double dx_cm = ctx.dx * 100;
    double D_e6 = ctx.D * 1e6;
    snprintf(str_info, 128, " D: %g (mm^2/s)   n: %d   dx: %g (cm)   dt: %g (ms)", D_e6, ctx.n, dx_cm, ctx.dt * 1e3);

    initScene(n, next_step, reset_application, (void *)&ctx, on_canvas_click, string(str_info));
    renderScene(&argc, argv);

    // double tm;
    // for (int i = 0; i < 10; i++) {
    //     next_step(i, (void *)&ctx, NULL, &tm, true);
    // }

    cout << "Fine" << endl;
    ctx.deallocate();

    return 0;
}