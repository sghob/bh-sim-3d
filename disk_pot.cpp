#include <iostream>
#include <cmath>

void r_init(int N, double L, double *r) {
    
    double h = L / (N - 1);
    for (int i = 0; i < N / 2 + 1; i++) {
        double x = i * h;
        r[N/2+i] = x;
        r[N/2-i] = -1.0 * x;
    }
    
    return;
}

void rho_init(int N, double L, double **rho, double A, double k) {
    
    double h = L / (N - 1);
    for (int i = 0; i < N / 2 + 1; i++) {
        double x = i * h;
        double density = A * exp(-1.0 * k * x);
        rho[N/2][N/2+i] = density;
        rho[N/2][N/2-i] = density;
    }

    return;
}

void relax(int N, double L, double **rho, double **phi, int max_iters, double eps) {
    double h = L / (N - 1);

    for (int i = 0; i < max_iters; i++) {
        double f_error = 0.0;

        for (int j = 1; j < N - 1; j++) {
            for (int k = 1; k < N - 1; k++) {
                double phi_new = 0.25 * (phi[j - 1][k] + phi[j + 1][k] + phi[j][k - 1] + phi[j][k + 1] - 4 * 3.14 * 6.67E-11 * h * h * rho[j][k]);
                double error = fabs(phi[j][k] - phi_new);
                if (error > f_error) {
                    f_error = error;
                }
                phi[j][k] = phi_new;
            }
        }
        if (f_error < eps) {
            printf("Converged at iteration: %d\n", i + 1);
            break;
        }
    }
}

int main() {
    int N = 200;
    double r[N];
    double density[N][N];
    double phi[N][N];
    double *a[N];
    double *b[N];

    for (int i = 0; i < N; i++) {
        r[i] = 0.0;
        for (int j = 0; j < N; j++) {
            phi[i][j] = 0.0;
            density[i][j] = 0.0;
        }
    }

    for (int i = 0; i < N; i++) {
        a[i] = density[i];
    }
    for (int k = 0; k < N; k++) {
        b[k] = phi[k];
    }

    double R_g = 2.0 * 3.086E19;
    double A = 3.0;
    double k = 1.0 / R_g;
    double L = 20.0 * R_g;
    double h = L / (N - 1);
    double pot0 = 4.0 * 6.67E-11 * 3.14 * 5E-19 * h * h;
    r_init(N, L, r);
    rho_init(N, L, a, 5E-19, k);
    relax(N, L, a, b, 100000, pot0 * 1E-4);

    for (int j = 0; j < N; j++) {
        printf("%E\n", r[j]);
    }
    printf("%E\n", phi[N/2][N/2]);
    //printf("\n%lf", phi[N/2][N/2]);

    return 1;
}

