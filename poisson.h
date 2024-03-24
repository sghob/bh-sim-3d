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

void z_init(int N, double L, double *z) {
    
    double h = L / (N - 1);
    for (int i = 0; i < N / 2 + 1; i++) {
        double x = i * h;
        z[N/2+i] = x;
        z[N/2-i] = -1.0 * x;
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
                double phi_new = 0.25 * (phi[j - 1][k] + phi[j + 1][k] + phi[j][k - 1] + phi[j][k + 1] - 4 * M_PI * 6.67E-11 * h * h * rho[j][k]);
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
