#include <vector>
#include <cmath>
#include "integration.h"

//double bulge_int = 16.1533; //in kpc^3
//double bulge_int = 0.0478; //in kpc^3
double bulge_int = 0.189;
//placeholder for now, needs to actually be integrated
double rhob_norm = 3.58E-19;
double kpc_m = 3.086E19;

void r_init(int N, double L, double *r) {
    
    double h = L / (N - 1);
    for (int i = 0; i < N / 2 + 1; i++) {
        double x = i * h;
        r[N/2+i] = x;
        r[N/2-i] = -1.0 * x;
    }

    return;
}

/*
void r2_init(int N1, int N2, double L1, double L2, double *r) {
    
    double h1 = L1 / (N1 - 1);
    double h2 = L2 / (N2 - 1);
    for (int i = 0; i < N1 / 2; i++) {
        double x = i * h1;
        for (int j = 0; j < N2 / 2; j++) {
            double y = j * h2;
            r[N1/2+i] = x;
            r[N1/2-i] = -1.0 * x;
            r[N2/2 + j] = y;
            r[N2/2 - j] = -1.0 * y;
        }
    }

    return;
}
*/
void z_init(int N, double L, double *z) {
    
    double h = L / (N - 1);
    for (int i = 0; i < N / 2 + 1; i++) {
        double x = i * h;
        z[N/2+i] = x;
        z[N/2-i] = -1.0 * x;
    }
    
    return;
}

double rho_stellar(double r, double z, const std::vector<double>& params) {
    //return 2 * params[1] * params[0] * exp(-1.0 * r / params[1]) * ((0.5 / (2*params[2])) * exp(-1.0 * abs(z) / params[2]) + (0.5 / (2*params[3])) * exp(-1.0 * abs(z) / params[3]));
    return params[0] * exp(-1.0 * r / params[1]) * ((0.5 / (2*params[2])) * exp(-1.0 * abs(z) / params[2]) + (0.5 / (2*params[3])) * exp(-1.0 * abs(z) / params[3]));

}

void rhos_init(int N, double L, double **rhos, const std::vector<double>& params) {
    double h = L / (N - 1);
    for (int i = 0; i < N/2; i++) {
        double r = i * h;
        for (int j = 0; j < N/2; j++) {
            double z = j * h;
            double density;
            density = rho_stellar(r, z, params);
            //printf("%E\n", density);
            rhos[N/2+j][N/2+i] = density;
            rhos[N/2-j][N/2+i] = density;
            rhos[N/2+j][N/2-i] = density;
            rhos[N/2-j][N/2-i] = density;
            
        }
    }

    return;
}

void rhos_init2(int N, int N2, double L, double L2, double **rhos, const std::vector<double>& params) {
    double h = L / (N - 1);
    double h2 = L2 / (N2 - 1);
    for (int i = 0; i < N/2; i++) {
        double r = i * h;
        for (int j = 0; j < N2/2; j++) {
            double z = j * h2;
            double density;
            density = rho_stellar(r, z, params);
            //printf("%E\n", density);
            rhos[N2/2+j][N/2+i] = density;
            rhos[N2/2-j][N/2+i] = density;
            rhos[N2/2+j][N/2-i] = density;
            rhos[N2/2-j][N/2-i] = density;
            
        }
    }

    return;
}


double rho_gas(double r, double z, double r_g, double z_g, double A) {
    return A * pow(M_E, -1.0 * r / r_g - abs(z) / z_g);
}

double rho_gas_f(double r, double z, const std::vector<double>& params) {
    //cannot integrate this function normally for some reason, need to use evenness symmetry
    return 2 * M_PI * r * params[0] * pow(M_E, (-1.0 * r / params[1]) - 1.0 * (z / params[2])); //note, no abs
}

void rhog_init(int N, double L, double **rhog, double r_g, double z_g, double A) {
    double h = L / (N - 1);
    for (int i = 0; i < N/2; i++) {
        double r = i * h;
        for (int j = 0; j < N/2; j++) {
            double z = j * h;
            double density;
            density = rho_gas(r, z, r_g, z_g, A);
            
            rhog[N/2+j][N/2+i] = density;
            rhog[N/2-j][N/2+i] = density;
            rhog[N/2+j][N/2-i] = density;
            rhog[N/2-j][N/2-i] = density;
            
        }
    }

    return;
}

void rhog2_init(int N, int N2, double L, double L2, double **rhog, double r_g, double z_g, double A) {
    double h = L / (N - 1);
    double h2 = L2 / (N2 - 1);
    for (int i = 0; i < N/2; i++) {
        double r = i * h;
        for (int j = 0; j < N2/2; j++) {
            double z = j * h2;
            double density;
            density = rho_gas(r, z, r_g, z_g, A);
            
            rhog[N2/2+j][N/2+i] = density;
            rhog[N2/2-j][N/2+i] = density;
            rhog[N2/2+j][N/2-i] = density;
            rhog[N2/2-j][N/2-i] = density;
            printf("%E\n", density);
        }
    }

    return;
}

void rho_init(int N, double L, double **rho, double A, double k) {
    
    double h = L / (N - 1);
    for (int i = 0; i < N / 2 + 1; i++) {
        double x = i * h;
        double density = A * exp(-1.0 * k * x); //NOTE THAT WE DON"T PUT "FLAT" CORRECTION
        rho[N/2][N/2+i] = density;
        rho[N/2][N/2-i] = density;
        /*
        rho[N/2 + 1][N/2+i] = density;
        rho[N/2 + 1][N/2-i] = density;
        rho[N/2 - 1][N/2+i] = density;
        rho[N/2 - 1][N/2-i] = density;
        rho[N/2 + 2][N/2+i] = density;
        rho[N/2 + 2][N/2-i] = density;
        rho[N/2 - 2][N/2+i] = density;
        rho[N/2 - 2][N/2-i] = density;
        rho[N/2 + 3][N/2+i] = density;
        rho[N/2 + 3][N/2-i] = density;
        rho[N/2 - 3][N/2+i] = density;
        rho[N/2 - 3][N/2-i] = density;
        rho[N/2 + 4][N/2+i] = density;
        rho[N/2 + 4][N/2-i] = density;
        rho[N/2 - 4][N/2+i] = density;
        rho[N/2 - 4][N/2-i] = density;
        */
    }

    return;
}


double rho_bulge(double r, double z, double rho0, double qb, double ab, double rb) {
    double m = pow(r * r + (z * z) / (qb * qb), 0.5);
    return rho0 * pow(m/ab, -1.8) * pow(M_E, -1.0 * m * m / (rb * rb));
}

double rho_bulge_f(double r, double z, const std::vector<double>& params) {
    double m = pow(r * r + (z * z) / (params[1] * params[1]), 0.5);
    return params[0] * pow(m / params[2], -1.8) * pow(M_E, -1.0 * m * m / (params[3] * params[3]));
}

double fs_v(double vc, double vs, double vesc) { //gamma=1.6
    if (vs > pow(2,0.5) * vc) {
        return 0.0;
    } else if (vs > 0 && vs < pow(2,0.5) * vc) {
        return (0.169 / pow(vc, 2*1.6)) * pow(2 * vc * vc - vs * vs, 1.6-1.5);
    } else {
        return 0.0;
    }
}

void bulge_init(int N, double L, double **rho_b, double rho0, double qb, double Rb) {
    double ab = 1.0 * kpc_m;//Rb / 4.0;               //EDIT THIS FOR VARIABLE SCALE R
    double rb = 1.9 * kpc_m;//* Rb / 4.0;
    //double ab = 1.0 * Rb / 4.0;
    //double rb = 1.0 * Rb / 4.0;
    double h = L / (N - 1);
    for (int i = 0; i < N / 2; i++) {
        double r = i * h;
        for (int j = 0; j < N / 2; j++) {
            double z = j * h;
            double density;
            if (r == 0 && z == 0) {
                density = 5.16E-16 * (rho0 != 0.0);
                //density = rho_bulge(r + h, z, rho0, qb, ab, rb);
            } else if (r > 3.086E19 || z > 3.086E19) {
                density = 0.0;
            } else {
                density = rho_bulge(r, z, rho0, qb, ab, rb);
            }

            rho_b[N/2+j][N/2+i] = density;
            rho_b[N/2-j][N/2+i] = density;
            rho_b[N/2+j][N/2-i] = density;
            rho_b[N/2-j][N/2-i] = density;
            
        }
    }

    return;
}

void bulge2_init(int N, int N2, double L, double L2, double **rho_b, double rho0, double qb, double Rb) {
    double ab = 1.0 * kpc_m;//Rb / 4.0;               //EDIT THIS FOR VARIABLE SCALE R
    double rb = 1.9 * kpc_m;//* Rb / 4.0;
    //double ab = 1.0 * Rb / 4.0;
    //double rb = 1.0 * Rb / 4.0;
    double h = L / (N - 1);
    double h2 = L2 / (N2 - 1);
    for (int i = 0; i < N / 2; i++) {
        double r = i * h;
        for (int j = 0; j < N2 / 2; j++) {
            double z = j * h2;
            double density;
            if (r == 0 && z == 0) {
                density = 8.048E-15 * (rho0 != 0.0); //DANGER DANGER DANGER DANGER, NEED TO ADJUST FOR OMEGA2
            } else if (r > 3.086E19 || z > 3.086E19) {
                density = 0.0;
            } else {
                density = rho_bulge(r, z, rho0, qb, ab, rb);
            }

            rho_b[N2/2+j][N/2+i] = density;
            rho_b[N2/2-j][N/2+i] = density;
            rho_b[N2/2+j][N/2-i] = density;
            rho_b[N2/2-j][N/2-i] = density;
            
        }
    }

    return;
}

double rho_halo(double r, double z, const std::vector<double>& params) {
    double m = pow(r * r + z * z, 0.5);
    return params[0] * (params[1] / m) * pow(1 + m / params[1], -2);
}

void halo_init(int N, double L, double **rhoh, const std::vector<double>& params) {
    double h = L / (N - 1);
    for (int i = 0; i < N/2; i++) {
        double r = i * h;
        for (int j = 0; j < N/2; j++) {
            double z = j * h;
            double density;
            if (r == 0.0 && z == 0.0) {
                continue;
            }
            density = rho_halo(r, z, params);
            //printf("%E\n", density);
            rhoh[N/2+j][N/2+i] = density;
            rhoh[N/2-j][N/2+i] = density;
            rhoh[N/2+j][N/2-i] = density;
            rhoh[N/2-j][N/2-i] = density;
            
        }
    }

    return;

}

void relax(int N, double L, double **rho, double **phi, int max_iters, double eps) {
    double h = L / (N - 1);

    for (int i = 0; i < max_iters; i++) {
        double f_error = 0.0;

        for (int j = 1; j < N - 1; j++) {
            for (int k = 1; k < N - 1; k++) {
                //double phi_new = 0.25 * (phi[j - 1][k] + phi[j + 1][k] + phi[j][k - 1] + phi[j][k + 1] - 4 * M_PI * 6.67E-11 * h * h * rho[j][k]);
                double r = abs(N / 2 - k) * h;
                double phi_new = (/*(phi[j][k + 1] - phi[j][k - 1]) / (2 * r * h) + */(1 / (h * h)) * (phi[j][k + 1] + phi[j][k - 1]) + (1 / (h * h)) * (phi[j+1][k] + phi[j-1][k]) - 4 * M_PI * 6.67E-11 * rho[j][k]) / (2 / (h * h) + 2/ (h * h));
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

void relax_amr(int N, int M, double L, double W, double **rho, double **phi, int max_iters, double eps) {
    double h = L / (N - 1);
    double hw = W / (M - 1);
    for (int i = 0; i < max_iters; i++) {
        double f_error = 0.0;

        for (int j = 1; j < M-1; j++) {
            for (int k = 1; k < N-1; k++) {
                //double phi_new = 0.25 * (phi[j - 1][k] + phi[j + 1][k] + phi[j][k - 1] + phi[j][k + 1] - 4 * M_PI * 6.67E-11 * h * hw * rho[j][k]);
                double phi_new = (/*(phi[j][k + 1] - phi[j][k - 1]) / (2 * r * h) + */(1 / (h * h)) * (phi[j][k + 1] + phi[j][k - 1]) + (1 / (hw * hw)) * (phi[j+1][k] + phi[j-1][k]) - 4 * M_PI * 6.67E-11 * rho[j][k]) / (2 / (h * h) + 2/ (hw * hw));
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
