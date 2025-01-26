#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <cmath> //need this for bessel and few other things
#include "poisson.h" //for computing galaxy potential
//#include "integration.h" //2D scaling integrals
#include <gsl/gsl_interp2d.h> //interpolate lookup table

//math constants
double PI = 3.14159265;
double e = 2.718281828;
//conversion factors
double yr_s = 3.154E7;
double sl_kg = 1.99E30;
double pc_m = 3.086E16;

//physics constants
double G = 6.67E-11;
double kb = 1.38E-23;
double m_p = 1.67E-27;
double c = 3.0E8;


double gama = 1.8;/*>0.6 */
double gama_e = 1.8;
double alpha = 4.0;
double alpha_b = 1.8;
double q_b = 0.6;

//parameters of sim
double M1 = 1E6; //solar mass
double q = 1.0 / 9.0;
double M2 = M1 * q;
double vg = -0.7; //think old plot can be restored if this is 0, breaks stdf tho

double R_sd = log10(M1 / pow(10,5)) * 1000 * pc_m;
double R_gd = 2.0 * R_sd;

double rho_bm = 1E8;
double M_s = 2.0E30;
double M_e = 5.97E24;
double n_gd = 300 * pow(100, 3);
double fg = 0.5;
//double z_thick = R_gd / 20.0;
double z_thick = 100.0 * pc_m;
//double M2 = 1E6;
double sigma_star = 1E5;

double rho_s0 = 5.01E-19; //subject to change, needs to follow gas frac definition involving integration

// Define the state of the system
struct State { //define phase space coord
    double r;
    double rdot;
    double theta;
    double thetadot;
    double z;
    double zdot;
};

//DECLARE DENSITY & POTENTIAL GRID
//Update 11/09/24: including AMR grids for Omega2 (disk boundary)
double *r_arr;
double *z_arr;
double *r2_arr;
double *z2_arr;

double **rho;
double **rho2;

double **rho_g;
double **rho_g2;

double **rho_s;
double **rho_s2;

double **rho_st;
double **rho_st2;

double **rho_bg;
double **rho_bg2;

double *rho_flat;
double *rho_flat2;

double *rho_g_flat;
double *rho_g_flat2;

double *rho_s_flat;
double *rho_s_flat2;

double *rho_st_flat;
double *rho_st_flat2;

//come back to dm halo later
double **rho_h;
double *rho_h_flat;

double **phi;
double **phi2;

double **phi_b;
double *phi_flat;
double *phi_flat2;

// Create an interpolation object
gsl_interp2d *interp;
gsl_interp_accel *xacc;
gsl_interp_accel *yacc;

gsl_interp2d *interp_r;
gsl_interp_accel *xacc_r;
gsl_interp_accel *yacc_r;

gsl_interp2d *interp_amr;
gsl_interp_accel *xacc_amr;
gsl_interp_accel *yacc_amr;

//Solver helper functions

double sgn(double x) {
    return (x > 0) - (x < 0);
}

double intStep(double a, double b, double (*f)(double)) { //numerical integration step using simpson's rule
    return ((b - a) / 6.0) * (f(a) + 4.0 * f((a + b) / 2.0) + f(b));
}

double numInt(double a, double b, double (*f) (double)) { //composite numerical integration with simpson's rule
    double s = 0.0;
    double n = 100.0;
    double x = a;
    double dx = (b - a) / n;
    for (int i = 0; i < n; i++) {
        s += intStep(x, x + dx, f);
        x += dx;
    }
    return s;
}

double grad(double z, double (*f) (double)) {
    double dz = z * 1E-8;
    return (f(z + dz) - f(z)) / dz;
}

/* Relevant potentials/forces (disk, bulge, stellar DF, gas DF, velocity distribution, etc. */


double rho_gd(double r) {
    return n_gd * m_p * pow(e, -1.0 * r / (R_gd));
}

double rho_rz(double r, double z) {
    return gsl_interp2d_eval(interp_r, r_arr, z_arr, rho_g_flat, r, z, xacc_r, yacc_r);
}


double rhos_rz(double r, double z) {
    return gsl_interp2d_eval(interp_r, r_arr, z_arr, rho_s_flat, r, z, xacc_r, yacc_r);
}

double rhoh_rz(double r, double z) {
    return gsl_interp2d_eval(interp_r, r_arr, z_arr, rho_h_flat, r, z, xacc_r, yacc_r);
}

//think might need to replace rho_st_flat with rho_b_flat?
double rhost_rz(double r, double z) {
    //printf("%E\n", gsl_interp2d_eval(interp_r, r_arr, z_arr, rho_st_flat, r, z, xacc_r, yacc_r));
    return gsl_interp2d_eval(interp_r, r_arr, z_arr, rho_st_flat, r, z, xacc_r, yacc_r);
}

double sigma_gd(double r) {
    return rho_gd(r) * (2 * z_thick);
}

double gasDFb(double r, double z) {
    //return 4 * PI * pow((G * M2 * sl_kg), 2) * rho_gd(r);
    //double rhot = gsl_interp2d_eval(interp_r, r_arr, z_arr, rho_flat, r, z, xacc_r, yacc_r);
    //printf("%E\t%E\n", rhot, rho_gd(r));
    //printf("%E\n", 4 * PI * pow((G * M2 * sl_kg), 2) * rho_gd(r));
    //return 4 * PI * pow((G * M2 * sl_kg), 2) * rho_gd(r);
    //printf("%E\n", rho_rz(r,z));
    return 4 * PI * pow((G * M2 * sl_kg), 2) * rho_rz(r, z);
    //return 4 * PI * pow((G * M2 * sl_kg), 2) * rhot; //DONT DIVIDE BY V^2 YET
}

//UPDATE ALL DFS WITH Z PARAM PASSING 4/20
double gasDFr(double r, double rdot, double tdot, double z, double cs, double vd) {
    
    double rt = r * tdot;
    double v = pow((rt + vd) * (rt + vd), 0.5);
    double mach = v / cs;
    if (mach < 1E-4) {
        //printf("zero vel");
        return 0;
    }
    
    double I_r;

    if (mach < 1.1) {
        I_r = pow(mach, 2) * pow(10, 3.51 * mach - 4.22);
    } else if (mach < 4.4) {
        I_r = 0.5 * log(9.33 * pow(mach, 2) * (pow(mach, 2) - 0.95));
    } else {
        I_r = 0.3 * pow(mach, 2);
    }
    //return 0.0;
    return (-1.0 * gasDFb(r, z) * I_r) / (v * v) ;
}

//INVESTIGATE KICK
//DISK IS THIN RN, BULGE IS OFF, Z0=0.5pc, THERE IS AN ISSUE WITH TERMINATION CONDITION (FIXED)
//9/21/24: Issue with rdot causing -NaN, but *super* close to lily's result (FIXED)
// also, investigate behavior in <100 pc regime where waviness occurs for some reason idfk
double gasDFt(double r, double rdot, double tdot, double z, double cs, double vd) {

    
    double rt = r * tdot;
    //printf("%E\n", rt + vd);
    double v = pow((rt + vd) * (rt + vd), 0.5);
    double mach = v / cs;
    if (mach < 1E-4) {
        //printf("zero vel");
        return 0;
    }
    
    double I_t = 0.0;
    double r_min = r / 10.0;
    if (mach < 1.0 && mach > 0.0) {
        I_t = 0.7706 * log((1 + mach) / (1.0004 - 0.9185 * mach)) - 1.4703 * mach;
        //printf("%E", I_t);
    } else if (mach < 4.4) {
        I_t = log(330 * 10.0 * pow(mach, -9.58) * pow(mach - 0.71, 5.72));
    } else if (mach > 4.4) {
        I_t = log(10 / (0.11 * mach + 1.65));
    } else {
        I_t = 0.0;
    }
    //printf("\n%E\t%E", I_t, mach);

    //return 0.0;
    //printf("%E\n", (-1.0 * gasDFb(r, z) * I_t) / (v * v));
    return (-1.0 * gasDFb(r, z) * I_t) / (v * v);
}

double gasDFz(double r, double z, double zdot, double cs) {

    double mach = abs(zdot) / cs;

    if (mach < 1E-4) {
        //printf("zero vel");
        //printf("%E\n", mach);
        return 0.0;
    }
    
    double I = 0.0;
    //printf("%E\n", mach);

    if (mach > 0.5 && mach < 0.99999) { //GOTCHA BITCH
        I = (0.5 * log((1-mach) / (1+mach)) - mach); //dfz "kicking" the BH
    } else if (mach > 1.0) {
        I = (0.5 * log(1 - 1 / (mach * mach)) + 10);
    } else if (mach < 0.1) {
        I = mach * mach * mach / 3.0;
    } else {
        I = 0.0;
    }
    //printf("%E\t%E\n", mach, (gasDFb(r, z) * I) / (zdot * zdot));
    return (gasDFb(r, z) * I) / (zdot * zdot);
    //return 0.0;
}

double oomDF(double r, double v, double cs) {
    return 4 * PI * G * G * M2 * M2 * n_gd * m_p / (cs * cs);
}

double df_coeff(double M2, double rho_s, double v) {
    //printf("v: %E\n", v);
    if (v < 1E-4) {
        return 0.0;
    }
    return -4 * M_PI * G * G * pow(M2 * sl_kg, 2) * 1.0 * rho_s / pow(v, 3.0);
}

double df_dm_test(double r, double z, double v) {
    if (v < 1E-4) {
        return 0.0;
    }
    return -4 * M_PI * G * G * pow(M2 * sl_kg, 2) * 1.0 * rhoh_rz(r,z) / pow(0.5 * v, 3.0);
}

double df_slow_r_disk(double v_g, double v, double r, double z, double rdot) {
    if (v > abs(v_g)) {
        //return 0.01 * rdot * df_coeff(M2, rho_rz(r, z), v) * 4 * M_PI * v_g * v_g * log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
        //return 0.0;
        return rdot * df_coeff(M2, rhos_rz(r, z), v) * log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
    } else {
        return 0.0;
    }
}

double df_slow_t_disk(double v_g, double v, double r, double z, double thetadot) {
    if (v > abs(v_g)) {
        //printf("%E\n", (M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
        //return 1.0 * (r * thetadot) * df_coeff(M2, 0.1 * rho_rz(r,z), v) * 4 * M_PI * v_g * v_g * log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
        return 1.0 * (r * thetadot) * df_coeff(M2, 0.1 * rhos_rz(r,z), v) * log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
    } else {
        return 0.0;
    }
}

double df_slow_z_disk(double v_g, double v, double r, double z, double zdot) {
    if (v > abs(v_g)) {
        //return zdot * df_coeff(M2, rho_rz(r,z), v) * 4 * M_PI * v_g * v_g * log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
        //return 0.0;
        return zdot * df_coeff(M2, rhos_rz(r,z), v) * log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
    } else {
        return 0.0;
    }
}

double df_fast_r_disk(double v_g, double v, double r, double z, double rdot) {
    if (v < abs(v_g)) {
        //return 0.01 * rdot * df_coeff(M2, rho_rz(r, z), v) * 4 * M_PI * v_g * v_g * log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
        return rdot * df_coeff(M2, rhos_rz(r, z), v) * (log((v_g + v) / (v_g - v)) - 2 * v / v_g);
        //return 0.0;
    } else {
        return 0.0;
    }
}

double df_fast_t_disk(double v_g, double v, double r, double z, double thetadot) {
    if (v < abs(v_g)) {
        //return 0.01 * rdot * df_coeff(M2, rho_rz(r, z), v) * 4 * M_PI * v_g * v_g * log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
        return r * thetadot * df_coeff(M2, rhos_rz(r, z), v) * (log((v_g + v) / (v_g - v)) - 2 * v / v_g);
    } else {
        return 0.0;
    }
}

double df_fast_z_disk(double v_g, double v, double r, double z, double zdot) {
    if (v < abs(v_g)) {
        //return 0.01 * rdot * df_coeff(M2, rho_rz(r, z), v) * 4 * M_PI * v_g * v_g * log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
        //return 0.0;
        return zdot * df_coeff(M2, rhos_rz(r, z), v) * (log((v_g + v) / (v_g - v)) - 2 * v / v_g);
    } else {
        return 0.0;
    }
}


//OLD
//Potential in the plane of the disk (analytic)
double testPhi(double r) { //test potential
    double a = 1.0 / R_gd;
    //printf("\n%E", a * r / 2.0);
    return PI * G * sigma_gd(0) * r * (std::cyl_bessel_i(1.0, a * r / 2.0) * std::cyl_bessel_k(0.0, a * r / 2) - std::cyl_bessel_i(0.0, a * r / 2.0) * std::cyl_bessel_k(1.0, a * r / 2.0));
}

double phi_rz(double r, double z) {
    
    return gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc);
}

double gF(double r) { //SPECIFIC FORCE
    double dr = r * 1E-8;
    return -1.0 * (testPhi(r + dr) - testPhi(r)) / dr;
    //return -1.0 * G * M1 * M2 * sl_kg * sl_kg / (r * r);
}

double g_fr(double r, double z) {
    double dr = r * 1E-8;
    //printf("%E\n", gsl_interp2d_eval_deriv_x(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc));
    //printf("%E\n", -1.0 * (testPhi(r + dr) - testPhi(r)) / dr);
    //return -1.0 * (testPhi(r + dr) - testPhi(r)) / dr;
    //printf("%E\n", r / pc_m);
    //printf("\n%E", gsl_interp2d_eval_deriv_x(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc));
    if (abs(z) < 0 * 0.04 * R_gd) {
        //printf("%E\n", (-1.0 * gsl_interp2d_eval_deriv_y(interp_amr, r2_arr, z2_arr, phi_flat2, r, z, xacc_amr, yacc_amr)) / (-1.0 * gsl_interp2d_eval_deriv_y(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc)) );
        //printf("o2\n");
        return -1.0 * gsl_interp2d_eval_deriv_x(interp_amr, r2_arr, z2_arr, phi_flat2, r, z, xacc_amr, yacc_amr);
    }
    return -1.0 * gsl_interp2d_eval_deriv_x(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc);
    //return 0.0;
}

double g_fz(double r, double z) {
    double dz = z * 1E-8;
    //printf("%E\n", z / pc_m);
    if (abs(z) < 0 * 0.04 * R_gd) {
        //printf("%E\n", (-1.0 * gsl_interp2d_eval_deriv_y(interp_amr, r2_arr, z2_arr, phi_flat2, r, z, xacc_amr, yacc_amr)) / (-1.0 * gsl_interp2d_eval_deriv_y(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc)) );
        return -1.0 * gsl_interp2d_eval_deriv_y(interp_amr, r2_arr, z2_arr, phi_flat2, r, z, xacc_amr, yacc_amr);
    }
    return -1.0 * gsl_interp2d_eval_deriv_y(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc);
    //return -1.0 * (disk_pot(r, z + dz) - disk_pot(r, z)) / dz;
}

double vc(double r, double z) { 

    return pow(r * -1.0 * g_fr(r, 0.0), 0.5);
}

//ISSUE WITH FS_V
//params[3] is v
double bulge_integrand(double vs, const std::vector<double>& params) {
    //double v_esc = pow(gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, params[0], params[1], xacc, yacc), 0.5);
    //printf("%E\t%E\t%E\n", fs_v(params[2], vs, params[5]), vs, params[5]);
    //if (params[3] > vs) {
    if ((1 / (q * params[4] * params[4])) * (params[3] * params[3] - vs * vs) > 1.0) {
        return 4 * M_PI * fs_v(params[2], vs, params[5]) * vs * vs * log((1 / (q * params[4] * params[4])) * (params[3] * params[3] - vs * vs));
    } else if (params[3] < vs && params[3] > 0.0) {
        return 4 * M_PI * fs_v(params[2], vs, params[5]) * vs * vs * (log((vs + params[3]) / (vs - params[3])) - 2 * params[3] / vs);
    } else {
        return 0.0;
    }
}

double df_bulge_r(double r, double z, double v, double rdot) {
    double v_esc = pow(-2.0 * gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc), 0.5);
    double df_slow = df_coeff(M2, rhost_rz(r, z), v) * integrate1D(bulge_integrand, 0, v, 20, {r, z, vc(r, z), v, sigma_star, v_esc});
    double df_fast = df_coeff(M2, rhost_rz(r, z), v) * integrate1D(bulge_integrand, v, v_esc, 20, {r, z, vc(r, z), v, sigma_star, v_esc});
    return rdot * (df_slow + df_fast);
    //return 0.0;
}
//ADD IN FAST DF BACK IN LATER (negligible for most purposes, but for completeness)
double df_bulge_t(double r, double z, double v, double tdot) {
    double v_esc = pow(-2.0 * gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc), 0.5);
    //printf("%E\t", integrate1D(bulge_integrand, 0, v, 20, {r, z, vc(r), v, sigma_star, v_esc}));
    //printf("%E\n", integrate1D(bulge_integrand, v, v_esc, 20, {r, z, vc(r), v, sigma_star, v_esc}));
    if (v > 0.0 * sigma_star) {
        double df_slow = 1.0 * df_coeff(M2, rhost_rz(r, z), v) * integrate1D(bulge_integrand, 0, v, 20, {r, z, vc(r, z), v, sigma_star, v_esc});
        double df_fast = 1.0 * df_coeff(M2, rhost_rz(r, z), v) * integrate1D(bulge_integrand, v, v_esc, 20, {r, z, vc(r, z), v, sigma_star, v_esc});
        //return r * tdot * (df_slow);
        //printf("%E\n", r * tdot * (df_slow + df_fast));
        //printf("%E\n", df_fast/df_slow);
        
        return r * tdot * (df_slow + df_fast);
    } else {
        return 0.0;
    }
}

double df_bulge_z(double r, double z, double v, double zdot) {
    double v_esc = pow(-2.0 * gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc), 0.5);
    double df_slow = df_coeff(M2, rhost_rz(r, z), v) * integrate1D(bulge_integrand, 0, v, 20, {r, z, vc(r, z), v, sigma_star, v_esc});
    double df_fast = df_coeff(M2, rhost_rz(r, z), v) * integrate1D(bulge_integrand, v, v_esc, 20, {r, z, vc(r, z), v, sigma_star, v_esc});
    return 1.0 * zdot * (df_slow + df_fast);
}

double df_dm_z(double r, double z, double v, double zdot) {
    double v_esc = pow(-2.0 * gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc), 0.5);
    double df_slow = df_coeff(M2, rhoh_rz(r, z), v) * integrate1D(bulge_integrand, 0, v, 20, {r, z, vc(r, z), v, sigma_star, v_esc});
    double df_fast = df_coeff(M2, rhoh_rz(r, z), v) * integrate1D(bulge_integrand, v, v_esc, 20, {r, z, vc(r, z), v, sigma_star, v_esc});
    return zdot * (df_slow + df_fast);
}

//total force governing EoM
double fr(double r, double rdot, double thetadot, double z, double zdot, double cs) {
    
    //double v = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0), 0.5);
    double v = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0) + pow(zdot, 2.0), 0.5);
    //printf("\n%E", (1.0 / (M2 * sl_kg)) * ((-1.0 * G * M1 * M2 * pow(sl_kg, 2)) / pow(r , 2.0) + (M2 * sl_kg) * r * pow(thetadot, 2) + gasDFr(r, v, cs)));
    //printf("\n%E\t%E\n", gasDFr(r, v, cs), (-1.0 * G * M1 * M2 * pow(sl_kg, 2)) / (pow(r , 2.0) + 1.0E-10));
    //printf("\n%E\t%E", (M2 * sl_kg) * r * pow(thetadot, 2), (-1.0 * G * M1 * M2 * pow(sl_kg, 2)) / (pow(r , 2.0)));
    //printf("vc: %E\nGravity: %E\nCentri: %E\nGasDF: %E\nStDF: %E\n", vc(r), M2 * sl_kg * g_fr(r,z), (M2 * sl_kg) * r * pow(thetadot, 2), gasDFr(r, rdot, thetadot, z, cs, vg * vc(r)), df_slow_r_disk(vg * vc(r), v, r, z, rdot));
    return (1.0 / (M2 * sl_kg)) * (M2 * sl_kg * g_fr(r,z) + (M2 * sl_kg) * r * pow(thetadot, 2) + 1.0 * gasDFr(r, rdot, thetadot, z, cs, vg * vc(r, z)) + 1.0 * df_slow_r_disk(vg * vc(r, z), v, r, z, rdot) + 1.0 * df_fast_r_disk(vg * vc(r, z), v, r, z, rdot) + 1.0 * df_bulge_r(r, z, v, rdot));

    //return (1.0 / (M2 * sl_kg)) * ((-1.0 * G * M1 * M2 * pow(sl_kg, 2)) / (pow(r , 2.0)) + (M2 * sl_kg) * r * pow(thetadot, 2) + gasDFr(r, r * thetadot, cs));

}

double ftheta(double r, double rdot, double theta, double thetadot, double z, double zdot, double cs) {//10/30 ADD ZDOT AS PARAM (UNRESOLVED)
    
    //double v = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0), 0.5);
    double v = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0) + pow(zdot, 2.0), 0.5);
    //printf("\n%E", (1.0 / (r * M2 * sl_kg)) * (gasDFt(r, r * thetadot, cs) + (-2.0) * rdot * thetadot * M2 * (sl_kg)));
    //printf("\n%E", gasDFt(r, v, cs));
    //printf("%E\n", df_slow_t_disk(vg * vc(r), v, r, z, thetadot) / gasDFt(r, rdot, thetadot, z, cs, vg * vc(r)));
    //printf("%E\n", r * thetadot / cs);
    //printf("%E\n", df_bulge_t(r, z, v, thetadot));
    return (1.0 / (r * M2 * sl_kg)) * ((-2.0) * rdot * thetadot * M2 * (sl_kg) + 1.0 * gasDFt(r, rdot, thetadot, z, cs, vg * vc(r, z)) + 1.0 * df_slow_t_disk(vg * vc(r, z), v, r, z, thetadot) + 1.0 * df_fast_t_disk(vg * vc(r, z), v, r, z, thetadot) + 1.0 * df_bulge_t(r, z, v, thetadot));
}

double fz(double r, double rdot, double thetadot, double z, double zdot, double cs) {
    //printf("%E\n", g_fz(r,z) / (-1.0 * sgn(zdot) * gasDFz(r, z, zdot, cs) / (M2 * sl_kg)));
    //printf("%f\n", zdot/cs);
    //printf("%E\n", (gasDFz(r, z, zdot, cs) / (M2 * sl_kg)) / g_fz(r, z));
    double v = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0) + pow(zdot, 2.0), 0.5);
    return g_fz(r, z) + -1.0 * sgn(zdot) * gasDFz(r, z, zdot, cs) / (M2 * sl_kg) + 1.0 * df_bulge_z(r, z, v, zdot) / (M2 * sl_kg);// + 0.0 * df_dm_z(r, z, v, zdot) / (M2 * sl_kg);
    //return 0.0;
}

// RK4 Solver function
State RK4Solver(State initial_state, double dt, double total_time) {
    State current_state = initial_state;
    double t = 0.0;
    char filename[20];
    char filename1[20];
    char filename2[20];
    char filename3[20];
    char filename_avgt[20];
    char filename_avgr[20];
    sprintf(filename, "bh_test2.txt");
    sprintf(filename1, "energy.txt");
    sprintf(filename2, "z.txt");
    sprintf(filename3, "tau.txt");
    sprintf(filename_avgr, "avgr.txt");
    FILE* fp = fopen(filename, "w");
    FILE* fm = fopen(filename1, "w");
    FILE* fv = fopen(filename2, "w");
    FILE* ft = fopen(filename3, "w");
    FILE* frad = fopen(filename_avgr, "w");
    
    double E0;
    double tau_z = 0.0;
    double t0_z = 0.0;
    bool peak = false;
    bool upz = current_state.zdot > 0.0;
    double dt0 = dt;

    double T = 1E4;
    double cs = pow((5 * kb * T) / (3 * m_p), 0.5);
    
    double mach = 0.0;
        
    while (t < total_time && pow(pow(current_state.r, 2.0) + pow(current_state.z, 2.0), 0.5) > 70.0 * pc_m && current_state.r < 1.0E4 * pc_m) {
        double r = current_state.r;
        //printf("%E\n", r);
        if (r < 0.0) {
            break;
        }
        double rdot = current_state.rdot;
        double theta = current_state.theta;
        double thetadot = current_state.thetadot;
        double z = current_state.z;
        double zdot = current_state.zdot;

        //printf("%E\n", r / pc_m);

        /* SPECIFIC ENERGY and ANGULAR MOMENTUM */
        double E;
        if (t == 0.0) {
            E0 = 0.5 * (rdot * rdot + r * r * thetadot * thetadot + zdot * zdot) + phi_rz(r, z);
            E = E0;
        } else {
            E = 0.5 * (rdot * rdot + r * r * thetadot * thetadot + zdot * zdot) + phi_rz(r, z);
        }

        double L_z = r * r * thetadot;
        double L_lat = r * zdot;
        double L_total = pow(L_z * L_z + L_lat * L_lat, 0.5);
        double v = pow(rdot * rdot + r * r * thetadot * thetadot + zdot * zdot, 0.5);

        //Introduce Toomre criterion
        cs = (G * M_PI * 2 * rho_rz(r, z) * z_thick) / ((2 * abs(vg) * vc(r, z)) / r);
        T = (3 * m_p * cs * cs) / (5 * kb);
        T += 1E4;
        cs = pow(5 * kb * T / (3 * m_p), 0.5);
        //cs = r * thetadot;
        mach = r * thetadot / cs;
        //double mach = zdot / cs;
        //printf("\n%E\t%E", mach, T);

        /*
        if (mach != mach) {
            break;
        }
        */

        //printf("%E\t%E\n", cs, T);

        //printf("\n%E", r / pc_m);
        //printf("\n%E", testPhi(r));
        
        //printf("\n%E", mach);
        
        //PERIODICITY TESTS FOR Z
        
        
        if (upz != zdot > 0.0) {
            upz = zdot > 0.0;
            if (peak) {
                tau_z = t - t0_z;
                t0_z = t;
                peak = false;
                //printf("%E\n", tau_z / (yr_s * 1E9));
                fprintf(ft, "%E\t%E\n", t, tau_z / (yr_s * 1E9));
                fprintf(frad, "%E\t%E\n", t, current_state.r);
            } else {
                peak = true;
            }
        }
        
        double k1r = rdot;
        double k1rd = fr(r, rdot, thetadot, z, zdot, cs);
        double k1t = thetadot;
        double k1td = ftheta(r, rdot, theta, thetadot, z, zdot, cs);
        double k1z = zdot;
        double k1zd = fz(r, rdot, thetadot, z, zdot, cs);
        
        //printf("\n%E", r);

        double k2r = rdot + 0.5 * k1rd * dt;
        double k2rd = fr(r + 0.5 * k1r * dt, rdot + 0.5 * k1rd * dt, thetadot + 0.5 * k1td * dt, z + 0.5 * k1z * dt, zdot + 0.5 * k1zd * dt, cs);
        double k2t = thetadot + 0.5 * k1td * dt;
        double k2td = ftheta(r + 0.5 * k1r * dt, rdot + 0.5 * k1rd * dt, theta + 0.5 * k1t * dt, thetadot + 0.5 * k1td * dt, z + 0.5 * k1z * dt, zdot + 0.5 * k1zd * dt, cs);
        double k2z = zdot + 0.5 * k1zd * dt;
        double k2zd = fz(r + 0.5 * k1r * dt, rdot + 0.5 * k1rd * dt, thetadot + 0.5 * k1td * dt, z + 0.5 * k1z * dt, zdot + 0.5 * k1zd * dt, cs);

        double k3r = rdot + 0.5 * k2rd * dt;
        double k3rd = fr(r + 0.5 * k2r * dt, rdot + 0.5 * k2rd * dt, thetadot + 0.5 * k2td * dt, z + 0.5 * k2z * dt, zdot + 0.5 * k2zd * dt, cs);
        double k3t = thetadot + 0.5 * k2td * dt;
        double k3td = ftheta(r + 0.5 * k2r * dt, rdot + 0.5 * k2rd * dt, theta + 0.5 * k2t * dt, thetadot + 0.5 * k2td * dt, z + 0.5 * k2z * dt, zdot + 0.5 * k2zd * dt, cs);
        double k3z = zdot + 0.5 * k2zd * dt;
        double k3zd = fz(r + 0.5 * k2r * dt, rdot + 0.5 * k2rd * dt, thetadot + 0.5 * k2td * dt, z + 0.5 * k2z * dt, zdot + 0.5 * k2zd * dt, cs);
        

        double k4r = rdot + k3rd * dt;
        double k4rd = fr(r + k3r * dt, rdot + k3rd * dt, thetadot + k3td * dt, z + k3z * dt, zdot + k3zd * dt, cs);
        double k4t = thetadot + k3td * dt;
        double k4td = ftheta(r + k3r * dt, rdot + k3rd * dt, theta + k3t * dt, thetadot + k3td * dt, z + k3z * dt, zdot + k3zd * dt, cs);
        double k4z = zdot + k3zd * dt;
        double k4zd = fz(r + k3r * dt, rdot + k3rd * dt, thetadot + k3td * dt, z + k3z * dt, zdot + k3zd * dt, cs);
        
        current_state.r += (k1r + 2.0 * k2r + 2.0 * k3r + k4r) * (dt / 6.0);
        current_state.rdot += (k1rd + 2.0 * k2rd + 2.0 * k3rd + k4rd) * (dt / 6.0);
        current_state.theta += (k1t + 2.0 * k2t + 2.0 * k3t + k4t) * (dt / 6.0);
        current_state.thetadot += (k1td + 2.0 * k2td + 2.0 * k3td + k4td) * (dt / 6.0);
        current_state.z += (k1z + 2.0 * k2z + 2.0 * k3z + k4z) * (dt / 6.0);
        current_state.zdot += (k1zd + 2.0 * k2zd + 2.0 * k3zd + k4zd) * (dt / 6.0);
        
        if (current_state.theta > 2.0 * PI) {
            current_state.theta -= 2.0 * PI;
        } else if (current_state.theta < 0.0) {
            current_state.theta += 2.0 * PI;
        }

        fprintf(fp, "%E\t%E\n", t, current_state.r);
        fprintf(fm, "%E\t%E\n", t, E);
        fprintf(fv, "%E\t%E\n", t, z);
        
        if (dt > (total_time / 1E8)) {
            dt = 0.001 * 2 * PI * current_state.r / abs(vc(r, z));
        } else {
            dt = dt0;
        }
        
        t += dt;
    }
    
    return current_state;
}

int main() {
    
    //Initialize galaxy
    int N = 200;
    int N2 = 200;
    r_arr = new double[N];
    z_arr = new double[N];
    r2_arr = new double[N];
    z2_arr = new double[N];

    double *a[N]; //rho
    double *a2[N2];
    double *b[N]; //phi
    double *b2[N2];
    double *c[N]; //rho_b
    double *c2[N2];
    double * rg[N];
    double * rg2[N2];
    double * rs[N];
    double * rs2[N2];
    double * rst[N];
    double * rst2[N2];
    double * rh[N];

    phi = new double*[N];
    phi2 = new double*[N2];
    rho = new double*[N];
    rho2 = new double*[N2];
    rho_s = new double*[N];
    rho_s2 = new double*[N2];
    rho_st = new double*[N];
    rho_st2 = new double*[N2];

    rho_g = new double*[N];
    rho_g2 = new double*[N2];

    rho_bg = new double*[N];
    rho_bg2 = new double*[N2];
    rho_h = new double*[N];
    phi_b = new double*[N];

    for (size_t i = 0; i < N; ++i) {
        phi[i] = new double[N];
        rho[i] = new double[N];
        rho_s[i] = new double[N];
        rho_st[i] = new double[N];
        rho_g[i] = new double[N];
        rho_bg[i] = new double[N];
        rho_h[i] = new double[N];
        phi_b[i] = new double[N];
    }
    for (size_t i = 0; i < N2; ++i) {
        phi2[i] = new double[N2];
        rho2[i] = new double[N2];
        rho_s2[i] = new double[N2];
        rho_st2[i] = new double[N2];
        rho_bg2[i] = new double[N2];
        rho_g2[i] = new double[N2];
    }
    for (int i = 0; i < N; i++) {
        r_arr[i] = 0.0;
        z_arr[i] = 0.0;
        for (int j = 0; j < N; j++) {
            phi[i][j] = -0.0 * 1E0;
            rho[i][j] = 0.0;
            rho_s[i][j] = 0.0;
            rho_st[i][j] = 0.0;
            rho_g[i][j] = 0.0;
            rho_h[i][j] = 0.0;
            phi_b[i][j] = 0.0;
            rho_bg[i][j] = 0.0;
        }
    }
    for (int i = 0; i < N2; i++) {
        r2_arr[i] = 0.0;
        z2_arr[i] = 0.0;
        for (int j = 0; j < N2; j++) {
            phi2[i][j] = -0.0 * 1E0;
            rho2[i][j] = 0.0;
            rho_s2[i][j] = 0.0;
            rho_st2[i][j] = 0.0;
            rho_g2[i][j] = 0.0;
            rho_bg2[i][j] = 0.0;
        }
    }

    for (int i = 0; i < N; i++) {
        a[i] = rho[i];
        rg[i] = rho_g[i];
        rs[i] = rho_s[i];
        rst[i] = rho_st[i];
        rh[i] = rho_h[i];
        b[i] = phi[i];
        c[i] = rho_bg[i];
    }
    for (int i = 0; i < N2; i++) {
        a2[i] = rho2[i];
        rg2[i] = rho_g2[i];
        rs2[i] = rho_s2[i];
        rst2[i] = rho_st2[i];
        b2[i] = phi2[i];
        c2[i] = rho_bg2[i];
    }
    /*
    for (int k = 0; k < N; k++) {
        b[k] = phi[k];
    }
    */
    double L = 2.0 * R_gd;
    double L2 = 0.1 * R_gd;
    double h = L / (N - 1);
    double h2 = L2 / (N2 - 1); //NOTE N not M

    //double rho0 = 10.0 * rho_gd(0.0); //FACTOR OF 10 TO ACCOUNT FOR FLATNESS, WILL BE REMOVED LATER WHEN THICK
    double rho0 = 1.0 * rho_gd(0.0);
    double rhob0 = 0.0 * 9.5E-21; //OFF
    double k = 1.0 / R_gd;

    double rho_sd = 0.25 * rho_gd(0.0); //to be calculated with integration later
    std::vector<double> params_std = {3.10448, R_sd, 0.1 * R_sd, R_sd / 3.0};
    std::vector<double> params_halo = {1.0 * 4.839E-20, 1.0 * 4.21E19};
    r_init(N, L, r_arr);
    z_init(N, L, z_arr);
    r_init(N, L, r2_arr);
    z_init(N2, L2, z2_arr);
    //rho_init(N, L, a, rho0, k);
    //printf("\nGotcha\n");
    //rhog_init(N, L, a, R_gd, z_thick, 0.0);
    rhos_init(N, L, rst, {0.0, R_sd, z_thick, 3.0 * z_thick});
    rhos_init(N, L, rs, params_std);
    rhog_init(N, L, rg, R_gd, z_thick, 1.0 * rho0);
    //halo_init(N, L, rh, params_halo);
    bulge_init(N, L, c, 1.0 * rhob0, 0.6, R_gd / 2.0);

    //std::vector<double> params_std2 = {3.10448, R_sd, 0.1 * R_sd, R_sd / 3.0};
    rhos_init2(N, N2, L, L2, rs2, params_std);
    rhog2_init(N, N2, L, L2, rg2, R_gd, z_thick, 1.0 * rho0);
    bulge2_init(N, N2, L, L2, c2, rhob0, 0.6, R_gd / 2.0);
    
    
    /*
    switching to flat conservative orbit:
    turn off DFs above,
    turn off bulge
    switch disk init
    */

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            rho_st[i][j] += (rho_bg[i][j] + rho_s[i][j]);
            rho[i][j] += (1.0 * rho_bg[i][j] + 1.0 * rho_g[i][j] + 1.0 * rho_s[i][j] + 0.0 * rho_h[i][j]); //HALO SWITCH
        }
    }
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N2; j++) {
            rho_st2[i][j] += (rho_bg2[i][j] + rho_s2[i][j]);
            rho2[i][j] += 1.0 * (1.0 * rho_bg2[i][j] + 1.0 * rho_g2[i][j] + 1.0 * rho_s2[i][j]);
            //printf("%E\n", rho2[i][j]);
        }
    }
    
    std::vector<double> dmass_params = {1.0 * rho_gd(0.0), R_gd, z_thick};
    double dmass = 2 * integrate2D(rho_gas_f, 0.0, 20.0 * R_sd, 0.0, 20.0 * z_thick, 200, 200, dmass_params);
    //Galaxy status
    printf("\nDMass: %E\n", dmass / (sl_kg));
    //printf("Energy: %E\t ")

/*
    char filename5[20];
    sprintf(filename5, "rho_bef.csv");
    FILE* fpot2 = fopen(filename5, "w");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(fpot2, "%E,", gsl_interp2d_eval(interp, r_arr, z_arr, rho_flat, (j - N / 2) * h, (i - N / 2) * h, xacc, yacc));
            //fprintf(fpot, "%E,", phi2[i][j]);
            //fprintf(fpot, "%E,", gsl_interp2d_eval(interp, z2_arr, r2_arr, rho2_flat, (j - N / 2) * dzh, (i - N2 / 2) * drh, xacc, yacc));
        }
        fprintf(fpot2, "\n");
    }
*/    
    relax(N, L, a, b, 50000, 4 * M_PI * G * rho0 * 1E-4 * h * h);
    
    // Create an interpolation object
    interp = gsl_interp2d_alloc(gsl_interp2d_bicubic, N, N);
    xacc = gsl_interp_accel_alloc();
    yacc = gsl_interp_accel_alloc();

    
    interp_r = gsl_interp2d_alloc(gsl_interp2d_bicubic, N, N);
    xacc_r = gsl_interp_accel_alloc();
    yacc_r = gsl_interp_accel_alloc();
    
    interp_amr = gsl_interp2d_alloc(gsl_interp2d_bicubic, N, N);
    xacc_amr = gsl_interp_accel_alloc();
    yacc_amr = gsl_interp_accel_alloc();
    
    double drh = L / (N + 1);
    double dzh = L2 / (N2 + 1);

    phi_flat = new double[N * N];
    rho_flat = new double[N * N];
    rho_g_flat = new double[N * N];
    rho_s_flat = new double[N * N];
    rho_st_flat = new double[N * N];
    rho_h_flat = new double[N * N];

    int i1 = 0;
    int i2 = 0;
    int i3 = 0;
    int i4 = 0;
    int i5 = 0;
    int i6 = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            rho_flat[i1++] = rho[i][j];
            phi_flat[i2++] = phi[i][j];
            rho_g_flat[i3++] = rho_g[i][j];
            rho_s_flat[i4++] = rho_s[i][j];
            rho_st_flat[i5++] = rho_st[i][j];
            rho_h_flat[i6++] = rho_h[i][j];
        }
    }

    // Initialize the interpolation object
    
    gsl_interp2d_init(interp, r_arr, z_arr, phi_flat, N, N);
    gsl_interp2d_init(interp_r, r_arr, z_arr, rho_flat, N, N);
    
    /*
    for (int j = 0; j < N/2; j++) {
        //printf("%E\n", phi[N/2 - 20][j]);
        //printf("%E\n", z_arr[j]);
        printf("%E\t%E\n", gsl_interp2d_eval_deriv_y(interp, r_arr, z_arr, phi_flat, 1000.0 * pc_m, -1.0 * j * h, xacc, yacc), gsl_interp2d_eval_deriv_y(interp, r_arr, z_arr, phi_flat, 1000.0 * pc_m, 1.0 * j * h, xacc, yacc));
    }
    */
    
    //sus ordering???
    for (int i = 0; i < N2; i++) {
        for (int j = 0; j < N; j++) {
            //if (i == 0.0 || j == 0.0 || i == N2 - 1 || j == N - 1) {
                phi2[i][j] = gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, (j - N / 2) * drh, (i - N / 2) * dzh, xacc_amr, yacc_amr);
            //}
            //CHECK LEFT AND RIGHT EDGES, THEY NEED TO BE 0
        }
    }
    
    char filenameb[20];
    sprintf(filenameb, "rho.csv");
    FILE* fpotb = fopen(filenameb, "w");
    for (int i = 0; i < N2; i++) {
        for (int j = 0; j < N; j++) {
            //fprintf(fpot, "%E,", gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, (j - N / 2) * h, (i - N / 2) * h, xacc, yacc));
            fprintf(fpotb, "%E,", rho2[i][j]);
            //printf("%E\n", rho2[i][j]);
            //fprintf(fpot2, "%E,", gsl_interp2d_eval(interp, r2_arr, z2_arr, phi_flat2, (i - N / 2) * drh, (j - N2 / 2) * dzh, xacc, yacc));
        }
        fprintf(fpotb, "\n");
    }
    
    relax_amr(N, N2, L, L2, a2, b2, 50000, 4 * M_PI * G * rho0 * 1E-4 * h * h2);

    phi_flat2 = new double[N * N2];
    rho_flat2 = new double[N * N2];
    rho_g_flat2 = new double[N * N2];
    rho_s_flat2 = new double[N * N2];
    rho_st_flat2 = new double[N * N2];


    i1 = 0;
    i2 = 0;
    i3 = 0;
    i4 = 0;
    i5 = 0;
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N2; j++) {
            rho_flat2[i1++] = rho2[i][j];
            phi_flat2[i2++] = phi2[i][j];
            rho_g_flat2[i3++] = rho_g2[i][j];
            rho_s_flat2[i4++] = rho_s2[i][j];
            rho_st_flat2[i5++] = rho_st2[i][j];
            //printf("%E\n", rho_g2[i][j]);
        }
    }  


    gsl_interp2d_init(interp_amr, r2_arr, z2_arr, phi_flat2, N, N);
    
    char filenamebc[20];
    sprintf(filenamebc, "inner.csv");
    FILE* fpot2 = fopen(filenamebc, "w");
    for (int i = 0; i < N2; i++) {
        for (int j = 0; j < N; j++) {
            //fprintf(fpot, "%E,", gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, (j - N / 2) * h, (i - N / 2) * h, xacc, yacc));
            fprintf(fpot2, "%E,", phi2[i][j]);
            //printf("%E\n", rho_g2[i][j]);
            //fprintf(fpot2, "%E,", gsl_interp2d_eval(interp, r2_arr, z2_arr, phi_flat2, (i - N / 2) * drh, (j - N2 / 2) * dzh, xacc, yacc));
        }
        fprintf(fpot2, "\n");
    }
    
    /*
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0 || j == 0 || i == N-1 || j == N-1) {
                rho[i][j] = 0.0;
            }
        }
        //fprintf(fpot2, "\n");
    }
    */
    char filename4[20];
    sprintf(filename4, "phi.csv");
    FILE* fpot = fopen(filename4, "w");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            //fprintf(fpot, "%E,", gsl_interp2d_eval(interp, r_arr, z_arr, rho_flat, (j - N / 2) * h, (i - N / 2) * h, xacc, yacc));
            fprintf(fpot, "%E,", phi[i][j]);
            //fprintf(fpot, "%E,", phi2[i][j]);
            //fprintf(fpot, "%E,", gsl_interp2d_eval(interp, z2_arr, r2_arr, rho2_flat, (j - N / 2) * dzh, (i - N2 / 2) * drh, xacc, yacc));
        }
        fprintf(fpot, "\n");
    }
    
    
    //printf("%E\n", gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, 1.0 * pc_m, -1.0 * pc_m, xacc, yacc));
    double cs0 = pow(5 * kb * 1E4 / (3 * m_p), 0.5);
    double r0 = 1.00 * 1000.0 * pc_m;
    double rd0 = 0.0;
    double t0 = 0.0;
    double z0 = 0.0 * pc_m; //standard of 0.5 pc
    //double z0 = 1.0 * z_thick;
    double zd0 = 0.0 * cs0;
    double f = 0.5;
    double v_c = vc(r0, z0);
    double td0 = 1.0 * v_c * pow(1 + f, 0.5) / r0;

    
    double v_0 = r0 * td0;
    double s_0 = pow(r0 * r0 + z0 * z0, 0.5);
    double theta_z = (30.0 / 180) * M_PI;
    r0 = s_0 * cos(theta_z);
    z0 = s_0 * sin(theta_z);
    //z0 = -0.0 * pc_m;
    td0 = 1.0 * v_0 / r0;

    printf("%E\t%E\n", r0 / (1000.0 * pc_m), z0 / (1000.0 * pc_m));
    double total_time = 1.2E10 * yr_s; //BE CAREFUL HERE
    double n = 1000000.0;
    double time_step = total_time / n;
    //printf("%lf\n", numInt(0.0, 3.0, testPhi));
    
    State initial_state{ r0, rd0, t0, td0, z0, zd0};
    State final_state = RK4Solver(initial_state, time_step, total_time);

    printf("\nTerminated at radius %E\n", final_state.r / pc_m);
    return 0;
}
