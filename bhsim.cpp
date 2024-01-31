#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <cmath> //need this for bessel and few other things

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

double R_sd = log10(M1 / pow(10,5));
double R_gd = 2.0 * R_sd;

double rho_bulge = 1E8;
double M_s = 2.0E30;
double M_e = 5.97E24;
double n_gd = 300 * pow(0.01, 3);
//double M2 = 1E6;

// Define the state of the system
struct State { //define phase space coord
    double r;
    double rdot;
    double theta;
    double thetadot;
};

//Solver helper functions
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
    double dz = 1E-8;
    return (f(z + dz) - f(z)) / dz;
}

/* Relevant potentials/forces (disk, bulge, stellar DF, gas DF, velocity distribution, etc. */


double rho_gd(double r) {
    return n_gd * m_p * pow(e, -1.0 * r / R_gd);
}

double gasDFb(double r, double v) {
    return 4 * PI * pow((G * M2 * sl_kg), 2) * rho_gd(r) / pow(v, 2);
}

double gasDFr(double r, double v, double cs) {
    double mach = v / cs;
    double I_r;

    if (mach < 1.1) {
        I_r = pow(mach, 2) * pow(10, 3.51 * mach - 4.22);
    } else if (mach < 4.4) {
        I_r = 0.5 * log(9.33 * pow(mach, 2) * (pow(mach, 2) - 0.95));
    } else {
        I_r = 0.3 * pow(mach, 2);
    }

    return -1.0 * gasDFb(r, v) * I_r;
}

double gasDFt(double r, double v, double cs) {
    double mach = v / cs;
    double I_t;
    double r_min = r / 10.0;
    if (mach < 1.0) {
        I_t = 0.7706 * log((1 + mach) / (1.0004 - 0.9185 * mach));
    } else if (mach < 4.4) {
        I_t = log(330 * (r / r_min) * pow(mach, -9.58) * pow(mach - 0.71, 5.72));
    } else {
        I_t = log((r / r_min) / (0.11 * mach + 1.65));
    }

    return -1.0 * gasDFb(r, v) * I_t;
}

double stelDF(double s, double velocity_c, double velocity, double p_max, double M2, double vg)
{
    if (abs(velocity) < pow(2.0, 0.5) * velocity_c) {
        /* g is the field star velocity distribution*/
        double g = (tgamma(gama + 1) / tgamma(gama - 0.5)) * pow((2 * pow(velocity_c, 2) - pow(s, 2)), (gama - 1.5)) / (pow(2, gama) * pow(PI, 1.5) * pow(velocity_c, (2 * gama)));
        /*double g = exp(-1*pow(s, 2) / pow(2E+5*pow((M1/(1.9E+8)),(1.0/5.1)),2)) / ( pow(PI, 1.5)* pow(2E+5*pow((M1/(1.9E+8)),(1.0/5.1)),3) ) ;*/
        double f = 4 * PI * g * pow(s, 2) * log((p_max / ((6.67E-11) * M2 * (2E+30))) * (pow(velocity, 2) - pow(s, 2)));
        return f;

    }
    else {
        double f = 0.0;
        return f;
    }
}

double sumintegral_10(double velocity_c, double velocity, double p_max, double M2, double vg)
{
    int n = 20;
    double lowbound = 0.0;
    double upbound = abs(velocity);
    double dx = (double)(upbound - lowbound) / n;
    double cumsum = 0;
    for (int i = 1; i < n; i++)
    {
        double xi = lowbound + i * dx;
        double function_value = stelDF(xi, velocity_c, abs(velocity), p_max, M2, vg);
        double rectangle_area = function_value * dx;
        cumsum += rectangle_area;
    }
    return cumsum;
}

double stelDF_func(double v) {
    return sumintegral_10(0.1, v, 1.0, 1.0, 0.5); //dummy vars for now
}
//total potential (no dissipative forces yet)
double testPhi(double z) { //test potential
    double mu = 1.0;//pow(2 * PI, 2);
    double k = 1.0;
    //return -1 * mu / (k + abs(z));
    return k * z * z;
}

//total force governing EoM
double fr(double r, double rdot, double thetadot, double cs) {
    //include dissipative forces later
    
    //double df1 = /* (4 * PI * pow(6.67E-11, 2) * pow(M2 * (2E+30), 2) * rho_bulge * (2E+30 / (pow(3.086E+19, 3))) */ (1.0 / pow(abs(v), 2)) * stelDF_func(abs(v));
    double v = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0), 0.5);
    return (1.0 / (M2 * sl_kg)) * ((-1.0 * G * M1 * M2 * pow(sl_kg, 2)) / pow(r , 2.0) + r * pow(thetadot, 2) + gasDFr(r, v, cs));
}

double ftheta(double r, double rdot, double theta, double thetadot, double cs) {
    
    double v = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0), 0.5);
    return (1.0 / (M2 * sl_kg)) * (gasDFt(r, v, cs) + (-2.0) * rdot * thetadot * M2 / (r)); //coriolis force to conserve angular momentum
}
//SUBSTITUTE L
// RK4 Solver function
State RK4Solver(State initial_state, double dt, double total_time) {
    State current_state = initial_state;
    double t = 0.0;
    char filename[20];
    sprintf(filename, "bh_test.txt");
    FILE* fp = fopen(filename, "w");

    while (t < total_time) {
        double r = current_state.r;
        double rdot = current_state.rdot;
        double theta = current_state.theta;
        double thetadot = current_state.thetadot;

        double T = 1E4;
        double cs = pow((5 * kb * T) / (3 * m_p), 0.5);
        // Calculate the four RK4 steps
        /*  
        double k1x = v;
        double k1v = customForce(x, v);

        double k2x = v + 0.5 * k1v * dt;
        double k2v = customForce(x + 0.5 * k1x * dt, v + 0.5 * k1v * dt);

        double k3x = v + 0.5 * k2v * dt;
        double k3v = customForce(x + 0.5 * k2x * dt, v + 0.5 * k2v * dt);

        double k4x = v + k3v * dt;
        double k4v = customForce(x + k3x * dt, v + k3v * dt);

        current_state.x += (k1x + 2.0 * k2x + 2.0 * k3x + k4x) * (dt / 6.0);
        current_state.v += (k1v + 2.0 * k2v + 2.0 * k3v + k4v) * (dt / 6.0);
        */

        double k1r = rdot;
        double k1rd = fr(r, rdot, thetadot, cs);
        double k1t = thetadot;
        double k1td = ftheta(r, rdot, theta, thetadot, cs);
        
        double k2r = rdot + 0.5 * k1rd * dt;
        double k2rd = fr(r + 0.5 * k1r * dt, rdot + 0.5 * k1rd * dt, thetadot + 0.5 * k1td * dt, cs);
        double k2t = thetadot + 0.5 * k1td * dt;
        double k2td = ftheta(r + 0.5 * k1r * dt, rdot + 0.5 * k1rd * dt, theta + 0.5 * k1t * dt, thetadot + 0.5 * k1td * dt, cs);
        
        double k3r = rdot + 0.5 * k2rd * dt;
        double k3rd = fr(r + 0.5 * k2r * dt, rdot + 0.5 * k2rd * dt, thetadot + 0.5 * k2td * dt, cs);
        double k3t = thetadot + 0.5 * k2td * dt;
        double k3td = ftheta(r + 0.5 * k2r * dt, rdot + 0.5 * k2rd * dt, theta + 0.5 * k2t * dt, thetadot + 0.5 * k2td * dt, cs);
        
        double k4r = rdot + k3rd * dt;
        double k4rd = fr(r + k3r * dt, rdot + k3rd * dt, thetadot + k3td * dt, cs);
        double k4t = thetadot + k3td * dt;
        double k4td = ftheta(r + k3r * dt, rdot + k3rd * dt, theta + k3t * dt, thetadot + k3td * dt, cs);
        
        current_state.r += (k1r + 2.0 * k2r + 2.0 * k3r + k4r) * (dt / 6.0);
        current_state.rdot += (k1rd + 2.0 * k2rd + 2.0 * k3rd + k4rd) * (dt / 6.0);
        current_state.theta += (k1t + 2.0 * k2t + 2.0 * k3t + k4t) * (dt / 6.0);
        current_state.thetadot += (k1td + 2.0 * k2td + 2.0 * k3td + k4td) * (dt / 6.0);
        
        if (current_state.theta > 2.0 * PI) {
            current_state.theta -= 2.0 * PI;
        } else if (current_state.theta < 0.0) {
            current_state.theta += 2.0 * PI;
        }

        fprintf(fp, "%E\t%f\n", t, current_state.r);
        t += dt;
    }
    
    return current_state;
}

int main() {
    //double initial_position = 0.0;
    //double initial_velocity = 5.0;

    double r0 = 3000 * pc_m;
    double rd0 = 0.0;
    double t0 = 0.0;
    double td0 = 1E-18;
    double total_time = 1.0E9 * yr_s;
    int n = 10000;
    double time_step = total_time / n;
    //printf("%lf\n", numInt(0.0, 3.0, testPhi));
    double l = pow(r0, 2) * td0;

    State initial_state{ r0, rd0, t0, td0};
    State final_state = RK4Solver(initial_state, time_step, total_time);

    return 0;
}
