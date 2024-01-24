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

//physics constants
double G = 6.67E-11;
double gama = 1.8;/*>0.6 */
double gama_e = 1.8;
double alpha = 4.0;
double alpha_b = 1.8;
double q_b = 0.6;

//parameters of sim
double rho_bulge = 1E8;
double M_s = 2.0E30;
double M_e = 5.97E24;

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
double fr(double r, double rdot, double tdot, double l) {
    //include dissipative forces later
    
    //double df1 = /* (4 * PI * pow(6.67E-11, 2) * pow(M2 * (2E+30), 2) * rho_bulge * (2E+30 / (pow(3.086E+19, 3))) */ (1.0 / pow(abs(v), 2)) * stelDF_func(abs(v));
    
    return (-1.0 * G * M_s) / pow(r , 2.0) + pow(l, 2) / pow(r, 3);
}

double ftheta(double theta, double thetadot) {
    
    return 0;
}
//SUBSTITUTE L DUMBASS
// RK4 Solver function
State RK4Solver(State initial_state, double dt, double total_time, double l) {
    State current_state = initial_state;
    double t = 0.0;
    char filename[20];
    sprintf(filename, "osc_sun.txt");
    FILE* fp = fopen(filename, "w");

    while (t < total_time) {
        double r = current_state.r;
        double rdot = current_state.rdot;
        double theta = current_state.theta;
        double thetadot = current_state.thetadot;

        
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
        double k1rd = fr(r, rdot, thetadot, l);
        double k1t = thetadot;
        double k1td = ftheta(theta, thetadot);
        
        double k2r = rdot + 0.5 * k1rd * dt;
        double k2rd = fr(r + 0.5 * k1r * dt, rdot + 0.5 * k1rd * dt, thetadot + 0.5 * k1td * dt, l);
        double k2t = thetadot + 0.5 * k1td * dt;
        double k2td = ftheta(theta + 0.5 * k1t * dt, thetadot + 0.5 * k1td * dt);
        
        double k3r = rdot + 0.5 * k2rd * dt;
        double k3rd = fr(r + 0.5 * k2r * dt, rdot + 0.5 * k2rd * dt, thetadot + 0.5 * k2td * dt, l);
        double k3t = thetadot + 0.5 * k2td * dt;
        double k3td = ftheta(theta + 0.5 * k2t * dt, thetadot + 0.5 * k2td * dt);
        
        double k4r = rdot + k3rd * dt;
        double k4rd = fr(r + k3r * dt, rdot + k3rd * dt, thetadot + k3td * dt, l);
        double k4t = thetadot + k3td * dt;
        double k4td = ftheta(theta + k3t * dt, thetadot + k3td * dt);
        
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

    double r0 = 1.47E11;
    double rd0 = 0.0;
    double t0 = 0.0;
    double td0 = 2.06E-7;
    double total_time = 3.2E7;
    int n = 100000;
    double time_step = total_time / n;
    //printf("%lf\n", numInt(0.0, 3.0, testPhi));
    double l = pow(r0, 2) * td0;

    State initial_state{ r0, rd0, t0, td0};
    State final_state = RK4Solver(initial_state, time_step, total_time, l);

    return 0;
}
