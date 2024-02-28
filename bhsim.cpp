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
double vg = 0.5;

double R_sd = log10(M1 / pow(10,5)) * 1000 * pc_m;
double R_gd = 2.0 * R_sd;

double rho_bulge = 1E8;
double M_s = 2.0E30;
double M_e = 5.97E24;
double n_gd = 300 * pow(100, 3);
double z_thick = R_gd / 10.0;
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
    double dz = z * 1E-8;
    return (f(z + dz) - f(z)) / dz;
}

/* Relevant potentials/forces (disk, bulge, stellar DF, gas DF, velocity distribution, etc. */


double rho_gd(double r) {
    return n_gd * m_p * pow(e, -1.0 * r / (R_gd));
}

double sigma_gd(double r) {
    return rho_gd(r) * (2 * z_thick);
}


double disk_pot(double r, double z) {
    int n1 = 20;
    int n2 = 20;

    double dk = 7.0 / n1;
    double dr = 2.0 * R_gd / n2;
    double cumsum = 0.0;
    
    double k = 0;

    for (int i = 0; i < n1; i++) {
        double rp = 0;
        for (int j = 0; j < n2; j++) {
            cumsum += exp(-1.0 * k * abs(z)) * std::cyl_bessel_j(0.0, r) * rp * std::cyl_bessel_j(0.0, k * rp) * sigma_gd(rp) * dr * dk;
            rp += dr;
        }
        k += dk;
    }

    return -2.0 * PI * G * cumsum;
}


double gasDFb(double r, double tdot) {
    return 4 * PI * pow((G * M2 * sl_kg), 2) * rho_gd(r); //DONT DIVIDE BY V^2 YET
}

double gasDFr(double r, double rdot, double tdot, double cs, double vd) {
    
    double rt = r * tdot;
    double v = pow(rdot * rdot + (rt - vd) * (rt - vd), 0.5);
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

    return (-1.0 * gasDFb(r, v) * I_r) / (v * v) ;
}

double gasDFt(double r, double rdot, double tdot, double cs, double vd) {

    
    double rt = r * tdot;
    double v = pow(rdot * rdot + (rt - vd) * (rt - vd), 0.5);
    double mach = v / cs;
    if (mach < 1E-4) {
        //printf("zero vel");
        return 0;
    }
    
    double I_t = 0.0;
    double r_min = r / 10.0;
    if (mach < 1.0 && mach > 0.0) {
        I_t = 0.7706 * log((1 + mach) / (1.0004 - 0.9185 * mach));
        //printf("%E", I_t);
    } else if (mach < 4.4) {
        I_t = log(330 * 10.0 * pow(mach, -9.58) * pow(mach - 0.71, 5.72));
    } else if (mach > 4.4) {
        I_t = log(10 / (0.11 * mach + 1.65));
    }
    //printf("\n%E", mach);
    return (-1.0 * gasDFb(r, v) * I_t) / (v * v);
}

double oomDF(double r, double v, double cs) {
    return 4 * PI * G * G * M2 * M2 * n_gd * m_p / (cs * cs);
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
double testPhi(double r) { //test potential
    double a = 1.0 / R_gd;
    //printf("\n%E", a * r / 2.0);
    return PI * G * sigma_gd(0) * r * (std::cyl_bessel_i(1.0, a * r / 2.0) * std::cyl_bessel_k(0.0, a * r / 2) - std::cyl_bessel_i(0.0, a * r / 2.0) * std::cyl_bessel_k(1.0, a * r / 2.0));
}

double gF(double r) { //SPECIFIC FORCE
    double dr = r * 1E-8;
    return -1.0 * (testPhi(r + dr) - testPhi(r)) / dr;
    //return -1.0 * G * M1 * M2 * sl_kg * sl_kg / (r * r);
}

double g_fr(double r, double z) {
    double dr = r * 1E-8;
    return -1.0 * (disk_pot(r + dr, z) - disk_pot(r, z)) / dr;
}

double g_fz(double r, double z) {
    double dz = z * 1E-8;
    return -1.0 * (disk_pot(r, z + dz) - disk_pot(r, z)) / dz;
}
double vc(double r) { 
    return pow(r * -1.0 * gF(r), 0.5);
}

//total force governing EoM
double fr(double r, double rdot, double thetadot, double cs) {
    //include dissipative forces later
    
    //double df1 = /* (4 * PI * pow(6.67E-11, 2) * pow(M2 * (2E+30), 2) * rho_bulge * (2E+30 / (pow(3.086E+19, 3))) */ (1.0 / pow(abs(v), 2)) * stelDF_func(abs(v));
    double v = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0), 0.5);
    //printf("\n%E", (1.0 / (M2 * sl_kg)) * ((-1.0 * G * M1 * M2 * pow(sl_kg, 2)) / pow(r , 2.0) + (M2 * sl_kg) * r * pow(thetadot, 2) + gasDFr(r, v, cs)));
    //printf("\n%E\t%E\n", gasDFr(r, v, cs), (-1.0 * G * M1 * M2 * pow(sl_kg, 2)) / (pow(r , 2.0) + 1.0E-10));
    //printf("\n%E\t%E", (M2 * sl_kg) * r * pow(thetadot, 2), (-1.0 * G * M1 * M2 * pow(sl_kg, 2)) / (pow(r , 2.0)));
    return (1.0 / (M2 * sl_kg)) * (M2 * sl_kg * gF(r) + (M2 * sl_kg) * r * pow(thetadot, 2) + gasDFr(r, rdot, thetadot, cs, vg * vc(r)));

    //return (1.0 / (M2 * sl_kg)) * ((-1.0 * G * M1 * M2 * pow(sl_kg, 2)) / (pow(r , 2.0)) + (M2 * sl_kg) * r * pow(thetadot, 2) + gasDFr(r, r * thetadot, cs));

}

double ftheta(double r, double rdot, double theta, double thetadot, double cs) {
    
    double v = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0), 0.5);
    //printf("\n%E", (1.0 / (r * M2 * sl_kg)) * (gasDFt(r, r * thetadot, cs) + (-2.0) * rdot * thetadot * M2 * (sl_kg)));
    //printf("\n%E", gasDFt(r, v, cs));
    return (1.0 / (r * M2 * sl_kg)) * ((-2.0) * rdot * thetadot * M2 * (sl_kg) + gasDFt(r, rdot, thetadot, cs, vg * vc(r))); //coriolis term to conserve angular momentum
}

double fz(double r, double rdot, double theta, double thetadot, double z, double zdot, double cs) {
    return g_fz(r, z);
}

// RK4 Solver function
State RK4Solver(State initial_state, double dt, double total_time) {
    State current_state = initial_state;
    double t = 0.0;
    char filename[20];
    char filename1[20];
    sprintf(filename, "bh_test2.txt");
    sprintf(filename1, "mach.txt");
    FILE* fp = fopen(filename, "w");
    FILE* fm = fopen(filename1, "w");
    while (t < total_time && current_state.r > 10.0 * pc_m && current_state.r < 1.0E4 * pc_m) {
        double r = current_state.r;

        if (r < 0.0) {
            break;
        }
        double rdot = current_state.rdot;
        double theta = current_state.theta;
        double thetadot = current_state.thetadot;
        //printf("\n%E\n", r);
        double T = 1E4;
        double cs = pow((5 * kb * T) / (3 * m_p), 0.5);
        
        double v = pow(rdot * rdot + r * r * thetadot * thetadot, 0.5);
        //printf("\n%E", r);
        double mach = r * thetadot / cs;
        if (mach != mach) {
            break;
        }
        printf("\n%E", r / pc_m);
        if (t < dt * 50) {
            //printf("\n%E", mach);
        }
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
        
        //printf("\n%E", r);

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

        fprintf(fp, "%E\t%E\n", t, current_state.r);
        fprintf(fm, "%E\t%E\n", t, mach);
        //dt = 0.01 * 2 * PI * current_state.r / vc(r);
        t += dt;
    }
    
    return current_state;
}

int main() {
    //double initial_position = 0.0;
    //double initial_velocity = 5.0;
    
    double r0 = 1000.0 * pc_m;
    double rd0 = 0.0;
    double t0 = 0.0;

    double f = 0.5;
    double v_c = vc(r0);
    double td0 = v_c * pow(1 + f, 0.5) / r0;
    //printf("%E", vc);
    double total_time = 1.0E10 * yr_s;
    double n = 1000000.0;
    double time_step = total_time / n;
    //printf("%lf\n", numInt(0.0, 3.0, testPhi));
    double l = pow(r0, 2) * td0;
    
    State initial_state{ r0, rd0, t0, td0};
    State final_state = RK4Solver(initial_state, time_step, total_time);

    printf("\ndone\n");
    return 0;
}
