/*
===============================================================================
  bhsim_2.0.cpp  --  Massive black hole binary (MBHB) inspiral simulation
===============================================================================

  WHAT THIS CODE DOES
  -------------------
  Simulates the orbital decay of a secondary supermassive black hole (SMBH,
  mass M2) as it spirals into a primary SMBH (mass M1) sitting at the center
  of a gas-rich host galaxy. The secondary is treated as a test particle whose
  orbit is integrated with a 4th-order Runge-Kutta (RK4) scheme in cylindrical
  coordinates (r, theta, z).

  The secondary feels four sources of acceleration:
    (1) Gravity from the galactic potential (gas disk + stellar disk + bulge),
        computed from a Poisson solve on a 2D axisymmetric grid.
    (2) Gravity from the primary SMBH (point-mass Keplerian term, added in
        analytically since a point mass cannot be resolved on the grid).
    (3) Gas dynamical friction (Ostriker 1999 + Kim & Kim 2007 fits, with
        Mach-number-dependent Coulomb logarithms).
    (4) Stellar dynamical friction, split into:
          - Disk component: anisotropic, "slow" vs "fast" regimes relative to
            the disk rotation speed (Antonini & Merritt 2012-style).
          - Bulge component: isotropic, integrated over an assumed stellar
            velocity distribution (see fs_v() in poisson.h).

  Simultaneously, the code tracks Bondi-Hoyle-Lyttleton (BHL) accretion onto
  both black holes, updates their masses in real time, and estimates jet
  power using a Blandford-Znajek-like prescription.

  The potential grid uses adaptive mesh refinement (AMR): a coarse outer grid
  covers the whole galaxy, and a finer "thin" grid resolves the disk midplane
  where vertical gradients are steep.

  -------------------------------------------------------------------------
  HOW TO RUN
  -------------------------------------------------------------------------
  Compile (needs GSL):
      g++ bhsim_2_0.cpp -o bhsim -lgsl -lgslcblas -lm

  Execute with one argument -- the initial orbital inclination in DEGREES:
      ./bhsim 30          # 30-degree inclined orbit
      ./bhsim 0           # coplanar (in-disk) orbit
      ./bhsim 90          # polar orbit

  The code expects a pre-computed potential file `pot_test.csv` (outer grid)
  and `pot_test2.csv` (AMR grid) in the working directory. To regenerate
  these from scratch, set `loaded = false` in main() -- this triggers the
  Poisson relaxation solver in poisson.h (slow, ~minutes).

  OUTPUTS (all written to the working directory):
      bh_test2.txt   (t, r)               main orbit trace
      orbit.txt      (t, r, theta, z)     full 3D position
      energy.txt     (t, E)               specific orbital energy
      z.txt          (t, z)               vertical position vs time
      mass.txt       (t, M2)              secondary mass growth
      mass1.txt      (t, M1)              primary mass growth
      mdot.txt       (t, Mdot2)           secondary accretion rate
      pbz.txt        (t, P_jet2)          secondary jet power
      pot_test.csv                        2D potential on outer grid
      ... and many more diagnostics (pr.txt, pt.txt, pz.txt, etc.)

  -------------------------------------------------------------------------
  SCENARIOS YOU CAN RUN (set in "SIMULATION PARAMETERS" block below)
  -------------------------------------------------------------------------
    * Primary mass:        change M1   (e.g. 1e7, 1e8, 1e9 solar masses)
    * Mass ratio:          change q    (M2 = q * M1; default 1/9)
    * Gas density:         change n_gd (central gas number density)
    * Gas fraction:        change fg   (gas vs stars in the disk)
    * Prograde/retrograde: change vg sign  (negative = prograde; DO NOT set 0)
    * Initial inclination: command-line argument in degrees
    * Halo on/off:         edit the "HALO SWITCH" line in main() (~line 1509)
    * Bulge/disk DF:       multiply the relevant df_*_disk or df_bulge_*
                           terms by 0 in fr/ftheta/fz to disable them
    * Load vs recompute potential:  toggle `loaded` in main()

  Units: SI throughout (meters, kilograms, seconds). Conversion factors
  (yr_s, sl_kg, pc_m) are defined below.

===============================================================================
*/

#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <algorithm>
#include <cmath>               // special math functions (bessel, etc.)
#include "poisson.h"           // galaxy density profiles + Poisson relaxation solver
                               // (also pulls in integration.h for 1D/2D Simpson's rule)
#include <gsl/gsl_interp2d.h>  // 2D bilinear interpolation of potential & density grids
#include <gsl/gsl_errno.h>     // GSL error handler (useful for debugging out-of-bounds)

/*
 GSL error handler. GSL aborts by default on errors, which is usually a
 nuisance during development. Install this with gsl_set_error_handler() to
 get a message instead. The most common error you'll see is the secondary
 flying off the edge of the interpolation grid.
*/
void my_gsl_error_handler(const char *reason, const char *file, int line, int gsl_errno) {
    fprintf(stderr, "GSL ERROR: %s\nIn file: %s at line %d\nError code: %d\n", reason, file, line, gsl_errno);
}


// ---------------------------------------------------------------------------
// Math constants (redundant with cmath's M_PI/M_E, but used throughout -- do
// not remove unless you plan to sweep the whole file).
// ---------------------------------------------------------------------------
double PI = 3.14159265;
double e = 2.718281828;


// ---------------------------------------------------------------------------
// Unit conversions. The code is in SI; these are used only when converting
// human-readable inputs (solar masses, parsecs) or when printing results.
// ---------------------------------------------------------------------------
double yr_s = 3.154E7;   // years -> seconds
double sl_kg = 1.99E30;  // solar masses -> kg
double pc_m = 3.086E16;  // parsecs -> meters

// ---------------------------------------------------------------------------
// Physical constants (SI)
// ---------------------------------------------------------------------------
double G    = 6.67E-11;   // gravitational constant        [m^3 kg^-1 s^-2]
double kb   = 1.38E-23;   // Boltzmann constant            [J/K]
double m_p  = 1.67E-27;   // proton mass                   [kg]
double c    = 3.0E8;      // speed of light                [m/s]
double mu_0 = 1.257E-6;   // vacuum magnetic permeability  [H/m]


// ===========================================================================
// SIMULATION PARAMETERS -- these are the main knobs you tune between runs.
// ===========================================================================
double M1    = 1.0E7;            // Primary BH mass [SOLAR MASSES, not kg]
double q     = 1.0 / 9.0;        // Mass ratio M2/M1 (1:9 in the first paper)
double M2    = M1 * q;           // Secondary BH mass [solar masses]
double n_gd  = 100 * pow(100,3); // Central gas number density [m^-3]
                                 // (multiply by m_p for mass density)
double fg    = 0.7;              // Gas fraction in the disk (rest is stars)
double vg    = -0.7;             // Disk rotation (fraction of v_c):
                                 //   vg < 0  -> PROGRADE orbit
                                 //   vg > 0  -> RETROGRADE orbit
                                 //   NEVER set vg = 0 (division by zero downstream)

// Disk scale radii, chosen to scale with M1 (rough empirical fit)
double R_sd  = log10(M1 / pow(10,5)) * 1000 * pc_m;  // stellar disk scale radius [m]
double R_gd  = 2.0 * R_sd;                            // gas disk scale radius      [m]

// Gas disk vertical scale height. Fixed at 100 pc in the first paper, but
// scales with mass in principle.
double z_thick    = 100.0 * pc_m;
double sigma_star = 6.3E4;        // 1D stellar velocity dispersion [m/s] (Antonini & Merritt)


// ---------------------------------------------------------------------------
// Phase-space state of the secondary in cylindrical coordinates.
// ---------------------------------------------------------------------------
struct State {
    double r;         // cylindrical radius           [m]
    double rdot;      // radial velocity              [m/s]
    double theta;     // azimuthal angle              [rad]
    double thetadot;  // angular velocity             [rad/s]
    double z;         // vertical position            [m]
    double zdot;      // vertical velocity            [m/s]
};

// ===========================================================================
// GLOBAL GRIDS (populated in main() before the integrator runs).
//
// Two grids are used: a coarse outer grid (N x N covering ~2 * R_gd) and a
// fine inner AMR grid (N x N2 with a much thinner vertical extent L2). The
// secondary's force is computed from whichever grid it lies in.
//
// For each density component we store both the 2D array (rho_*) and a
// flattened copy (rho_*_flat) that GSL's 2D interpolator needs as a 1D array.
// ===========================================================================

// Coordinate axes for the two grids
double *r_arr;    double *z_arr;    // outer grid axes
double *r2_arr;   double *z2_arr;   // AMR grid axes

// Sound-speed lookup (unused in main loop -- kept for diagnostics)
double **cs_arr;

// Total density = gas + stars + bulge (+ halo if enabled)
double **rho;     double **rho2;

// Gas disk only (used for gas dynamical friction and BHL accretion)
double **rho_g;   double **rho_g2;

// Stellar disk only (used for disk-DF)
double **rho_s;   double **rho_s2;

// Total stellar = disk + bulge (used for bulge-DF integrand)
double **rho_st;  double **rho_st2;

// Bulge only (added into rho_st and rho in main())
double **rho_bg;  double **rho_bg2;

// Flattened copies for GSL interpolation
double *rho_flat;    double *rho_flat2;
double *rho_g_flat;  double *rho_g_flat2;
double *rho_s_flat;  double *rho_s_flat2;
double *rho_st_flat; double *rho_st_flat2;

// Dark matter halo (NOT used in the first paper). A proper halo would
// require a third, even larger grid: the halo scale radius greatly exceeds
// the disk, and a small grid introduces spurious potential gradients since
// vacuum boundary conditions break down. Left in as a stub.
double **rho_h;
double *rho_h_flat;

// Gravitational potential grids (self-gravity of the galaxy -- does NOT
// include the primary BH point mass, which is added analytically later)
double **phi;     double **phi2;
double **phi_b;                       // (unused scratch array)
double *phi_flat; double *phi_flat2;

// GSL interpolation objects, one per grid/quantity
gsl_interp2d *interp;       gsl_interp_accel *xacc;     gsl_interp_accel *yacc;      // potential (outer)
gsl_interp2d *interp_r;     gsl_interp_accel *xacc_r;   gsl_interp_accel *yacc_r;    // density (outer)
gsl_interp2d *interp_amr;   gsl_interp_accel *xacc_amr; gsl_interp_accel *yacc_amr;  // potential (AMR)


// ===========================================================================
// SOLVER HELPER FUNCTIONS
// ===========================================================================

/*
 Load a previously-generated potential from a CSV file into phi[][].
 This is a big speedup: the Poisson relaxation solver is often the slowest
 part of a run, so once you've generated a potential for a given galaxy
 configuration you can reuse it across many orbital runs.

 Format: one row per grid row, comma-separated floats.
*/
int load_potential_csv_ptrptr(const char *filename,
                              double **phi, int max_rows, int max_cols,
                              int *out_rows, int *out_cols)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) { perror("fopen"); return 0; }

    char line[65536];
    int rows = 0, cols = -1;

    while (fgets(line, sizeof line, fp)) {
        int col = 0;
        char *p = line;

        while (*p) {
            while (*p == ' ' || *p == '\t' || *p == ',') p++;
            if (*p == '\0' || *p == '\n' || *p == '\r') break;

            char *endp;
            errno = 0;
            double v = strtod(p, &endp);
            if (p == endp) break;
            if (errno) { fclose(fp); return 0; }

            if (rows >= max_rows || col >= max_cols) { fclose(fp); return 0; }
            phi[rows][col] = v;
            col++;
            p = endp;
        }

        if (col == 0) continue;
        if (cols == -1) cols = col;
        else if (col != cols) { fclose(fp); return 0; }
        rows++;
        if (rows > max_rows) { fclose(fp); return 0; }
    }

    fclose(fp);
    if (rows == 0 || cols <= 0) return 0;

    if (out_rows) *out_rows = rows;
    if (out_cols) *out_cols = cols;
    return 1;
}

// Sign function: sgn(-3) = -1, sgn(0) = 0, sgn(+7) = +1
double sgn(double x) {
    return (x > 0) - (x < 0);
}

// ===========================================================================
// DENSITY AND POTENTIAL QUERIES
//
// Thin wrappers around GSL's 2D interpolator that look up the pre-computed
// density and potential grids at an arbitrary (r, z) point. All return SI
// mass density (kg/m^3) or potential (J/kg = m^2/s^2).
// ===========================================================================

// Analytic gas density on the midplane (z = 0), used as a reference value
double rho_gd(double r) {
    return n_gd * m_p * pow(e, -1.0 * r / (R_gd));
}

// Gas density at (r, z) -- bilinear interpolation from the grid
double rho_rz(double r, double z) {
    return gsl_interp2d_eval(interp_r, r_arr, z_arr, rho_g_flat, r, z, xacc_r, yacc_r);
}

// Stellar-disk density at (r, z)
double rhos_rz(double r, double z) {
    return gsl_interp2d_eval(interp_r, r_arr, z_arr, rho_s_flat, r, z, xacc_r, yacc_r);
}

// Dark matter halo density at (r, z) -- stub, only meaningful if halo is enabled
double rhoh_rz(double r, double z) {
    return gsl_interp2d_eval(interp_r, r_arr, z_arr, rho_h_flat, r, z, xacc_r, yacc_r);
}

// TOTAL stellar density (disk + bulge) at (r, z) -- feeds the bulge DF integrand
double rhost_rz(double r, double z) {
    return gsl_interp2d_eval(interp_r, r_arr, z_arr, rho_st_flat, r, z, xacc_r, yacc_r);
}

// Gas surface density (integrated vertically assuming thickness 2*z_thick)
double sigma_gd(double r) {
    return rho_gd(r) * (2 * z_thick);
}

// ===========================================================================
// DYNAMICAL FRICTION (DF)
//
// Dynamical friction is the drag force a massive perturber (the secondary BH)
// feels as it moves through a medium of lighter particles (gas molecules,
// stars). The medium's self-gravity pulls the perturber back, draining its
// orbital energy. This is the main mechanism that drives the inspiral in the
// hundreds-of-parsecs regime simulated here.
//
// Two flavors are implemented:
//   - Gas DF (Ostriker 1999; Kim & Kim 2007 fits): depends on the Mach number
//     v / c_s. Piecewise fits are used for subsonic, transonic, and supersonic.
//   - Stellar DF (Chandrasekhar; Antonini & Merritt 2012 for the disk): the
//     disk version is anisotropic and splits into "slow" (v > v_disk) and
//     "fast" (v < v_disk) regimes because stars in a rotating disk carry net
//     azimuthal velocity.
// ===========================================================================

// Common prefactor for gas DF: 4*pi*(G*M2)^2 * rho_gas(r,z).
// Appears in all three (r, theta, z) components below.
double gasDFb(double r, double z) {
    // Guard against a rare interpolation artifact where the bilinear interp
    // produces a very slightly negative density; fall back to the midplane value.
    if (rho_rz(r, z) < 0.0) {
        return 4 * PI * pow((G * M2 * sl_kg), 2) * rho_gd(0.0);
    }
    return 4 * PI * pow((G * M2 * sl_kg), 2) * rho_rz(r, z);
}

// Radial gas DF. Uses the radial component of the secondary's velocity
// relative to the (rotating) gas, and the Kim & Kim 2007 I_r fit.
double gasDFr(double r, double rdot, double tdot, double z, double zdot, double cs, double vd) {
    double rt   = r * tdot;
    double v    = pow((rt + vd) * (rt + vd), 0.5);
    double mach = abs(v) / cs;
    if (mach < 1E-4) return 0;

    double I_r;
    if      (mach < 1.1) I_r = pow(mach, 2) * pow(10, 3.51 * mach - 4.22);
    else if (mach < 4.4) I_r = 0.5 * log(9.33 * pow(mach, 2) * (pow(mach, 2) - 0.95));
    else                 I_r = 0.3 * pow(mach, 2);

    return (-1.0 * gasDFb(r, z) * I_r) / (v * v);
}

// Azimuthal gas DF. This is usually the DOMINANT orbital energy loss channel
// for prograde in-disk motion, and is where most numerical issues tend to
// originate -- if you're chasing a bug, check here first.
double gasDFt(double r, double rdot, double tdot, double z, double zdot, double cs, double vd) {
    double rt = r * tdot;
    // Velocity relative to the rotating disk (vd = vg * v_c):
    double v  = pow((rt + vd) * (rt + vd), 0.5);
    double mach = v / cs;
    if (v < 0.01) return 0;

    double I_t = 0.0;
    if      (mach < 1.0 && mach > 0.05) I_t = 0.7706 * log((1 + mach) / (1.0004 - 0.9185 * mach)) - 1.4703 * mach;
    else if (mach < 4.4 && mach > 1.0)  I_t = log(330 * 10.0 * pow(mach, -9.58) * pow(mach - 0.71, 5.72));
    else if (mach > 4.4)                I_t = log(10 / (0.11 * mach + 1.65));

    if (I_t < 0.0) return 0.0;

    return (-1.0 * gasDFb(r, z) * I_t) / (v * v);
}

// Vertical (z) gas DF. Uses only |zdot|/cs as the relevant Mach number since
// the vertical motion is roughly orthogonal to the disk flow.
double gasDFz(double r, double rdot, double thetadot, double z, double zdot, double cs, double vd) {
    double rt   = r * thetadot;
    double v    = pow((rt + vd) * (rt + vd) + rdot * rdot + zdot * zdot, 0.5);
    double mach = abs(zdot) / cs;
    if (mach < 1E-4) return 0.0;

    double I = 0.0;
    if      (mach > 0.1 && mach < 0.999) I = 0.5 * log((1+mach) / (1-mach)) - mach;
    else if (mach > 1.0)                 I = 0.5 * log(1 - 1 / (mach * mach)) + 10;
    else if (mach < 0.1)                 I = mach * mach * mach / 3.0;

    return (gasDFb(r, z) * I) / (zdot * zdot);
}

// ---------------------------------------------------------------------------
// Stellar dynamical friction (Chandrasekhar form).
// ---------------------------------------------------------------------------

// Common prefactor: -4*pi*G^2*M2^2 * rho_star / v^3
double df_coeff(double M2, double rho_s, double v) {
    if (v < 0.01 * sigma_star) return 0.0;
    return -4 * M_PI * G * G * pow(M2 * sl_kg, 2) * 1.0 * rho_s / pow(v, 3.0);
}

// Dark matter DF (only active if halo is enabled). Kept for completeness.
double df_dm_test(double r, double z, double v) {
    if (v < 1E-4) return 0.0;
    return -4 * M_PI * G * G * pow(M2 * sl_kg, 2) * 1.0 * rhoh_rz(r,z) / pow(0.5 * v, 3.0);
}

// ---------------------------------------------------------------------------
// Disk stellar DF -- "slow" branch (perturber moves faster than disk).
// Anisotropic Chandrasekhar, see Antonini & Merritt 2012 for derivation.
// ---------------------------------------------------------------------------
double df_slow_r_disk(double v_g, double v, double r, double z, double rdot) {
    if (v > abs(v_g)) {
        if ((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g) < 1.0) return 0.0;
        return rdot * df_coeff(M2, rhos_rz(r, z), v) * log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
    }
    return 0.0;
}

double df_slow_t_disk(double v_g, double v, double r, double z, double thetadot) {
    if (v > abs(v_g)) {
        if (log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g)) < 0.0 || rhos_rz(r, z) < 0.0) return 0.0;
        return 1.0 * (r * thetadot) * df_coeff(M2, rhos_rz(r,z), v) * log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
    }
    return 0.0;
}

double df_slow_z_disk(double v_g, double v, double r, double z, double zdot) {
    if (v > abs(v_g) && log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g)) > 0.0) {
        return zdot * df_coeff(M2, rhos_rz(r,z), v) * log((M1 * sl_kg / (v_g * v_g * M2 * sl_kg)) * (v * v - v_g * v_g));
    }
    return 0.0;
}

// ---------------------------------------------------------------------------
// Disk stellar DF -- "fast" branch (perturber moves SLOWER than disk).
// This regime uses a different Coulomb log derived from resonant-scattering
// arguments, hence the log((v_g+v)/(v_g-v)) form.
// ---------------------------------------------------------------------------
double df_fast_r_disk(double v_g, double v, double r, double z, double rdot) {
    if (v < abs(v_g)) {
        if ((v_g + v) / (v_g - v) < 1.0) return 0.0;
        return rdot * df_coeff(M2, rhos_rz(r, z), v) * (log((v_g + v) / (v_g - v)) - 2 * v / v_g);
    }
    return 0.0;
}

double df_fast_t_disk(double v_g, double v, double r, double z, double thetadot) {
    if (v < abs(v_g)) {
        if ((v_g + v) / (v_g - v) < 1.0) return 0.0;
        return r * thetadot * df_coeff(M2, rhos_rz(r, z), v) * (log((v_g + v) / (v_g - v)) - 2 * v / v_g);
    }
    return 0.0;
}

double df_fast_z_disk(double v_g, double v, double r, double z, double zdot) {
    if (v < abs(v_g) && log((v_g + v) / (v_g - v)) > 0.0) {
        return zdot * df_coeff(M2, rhos_rz(r, z), v) * (log((v_g + v) / (v_g - v)) - 2 * v / v_g);
    }
    return 0.0;
}

// Bulge DF is defined farther down because it needs the escape velocity
// (and hence the total potential), which requires phi_rz() below.

// ===========================================================================
// GRAVITATIONAL POTENTIAL AND ACCELERATION
//
// The galaxy's self-gravity is stored in phi[][] (outer grid) and phi2[][]
// (AMR grid). The primary SMBH is added as a point-mass term so we don't
// have to resolve its potential on the grid.
// ===========================================================================

// Total potential at (r, z) = galaxy potential (interpolated) + primary BH point mass.
double phi_rz(double r, double z) {
    return gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc)
         - (G * M1 * sl_kg) / (pow(r * r + z * z, 0.5));
}

// Radial gravitational acceleration, g_r = -d(phi)/dr.
// The AMR switch below decides which grid to take the derivative from:
// when |z| < 0.049 * R_gd the secondary is near the midplane and we use
// the fine z-resolution AMR grid (phi_flat2); otherwise we use the outer
// grid. The primary BH Keplerian term is always added analytically.
double g_fr(double r, double z) {
    double dr = r * 1E-8;

    // AMR transition criterion: near the midplane, switch to the fine grid.
    if (abs(z) < 0.049 * R_gd) {
        return -1.0 * gsl_interp2d_eval_deriv_x(interp_amr, r2_arr, z2_arr, phi_flat2, r, z, xacc_amr, yacc_amr)
             - (G * M1 * sl_kg * r) / (pow(r * r + z * z, 1.5));
    }

    // Guard: if the grid derivative is positive (an artifact near the outer
    // boundary where phi flattens out), zero it out so we don't get unphysical
    // outward gravity.
    if (-1.0 * gsl_interp2d_eval_deriv_x(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc) > 0.0) {
        return 0.0;
    }

    return -1.0 * gsl_interp2d_eval_deriv_x(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc)
         - (G * M1 * sl_kg * r) / (pow(r * r + z * z, 1.5));
}

// Vertical gravitational acceleration, g_z = -d(phi)/dz. Same AMR switch.
double g_fz(double r, double z) {
    if (abs(z) < 0.049 * R_gd) {
        return -1.0 * gsl_interp2d_eval_deriv_y(interp_amr, r2_arr, z2_arr, phi_flat2, r, z, xacc_amr, yacc_amr)
             - (G * M1 * sl_kg * z) / (pow(r * r + z * z, 1.5));
    }
    return -1.0 * gsl_interp2d_eval_deriv_y(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc)
         - (G * M1 * sl_kg * z) / (pow(r * r + z * z, 1.5));
}

// Circular (rotational) speed at cylindrical radius r, computed from the
// midplane radial gravity:  v_c = sqrt(r * |g_r(r, 0)|).
double vc(double r, double z) {
    return pow(r * -1.0 * g_fr(r, 0.0), 0.5);
}

// Escape speed from the BH at (r, z) (used only for diagnostics).
double v_esc_bin(double r, double z) {
    return pow((2.0 * G * M1 * sl_kg) / pow(r * r + z * z, 0.5), 0.5);
}

// ---------------------------------------------------------------------------
// BULGE DYNAMICAL FRICTION
//
// The bulge has an ISOTROPIC stellar velocity distribution fs_v() (defined
// in poisson.h), so the Chandrasekhar integral must be done numerically over
// stellar speeds. Below, params = {r, z, v_c, v_perturber, sigma_star, v_esc}.
// ---------------------------------------------------------------------------

// Integrand of the bulge DF integral. Has two branches depending on whether
// the secondary is faster or slower than a given background star (vs):
//   - If the perturber is faster (but not too much): standard Chandrasekhar log
//   - If the perturber is slower:                    resonant-scattering form
double bulge_integrand(double vs, const std::vector<double>& params) {
    if ((1 / (q * params[4] * params[4])) * (params[3] * params[3] - vs * vs) > 1.0) {
        return 4 * M_PI * fs_v(params[2], vs, params[5]) * vs * vs
             * log((1 / (q * params[4] * params[4])) * (params[3] * params[3] - vs * vs));
    } else if (params[3] < vs && params[3] > 0.0) {
        return 4 * M_PI * fs_v(params[2], vs, params[5]) * vs * vs
             * (log((vs + params[3]) / (vs - params[3])) - 2 * params[3] / vs);
    }
    return 0.0;
}

// Radial, azimuthal, and vertical components of the bulge DF force.
// Each integrates the "slow" contribution (stars slower than the perturber)
// and "fast" contribution (stars faster) and sums them.
double df_bulge_r(double r, double z, double v, double rdot) {
    double v_esc = pow(-2.0 * gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc), 0.5);
    double df_slow = df_coeff(M2, rhost_rz(r, z), v) * integrate1D(bulge_integrand, 0, v,     20, {r, z, vc(r, z), v, sigma_star, v_esc});
    double df_fast = df_coeff(M2, rhost_rz(r, z), v) * integrate1D(bulge_integrand, v, v_esc, 20, {r, z, vc(r, z), v, sigma_star, v_esc});
    return rdot * (df_slow + df_fast);
}

double df_bulge_t(double r, double z, double v, double tdot) {
    double v_esc = pow(-2.0 * gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc), 0.5);
    if (v > 0.0 * sigma_star) {
        double df_slow = 1.0 * df_coeff(M2, rhost_rz(r, z), v) * integrate1D(bulge_integrand, 0, v,     20, {r, z, vc(r, z), v, sigma_star, v_esc});
        double df_fast = 1.0 * df_coeff(M2, rhost_rz(r, z), v) * integrate1D(bulge_integrand, v, v_esc, 20, {r, z, vc(r, z), v, sigma_star, v_esc});
        return r * tdot * (df_slow + df_fast);
    }
    return 0.0;
}

double df_bulge_z(double r, double z, double v, double zdot) {
    double v_esc = pow(-2.0 * gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc), 0.5);
    double df_slow = df_coeff(M2, rhost_rz(r, z), v) * integrate1D(bulge_integrand, 0, v,     20, {r, z, vc(r, z), v, sigma_star, v_esc});
    double df_fast = df_coeff(M2, rhost_rz(r, z), v) * integrate1D(bulge_integrand, v, v_esc, 20, {r, z, vc(r, z), v, sigma_star, v_esc});
    return 1.0 * zdot * (df_slow + df_fast);
}

// Dark matter halo DF z-component (stub -- disabled unless halo is enabled).
double df_dm_z(double r, double z, double v, double zdot) {
    double v_esc = pow(-2.0 * gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat, r, z, xacc, yacc), 0.5);
    double df_slow = df_coeff(M2, rhoh_rz(r, z), v) * integrate1D(bulge_integrand, 0, v,     20, {r, z, vc(r, z), v, sigma_star, v_esc});
    double df_fast = df_coeff(M2, rhoh_rz(r, z), v) * integrate1D(bulge_integrand, v, v_esc, 20, {r, z, vc(r, z), v, sigma_star, v_esc});
    return zdot * (df_slow + df_fast);
}

// ===========================================================================
// EQUATIONS OF MOTION (RHS of the RK4)
//
// These compute the accelerations on the secondary in cylindrical coordinates.
// Each is the sum of:
//    (1) Galactic gravity (self-gravity of gas/stars/bulge + primary BH),
//    (2) Centrifugal / Coriolis terms (for r and theta equations),
//    (3) Gas dynamical friction,
//    (4) Stellar disk DF (slow + fast),
//    (5) Bulge DF.
//
// Each acceleration is divided by M2 * sl_kg to get m/s^2 (the DF forces are
// written above as forces, not accelerations).
//
// The `eps` offset near r < 0.5 pc is a crude regularization to avoid the
// 1/r singularity in the theta equation at the very center.
// ===========================================================================

// Radial acceleration (centrifugal term included):
//   r_dotdot = g_r + r*thetadot^2 + DF_r
double fr(double r, double rdot, double thetadot, double z, double zdot, double cs) {
    double v   = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0) + pow(zdot, 2.0), 0.5);
    double eps = 1.0 * pc_m;  // regularization offset near r = 0

    if (r < 0.5 * pc_m) {
        return (1.0 / (M2 * sl_kg)) * (
              M2 * sl_kg * g_fr(r + eps, z)
            + (M2 * sl_kg) * (r + eps) * pow(thetadot, 2)
            + 1.0 * gasDFr(r + eps, rdot, thetadot, z, zdot, cs, vg * vc(r + eps, z))
            + 1.0 * df_slow_r_disk(vg * vc(r + eps, z), v, r + eps, z, rdot)
            + 1.0 * df_fast_r_disk(vg * vc(r + eps, z), v, r + eps, z, rdot)
            + 1.0 * df_bulge_r(r + eps, z, v, rdot));
    }
    return (1.0 / (M2 * sl_kg)) * (
          M2 * sl_kg * g_fr(r, z)
        + (M2 * sl_kg) * r * pow(thetadot, 2)
        + 1.0 * gasDFr(r, rdot, thetadot, z, zdot, cs, vg * vc(r, z))
        + 1.0 * df_slow_r_disk(vg * vc(r, z), v, r, z, rdot)
        + 1.0 * df_fast_r_disk(vg * vc(r, z), v, r, z, rdot)
        + 1.0 * df_bulge_r(r, z, v, rdot));
}

// Angular acceleration (Coriolis term included):
//   thetadot_dot = (1/r) * [ -2*rdot*thetadot + DF_theta / M2 ]
double ftheta(double r, double rdot, double theta, double thetadot, double z, double zdot, double cs) {
    double v   = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0) + pow(zdot, 2.0), 0.5);
    double eps = 1.0 * pc_m;

    if (r < 0.5 * pc_m) {
        return (1.0 / ((abs(r + eps)) * M2 * sl_kg)) * (
              (-2.0) * rdot * thetadot * M2 * sl_kg
            + 1.0 * gasDFt(r + eps, rdot, thetadot, z, zdot, cs, vg * vc(r + eps, z))
            + 1.0 * df_slow_t_disk(vg * vc(r + eps, z), v, r + eps, z, thetadot)
            + 1.0 * df_fast_t_disk(vg * vc(r + eps, z), v, r + eps, z, thetadot)
            + 1.0 * df_bulge_t(r + eps, z, v, thetadot));
    }
    return (1.0 / (abs(r) * M2 * sl_kg)) * (
          (-2.0) * rdot * thetadot * M2 * sl_kg
        + 1.0 * gasDFt(r, rdot, thetadot, z, zdot, cs, vg * vc(r, z))
        + 1.0 * df_slow_t_disk(vg * vc(r, z), v, r, z, thetadot)
        + 1.0 * df_fast_t_disk(vg * vc(r, z), v, r, z, thetadot)
        + 1.0 * df_bulge_t(r, z, v, thetadot));
}

// Vertical acceleration:
//   z_dotdot = g_z + DF_z / M2
double fz(double r, double rdot, double thetadot, double z, double zdot, double cs) {
    double v   = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0) + pow(zdot, 2.0), 0.5);
    double eps = 1.0 * pc_m;

    if (r < 0.5 * pc_m) {
        return g_fz(r + eps, z)
             + -1.0 * sgn(zdot) * gasDFz(r + eps, rdot, thetadot, z, zdot, cs, vg * vc(r, z)) / (M2 * sl_kg)
             + 1.0 * df_bulge_z(r + eps, z, v, zdot) / (M2 * sl_kg)
             + 1.0 * df_slow_z_disk(vg * vc(r + eps, z), v, r + eps, z, zdot) / (M2 * sl_kg)
             + 1.0 * df_fast_z_disk(vg * vc(r + eps, z), v, r + eps, z, zdot) / (M2 * sl_kg);
    }
    return g_fz(r, z)
         + -1.0 * sgn(zdot) * gasDFz(r, rdot, thetadot, z, zdot, cs, vg * vc(r, z)) / (M2 * sl_kg)
         + 1.0 * df_bulge_z(r, z, v, zdot) / (M2 * sl_kg)
         + 1.0 * df_slow_z_disk(vg * vc(r, z), v, r, z, zdot) / (M2 * sl_kg)
         + 1.0 * df_fast_z_disk(vg * vc(r, z), v, r, z, zdot) / (M2 * sl_kg);
}

// ---------------------------------------------------------------------------
// Cartesian versions of the force, used by the optional Cartesian RK4 branch
// in RK4Solver() (currently disabled by the `if (1E30 < 100*pc)` guard, but
// kept for testing behavior through r = 0 -- the cylindrical coords are
// singular at the axis).
// ---------------------------------------------------------------------------
double force_r(double r, double rdot, double thetadot, double z, double zdot, double cs) {
    double v = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0) + pow(zdot, 2.0), 0.5);
    return (1.0 / (M2 * sl_kg)) * (
          M2 * sl_kg * g_fr(r, z)
        + 1.0 * gasDFr(r, rdot, thetadot, z, zdot, cs, vg * vc(r, z))
        + 1.0 * df_slow_r_disk(vg * vc(r, z), v, r, z, rdot)
        + 1.0 * df_fast_r_disk(vg * vc(r, z), v, r, z, rdot)
        + 1.0 * df_bulge_r(r, z, v, rdot));
}

double force_theta(double r, double rdot, double theta, double thetadot, double z, double zdot, double cs) {
    double v = pow(pow(rdot, 2.0) + pow(r * thetadot, 2.0) + pow(zdot, 2.0), 0.5);
    return (1.0 / (M2 * sl_kg)) * (
          1.0 * gasDFt(r, rdot, thetadot, z, zdot, cs, vg * vc(r, z))
        + 1.0 * df_slow_t_disk(vg * vc(r, z), v, r, z, thetadot)
        + 1.0 * df_fast_t_disk(vg * vc(r, z), v, r, z, thetadot)
        + 1.0 * df_bulge_t(r, z, v, thetadot));
}

// Cartesian accelerations (x, y, z). Converts back to cylindrical internally
// to call the force_r / force_theta helpers above.
double a_x(double x, double xdot, double y, double ydot, double z, double zdot, double cs) {
    double r = pow(x * x + y * y, 0.5);
    double theta = atan(y / x);
    double rdot = xdot * cos(theta) + ydot * sin(theta);
    double thetadot = (cos(theta) * ydot - xdot * sin(theta)) / r;
    if (r * thetadot > 3E8) return 0.0;  // guard against superluminal junk
    return force_r(r, rdot, thetadot, z, zdot, cs) * cos(theta)
         - force_theta(r, rdot, theta, thetadot, z, zdot, cs) * sin(theta);
}

double a_y(double x, double xdot, double y, double ydot, double z, double zdot, double cs) {
    double r = pow(x * x + y * y, 0.5);
    double theta = atan(y / x);
    double rdot = xdot * cos(theta) + ydot * sin(theta);
    double thetadot = (cos(theta) * ydot - xdot * sin(theta)) / r;
    if (r * thetadot > 3E8) return 0.0;
    return force_r(r, rdot, thetadot, z, zdot, cs) * sin(theta)
         + force_theta(r, rdot, theta, thetadot, z, zdot, cs) * cos(theta);
}

double a_z(double x, double xdot, double y, double ydot, double z, double zdot, double cs) {
    double r = pow(x * x + y * y, 0.5);
    double theta = atan(y / x);
    double rdot = xdot * cos(theta) + ydot * sin(theta);
    double thetadot = (cos(theta) * ydot - xdot * sin(theta)) / r;
    if (r * thetadot > 3E8) return 0.0;
    return fz(r, rdot, thetadot, z, zdot, cs);
}

// ===========================================================================
// RK4 ORBIT INTEGRATOR
//
// Integrates the secondary's equations of motion (fr, ftheta, fz above) from
// the supplied initial State to `total_time`, using classical 4th-order
// Runge-Kutta with a user-supplied time step `dt` (adaptive refinement kicks
// in near the center).
//
// At each step, in addition to advancing the orbit, the routine:
//   - Tracks the specific energy E and semi-major axis `ak`
//   - Detects orbital peaks in z (for oscillation-period diagnostics)
//   - Computes BHL accretion rates onto both black holes
//   - Updates M1, M2 via accretion
//   - Estimates jet power via a Blandford-Znajek-like scaling
//   - Decides whether the run was a merger (flyby=1) or a flyby (flyby=2)
//
// The final exit condition is:
//   - t exceeds total_time (likely a flyby that didn't merge), OR
//   - the secondary reaches within 10 pc of the primary AND is gravitationally
//     bound (-> merger), OR
//   - the secondary flies out past 10 kpc in r (numerical escape).
//
// Every diagnostic is logged to one of many text files listed below.
// ===========================================================================
State RK4Solver(State initial_state, double dt, double total_time, double angle) {
    State current_state = initial_state;
    double t = 0.0;

    // -------------------------------------------------------------------
    // Output filenames (one file per diagnostic). Opening them up front
    // keeps the main loop clean.
    // -------------------------------------------------------------------
    // Orbital state & energetics
    FILE* fp       = fopen("bh_test2.txt", "w");    // (t, r)
    FILE* fm       = fopen("energy.txt",   "w");    // (t, E)
    FILE* fv       = fopen("z.txt",        "w");    // (t, z)
    FILE* ft       = fopen("tau.txt",      "w");    // (t, period_z)  -- vertical oscillation period
    FILE* flz      = fopen("lz.txt",       "w");    // (t, L_z)
    FILE* forb     = fopen("orbit.txt",    "w");    // (t, r, theta, z)

    // Peak-of-z diagnostics
    FILE* frad     = fopen("avgr.txt",     "w");
    FILE* favgz    = fopen("avgz.txt",     "w");
    FILE* finc     = fopen("inc.txt",      "w");
    FILE* flinc    = fopen("l_inc.txt",    "w");
    FILE* fany     = fopen("firstfew.txt", "w");    // early-time trajectory (r, z) points
    FILE* fphi     = fopen("phi_period.txt","w");   // azimuthal orbital period

    // Mach numbers and energy loss rates
    FILE* fmach    = fopen("mach1.txt",    "w");
    FILE* fpr      = fopen("pr.txt",       "w");    // dE_r/dt
    FILE* fpt      = fopen("pt.txt",       "w");    // dE_theta/dt
    FILE* fpz      = fopen("pz.txt",       "w");    // (t, z, zdot)

    // Accretion, jets, mass growth
    FILE* fmdot    = fopen("mdot.txt",     "w");
    FILE* fmdote   = fopen("mdot_edd.txt", "w");
    FILE* fflux    = fopen("flux.txt",     "w");
    FILE* fpbz     = fopen("pbz.txt",      "w");    // secondary jet power
    FILE* fpbz1    = fopen("pbz1.txt",     "w");    // primary jet power
    FILE* fpbzavg  = fopen("pbz_avg.txt",  "w");
    FILE* ffluxd   = fopen("fluxd.txt",    "w");
    FILE* fmass    = fopen("mass.txt",     "w");
    FILE* fmass1   = fopen("mass1.txt",    "w");
    FILE* fq       = fopen("q.txt",        "w");
    FILE* fmdot1   = fopen("mdot1.txt",    "w");
    FILE* faccavg  = fopen("accavg.txt",   "w");
    FILE* feacc    = fopen("e_acc.txt",    "w");
    FILE* fejet    = fopen("e_jet.txt",    "w");

    // Run-summary files (APPEND mode -- persist across multiple invocations)
    FILE* ftime       = fopen("decay_time_params_frfrfr.txt", "a");
    FILE* flog        = fopen("tau_log.txt",  "a");
    FILE* fmachavg    = fopen("time_machavg.txt", "a");
    FILE* facc_log0   = fopen("acc_log0.txt",   "a");
    FILE* facc_log1   = fopen("acc_log1.txt",   "a");
    FILE* facc_logmax = fopen("acc_logmax.txt", "a");

    // -------------------------------------------------------------------
    // Scratch / state variables used inside the main integration loop.
    // -------------------------------------------------------------------

    // Orbital invariants
    double E0;                   // initial specific energy (kept for drift diagnostics)
    double ak;                   // semi-major axis (instantaneous, from vis-viva)
    double dt0 = dt;             // save the nominal time step (gets varied near center)
    int    flyby = 0;            // outcome flag: 0 = unfinished, 1 = merger, 2 = flyby

    // Period-detection bookkeeping
    double tau_z  = 0.0;         // time between vertical peaks
    double t0_z   = 0.0;         // last time a z-peak was detected
    double tau_p  = 0.0;         // time between azimuthal crossings (theta wrap)
    double t_phi  = 0.0;         // last time theta wrapped
    double phi_0  = initial_state.theta;
    bool   peak   = false;
    bool   upz    = current_state.zdot > 0.0;

    // Gas thermodynamics (updated every step from the Toomre-based prescription)
    double T  = 1E4;
    double cs = pow((5 * kb * T) / (3 * m_p), 0.5);

    // Mach numbers & instantaneous energy-loss rates (for diagnostics)
    double mach    = 0.0;   // in-plane orbital Mach number (r*thetadot / cs)
    double macheff = 0.0;   // effective Mach relative to rotating gas
    double machv   = 0.0;   // vertical Mach number (zdot / cs)
    double dErdt   = 0.0;   // radial power loss
    double dEtdt   = 0.0;   // azimuthal power loss
    double dEzdt   = 0.0;   // vertical power loss
    double machavg = 0.0;   // time-averaged |machv|
    double tauavg  = 0.0;   // time-averaged vertical period

    // Magnetic-flux / jet diagnostics
    double b_inf    = 0.0;  // asymptotic B field on the BH
    double flux     = 0.0;  // dimensionful magnetic flux
    double flux_dim = 0.0;  // dimensionless magnetic flux (phi_BH)
    double pbz      = 0.0;  // secondary jet power
    double pbz1     = 0.0;  // primary jet power

    // Per-BH geometric / accretion quantities
    double a = 0.9;                                  // dimensionless BH spin
    double r_g  = G * M2 * sl_kg / (c * c);          // gravitational radius of secondary
    double r_g1;                                     // gravitational radius of primary (set in loop)
    double r_a  = 0.0;                               // accretion radius
    double r_a1 = 0.0;
    double r_b1 = 0.0;                               // Bondi radius (primary)
    double r_b2 = 0.0;                               // Bondi radius (secondary)
    double m_bhl = 0.0;                              // secondary BHL accretion rate
    double m_edd = 0.0;                              // secondary Eddington rate
    double m1_b  = 0.0;                              // primary Bondi accretion rate
    double m1_edd = 0.0;                             // primary Eddington rate
    double q     = M2 / M1;                          // (shadows global q -- re-evaluated as masses grow)
    double cs00  = 0.0;                              // sound speed at primary's Bondi radius
    double T0    = 0.0;                              // gas temperature at the primary
    double M10   = M1;                               // initial primary mass (for accretion bookkeeping)
    double m1_max = 0.0;
    double flux1  = 0.0;
    double ra_max1 = -1.0 * G * M1 * sl_kg / (0.1 * phi_rz(1.0 * pc_m, 0.0));
    double p_avg   = 0.0;                            // time-averaged jet power
    double acc_avg = 0.0;                            // time-averaged accretion rate

    // Alpha-disk (Shakura-Sunyaev) reference scales. These set the thin-disk
    // accretion rate, used in the viscous-timescale prescription below.
    double m_ad = 0.01 * M1 * sl_kg;            // disk mass (primary)
    double r_ad = 1.0 * (M1 / 1E8) * pc_m;      // disk outer radius (primary)
    double h_ad = 0.01 * r_ad;                  // disk scale height (primary)
    double m_ad2 = 0.01 * M2 * sl_kg;
    double r_ad2 = 1.0 * (M2 / (1E8 / 9)) * pc_m;
    double h_ad2 = 0.01 * r_ad2;
    double alpha = 0.1;                         // Shakura-Sunyaev viscosity parameter

    double t_nu  = (r_ad  * r_ad  / (alpha * h_ad  * h_ad )) * pow(pow(r_ad,  3.0) / (G * M1 * sl_kg), 0.5);
    double t_nu2 = (r_ad2 * r_ad2 / (alpha * h_ad2 * h_ad2)) * pow(pow(r_ad2, 3.0) / (G * M2 * sl_kg), 0.5);
    double m_dot_thin  = m_ad  / t_nu;
    double m_dot_thin2 = m_ad2 / t_nu2;
    double acc_max2 = 0.0, acc_max2_edd = 0.0;

    // Radiative / jet efficiencies (updated each step based on Eddington ratio)
    double e_acc2 = 1.0, e_rad2 = 0.0, e_jet2 = 0.0;
    double e_acc1 = 1.0, e_rad1 = 0.0, e_jet1 = 0.0;
    double rtr_0 = 1E4;     // transition radius normalization (thin to thick disk)
    double rth1 = 0.0, rth2 = 0.0;

    // BH spin geometry: ISCO radius, horizon radius, horizon angular velocity.
    double z_1 = 1 + pow((1 - a*a), 1.0/3.0) * (pow(1 + a, 1.0/3.0) + pow(1 - a, 1.0/3.0));
    double z_2 = pow(3 * a * a + z_1 * z_1, 0.5);
    double r_isco   = 3 + z_2 - pow((3.0 - z_1) * (3 + z_1 + 2.0 * z_2), 0.5);
    double r_h      = 1 + pow(1 - a * a, 0.5);
    double omega_h  = a / (2 * r_h);
    // Tchekhovskoy et al. fit for dimensionless magnetic flux at saturation
    double flux_d0 = -20.2 * pow(a, 3.0) - 14.9 * pow(a, 2.0) + 34 * a + 52.6;
    double flux_d1 = 0.0, flux_d2 = 0.0;
    double kappa = 0.05;                          // jet-power normalization

    // Eddington ratios (f_Edd) for both BHs: d = disk (geometric/BHL), h = hole (post-efficiency)
    double fd_edd1 = 0.0, fh_edd1 = 0.0;
    double fd_edd2 = 0.0, fh_edd2 = 0.0;

    double mdot_h1 = 0.0, mdot_h2 = 0.0;     // accretion rates onto the holes after efficiency cuts
    double mdot_bh1 = 0.0, mdot_bh2 = 0.0;   // mass actually added (accounting for rad + jet losses)

    // Vertical peak detection: watch for zdot changing from + to -. The noise
    // threshold `epsilon` is a tiny fraction of the in-plane velocity.
    double prev_zdot = current_state.zdot;
    const double epsilon = 1e-4 * current_state.r * current_state.thetadot;

    // ===================================================================
    // MAIN INTEGRATION LOOP
    //
    // Exit when:  time runs out,  OR  secondary hits center (r < 0),  OR
    //             secondary escapes past 10 kpc (r > 10^4 pc).
    // ===================================================================
    while (t < total_time && pow(pow(current_state.r, 2.0) + pow(current_state.z, 2.0), 0.5) > 0.0 * pc_m && current_state.r < 1.0E4 * pc_m) {
        double r = current_state.r;
        if (r < 0.0) break;
        double rdot     = current_state.rdot;
        double theta    = current_state.theta;
        double thetadot = current_state.thetadot;
        double z        = current_state.z;
        double zdot     = current_state.zdot;

        // ---------------------------------------------------------------
        // Specific energy E and semi-major axis ak (Keplerian vis-viva).
        // Logged every step so we can monitor energy drift (a good sanity
        // check on the integrator).
        // ---------------------------------------------------------------
        double E;
        if (t == 0.0) {
            E0 = 0.5 * (rdot * rdot + r * r * thetadot * thetadot + zdot * zdot) + phi_rz(r, z);
            ak = (-1.0 * G * (M1 + M2) * 2E30 * pow(r * r + z * z, 0.5)) / (pow(r * r + z * z, 0.5) * (rdot * rdot + r * r * thetadot * thetadot + zdot * zdot) - 2 * G * (M1 + M2) * 2E30);
            E = E0;
        } else {
            E = 0.5 * (rdot * rdot + r * r * thetadot * thetadot + zdot * zdot) + phi_rz(r, z);
            ak = (-1.0 * G * (M1 + M2) * 2E30 * pow(r * r + z * z, 0.5)) / (pow(r * r + z * z, 0.5) * (rdot * rdot + r * r * thetadot * thetadot + zdot * zdot) - 2 * G * (M1 + M2) * 2E30);
        }

        // Angular momentum components and inclination
        double L_z     = r * r * thetadot;                 // z-component (conserved if axisymmetric)
        double L_lat   = r * zdot;                          // "lateral" component (z motion)
        double L_total = pow(L_z * L_z + L_lat * L_lat, 0.5);
        double L_inc   = atan(L_z / L_lat);                 // orbital inclination from L

        double v = pow(rdot * rdot + r * r * thetadot * thetadot + zdot * zdot, 0.5);

        // ---------------------------------------------------------------
        // Merger / flyby detection. If the secondary gets within 50 pc,
        // check whether it's gravitationally bound enough to stay:
        //   - r < 10 pc OR   (E deeper than phi at 10 pc along r and z)
        //       -> treat as MERGER (flyby = 1, despite the misleading name)
        //   - otherwise it's a close pass without capture (flyby = 2)
        // ---------------------------------------------------------------
        if (pow(r * r + z * z, 0.5) < 50.0 * pc_m) {
            if (pow(r * r + z * z, 0.5) < 10.0 * pc_m || (E < phi_rz(10.0 * pc_m, 0.0) && E < phi_rz(0.0, 10.0 * pc_m))) {
                flyby = 1;  // merger (name is historical; keep for output compatibility)
                break;
            } else {
                flyby = 2;  // real flyby (no merger)
            }
        }

        // ---------------------------------------------------------------
        // Gas sound speed from the Toomre criterion: a marginally-stable
        // disk (Q ~ 1) gives  cs = pi*G*Sigma / kappa ~ G*pi*Sigma / Omega.
        // We add a 1e4 K floor so the gas doesn't get unphysically cold.
        // ---------------------------------------------------------------
        cs = (G * M_PI * 2 * rho_rz(r, z) * z_thick) / ((2 * abs(vg) * vc(r, z)) / r);
        T  = (3 * m_p * cs * cs) / (5 * kb);
        T += 1E4;                                    // warm-gas floor
        cs = pow(5 * kb * T / (3 * m_p), 0.5);

        // Mach numbers: in-plane (mach), against rotating gas (macheff), vertical (machv)
        mach    = r * thetadot / cs;
        macheff = abs(r * thetadot + vg * vc(r, z)) / cs;
        machv   = zdot / cs;
        if (t == 0.0) machavg = 0.0;
        else          machavg += dt * abs(machv);

        // Sound speed at the primary BH's location (for its own Bondi rate)
        cs00 = (G * M_PI * 2 * rho_rz(0.0 * pc_m, 0) * z_thick) / ((2 * abs(vg) * vc(10.0 * pc_m, 0)) / (10.0 * pc_m));
        T0   = (3 * m_p * cs00 * cs00) / (5 * kb);
        T0  += 1E4;
        cs00 = 100.0 * pow(5 * kb * T0 / (3 * m_p), 0.5);

        // ---------------------------------------------------------------
        // Azimuthal-period detection: theta wraps past 0, so a jump down
        // in theta means one full orbit has elapsed.
        // ---------------------------------------------------------------
        if (current_state.theta < phi_0) {
            tau_p = t - t_phi;
            t_phi = t;
            fprintf(fphi, "%E\t%E\n", t, tau_p / (yr_s * 1E9));
        }
        phi_0 = current_state.theta;

        // ---------------------------------------------------------------
        // Vertical-peak detection: if zdot flips from + to -, we just
        // passed apoapsis in z. Record the period, radius, and inclination.
        // ---------------------------------------------------------------
        if (prev_zdot > epsilon && zdot < -epsilon) {
            tau_z = t - t0_z;
            t0_z  = t;
            fprintf(ft,    "%E\t%E\n", t, tau_z / (yr_s * 1E9));
            fprintf(frad,  "%E\t%E\n", t, current_state.r);
            fprintf(favgz, "%E\t%E\n", t, current_state.z);
            fprintf(finc,  "%E\t%E\n", t, atan(current_state.z / current_state.r));
        }
        tauavg += tau_z * dt;

        if (t < 2E8 * yr_s) {
            // Early-time (r, z) dump: useful for plotting the first few orbits
            fprintf(fany, "%E\t%E\n", current_state.r, current_state.z);
        }

        prev_zdot = zdot;  // save for next step's peak check

        // ---------------------------------------------------------------
        // Instantaneous energy-loss rates along each axis (v dot F).
        // These tell you which DF channel is dominant at each moment.
        // ---------------------------------------------------------------
        dErdt = rdot * force_r(r, rdot, thetadot, z, zdot, cs);
        dEtdt = r * thetadot * force_theta(r, rdot, theta, thetadot, z, zdot, cs);
        dEzdt = zdot * fz(r, rdot, thetadot, z, zdot, cs);
        

        // /* Magnetic stuff */
        // b_inf = pow(2.0 * mu_0 * rho_rz(r, z) * cs * cs, 0.5) * 1E4; //GAUSS
        // if (r * thetadot + vg * vc(r,z) > 0.0) {
        //     r_a = 2.0 * G * (M2 * sl_kg) / pow(r * thetadot + vg * vc(r,z), 2.0); //meter
        //     //flux can't be accreted that fast
        // } else {
        //     r_a = 0.0;
        // }

        r_g = G * M2 * (sl_kg) / (c * c);
        r_g1 = G * M1 * sl_kg / (c * c);
        // ===============================================================
        // BHL ACCRETION, RADIATIVE + JET EFFICIENCIES, MASS UPDATE
        //
        // Each BH accretes gas via the Bondi-Hoyle-Lyttleton formula:
        //     Mdot_BHL = 2*pi * (G*M)^2 * rho / (v^2 + cs^2)^(3/2)
        // The Eddington rate m_edd sets the upper envelope. The ratio
        // f_Edd = Mdot_BHL / m_edd selects the accretion regime:
        //   f_Edd < 0.01  -> radiatively inefficient (RIAF):
        //                    e_acc < 1, e_rad = 0   (mass flows outward)
        //   0.01-1        -> thin disk (Shakura-Sunyaev):
        //                    e_acc = 1, e_rad = NT efficiency at r_isco
        //   > 1           -> super-Eddington / slim disk:
        //                    e_acc = 0.01, e_rad = 0
        // Jet power is added via a Blandford-Znajek-like scaling with the
        // saturated dimensionless flux flux_d (Tchekhovskoy-style fit).
        // ===============================================================

        // Secondary BH: BHL rate, Eddington rate, Bondi radius
        m_bhl  = 2.0 * M_PI * pow(G * M2 * sl_kg, 2.0) * rho_rz(r, z)
               / pow(pow(r * thetadot + vg * vc(r, z), 2.0) + cs * cs, 1.5);   // [kg/s]
        m_edd  = 2.2E-8 * M2 * sl_kg / yr_s;                                    // [kg/s]
        r_b2   = G * M2 * sl_kg / (pow(r * thetadot + vg * vc(r, z), 2.0) + cs * cs);

        // Primary BH: Bondi rate using the central gas density + cs00
        m1_b   = 2.0 * M_PI * pow(G * M1 * sl_kg, 2.0) * rho_rz(0.0, 0.0) / pow(cs00 * cs00, 1.5);
        m1_edd = 2.2E-8 * M1 * sl_kg / yr_s;
        r_b1   = G * M1 * sl_kg / (cs00 * cs00);

        // Disk Eddington ratios (before efficiency factors)
        fd_edd2 = m_bhl / m_edd;
        fd_edd1 = m1_b / m1_edd;

        // Thin-to-thick truncation radius (min of Bondi/r_g and f_Edd-based cutoff)
        rth2 = std::min(r_b2 / r_g,  rtr_0 * pow(0.01 / fd_edd2, 2.0));
        rth1 = std::min(r_b1 / r_g1, rtr_0 * pow(0.01 / fd_edd1, 2.0));

        // Secondary BH efficiency branches (RIAF / thin / super-Edd)
        if (fd_edd2 < 0.01) {
            flux_d2 = flux_d0;
            e_acc2  = pow(10 / rth2, 0.5);
            e_rad2  = 0.0;
        } else if (fd_edd2 < 1.0) {
            flux_d2 = flux_d0 * pow(fh_edd2 / 1.88, 1.29) / (1 + pow(fh_edd2 / 1.88, 1.29));
            e_acc2  = 1.0;
            e_rad2  = 1 - pow(1.0 - 2.0 / (3.0 * r_isco), 0.5);   // Novikov-Thorne at ISCO
        } else {
            flux_d2 = flux_d0 * pow(fh_edd2 / 1.88, 1.29) / (1 + pow(fh_edd2 / 1.88, 1.29));
            e_acc2  = 0.01;
            e_rad2  = 0.0;
        }
        // Blandford-Znajek jet efficiency
        e_jet2 = (kappa / (4.0 * M_PI)) * pow(flux_d2, 2.0) * pow(omega_h, 2.0)
               * (1.0 + 1.38 * pow(omega_h, 2.0) - 9.2 * pow(omega_h, 4.0));

        // Primary BH: same branching
        if (fd_edd1 < 0.01) {
            flux_d1 = flux_d0;
            e_acc1  = pow(10 / rth1, 0.5);
            e_rad1  = 0.0;
        } else if (fd_edd1 < 1.0) {
            flux_d1 = flux_d0 * pow(fh_edd1 / 1.88, 1.29) / (1 + pow(fh_edd1 / 1.88, 1.29));
            e_acc1  = 1.0;
            e_rad1  = 1 - pow(1.0 - 2.0 / (3.0 * r_isco), 0.5);
        } else {
            flux_d1 = flux_d0 * pow(fh_edd1 / 1.88, 1.29) / (1 + pow(fh_edd1 / 1.88, 1.29));
            e_acc1  = 0.01;
            e_rad1  = 0.0;
        }
        e_jet1 = (kappa / (4.0 * M_PI)) * pow(flux_d1, 2.0) * pow(omega_h, 2.0)
               * (1.0 + 1.38 * pow(omega_h, 2.0) - 9.2 * pow(omega_h, 4.0));

        // Mass that actually crosses the horizon each step:
        //   mdot_h   = e_acc * m_BHL                           (flow into disk)
        //   mdot_bh  = (1 - e_rad - e_jet) * e_acc * m_BHL     (growth of M)
        mdot_h1  = e_acc1 * m1_b;
        mdot_h2  = e_acc2 * m_bhl;
        fh_edd1  = mdot_h1 / m1_edd;
        fh_edd2  = mdot_h2 / m_edd;
        mdot_bh1 = (1.0 - e_rad1 - e_jet1) * e_acc1 * m1_b;
        mdot_bh2 = (1.0 - e_rad2 - e_jet2) * e_acc2 * m_bhl;

        // Update both BH masses (remember: M1, M2 are stored in solar masses)
        M1 += (mdot_bh1 * dt) / sl_kg;
        M2 += (mdot_bh2 * dt) / sl_kg;

        // Jet powers in cgs-like units (multiplied by c^2)
        pbz  = e_jet2 * e_acc2 * m_bhl * c * c;
        pbz1 = e_jet1 * e_acc1 * m1_b  * c * c;

        // ===============================================================
        // RK4 STEP
        //
        // Two branches are available:
        //   (A) A Cartesian (x, y, z) RK4, used to handle trajectories
        //       that pass through the cylindrical-axis singularity.
        //   (B) The default cylindrical (r, theta, z) RK4.
        //
        // The Cartesian branch is GATED OFF by the `if (1E30 < 100*pc_m)`
        // dummy condition -- flip that to `if (r < <threshold>)` if you
        // ever need to enable the coordinate switch.
        // ===============================================================

        // Cartesian placeholders (consumed inside branch A only)
        double x    = r * cos(theta);
        double y    = r * sin(theta);
        double xdot = rdot * cos(theta) - r * thetadot * sin(theta);
        double ydot = rdot * sin(theta) + r * thetadot * cos(theta);
        double zsub  = z;
        double zdsub = zdot;

        if (1E30 < 100.0 * pc_m) {  // always false -- selects branch (A)
            // ---- (A) CARTESIAN RK4: disabled by default ----
            double k1r  = xdot;
            double k1rd = a_x(x, xdot, y, ydot, z, zdot, cs);
            double k1t  = ydot;
            double k1td = a_y(x, xdot, y, ydot, z, zdot, cs);
            double k1z  = zdot;
            double k1zd = a_z(x, xdot, y, ydot, z, zdot, cs);

            double k2r  = xdot + 0.5 * k1rd * dt;
            double k2rd = a_x(x + 0.5 * k1r * dt, xdot + 0.5 * k1rd * dt, y + 0.5 * k1t * dt, ydot + 0.5 * k1td * dt, z + 0.5 * k1z * dt, zdot + 0.5 * k1zd * dt, cs);
            double k2t  = ydot + 0.5 * k1td * dt;
            double k2td = a_y(x + 0.5 * k1r * dt, xdot + 0.5 * k1rd * dt, y + 0.5 * k1t * dt, ydot + 0.5 * k1td * dt, z + 0.5 * k1z * dt, zdot + 0.5 * k1zd * dt, cs);
            double k2z  = zdot + 0.5 * k1zd * dt;
            double k2zd = a_z(x + 0.5 * k1r * dt, xdot + 0.5 * k1rd * dt, y + 0.5 * k1t * dt, ydot + 0.5 * k1td * dt, z + 0.5 * k1z * dt, zdot + 0.5 * k1zd * dt, cs);

            double k3r  = xdot + 0.5 * k2rd * dt;
            double k3rd = a_x(x + 0.5 * k2r * dt, xdot + 0.5 * k2rd * dt, y + 0.5 * k2t * dt, ydot + 0.5 * k2td * dt, z + 0.5 * k2z * dt, zdot + 0.5 * k2zd * dt, cs);
            double k3t  = ydot + 0.5 * k2td * dt;
            double k3td = a_y(x + 0.5 * k2r * dt, xdot + 0.5 * k2rd * dt, y + 0.5 * k2t * dt, ydot + 0.5 * k2td * dt, z + 0.5 * k2z * dt, zdot + 0.5 * k2zd * dt, cs);
            double k3z  = zdot + 0.5 * k2zd * dt;
            double k3zd = a_z(x + 0.5 * k2r * dt, xdot + 0.5 * k2rd * dt, y + 0.5 * k2t * dt, ydot + 0.5 * k2td * dt, z + 0.5 * k2z * dt, zdot + 0.5 * k2zd * dt, cs);

            double k4r  = xdot + k3rd * dt;
            double k4rd = a_x(x + k3r * dt, xdot + k3rd * dt, y + k3t * dt, ydot + k3td * dt, z + k3z * dt, zdot + k3zd * dt, cs);
            double k4t  = ydot + k3td * dt;
            double k4td = a_y(x + k3r * dt, xdot + k3rd * dt, y + k3t * dt, ydot + k3td * dt, z + k3z * dt, zdot + k3zd * dt, cs);
            double k4z  = zdot + k3zd * dt;
            double k4zd = a_z(x + k3r * dt, xdot + k3rd * dt, y + k3t * dt, ydot + k3td * dt, z + k3z * dt, zdot + k3zd * dt, cs);

            // Standard RK4 update
            x     += (k1r  + 2.0 * k2r  + 2.0 * k3r  + k4r ) * (dt / 6.0);
            xdot  += (k1rd + 2.0 * k2rd + 2.0 * k3rd + k4rd) * (dt / 6.0);
            y     += (k1t  + 2.0 * k2t  + 2.0 * k3t  + k4t ) * (dt / 6.0);
            ydot  += (k1td + 2.0 * k2td + 2.0 * k3td + k4td) * (dt / 6.0);
            zsub  += (k1z  + 2.0 * k2z  + 2.0 * k3z  + k4z ) * (dt / 6.0);
            zdsub += (k1zd + 2.0 * k2zd + 2.0 * k3zd + k4zd) * (dt / 6.0);

            // Convert Cartesian result back to cylindrical for the state struct
            current_state.r        = pow(x * x + y * y, 0.5);
            current_state.theta    = atan(y / x);
            current_state.z        = zsub;
            current_state.rdot     = xdot * cos(current_state.theta) + y * sin(current_state.theta);
            current_state.zdot     = zdsub;
            current_state.thetadot = (ydot * cos(current_state.theta) - xdot * sin(current_state.theta)) / current_state.r;

        } else {
            // ---- (B) CYLINDRICAL RK4: the default path ----
            // k1 = derivatives at the current step
            double k1r  = rdot;
            double k1rd = fr(r, rdot, thetadot, z, zdot, cs);
            double k1t  = thetadot;
            double k1td = ftheta(r, rdot, theta, thetadot, z, zdot, cs);
            double k1z  = zdot;
            double k1zd = fz(r, rdot, thetadot, z, zdot, cs);

            // k2 = derivatives at midpoint using k1
            double k2r  = rdot + 0.5 * k1rd * dt;
            double k2rd = fr(r + 0.5 * k1r * dt, rdot + 0.5 * k1rd * dt, thetadot + 0.5 * k1td * dt, z + 0.5 * k1z * dt, zdot + 0.5 * k1zd * dt, cs);
            double k2t  = thetadot + 0.5 * k1td * dt;
            double k2td = ftheta(r + 0.5 * k1r * dt, rdot + 0.5 * k1rd * dt, theta + 0.5 * k1t * dt, thetadot + 0.5 * k1td * dt, z + 0.5 * k1z * dt, zdot + 0.5 * k1zd * dt, cs);
            double k2z  = zdot + 0.5 * k1zd * dt;
            double k2zd = fz(r + 0.5 * k1r * dt, rdot + 0.5 * k1rd * dt, thetadot + 0.5 * k1td * dt, z + 0.5 * k1z * dt, zdot + 0.5 * k1zd * dt, cs);

            // k3 = derivatives at midpoint using k2
            double k3r  = rdot + 0.5 * k2rd * dt;
            double k3rd = fr(r + 0.5 * k2r * dt, rdot + 0.5 * k2rd * dt, thetadot + 0.5 * k2td * dt, z + 0.5 * k2z * dt, zdot + 0.5 * k2zd * dt, cs);
            double k3t  = thetadot + 0.5 * k2td * dt;
            double k3td = ftheta(r + 0.5 * k2r * dt, rdot + 0.5 * k2rd * dt, theta + 0.5 * k2t * dt, thetadot + 0.5 * k2td * dt, z + 0.5 * k2z * dt, zdot + 0.5 * k2zd * dt, cs);
            double k3z  = zdot + 0.5 * k2zd * dt;
            double k3zd = fz(r + 0.5 * k2r * dt, rdot + 0.5 * k2rd * dt, thetadot + 0.5 * k2td * dt, z + 0.5 * k2z * dt, zdot + 0.5 * k2zd * dt, cs);

            // k4 = derivatives at the end-point
            double k4r  = rdot + k3rd * dt;
            double k4rd = fr(r + k3r * dt, rdot + k3rd * dt, thetadot + k3td * dt, z + k3z * dt, zdot + k3zd * dt, cs);
            double k4t  = thetadot + k3td * dt;
            double k4td = ftheta(r + k3r * dt, rdot + k3rd * dt, theta + k3t * dt, thetadot + k3td * dt, z + k3z * dt, zdot + k3zd * dt, cs);
            double k4z  = zdot + k3zd * dt;
            double k4zd = fz(r + k3r * dt, rdot + k3rd * dt, thetadot + k3td * dt, z + k3z * dt, zdot + k3zd * dt, cs);

            // Standard RK4 weighted sum: (k1 + 2*k2 + 2*k3 + k4) * dt/6
            current_state.r        += (k1r  + 2.0 * k2r  + 2.0 * k3r  + k4r ) * (dt / 6.0);
            current_state.rdot     += (k1rd + 2.0 * k2rd + 2.0 * k3rd + k4rd) * (dt / 6.0);
            current_state.theta    += (k1t  + 2.0 * k2t  + 2.0 * k3t  + k4t ) * (dt / 6.0);
            current_state.thetadot += (k1td + 2.0 * k2td + 2.0 * k3td + k4td) * (dt / 6.0);
            current_state.z        += (k1z  + 2.0 * k2z  + 2.0 * k3z  + k4z ) * (dt / 6.0);
            current_state.zdot     += (k1zd + 2.0 * k2zd + 2.0 * k3zd + k4zd) * (dt / 6.0);
        }

        // Keep theta in [0, 2*pi) for cleaner plotting
        if      (current_state.theta > 2.0 * PI) current_state.theta -= 2.0 * PI;
        else if (current_state.theta < 0.0)      current_state.theta += 2.0 * PI;

        // Running time-weighted averages of jet power and accretion rate
        p_avg   = (p_avg   * t + pbz   * dt) / (t + dt);
        acc_avg = (acc_avg * t + m_bhl * dt) / (t + dt);

        // ---------------------------------------------------------------
        // Log every quantity of interest for this step.
        // ---------------------------------------------------------------
        fprintf(fp,     "%E\t%E\n",             t, current_state.r);
        fprintf(fm,     "%E\t%E\n",             t, E);
        fprintf(fv,     "%E\t%E\n",             t, z);
        fprintf(flz,    "%E\t%E\n",             t, L_z);
        fprintf(fpr,    "%E\t%E\n",             t, dErdt);
        fprintf(fpt,    "%E\t%E\n",             t, dEtdt);
        fprintf(fpz,    "%E\t%E\t%E\n",         t, current_state.z, current_state.zdot);
        fprintf(fmach,  "%E\t%E\n",             t, machv);
        fprintf(forb,   "%E\t%E\t%E\t%E\n",     t, current_state.r, current_state.theta, current_state.z);
        fprintf(flinc,  "%E\t%E\n",             current_state.r, current_state.z);
        fprintf(fmdot,  "%E\t%E\n",             t, m_bhl * (3.15E7) / (2E30));  // solar masses / yr
        fprintf(fmdote, "%E\t%E\n",             t, mdot_bh2 / m_edd);
        fprintf(fflux,  "%E\t%E\n",             t, flux);
        fprintf(fpbz,   "%E\t%E\n",             t, pbz * 1E7);                   // erg/s (CGS)
        fprintf(ffluxd, "%E\t%E\n",             t, flux_dim);
        fprintf(fmass,  "%E\t%E\n",             t, M2);
        fprintf(fmass1, "%E\t%E\n",             t, M1);
        fprintf(fq,     "%E\t%E\n",             t, M2 / M1);
        fprintf(fpbz1,  "%E\t%E\n",             t, pbz1 * 1E7);
        fprintf(fpbzavg,"%E\t%E\n",             t, p_avg);
        fprintf(fmdot1, "%E\t%E\n",             t, m1_b / m1_edd);
        fprintf(faccavg,"%E\t%E\n",             t, acc_avg / m_edd);
        fprintf(feacc,  "%E\t%E\n",             t, e_acc2);
        fprintf(fejet,  "%E\t%E\n",             t, e_jet2);

        // ---------------------------------------------------------------
        // Adaptive time step: shrink dt as the orbit gets tight, so we
        // resolve each period with at least ~100-1000 steps. The pure
        // circular-orbit period is 2*pi*r/v_c.
        // ---------------------------------------------------------------
        if (current_state.r < 1.0 * pc_m && dt > (total_time / 1E8)) {
            dt = 0.001 * 2 * M_PI * current_state.r / abs(vc(r, z));
        } else if (dt > (total_time / 1E8)) {
            dt = 0.01 * 2 * PI * current_state.r / abs(vc(r, z));
        } else {
            dt = dt0;
        }

        t += dt;
    }

    // -------------------------------------------------------------------
    // Post-loop bookkeeping: normalize running averages, catch numerical
    // escapes (NaN), and write the one-line run summary.
    // -------------------------------------------------------------------
    machavg = machavg / t;
    tauavg  = tauavg  / t;
    printf("Average Z Mach: %E\n Time: %E\n", machavg, t / (1E9 * yr_s));

    if (current_state.r < 0.0) flyby = 1;
    if (flyby == 2)      t = 14.0 * (1E9 * yr_s);   // flyby: set t to a Hubble time
    if (t != t)          t = 14.0 * (1E9 * yr_s);   // NaN guard

    // Run-summary line:  M1  fg  vg  n_gd[cm^-3]  inclination  t_decay[Gyr]  outcome_flag
    fprintf(ftime, "%E\t%E\t%E\t%E\t%E\t%E\t%d\n", M1, fg, vg, n_gd / (1E6), angle, t / (1E9 * yr_s), flyby);
    fprintf(fmachavg, "%E\t%E\t%E\n", angle, t, machavg);
    fprintf(flog, "%E\t%E\n", initial_state.z / z_thick, tauavg / (1E9 * yr_s));
    fprintf(facc_logmax, "%E\t%E\n", acc_max2 * (yr_s / sl_kg), acc_max2_edd);
    return current_state;
}

// ===========================================================================
// MAIN: set up the galaxy, build the potential, kick off the orbit integrator.
//
// Command-line argument: initial orbital inclination in degrees.
// The secondary is placed at a starting radius r0 on a circular orbit
// rotated up by that inclination angle.
// ===========================================================================
int main(int argc, char *argv[]) {

    // Check usage: exactly one argument (inclination in degrees)
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <inclination_in_degrees>\n";
        return 1;
    }

    // `loaded = true`  -> read the pre-computed potential from pot_test.csv
    //                     and pot_test2.csv (fast, requires those files exist)
    // `loaded = false` -> solve Poisson from scratch via the relaxation
    //                     method in poisson.h (slow; use the first time you
    //                     change the galaxy parameters)
    bool loaded = true;

    float input = std::stof(argv[1]);
    printf("Inclination: %E\n", input);

    // Turn off GSL's abort-on-error behavior; guards inside the code handle
    // out-of-bounds interpolation gracefully.
    gsl_set_error_handler_off();

    // -------------------------------------------------------------------
    // GALAXY GRID ALLOCATION
    //
    // Outer grid: N x N, covering an L x L square centered on the galaxy.
    // AMR grid:   N x N2, covering L x L2 with L2 << L (disk midplane zoom-in).
    // -------------------------------------------------------------------
    int N  = 200;
    int N2 = 200;
    r_arr  = new double[N];   z_arr  = new double[N];
    r2_arr = new double[N];   z2_arr = new double[N];

    // Convenience pointers that alias the global 2D arrays (used when passing
    // them to rhos_init / rhog_init / etc. below).
    double *a[N];    double *a2[N2];   // rho
    double *b[N];    double *b2[N2];   // phi
    double *c[N];    double *c2[N2];   // rho_bg (bulge)
    double *rg[N];   double *rg2[N2];  // rho_g  (gas)
    double *rs[N];   double *rs2[N2];  // rho_s  (stellar disk)
    double *rst[N];  double *rst2[N2]; // rho_st (stellar total)
    double *rh[N];                      // rho_h  (halo; outer only)
    double *csa[N];                     // cs_arr (sound-speed diagnostic)

    // Allocate the global 2D grids
    phi     = new double*[N];   phi2     = new double*[N2];
    rho     = new double*[N];   rho2     = new double*[N2];
    rho_s   = new double*[N];   rho_s2   = new double*[N2];
    rho_st  = new double*[N];   rho_st2  = new double*[N2];
    rho_g   = new double*[N];   rho_g2   = new double*[N2];
    rho_bg  = new double*[N];   rho_bg2  = new double*[N2];
    rho_h   = new double*[N];
    phi_b   = new double*[N];
    cs_arr  = new double*[N];

    for (size_t i = 0; i < N; ++i) {
        phi[i]    = new double[N];   rho[i]    = new double[N];
        rho_s[i]  = new double[N];   rho_st[i] = new double[N];
        rho_g[i]  = new double[N];   rho_bg[i] = new double[N];
        rho_h[i]  = new double[N];   phi_b[i]  = new double[N];
        cs_arr[i] = new double[N];
    }
    for (size_t i = 0; i < N2; ++i) {
        phi2[i]    = new double[N2];  rho2[i]   = new double[N2];
        rho_s2[i]  = new double[N2];  rho_st2[i] = new double[N2];
        rho_bg2[i] = new double[N2];  rho_g2[i]  = new double[N2];
    }

    // Zero-initialize everything
    for (int i = 0; i < N; i++) {
        r_arr[i] = 0.0; z_arr[i] = 0.0;
        for (int j = 0; j < N; j++) {
            phi[i][j] = 0.0; rho[i][j] = 0.0; rho_s[i][j] = 0.0; rho_st[i][j] = 0.0;
            rho_g[i][j] = 0.0; rho_h[i][j] = 0.0; phi_b[i][j] = 0.0;
            rho_bg[i][j] = 0.0; cs_arr[i][j] = 0.0;
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

    // Alias globals into the local pointer arrays that the init routines take
    for (int i = 0; i < N; i++) {
        a[i]  = rho[i];     rg[i]  = rho_g[i];   rs[i]  = rho_s[i];
        rst[i]= rho_st[i];  rh[i]  = rho_h[i];   b[i]   = phi[i];
        c[i]  = rho_bg[i];  csa[i] = cs_arr[i];
    }
    for (int i = 0; i < N2; i++) {
        a2[i]  = rho2[i];   rg2[i] = rho_g2[i];  rs2[i] = rho_s2[i];
        rst2[i]= rho_st2[i];b2[i]  = phi2[i];    c2[i]  = rho_bg2[i];
    }

    // -------------------------------------------------------------------
    // Grid physical extents. L is the outer-grid box side; L2 is the
    // AMR-grid vertical extent (much smaller -> resolves disk midplane).
    // -------------------------------------------------------------------
    double L  = 2.0 * R_gd;    // outer box:  -R_gd .. +R_gd (and similarly in z)
    double L2 = 0.1 * R_gd;    // AMR vertical half-extent (thin slab)
    double h  = L  / (N  - 1); // outer grid spacing
    double h2 = L2 / (N2 - 1); // AMR z spacing

    // -------------------------------------------------------------------
    // Compute galaxy-wide mass normalizations by integrating the analytic
    // density profiles (see poisson.h). We do this so that the on-grid
    // densities correctly reproduce the desired gas mass, stellar mass,
    // bulge mass, etc.
    //
    //   dmass    = total gas disk mass
    //   dst_mass = total stellar disk mass  (from dmass & gas fraction fg)
    //   sigma_d  = stellar-disk central surface-density normalization
    //   rhob0    = bulge central density  (so that M_bulge ~ 1000 * M1 [solar])
    // -------------------------------------------------------------------
    std::vector<double> dmass_params = {1.0 * rho_gd(0.0), R_gd, z_thick};
    double dmass    = 2 * integrate2D(rho_gas_f, 0.0, 10.0 * R_gd, 0.0, 10.0 * R_gd, 200, 200, dmass_params);
    double bmassco  = 2 * integrate2D(rho_bulge_cof, 0.0, 2.0, 0.0, 2.0, 200, 200, {h / (1000.0 * pc_m), 1.0});
    double dst_mass = (1.0 / fg) * (dmass - fg * dmass);

    double dst_co = 2 * integrate2D(rho_stellar_f, 0.0, 10.0 * R_sd, 0.0, 10.0 * R_sd, 200, 200, {1.0, R_sd, 0.1 * R_sd, 0.3333 * R_sd});
    double sigma_d = dst_mass / dst_co;
    double sdmass  = 2.0 * integrate2D(rho_stellar_f, 0.0, 10.0 * R_sd, 0.0, 10.0 * R_sd, 200, 200, {sigma_d, R_sd, 0.1 * R_sd, 0.3333 * R_sd});

    double rho0  = 1.0 * rho_gd(0.0);                                    // central gas density
    double rhob0 = 1.0 * (1000.0 * M1 / bmassco) * sl_kg / pow(1000.0 * pc_m, 3);  // bulge central density
    double k     = 1.0 / R_gd;
    double f_tot = dmass / (dst_mass + 1000 * M1 * sl_kg + dmass);

    // Stellar-disk profile parameters: {sigma_d, R_sd, z_thin, z_thick}
    // The double-exponential profile is a sum of a thin + thicker component.
    std::vector<double> params_std  = {sigma_d, R_sd, 0.1 * R_sd, R_sd / 3.0};
    std::vector<double> params_halo = {1.0 * 4.839E-20, 1.0 * 4.21E19};

    // Populate the coordinate axes
    r_init(N, L,  r_arr);   z_init(N,  L,  z_arr);
    r_init(N, L, r2_arr);   z_init(N2, L2, z2_arr);

    // -------------------------------------------------------------------
    // Populate each component density grid:
    //   rhos_init   -> stellar disk (two-component exp)
    //   rhog_init   -> exponential gas disk
    //   bulge_init  -> cusp-cutoff bulge
    // Uncomment `halo_init` to include a dark matter halo.
    // -------------------------------------------------------------------
    rhos_init(N, L, rst, {0.0, R_sd, z_thick, 3.0 * z_thick});  // placeholder total-star (gets overwritten below)
    rhos_init(N, L, rs,  params_std);                            // stellar disk
    rhog_init(N, L, rg,  R_gd, 1.0 * z_thick, 1.0 * rho0);       // gas disk
    //halo_init(N, L, rh, params_halo);                          // (disabled) dark matter halo
    bulge_init(N, L, c,  1.0 * rhob0, 0.6, R_sd);                // bulge

    // Same on the AMR grid
    rhos_init2(N, N2, L, L2, rs2, params_std);
    rhog2_init(N, N2, L, L2, rg2, R_gd, 1.0 * z_thick, 1.0 * rho0);
    bulge2_init(N, N2, L, L2, c2, rhob0, 0.6, R_sd);

    // -------------------------------------------------------------------
    // Assemble the TOTAL density grids from the components.
    // rho_st = stellar (disk + bulge) -- used in bulge DF
    // rho    = gas + stars + bulge (+ halo, if enabled) -- source of phi
    //
    // *** HALO / COMPONENT SWITCH ***  Change the multipliers below to turn
    // individual components on/off in the Poisson source term. The `0.0 *
    // rho_h[i][j]` below is the HALO SWITCH: set to 1.0 to include the halo.
    // -------------------------------------------------------------------
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0 || j == 0 || i == N-1 || j == N-1) {
                // Zero out boundary rows/cols (required for vacuum BCs in Poisson solve)
                rho[i][j] = 0.0; rho_g[i][j] = 0.0; rho_s[i][j] = 0.0;
                rho_bg[i][j] = 0.0; rho_st[i][j] = 0.0; rho_h[i][j] = 0.0;
            } else {
                rho_st[i][j] += (rho_bg[i][j] + rho_s[i][j]);
                // HALO SWITCH here --------------------------------v
                rho[i][j] += (1.0 * rho_bg[i][j] + 1.0 * rho_g[i][j] + 1.0 * rho_s[i][j] + 0.0 * rho_h[i][j]);
            }
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N2; j++) {
            if (i == 0 || j == 0 || i == N-1 || j == N2-1) {
                rho2[i][j] = 0.0; rho_g2[i][j] = 0.0; rho_s2[i][j] = 0.0;
                rho_bg2[i][j] = 0.0; rho_st2[i][j] = 0.0;
            } else {
                rho_st2[i][j] += (rho_bg2[i][j] + rho_s2[i][j]);
                rho2[i][j]   += 1.0 * (1.0 * rho_bg2[i][j] + 1.0 * rho_g2[i][j] + 1.0 * rho_s2[i][j]);
            }
        }
    }
    
    
    // Total bulge mass in solar masses -- used in the status printout below
    double bmass = (rhob0 * pow(pc_m * 1000, 3) / sl_kg) * 2 * integrate2D(rho_bulge_cof, 0.0, 2.0, 0.0, 2.0, 200, 200, {h / (pc_m * 1000.0), 1.0});

    // Galaxy status report -- printed so you can sanity-check the setup
    printf("\nDGMass: %E\t DSMass: %E\t BGMass: %E\t SigmaD: %E\n",
           dmass / sl_kg, dst_mass / sl_kg, bmass, sigma_d);

    // ===================================================================
    // POISSON SOLVE (outer grid)
    //
    // Either load a pre-computed potential from pot_test.csv (fast), or
    // run the iterative Gauss-Seidel-like relaxation solver in poisson.h
    // (slow, ~50000 iterations). The convergence tolerance is proportional
    // to 4*pi*G*rho0*h^2, the scale of a single grid-cell's source term.
    // ===================================================================
    if (loaded) {
        load_potential_csv_ptrptr("pot_test.csv", phi, N, N, &N, &N);
    } else {
        relax(N, L, a, b, 50000, 4 * M_PI * G * rho0 * 1E-4 * h * h);
    }

    // -------------------------------------------------------------------
    // Allocate and initialize the GSL bilinear interpolation objects.
    // We use bilinear (not bicubic) specifically because bicubic can
    // produce the rare negative-density artifacts we guard against in
    // gasDFb() and friends.
    // -------------------------------------------------------------------
    interp     = gsl_interp2d_alloc(gsl_interp2d_bilinear, N, N);
    xacc       = gsl_interp_accel_alloc();
    yacc       = gsl_interp_accel_alloc();

    interp_r   = gsl_interp2d_alloc(gsl_interp2d_bilinear, N, N);
    xacc_r     = gsl_interp_accel_alloc();
    yacc_r     = gsl_interp_accel_alloc();

    interp_amr = gsl_interp2d_alloc(gsl_interp2d_bilinear, N, N);
    xacc_amr   = gsl_interp_accel_alloc();
    yacc_amr   = gsl_interp_accel_alloc();

    double drh = L  / (N  + 1);   // step used when resampling outer -> AMR grid
    double dzh = L2 / (N2 + 1);

    // Flatten each 2D array into 1D for GSL. Boundary cells are zeroed so
    // the potential / density go to zero at infinity (vacuum BC).
    phi_flat    = new double[N * N];
    rho_flat    = new double[N * N];
    rho_g_flat  = new double[N * N];
    rho_s_flat  = new double[N * N];
    rho_st_flat = new double[N * N];
    rho_h_flat  = new double[N * N];

    int i1 = 0, i2 = 0, i3 = 0, i4 = 0, i5 = 0, i6 = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0 || j == 0 || i == N-1 || j == N-1) {
                rho_flat[i1++]   = 0.0;  phi_flat[i2++]    = 0.0;
                rho_g_flat[i3++] = 0.0;  rho_s_flat[i4++]  = 0.0;
                rho_st_flat[i5++]= 0.0;  rho_h_flat[i6++]  = 0.0;
            } else {
                rho_flat[i1++]    = rho[i][j];
                phi_flat[i2++]    = phi[i][j];
                rho_g_flat[i3++]  = rho_g[i][j];
                rho_s_flat[i4++]  = rho_s[i][j];
                rho_st_flat[i5++] = rho_st[i][j];
                rho_h_flat[i6++]  = rho_h[i][j];
            }
        }
    }

    gsl_interp2d_init(interp,   r_arr, z_arr, phi_flat, N, N);
    gsl_interp2d_init(interp_r, r_arr, z_arr, rho_flat, N, N);

    // -------------------------------------------------------------------
    // Seed the AMR grid's potential by evaluating the (just-computed) outer
    // potential at the AMR grid nodes. This gives the AMR relaxation a warm
    // start so it converges much faster.
    // -------------------------------------------------------------------
    for (int i = 0; i < N2; i++) {
        for (int j = 0; j < N; j++) {
            phi2[i][j] = gsl_interp2d_eval(interp, r_arr, z_arr, phi_flat,
                                           (j - N / 2) * drh, (i - N / 2) * dzh,
                                           xacc_amr, yacc_amr);
        }
    }

    // Write out the bulge density for visual inspection (plotted separately)
    FILE* fpotb = fopen("rho.csv", "w");
    for (int i = 0; i < N2; i++) {
        for (int j = 0; j < N; j++) fprintf(fpotb, "%E,", rho_bg[i][j]);
        fprintf(fpotb, "\n");
    }

    // ===================================================================
    // POISSON SOLVE (AMR grid). Same branching as the outer solve.
    // ===================================================================
    if (loaded) {
        load_potential_csv_ptrptr("pot_test2.csv", phi2, N2, N2, &N2, &N2);
    } else {
        relax_amr(N, N2, L, L2, a2, b2, 50000, 4 * M_PI * G * rho0 * 1E-4 * h * h2);
    }

    // Flatten AMR grid for GSL
    phi_flat2    = new double[N * N2];
    rho_flat2    = new double[N * N2];
    rho_g_flat2  = new double[N * N2];
    rho_s_flat2  = new double[N * N2];
    rho_st_flat2 = new double[N * N2];

    i1 = 0; i2 = 0; i3 = 0; i4 = 0; i5 = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N2; j++) {
            if (i == 0 || j == 0 || i == N-1 || j == N2-1) {
                rho_flat2[i1++]   = 0.0;  phi_flat2[i2++]    = 0.0;
                rho_g_flat2[i3++] = 0.0;  rho_s_flat2[i4++]  = 0.0;
                rho_st_flat2[i5++]= 0.0;
            } else {
                rho_flat2[i1++]    = rho2[i][j];
                phi_flat2[i2++]    = phi2[i][j];
                rho_g_flat2[i3++]  = rho_g2[i][j];
                rho_s_flat2[i4++]  = rho_s2[i][j];
                rho_st_flat2[i5++] = rho_st2[i][j];
            }
        }
    }

    gsl_interp2d_init(interp_amr, r2_arr, z2_arr, phi_flat2, N, N);

    // Dump the AMR potential to pot_test2.csv so future runs can load it
    FILE* fpot2 = fopen("pot_test2.csv", "w");
    for (int i = 0; i < N2; i++) {
        for (int j = 0; j < N; j++) fprintf(fpot2, "%E,", phi2[i][j]);
        fprintf(fpot2, "\n");
    }

    // -------------------------------------------------------------------
    // Dump the equipartition / asymptotic magnetic field strength to b.csv
    // (optional diagnostic; assumes B^2/(8*pi) = rho * cs^2 in equipartition
    // with the gas pressure; scaled to Gauss).
    // -------------------------------------------------------------------
    FILE* fpotmg = fopen("b.csv", "w");
    for (int i = 0; i < N2; i++) {
        for (int j = 0; j < N; j++) {
            if (j - N / 2 == 0) {
                fprintf(fpotmg, "%E,", 0.0);
            } else {
                double b_inf = -1.0 * G * M_PI * z_thick * abs((j - N / 2) * h)
                             * rho_rz((j - N / 2) * h, (i - N / 2) * h)
                             * pow(2 * 1.26E-6 * rho_rz((j - N / 2) * h, (i - N / 2) * h), 0.5)
                             / (vg * vc((j - N / 2) * h, (i - N / 2) * h));
                fprintf(fpotmg, "%E,", b_inf * 1E4);   // x 1e4 converts T -> Gauss
            }
        }
        fprintf(fpotmg, "\n");
    }

    // -------------------------------------------------------------------
    // Pre-compute a sound-speed grid cs_arr[][] across the whole domain.
    // (Not currently consulted in the RK4 loop -- retained as a diagnostic
    // and for any future use that wants a full cs map.)
    // -------------------------------------------------------------------
    double csm = 0.0, Tm = 0.0;
    for (int i = 0; i < N / 2; i++) {
        double rt = i * h;
        for (int j = 0; j < N / 2; j++) {
            double zt = j * h;
            csm = (G * M_PI * 2 * rho_rz(rt, zt) * z_thick) / ((2 * abs(vg) * vc(rt, zt)) / rt);
            Tm  = (3 * m_p * csm * csm) / (5 * kb);
            Tm += 1E4;
            csm = pow(5 * kb * Tm / (3 * m_p), 0.5);
            // Mirror into all four quadrants (the profile is reflection-symmetric)
            cs_arr[N/2+j][N/2+i] = csm;   cs_arr[N/2-j][N/2+i] = csm;
            cs_arr[N/2+j][N/2-i] = csm;   cs_arr[N/2-j][N/2-i] = csm;
        }
    }

    // Dump the outer-grid potential to pot_test.csv for future loading
    FILE* fpot = fopen("pot_test.csv", "w");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) fprintf(fpot, "%E,", phi[i][j]);
        fprintf(fpot, "\n");
    }

    // ===================================================================
    // INITIAL CONDITIONS FOR THE ORBIT
    //
    // We place the secondary at 1 kpc on a near-circular orbit. The "f"
    // factor slightly boosts the angular velocity above v_c/r (giving an
    // apocentric orbit at r0) to avoid zero-eccentricity edge cases.
    //
    // The command-line inclination angle then rotates the orbit out of
    // the midplane: r0 -> r0*cos(theta),  z0 -> r0*sin(theta).
    // ===================================================================
    double cs0 = pow(5 * kb * 1E4 / (3 * m_p), 0.5);    // reference gas cs at 1e4 K
    double r0  = 1.00 * 1000.0 * pc_m;                   // initial cylindrical radius
    double rd0 = 0.0;                                    // initial rdot
    double t0  = 0.0;                                    // initial theta
    double z0  = 0.0 * pc_m;                             // initial z
    double zd0 = 0.0 * cs0;                              // initial zdot
    double f   = 0.5;                                    // eccentricity-like factor
    double v_c = vc(r0, z0);                             // local circular speed
    double td0 = 1.0 * v_c * pow(1 + f, 0.5) / r0;       // initial thetadot

    // Rotate the initial position upward by the user-supplied inclination angle
    double v_0     = r0 * td0;                           // initial tangential speed
    double s_0     = pow(r0 * r0 + z0 * z0, 0.5);        // total initial radius
    double theta_z = (input / 180) * M_PI;               // inclination in radians
    r0  = s_0 * cos(theta_z);
    z0  = s_0 * sin(theta_z);
    td0 = 1.0 * v_0 / r0;                                // keep tangential speed constant

    printf("%E\t%E\n", r0 / (1000.0 * pc_m), z0 / (1000.0 * pc_m));

    // Total simulation time and base time step
    double total_time = 1.4E10 * yr_s;   // 14 Gyr ~ Hubble time
    double n          = 1000000.0;        // target number of steps
    double time_step  = total_time / n;   // nominal dt (adapted near center)

    // Launch the integrator
    State initial_state{ r0, rd0, t0, td0, z0, zd0 };
    State final_state  = RK4Solver(initial_state, time_step, total_time,
                                   180.0 * theta_z / M_PI);

    printf("\nTerminated at radius %E\n",
           pow(final_state.r * final_state.r + final_state.z * final_state.z, 0.5) / pc_m);

    return 0;
}
