#include <iostream>
#include <gsl/gsl_interp2d.h>

int main() {
    // Define the 2D array (grid data)
    const size_t nx = 4; // number of points in x-direction
    const size_t ny = 4; // number of points in y-direction

    double x_vals[nx] = {0.0, 1.0, 2.0, 3.0}; // x-coordinates
    double y_vals[ny] = {0.0, 1.0, 2.0, 3.0}; // y-coordinates
    /*
    double f_vals[nx * ny] = {
        0.0, 1.0, 2.0, 3.0,
        4.0, 5.0, 6.0, 7.0,
        8.0, 9.0, 10.0, 11.0,
        12.0, 13.0, 14.0, 15.0
    }; // function values at grid points
    */

    double f_vals[nx][ny] = 
    {{1,2,3,4}, {2,4,6,8}, {3,6,9,12}, {4,8,12,16}};

    // Create an interpolation object
    gsl_interp2d *interp = gsl_interp2d_alloc(gsl_interp2d_bicubic, nx, ny);
    gsl_interp_accel *xacc = gsl_interp_accel_alloc();
    gsl_interp_accel *yacc = gsl_interp_accel_alloc();

    // Initialize the interpolation object
    gsl_interp2d_init(interp, x_vals, y_vals, reinterpret_cast<double*>(f_vals), nx, ny);

    // Interpolate at a specific point (e.g., (1.5, 1.5))
    double x = 1.5;
    double y = 0.0;
    double z;
    z = gsl_interp2d_eval(interp, x_vals, y_vals, reinterpret_cast<double*>(f_vals), x, y, xacc, yacc);

    // Output the interpolated value
    std::cout << "Interpolated value at (" << x << ", " << y << "): " << z << std::endl;

    // Free allocated memory
    gsl_interp2d_free(interp);
    gsl_interp_accel_free(xacc);
    gsl_interp_accel_free(yacc);

    return 0;
}
