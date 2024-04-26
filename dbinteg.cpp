#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include "integration.h"
#include "poisson.h"

// Function to integrate
double function(double x, double y, const std::vector<double>& params) {
    // Example function: f(x, y) = a*x^2 + b*y^2
    double result = 0.0;
    if (params.size() >= 2) {
        result = params[0] * x * x + params[1] * y * y;
    }
    return result;
}

double function1(double r, double z, const std::vector<double>& params) {
    double fx = 0.0;
    double m = 0.0;
    if (params.size() >= 2) {
        m = pow(r * r + z * z / (params[0] * params[0]), 0.5);
        fx = 2*M_PI*r*pow(M_E, -1.0 * m * m / (params[1] * params[1]))*pow(m/params[2], -1.8);
    }
    return fx;
}

int main() {
    //double kpc_m = 3.086E19;
    double ax = 0.00, bx = 1.0; // Integration limits for x
    double ay = -10.0, by = 10.0; // Integration limits for y
    int nx = 100000; // Number of subdivisions for x
    int ny = 100; // Number of subdivisions for y
    double rgas = 300.0 * 1E6 * 1.67E-27;
/*
    // Parameters for the function
    std::vector<double> params = {0.6, 1.9, 1.0}; // Example coefficients a and b

    // Perform double integration
    double result = integrate2D(function1, ax, bx, ay, by, nx, ny, params);
    std::cout << "Result of double integration: " << result << std::endl;
*/
    std::vector<double> params = {1.0, 2.0, 0.1};
    double result = integrate2D(rho_gas_f, 0.0, 100.0, 0.0, 100.0, 1000, 1000, params);
    printf("%f\n", result);
    return 0;
}
